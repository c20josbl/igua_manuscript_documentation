import argparse
import csv
import os
import io
import itertools
import tarfile
import glob
import re
import multiprocessing.pool
import sys

import Bio.SeqIO
import gecco.model
import gecco.hmmer
import gecco.crf
import pandas
import numpy
import scipy.sparse
import rich
import rich.progress
import tqdm


class TrackedFile(io.RawIOBase):

    def __init__(self, handle, progress, task, scale=0):
        self.handle = handle
        self.progress = progress
        self.task = task
        self.scale = scale

    def __enter__(self):
        self.handle.__enter__()
        return self

    def __exit__(self, exc_val, exc_ty, tb):
        self.handle.__exit__(exc_val, exc_ty, tb)
        return False

    def readable(self):
        return True

    def read(self, size=-1):
        block = self.handle.read(size)
        self.progress.update(self.task, advance=len(block) / (1024 ** self.scale))
        return block

    def close(self):
        self.handle.close()



GBK_RX = re.compile(r".*region(\d{3})\.gbk$")

RECORDS = [
    "proG3_total_final.gbk",
]


parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--pfam", required=True)
args = parser.parse_args()

pfam = gecco.hmmer.HMM("Pfam", version="35.0", url="", path=args.pfam, size=19632, relabel_with=r"s/(PF\d+).\d+/\1/")

HMMS = [ pfam ]


with rich.progress.Progress(
        rich.progress.SpinnerColumn(finished_text="[green]:heavy_check_mark:[/]"),
        "[progress.description]{task.description}",
        rich.progress.BarColumn(bar_width=60),
        rich.progress.DownloadColumn(),
        "[progress.percentage]{task.percentage:>3.0f}%",
        rich.progress.TimeElapsedColumn(),
        rich.progress.TimeRemainingColumn(),
) as progress:

    progress.console.print("[bold green]{:>12}[/] GECCO internal CRF".format("Loading"))
    crf = gecco.crf.ClusterCRF.trained()

    progress.console.print("[bold green]{:>12}[/] GCF representatives from {!r}".format("Loading", args.input))
    gcfs = pandas.read_table(args.input, comment="#")
    representative_ids = set(gcfs.gcf_representative)
    bgc_to_gcf_index = {row.gcf_representative:row.gcf_id for row in gcfs.itertuples()}

    # Extract clusters from antiSMASH results

    def region_to_cluster(record, bgc_id):
        genes = []
        for cds in filter(lambda f: f.type == "CDS", record.features):
            qualifiers = cds.qualifiers
            if "locus_tag" in qualifiers:
                id_ = qualifiers["locus_tag"][0]
            elif "protein_id" in qualifiers:
                id_ = qualifiers["protein_id"][0]
            elif "gene" in qualifiers:
                id_ = qualifiers["gene"][0]
            else:
                progress.print(cds); exit(1)
            seq = Bio.Seq.Seq(qualifiers['translation'][0])
            protein = gecco.model.Protein(id_, seq)
            gene = gecco.model.Gene(source=record, start=cds.location.start, end=cds.location.end, protein=protein, strand=gecco.model.Strand.Coding)
            genes.append(gene)
        return gecco.model.Cluster(bgc_id, genes)

    def record_to_cluster(record):
        genes = []
        for cds in filter(lambda f: f.type == "CDS", record.features):
            qualifiers = cds.qualifiers
            protein = gecco.model.Protein(record.name, Bio.Seq.Seq(qualifiers['translation'][0]))
            gene = gecco.model.Gene(record, start=cds.location.start, end=cds.location.end, protein=protein, strand=gecco.model.Strand.Coding)
            genes.append(gene)
        return gecco.model.Cluster(record.name, genes)

    clusters = []
    for filename in RECORDS:
        progress.console.print("[bold green]{:>12}[/] representatives from {!r}".format("Extracting", filename))
        task = progress.add_task("Reading")
        with progress.open(filename, "rb", task_id=task) as f:
            for record in Bio.SeqIO.parse(filename, "genbank"):
                if record.name.startswith("BGC"):
                    id_ = record.name.rsplit(".", 1)[0]
                    if id_ in representative_ids:
                        clusters.append(region_to_cluster(record, id_))
                elif record.name in representative_ids:
                    clusters.append(record_to_cluster(record))

with rich.progress.Progress(
        rich.progress.SpinnerColumn(finished_text="[green]:heavy_check_mark:[/]"),
        "[progress.description]{task.description}",
        rich.progress.BarColumn(bar_width=60),
        "[progress.completed]{task.completed:.1f}/{task.total:.1f}",
        "[progress.completed]{task.fields[unit]}",
        "[progress.percentage]{task.percentage:>3.0f}%",
        rich.progress.TimeElapsedColumn(),
        rich.progress.TimeRemainingColumn(),
) as progress:


    # Annotate with Pfam

    progress.console.print("[bold green]{:>12}[/] genes of {} BGCs".format("Annotating", len(clusters)))
    genes = [gene for cluster in clusters for gene in cluster.genes]
    for hmm in HMMS:
        hmmer = gecco.hmmer.PyHMMER(hmm, cpus=28)
        task = progress.add_task(hmm.id, total=hmm.size, unit="HMM")
        hmmer.run(genes, progress=lambda h,i: progress.update(task, advance=1))

    progress.console.print("[bold green]{:>12}[/] domains by p-value under 1e-5".format("Filtering"))
    genes = [
        gene.with_protein(gene.protein.with_domains([d for d in gene.protein.domains if d.pvalue < 1e-5]))
        for gene in genes
    ]

    progress.console.print("[bold green]{:>12}[/] genes by source BGC".format("Sorting"))
    genes.sort(key=lambda g: g.source.name)

    progress.console.print("[bold green]{:>12}[/] genes into source BGCs".format("Grouping"))
    all_clusters = []
    for bgc_id, cluster_genes in itertools.groupby(genes, lambda g: g.source.name):
        id_ = bgc_id.rsplit(".", 1)[0] if bgc_id.startswith("BGC") else bgc_id
        cluster = gecco.model.Cluster(id=id_, genes=list(cluster_genes))
        all_clusters.append(cluster)


    # --- Finalize

    progress.console.print("[bold green]{:>12}[/] {} clusters by GCF ID".format("Sorting", len(all_clusters)))
    all_clusters.sort(key=lambda cluster: bgc_to_gcf_index[cluster.id])

    progress.console.print("[bold green]{:>12}[/] all possible domains from clusters".format("Extracting"))
    all_possible = sorted({domain.name for cluster in all_clusters for gene in cluster.genes for domain in gene.protein.domains})

    os.makedirs(args.output, exist_ok=True)
    progress.console.print("[bold green]{:>12}[/] labels and domain compositions".format("Writing"))
    with open(os.path.join(args.output, "labels.tsv"), "w") as f:
        for cluster in all_clusters:
            f.write(cluster.id)
            f.write("\t")
            f.write(bgc_to_gcf_index[cluster.id])
            f.write("\n")

    with open(os.path.join(args.output, "domains.tsv"), "w") as f:
        for domain in all_possible:
            f.write(domain)
            f.write("\n")

    comp = numpy.array([c.domain_composition(all_possible) for c in all_clusters])
    scipy.sparse.save_npz(os.path.join(args.output, "compositions.npz"), scipy.sparse.coo_matrix(comp))

    # --- Check missing
    all_clusters_set = { cluster.id for cluster in all_clusters }
    progress.console.print("[bold red]{:>12}[/] compositions for {} BGCs".format("Missing", len(representative_ids.difference(all_clusters_set))))
    #for bgc_id in representative_ids.difference(all_clusters_set):
    #    progress.console.print("[bold red]{:>12}[/] compositions for {}".format("Missing", bgc_id))

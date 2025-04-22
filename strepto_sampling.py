from Bio import SeqIO
import random

# Function to perform random sampling from GenBank file
def random_sample(input_file, output_file, num_records, sample_size):
    with open(input_file, "r") as handle:
        records = list(SeqIO.parse(handle, "genbank"))

    sampled_records = random.sample(records, k=num_records)

    with open(output_file, "w") as out_handle:
        SeqIO.write(sampled_records, out_handle, "genbank")

# Input and output file paths
input_file = "cat_strep_mibig.gbk"
output_dir = "./sampled_genbank_files/"

# Number of records to sample
num_records_list = [2500, 2500,2500]

# Sample sizes for each run
sample_sizes = [2500]*3

# Perform random sampling for each combination
for i, (num_records, sample_size) in enumerate(zip(num_records_list, sample_sizes)):
    output_file = f"{output_dir}sample_{i+10}_size_{sample_size}.gb"
    random_sample(input_file, output_file, num_records, sample_size)

print("Sampling completed successfully.")

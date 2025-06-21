
import sys, subprocess
import pandas as pd

genome_info_file = sys.argv[1]
otu_range_file = sys.argv[2]

# check
otu_range = pd.read_csv(otu_range_file, sep="\t", header=None)
num_columns = otu_range.shape[1]
if num_columns != 4:
    genome_info = pd.read_csv(genome_info_file, sep="\t", dtype=object)

    pangenome_species_eq1 = genome_info["species_taxid"][genome_info["species_taxid"].map(genome_info["species_taxid"].value_counts()) == 1].tolist()
    otu_range_file_tokens = otu_range_file.strip().split("/")
    old_otu_range_file_name = otu_range_file_tokens[-1].split(".")[0] + "_old.txt"
    old_otu_range_file_tokens = otu_range_file_tokens[:-1] + [old_otu_range_file_name]
    old_otu_range_file = "/".join(old_otu_range_file_tokens)
    subprocess.run(f"mv {otu_range_file} {old_otu_range_file}", shell=True)

    with open(old_otu_range_file, "r") as f_in, open(otu_range_file, "w") as f_out:
        for line in f_in:
            tokens = line.strip().split("\t")
            if tokens[0] in pangenome_species_eq1:
                tokens.append("0")
            else:
                tokens.append("1")
            f_out.write("\t".join(tokens) + "\n")
# else:
#     print("This format don't need to modify")




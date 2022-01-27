import pickle
from pathlib import Path

from acc2homologsDict import acc2sequence

p = Path("./cache")
filtered_dict = p / "homologs_dict_filtered.pkl"
fasta_tmp = p / "homologs.fasta"

with open(filtered_dict, mode="rb") as f:
    homolog_groups = pickle.load(f)

    for group in homolog_groups:

        pident = group["pident"]
        data = group["data"]

        indiv_fastas = [acc2sequence(homolog["accession"]) for homolog in data]

        group_fasta = "".join(indiv_fastas)

        with open(f"./cache/{str(pident)}_percent_homologs.fasta", mode="w+") as out:
            out.write(group_fasta)
        print("fasta file for "+str(pident)+" percent homologs cached")


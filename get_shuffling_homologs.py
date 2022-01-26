from Bio.Blast.NCBIWWW import qblast
from Bio.Blast.NCBIXML import read, parse

import requests
import time
from pathlib import Path

p = Path("./cache")
blast_cache = p / "blast_tmp.xml"


def acc2sequence(accID):
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/?db=protein&id="+accID+"&rettype=fasta"
    response = requests.get(URL)
    if response.ok:
        return response.text
    else:
        print("bad request")
        print(response.status_code)



def blast(sequence):
    
    with open('cache/blast_tmp.xml', mode="w+") as f:
        print('entering blast function')

        blastStart = time.time()
        blast_results = qblast("blastp", "nr", sequence)  
        blastEnd = time.time()

        print('finished blast. Took '+str(blastEnd - blastStart)+" seconds.")
        f.write(blast_results.read())
        print('cached blast result')   


def check_blast():
        with open(blast_cache, mode="r") as f:

            records = [record for record in parse(f)]
            for record in records:
                homologs = [
                    {"accession":alignment.accession, 
                    "identity":round((hsp.identities/hsp.align_length)*100, 2),
                    "coverage": round((hsp.align_length/record.query_length)*100, 2)
                    }
                    for alignment in record.alignments
                        for hsp in alignment.hsps
                ]
                               
            print(homologs)


if __name__ == "__main__":

    acc = "WP_000113609.1"

    #seq = acc2sequence(acc)

    #blast(seq)

    check_blast()
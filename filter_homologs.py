import pickle
from pathlib import Path

p = Path("./cache")
raw_dict = p / "homologs_dict_raw.pkl"
filtered_dict = p / "homologs_dict_filtered.pkl"



def filter_homologs():
        with open(raw_dict, mode="rb") as f:

            homologs = pickle.load(f)
                               
                # Remove homologs with a coverage under 90%
            homologs = [i for i in homologs if i['coverage'] > 95]

                # Create homolog subsets based on percent identities
            homologs_90 = [i for i in homologs if i['identity'] < 95 and i['identity'] > 90]    # ~ 15 mutations
            print("95 group size: "+str(len(homologs_90)))

            homologs_80 = [i for i in homologs if i['identity'] < 90 and i['identity'] > 80]    # ~ 30 mutations
            print("85 group size: "+str(len(homologs_80)))

            homologs_60 = [i for i in homologs if i['identity'] < 80 and i['identity'] > 60]    # ~ 50 mutations
            print("70 group size: "+str(len(homologs_60)))

            homolog_groups = [
                    {"pident":"90", "data": homologs_90},\
                    {"pident":"80", "data": homologs_80},\
                    {"pident":"60", "data": homologs_60}
             ]

            with open(filtered_dict, mode="wb") as f:
                pickle.dump(homolog_groups, f)
                print("cached filtered homolog groups")


if __name__ == "__main__":

    filter_homologs()
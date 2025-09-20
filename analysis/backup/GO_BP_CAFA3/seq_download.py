import requests

# download Training fasta

ids = [line.strip() for line in open("Training_PosEntries_Euka_BP.tab")]

outfile = "Training_PosEntries_Euka_BP.tab.fasta"

with open(outfile, "w") as f:
    for uid in ids:
        print(f"Downloading: {uid}")
        url = f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"
        resp = requests.get(url)
        if resp.status_code == 200:
            f.write(resp.text + "\n")
        else:
            print(f"Warning: {uid} not found")

# download Bench fasta
ids = [line.strip() for line in open("Bench_PosEntries_Euka_BP.tab")]

outfile = "Bench_PosEntries_Euka_BP.tab.fasta"

with open(outfile, "w") as f:
    for uid in ids:
        print(f"Downloading: {uid}")
        url = f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"
        resp = requests.get(url)
        if resp.status_code == 200:
            f.write(resp.text + "\n")
        else:
            print(f"Warning: {uid} not found")

# python seq_download.py >Bench_download.log 
## 8% seqs were not downloaded
# python seq_download.py >Training_download.log
## 0% seqs were not downloaded

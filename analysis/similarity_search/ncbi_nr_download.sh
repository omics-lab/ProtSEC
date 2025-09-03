## How to download fasta: https://ncbiinsights.ncbi.nlm.nih.gov/2024/01/25/blast-fasta-unavailable-on-ftp/

# install ncbi-blast
sudo apt install ncbi-blast+
conda install -c bioconda blast

# download temp.sh from https://ftp.ncbi.nlm.nih.gov/blast/db/nr-prot-metadata.json
cd /mnt/c/nr_db
less temp.sh | grep ftp | tr -d '"' | tr -d ',' | awk '{print "wget -c", $1}' |  awk -F'\t' '{print}' >dnld.gz.sh
less temp.sh | grep ftp | tr -d '"' | tr -d ',' | awk '{print "wget -c", $1".md5"}' |  awk -F'\t' '{print}' >dnld.gz.md5.sh
wget -c ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr-prot-metadata.json 

# extract db files
for f in nr.*.tar.gz; do echo $f; tar -xvzf $f; done
for f in env*.tar.gz; do echo $f; tar -xvzf $f; done
rm *tar.gz 

# integrity check
blastdbcmd -db nr -info

# Download alternative blast db
update_blastdb.pl --decompress nr
update_blastdb.pl --decompress taxdb

# extract bacteria from blastdb by taxids
blastdbcmd -db nr -dbtype prot -taxids 2 -target_only -out nr.bacteria.fsa 

# extract archea from blastdb by taxids
blastdbcmd -db nr -dbtype prot -taxids 2157 -target_only -out nr.archaea.fsa 

# extract fungi from blastdb by taxids
blastdbcmd -db nr -dbtype prot -taxids 4751 -target_only -out nr.fungi.fsa 

# extract viruses from blastdb by taxids
blastdbcmd -db nr -dbtype prot -taxids 10239 -target_only -out nr.viruses.fsa 


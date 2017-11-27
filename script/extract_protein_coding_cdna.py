import sys, gzip, re
from Bio import SeqIO

in_file = gzip.open(sys.argv[1], 'rb')
out_file = gzip.open(sys.argv[2], 'wb')

records = []

for rec in SeqIO.parse(in_file, 'fasta'):
	sid = rec.id
	rec.id = sid.split('.')[0]
	if re.search('transcript_biotype:protein_coding', rec.description):
		records.append(rec)

SeqIO.write(records, out_file, 'fasta')
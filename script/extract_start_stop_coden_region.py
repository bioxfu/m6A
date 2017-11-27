import sys, gzip, re
from Bio import SeqIO

fa_file = gzip.open(sys.argv[1], 'rb')
cdna = SeqIO.to_dict(SeqIO.parse(fa_file, 'fasta'))

out_file = open(sys.argv[3], 'w')

with open(sys.argv[2]) as file:
	for line in file:
		lst = line.strip().split('\t')
		gene_id = lst[0]
		trans_id = lst[1]
		lst[4] = lst[4].replace('NA,', '')
		lst[4] = lst[4].replace(',NA', '')
		lst[5] = lst[5].replace('NA,', '')
		lst[5] = lst[5].replace(',NA', '')
		cds_start = [int(i) for i in lst[4].replace('NA,', '').split(',')][0] - 1
		cds_end = [int(i) for i in lst[5].replace('NA,', '').split(',')][-1]
		seq = cdna[trans_id].seq
		start_codon = seq[cds_start:(cds_start+3)]
		stop_codon = seq[(cds_end-3):cds_end]
		if start_codon == 'ATG' and (stop_codon == 'TAA' or stop_codon == 'TAG' or stop_codon == 'TGA'):
			seq_ext = 'N'*400 + seq + 'N'*400
			cds_start += 400
			cds_end += 400
			start_region = seq_ext[cds_start - 399 : cds_start + 400]
			stop_region = seq_ext[cds_end - 400 : cds_end + 399]
			out_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (gene_id, trans_id, start_codon, stop_codon, start_region, stop_region))

out_file.close()
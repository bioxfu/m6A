import sys, re

#motif = 'GAC[AT]'
#motif_len = 4
motif = sys.argv[1]
motif_len = int(sys.argv[2])

start_out_file = open(sys.argv[4], 'w')
stop_out_file = open(sys.argv[5], 'w')

with open(sys.argv[3]) as file:
	for line in file:
		lst = line.strip().split('\t')
		gene_id, trans_id, start, stop = lst[0], lst[1], lst[4], lst[5]
		start_hit = [gene_id, trans_id]
		stop_hit = [gene_id, trans_id]
		for i in range(len(start)-motif_len):
			win = start[i:(i+motif_len)]
			if re.search(motif, win) is not None:
				#print('start:%s\t%s\t%s' % (motif, win, i))
				start_hit.append('1')
			else:
				start_hit.append('0')

		for i in range(len(stop)-motif_len):
			win = stop[i:(i+motif_len)]
			if re.search(motif, win) is not None:
				#print('stop:%s\t%s\t%s' % (motif, win, i))
				stop_hit.append('1')
			else:
				stop_hit.append('0')

		start_out_file.write('%s\n' % '\t'.join(start_hit))
		stop_out_file.write('%s\n' % '\t'.join(stop_hit))


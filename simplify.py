import sys

inFile = open(sys.argv[1],'r')
inFile.next()

def getAln(stream):
	block = []
	for line in stream:
		if line.strip() == '':
			break
		block.append(line.strip().split())
	return block

def seq_metric(trp, ref):
	matches = 0
	mismatches = 0
	ref_deletion = 0
	ref_insertion = 0
	ns = 0
	for i in range(len(trp)):
		a = trp[i]
		b = ref[i]

		if a == b:
			matches += 1
		else:
			if b == '-':
				ref_deletion += 1
			elif a == '-':
				ref_insertion += 1
			elif b == 'n':
				ns += 1
			else:
				mismatches += 1
	return (matches, mismatches, ref_deletion, ref_insertion, ns)


while 1:
	aln = getAln(inFile)
	if len(aln) == 0:
		break

	ref_t, ref_name, ref_start, alnLength, ref_strand, ref_size, ref_seq = aln[1]
	trp_t, trp_name, trp_start, alnLength, trp_strand, trp_size, trp_seq = aln[2]

	ref_coord = ref_start + ',' + str(int(ref_start) + int(alnLength)) + ',' + ref_strand
	trp_coord = trp_start + ',' + str(int(trp_start) + int(alnLength)) + ',' + trp_strand

	m, mm, r_del, r_ins, ns = seq_metric(trp_seq, ref_seq)

	indel = r_del + r_ins
	alnLength_noIndel = int(alnLength) - r_del - r_ins

	print '\t'.join(map(str,[ref_name, ref_size, ref_coord, trp_name, trp_size, trp_coord, alnLength, ref_seq, trp_seq, m, mm, r_del, r_ins, ns]))


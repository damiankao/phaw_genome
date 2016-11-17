import sys
from collections import defaultdict

inFile = open(sys.argv[1],'r')

def processStack(header, stack):
	if len(stack) > 1:
		refStack = stack[1:]
		refStack = zip(*refStack)

		typecodes = []

		cov = 0

		tranSeq = stack[0]
		for c, col in enumerate(refStack):
			aligned = [b for b in col if b != ' ']
			if len(aligned) < 5:
				if len(aligned) > 0:
					cov += 1
				
					het = '0'
					vari = '1'
					indel = '0'

					tranBase = tranSeq[c]
					primary = aligned[0]

					if len(aligned) > 1:
						secondary = aligned[1]
					else:
						secondary = primary

					if primary != secondary:
						het = '1'

					if tranBase == primary or tranBase == secondary:
						vari = '0'

					if primary == '-' or secondary == '-':
						indel = '1'

					tc = het + vari + indel
					if tc != '000':
						typecodes.append(str(c) + ":" + tc)

		tidLength = len(stack[0])

		print header.strip().split()[0][1:] + '\t' + str(tidLength) + '\t' + str(cov) + '\t' + ','.join(typecodes)

header = inFile.next().strip()
stack = []
for line in inFile:
	if line[0] == '>':
		processStack(header,stack)

		header = line.strip()
		stack = []
	else:
		stack.append([x for x in line[:-1]])

processStack(header,stack)


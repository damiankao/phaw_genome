import sys
from Bio import SeqIO
from collections import defaultdict

faFile = open(sys.argv[1],'r')
inFile = open(sys.argv[2],'r')

seqs = {}
for record in SeqIO.parse(faFile,'fasta'):
	seqs[record.id] = str(record.seq)

def revComp(seq):
	alph = {'a':'t','t':'a','g':'c','c':'g','n':'n','-':'-'}
	return ''.join([alph[c] for c in seq][::-1])

def processBlock(block):
	blockRef = defaultdict(list)

	for data in block:
		blockRef[data[0]].append(data)

	blockRef = blockRef.values()
	blockRef = [x for x in blockRef if sum([int(y[9]) for y in x]) >= 30]
	blockRef.sort(key = lambda x : sum([int(y[9]) for y in x]), reverse = True)

	trp = block[0][3]
	trp_fullSeq = seqs[block[0][3]].lower()

	stack = [[x] for x in trp_fullSeq]

	stackRows = []
	stackRefs = defaultdict(lambda : {'trp_ins':0,'match':0,'mismatch':0,'trp_del':0,'ns':0})
	for refAln in blockRef:
		for col in stack:
			col.append(' ')

		totalTrpIns = 0
		for aln in refAln:
			ref, refSize, refCoord, trp, trpSize, trpCoord, aln, refSeq, trpSeq, m, mm, rDel, rIns, ns = aln
			refStart, refEnd, refStrand = refCoord.split(',')
			refStart = int(refStart)
			refEnd = int(refEnd)

			trpStart, trpEnd, trpStrand = trpCoord.split(',')
			trpStart = int(trpStart)
			trpEnd = int(trpEnd)

			if trpStrand == '-':
				trpStart = int(trpSize) - trpEnd
				trpEnd = int(trpSize) - trpStart
				trpSeq = revComp(trpSeq)
				refSeq = revComp(refSeq)

			insCount = 0
			for i, refC in enumerate(refSeq):
				trpC = trpSeq[i]
				if trpC != '-':
					stack[trpStart + i - insCount][-1] = refC
				else:
					insCount += 1
			totalTrpIns += insCount

		stackRows.append(refAln[0][0])
		stackRefs[refAln[0][0]]['trp_ins'] = totalTrpIns

	for col in stack:
		trpCol = col[0]
		for i, refName in enumerate(stackRows):
			refCol = col[i + 1]
			if refCol != ' ':
				if refCol == trpCol:
					stackRefs[refName]['match'] += 1
				else:
					if refCol == 'n':
						stackRefs[refName]['ns'] += 1
					elif refCol == '-':
						stackRefs[refName]['trp_del'] += 1
					else:
						stackRefs[refName]['mismatch'] += 1

	rowMeta = []
	for refName in stackRows:
		refMeta = stackRefs[refName]
		rowMeta.append(refName + ':' + str(refMeta['match']) + ':' + str(refMeta['mismatch']) + ':' + str(refMeta['trp_del']) + ':' + str(refMeta['trp_ins']) + ':' + str(refMeta['ns']))

	print '>' + trp + '\t' + ','.join(rowMeta)
	print '\n'.join(''.join(x) for x in zip(*stack))

	return stack

currentTrp = ''
block = []
for line in inFile:
	meta = line.strip().split()

	if currentTrp != meta[3]:
		if len(block) != 0:
			processBlock(block)

		block = []
		currentTrp = meta[3]
		block.append(meta)
	else:
		block.append(meta)

processBlock(block)


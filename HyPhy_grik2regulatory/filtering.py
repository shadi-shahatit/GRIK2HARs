import re  # for Regular expressions
import sys
from Bio import AlignIO

# infile = sys.argv[1]

infile = 'all.prunned.list'

with open(infile) as f:
    reflist = f.read().splitlines()

missing = r'[(+*)]'
ambiguous = r'[N]'

badsimiosNs = []
badsimiosAsk = []
goodsimios = []

for i in reflist:
    with open(i, 'r') as myfile:
        mydata = myfile.read()
    fasta = AlignIO.read(i, "fasta")
    
    if fasta.get_alignment_length() > 350:
        maxNs = len(re.findall(ambiguous, mydata, re.I))
        maxAsk = len(re.findall(missing, mydata, re.I))
        
        if maxNs > 2:
            print(')-B')
            badsimiosNs.append(i)
        elif maxAsk > 1:
            print(')-;')
            badsimiosAsk.append(i)
        else:
            print('lol')
            goodsimios.append(i)

# Write good alignments
with open('goodalignments.txt', 'w') as thefile3:
    for item in goodsimios:
        thefile3.write(f"{item}\n")

# Write ambiguous sequences
with open('ambiguous.txt', 'w') as thefile1:
    for item in badsimiosNs:
        thefile1.write(f"{item}\n")

# Write asterisks sequences
with open('asterisks.txt', 'w') as thefile2:
    for item in badsimiosAsk:
        thefile2.write(f"{item}\n")


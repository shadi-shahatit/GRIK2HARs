import re  # for Regular expressions
import sys

# No need to import 'sets', using the built-in 'set' in Python 3

infile = 'all.list'

# Open the file containing the list of fasta files
with open(infile) as f:
    reflist = f.read().splitlines()

for i in reflist:
    with open(i, 'r') as myfile:
        data = myfile.read().splitlines()

    # Create output file for pruned sequences
    with open('%s.prunned' % i, 'w') as output:
        results = {}
        seq = ""
        compName = ""

        # Process the data
        for line in data:
            if line.startswith(">"):
                if compName:
                    results[compName] = seq
                compName = line[1:]
                seq = ""
            else:
                seq += line

        results[compName] = seq

        # Initialize variables
        nsize = len(results[list(results.keys())[0]])  # First sequence length
        ns = len(results.keys())
        indels = set()  # Replacing 'sets.Set()' with built-in 'set'
        prunned = {}

        # Identify positions with indels (non-ACGT characters)
        for spec in results.keys():
            for i in range(nsize):
                base = results[spec][i]
                if base.upper() not in ["A", "C", "G", "T"]:
                    indels.add(i)

        # Prune the sequences based on indel positions
        nsize2 = nsize - len(indels)
        for spec in results.keys():
            seq = ""
            for i in range(nsize):
                if i not in indels:
                    seq += results[spec][i]
            prunned[spec] = seq

        # Write pruned sequences to the output file
        for f, spec in enumerate(prunned.keys()):
            if f > 0:
                output.write("\n")
            output.write(f">{spec}\n")
            output.write(prunned[spec])

print("Pruning complete.")

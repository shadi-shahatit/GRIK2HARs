import random
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

keys = "queries.list"
values = "neutral.list"

# Read query and neutral sequence lists
with open(keys) as f:
    keylist = f.read().splitlines()

with open(values) as f:
    valuelist = f.read().splitlines()

# Create dictionary with queries and random neutral sequences
dictionary = {}
for key in keylist:
    random.shuffle(valuelist)
    dictionary[key] = valuelist[:10]  # Get the first 10 neutral sequences

# Write the dictionary to a file
with open("global.dict", "w") as fielddict_file:
    for key, value in dictionary.items():
        fielddict_file.write(f"{key}: {value}\n")

# Process the alignments
reference = []
for i, j in dictionary.items():  # Use items() for Python 3
    n = 0
    # Initialize the combined alignment with empty sequences for each species
    combined_seq = MultipleSeqAlignment([
        SeqRecord(Seq('', generic_dna), id="hg19"),
        SeqRecord(Seq('', generic_dna), id="panTro4"),
        SeqRecord(Seq('', generic_dna), id="gorGor3"),
        SeqRecord(Seq('', generic_dna), id="rheMac3"),
        SeqRecord(Seq('', generic_dna), id="ponAbe2")
    ])
    
    combined_seq.sort()
    
    for ref in j:
        n += 1
        seq_records = AlignIO.read(ref, 'fasta')
        seq_records.description = ""
        seq_records.sort()
        combined_seq += seq_records
        combined_seq.description = ""
    
    # Write the combined alignment to a file
    with open(f'{i}.ref', 'w') as write_file:
        AlignIO.write(combined_seq, write_file, 'fasta')
    
    # Update the reference list
    with open('reference.list', 'a') as referencelist:
        referencelist.write(f'{i}\t{n}\n')

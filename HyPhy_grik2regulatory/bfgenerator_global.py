import sys
import random

# Get the query file from command-line arguments
query = sys.argv[1]

# Read the query list from the file
with open(query, 'r') as f:
    querylist = f.read().splitlines()

# Define models and branches
models = ['null', 'alt']
branches = ['hg38', 'panTro5']

# Generate batch files for each combination of model, branch, and query
for model in models:
    for branch in branches:
        for i in querylist:
            # Generate a random integer
            k = random.randint(1, 1000)
            
            # Define the filename for the batch file
            filename = f'{i}.{branch}.{model}.bf'
            
            # Write to the batch file
            with open(filename, 'w') as f:
                f.write(f'random_seed={k};\n')
                f.write(f'quer_seq_file="query/{i}.fa.prunned";\n')
                f.write(f'ref_seq_file="ref/{i}.ref.fa";\n')
                f.write('fit_repl_count=20;\n')
                f.write('tree="((((hg38,panTro5),gorGor5),ponAbe2),rheMac8)";\n')
                f.write(f'fgrnd_branch_name="{branch}";\n')
                f.write(f'res_file="res/{i}.{branch}.{model}.res";\n')
                f.write(f'#include "{model}4-fgrnd_spec.bf";\n')


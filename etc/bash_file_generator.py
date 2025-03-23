maf_files = [
    'maf_LR537121.csv', 'maf_LR537122.csv', 'maf_LR537123.csv', 'maf_LR537124.csv',
    'maf_LR537125.csv', 'maf_LR537126.csv', 'maf_LR537127.csv', 'maf_LR537128.csv',
    'maf_LR537129.csv', 'maf_LR537130.csv', 'maf_LR537131.csv', 'maf_LR537132.csv',
    'maf_LR537133.csv', 'maf_LR537134.csv', 'maf_LR537135.csv', 'maf_LR537136.csv',
    'maf_LR537137.csv', 'maf_LR537138.csv', 'maf_LR537139.csv', 'maf_LR537140.csv',
    'maf_LR537141.csv', 'maf_LR537142.csv', 'maf_LR537143.csv', 'maf_LR537144.csv'
]

template = """#!/bin/bash
#SBATCH --time=12:00:00             # Run time (days-hh:mm:ss) - (max 7days)
#SBATCH --job-name=fet_{job_name}
#SBATCH --cpus-per-task=12 
#SBATCH --mail-user=amitsb@bio.auth.gr
#SBATCH --mail-type=BEGIN,END,FAIL

perl fisher-test.pl --input ./{input_file} --output ./sync_{output_file} --min-count 1 --min-coverage 1 --max-coverage 10000 --window-size=1 --step-size=1
"""

fst_template = """#!/bin/bash
#SBATCH --time=08:00:00             # Run time (days-hh:mm:ss) - (max 7days)
#SBATCH --job-name=fst_{job_name}
#SBATCH --cpus-per-task=12 
#SBATCH --mail-user=amitsb@bio.auth.gr
#SBATCH --mail-type=BEGIN,END,FAIL

perl fst-sliding.pl --input ./{input_file} --output ./sync_{output_file} --min-count 1 --min-coverage 1 --max-coverage 100000 --window-size=1 --step-size=1 --pool-size=50:50:50:50:50:50:50:50:14:13:50:50:50:50:50:50:50:50:50:50
"""

snp_freq_template = """#!/bin/bash
#SBATCH --time=08:00:00             # Run time (days-hh:mm:ss) - (max 7days)
#SBATCH --job-name=snp_freq_{job_name}
#SBATCH --cpus-per-task=12 
#SBATCH --mail-user=amitsb@bio.auth.gr
#SBATCH --mail-type=BEGIN,END,FAIL

perl snp-frequency-diff.pl --input ./{input_file} --output-prefix sync_{output_prefix} --min-coverage=1 --max-coverage=10000
"""


for maf_file in maf_files:
    job_name = maf_file.split('_')[1].split('.')[0]
    script_name = f"Sp_aurata_{job_name}_fet.sh"
    input_file = maf_file.replace('.csv', '_popoolation.txt')
    output_file = f"Sp_aurata_{job_name}_popoolation.fet"
    
    with open(script_name, 'wb') as script_file:
        script_file.write(template.format(job_name=job_name, input_file=input_file, output_file=output_file))

print("Bash scripts generated successfully!")

for maf_file in maf_files:
    job_name = maf_file.split('_')[1].split('.')[0]
    input_file = maf_file.replace('.csv', '_popoolation.txt')
    output_file_fst = f"Sp_aurata_{job_name}_popoolation.fst"
    output_prefix_snp = f"Sp_aurata_{job_name}_popoolation"
    
    # Generate fst-sliding script
    fst_script_name = f"Sp_aurata_{job_name}_fst-sliding.sh"
    with open(fst_script_name, 'wb') as fst_script:
        fst_script.write(fst_template.format(job_name=job_name, input_file=input_file, output_file=output_file_fst))
    
    # Generate snp-frequency-diff script
    snp_freq_script_name = f"Sp_aurata_{job_name}_snp_freq.sh"
    with open(snp_freq_script_name, 'wb') as snp_freq_script:
        snp_freq_script.write(snp_freq_template.format(job_name=job_name, input_file=input_file, output_prefix=output_prefix_snp))

print("fst-sliding and snp-frequency-diff bash scripts generated successfully!")

# Dynamic_checker
A Dynamic Programming approach for verifying results of CRISPRme

First, create & activate the conda env:

```bash
conda create -n crisp_env python=3.10 numpy biopython pandas matplotlib matplotlib-venn -c conda-forge
conda activate crisp_env
```

To Run the Checker run the following in the command line:
```bash
python checker_command_line.py \
  --fasta path/to/enriched.fa \
  --rna "CTAACAGTTGCTTTTATCACNGG" \
  --max-mismatches 4 \
  --max-dna-gaps 1 \
  --max-rna-gaps 1 \
  --chr 22 \
  --output results.csv 
```
After Running the Dynamic script, run comparison.py:
```bash
python comparison.py \
  path/to/crisprme_main.tsv \
  path/to/crisprme_alt.tsv \
  results.csv \
  output_folder
```

After running comparison.py, your output_folder will contain:

CSVs:

crisprme_only.csv

dynamic_only.csv

both_matched.csv

combined.csv

Plots:

alignments_found.png

alignments_venn.png

unique_start_values.png

start_positions_venn.png

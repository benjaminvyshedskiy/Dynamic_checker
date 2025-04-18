import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import os


def main():
    parser = argparse.ArgumentParser(
        description="Compare CRISPRme and dynamic programming approaches."
    )
    parser.add_argument(
        "Allignments_Main",
        help="Path to the CRISPRme Main Allignments"
    )
    parser.add_argument(
        "Allignments_Alt",
        help="Path to the CRISPRme Alternate Allignments"
    )
    parser.add_argument(
        "dynamic_programing_file",
        help="Path to dynamic approach csv"
    )
    parser.add_argument(
        "output_folder",
        help="path to output folder"
    )
    parser.add_argument(
        "Filter_Pam",
        default=True,
        help="Remove Pams with Gaps and missmatches"
    )

    args = parser.parse_args()

    Crispr_Main_path = args.Allignments_Main
    Crispr_Alt_path = args.Allignments_Alt
    Dynamic_path = args.dynamic_programing_file
    Output_folder = args.output_folder
    Filter_pam = args.Filter_Pam

    table = pd.read_csv(Crispr_Alt_path, sep='\t')
    table2 = pd.read_csv(Crispr_Main_path, sep='\t')
    table = pd.concat([table, table2], axis=0, ignore_index=True)
    cleaned_crisprme_allignments=pd.DataFrame()
    cleaned_crisprme_allignments["CHR"] = table["Chromosome"]
    cleaned_crisprme_allignments["Start"] = table["Start_coordinate_(highest_CFD)"]
    cleaned_crisprme_allignments["RNA"] = table["Aligned_spacer+PAM_(highest_CFD)"]
    cleaned_crisprme_allignments["DNA"] = table["Aligned_protospacer+PAM_REF_(highest_CFD)"]
    cleaned_crisprme_allignments["Strand"] = table["Strand_(highest_CFD)"]
    cleaned_crisprme_allignments["mismatches"] = table["Mismatches_(highest_CFD)"]
    cleaned_crisprme_allignments["gaps"] = table["Bulges_(highest_CFD)"]
    cleaned_crisprme_allignments["gap_location"] = table["Bulge_type_(highest_CFD)"]
    cleaned_crisprme_allignments = cleaned_crisprme_allignments.sort_values(by="Start")

    crisprme_data = cleaned_crisprme_allignments
    dynamic_data = pd.read_csv(Dynamic_path)


    iupac_dict = {
        "A": {"A"},
        "C": {"C"},
        "G": {"G"},
        "T": {"T", "U"},
        "U": {"T", "U"},
        "R": {"A", "G"},
        "Y": {"C", "T"},
        "S": {"G", "C"},
        "W": {"A", "T"},
        "K": {"G", "T"},
        "M": {"A", "C"},
        "B": {"C", "G", "T"},
        "D": {"A", "G", "T"},
        "H": {"A", "C", "T"},
        "V": {"A", "C", "G"},
        "N": {"A", "C", "G", "T"},
        "-": {"-"}
    }

    def iupac_match(a, b):
        """Check if two nucleotides (possibly ambiguous per IUPAC) overlap."""
        return bool(iupac_dict.get(a.upper(), set()) & iupac_dict.get(b.upper(), set()))

    def exact_iupac_match(seq1, seq2):
        if len(seq1) != len(seq2):
            return False
        output = True
        for i in range(len(seq1)):
            output = output and iupac_match(seq1[i], seq2[i])
        return output

    # --- Existing function to check if a crispr row is in dynamic data ---
    def is_sequence_in_reference(chr_val, strand, start, RNA, DNA):
        # Filter dynamic_data for candidate rows.
        DF_potential_search = dynamic_data[(dynamic_data["CHR"] == chr_val) &
                                            (dynamic_data["Strand"] == strand) &
                                            (dynamic_data["Start"] == start)]
        RNA_MATCH = any(exact_iupac_match(rna_seq, RNA) for rna_seq in DF_potential_search["RNA"].dropna())
        DNA_MATCH = any(exact_iupac_match(dna_seq, DNA) for dna_seq in DF_potential_search["DNA"].dropna())
        return RNA_MATCH and DNA_MATCH

    # --- New helper: check if a dynamic row is found in the crisprme data ---
    def is_sequence_in_crisprme_data(chr_val, strand, start, RNA, DNA):
        DF_potential_search = filtered_crisprme_data[(filtered_crisprme_data["CHR"] == chr_val) &
                                                    (filtered_crisprme_data["Strand"] == strand) &
                                                    (filtered_crisprme_data["Start"] == start)]
        RNA_MATCH = any(exact_iupac_match(rna_seq, RNA) for rna_seq in DF_potential_search["RNA"].dropna())
        DNA_MATCH = any(exact_iupac_match(dna_seq, DNA) for dna_seq in DF_potential_search["DNA"].dropna())
        return RNA_MATCH and DNA_MATCH

    # --- Marking matches in filtered crisprme data ---
    # Here, note that in your working code you set chr_val manually; adjust if needed.
    def mark_in_dynamic(row):
        return is_sequence_in_reference(row['CHR'], row['Strand'], row['Start'], row['RNA'], row['DNA'])
    
    def valid_pam(sequence):
        last_chr = sequence[-1]
        second_last = sequence[-2]

        return(iupac_match(last_chr,"G")&iupac_match(second_last,"G"))


    crisprme_data['total_mistakes'] = crisprme_data['mismatches'] + crisprme_data['gaps']

    grouped_best = crisprme_data.groupby(['Start', 'RNA'])['total_mistakes'].transform('min')
    filtered_crisprme_data = crisprme_data[crisprme_data['total_mistakes'] == grouped_best]

    filtered_crisprme_data = filtered_crisprme_data.reset_index(drop=True)
    if Filter_pam:dynamic_data = dynamic_data[dynamic_data["DNA"].apply(valid_pam)]


    # Apply function to create a new boolean column "in_dynamic"
    filtered_crisprme_data['in_dynamic'] = filtered_crisprme_data.apply(mark_in_dynamic, axis=1)
    # --- Marking matches in dynamic data ---
    def mark_in_crispr(row):
        return is_sequence_in_crisprme_data(row['CHR'], row['Strand'], row['Start'], row['RNA'], row['DNA'])

    # Apply function to dynamic_data, adding column "in_crisprme"
    dynamic_data['in_crisprme'] = dynamic_data.apply(mark_in_crispr, axis=1)

    # --- Create the three DataFrames ---
    # 1. Rows found only in crisprme_filtered (no match in dynamic data)
    df_crisprme_only = filtered_crisprme_data[~filtered_crisprme_data['in_dynamic']].copy()

    # 2. Rows found only in dynamic filtered (no match in crisprme data)
    df_dynamic_only = dynamic_data[~dynamic_data['in_crisprme']].copy()

    # 3. Rows found in both.
    # To capture "both", we merge the matched rows on the common keys.
    df_both_crispr = filtered_crisprme_data[filtered_crisprme_data['in_dynamic']].copy()
    df_both_dynamic = dynamic_data[dynamic_data['in_crisprme']].copy()
    df_both = pd.merge(df_both_crispr, df_both_dynamic, 
                    left_on=['CHR', 'Start', 'Strand'], 
                    right_on=['CHR', 'Start', 'Strand'],
                    suffixes=('_crispr', '_dyn'))

    # --- Save each DataFrame as a CSV file ---
    df_crisprme_only.to_csv(f'{Output_folder}/crisprme_only.csv', index=False)
    df_dynamic_only.to_csv(f'{Output_folder}/dynamic_only.csv', index=False)
    df_both = df_both[df_both.apply(lambda row: exact_iupac_match(row['DNA_crispr'], row['DNA_dyn']), axis=1)]
    df_both = df_both[df_both.apply(lambda row: exact_iupac_match(row['RNA_crispr'], row['RNA_dyn']), axis=1)]
    df_both.to_csv(f'{Output_folder}/both_matched.csv', index=False)

    #Create and save graphs 
    categories = ['Dynamic only', 'Crisprme only', 'Both']
    counts = [len(df_dynamic_only), len(df_crisprme_only), len(df_both)]

    plt.figure(figsize=(8, 6))
    plt.bar(categories, counts)
    plt.title("Amount of possible alignments found")
    plt.ylabel("Number of Alignments")
    plt.savefig(f"{Output_folder}/alignments_found.png")
    plt.close()

    plt.figure(figsize=(6, 6))
    venn2(
        subsets=(len(df_dynamic_only), len(df_crisprme_only), len(df_both)),
        set_labels=('Dynamic', 'CRISPRme')
    )
    plt.title("Overlap of Alignments")
    plt.savefig(f"{Output_folder}/alignments_venn.png")
    plt.close()

    


    # your existing unique‑counts bar‑chart …
    unique_counts = [
        df_dynamic_only['Start'].nunique(),
        df_crisprme_only['Start'].nunique(),
        df_both['Start'].nunique()
    ]

    plt.figure(figsize=(8, 6))
    plt.bar(categories, unique_counts, color="skyblue")
    plt.title("Unique Start Values")
    plt.ylabel("Unique Value Count")
    plt.savefig(f"{Output_folder}/unique_start_values.png")
    plt.close()


    starts_dyn    = set(dynamic_data['Start'])
    starts_crispr = set(filtered_crisprme_data['Start'])
    only_dyn      = starts_dyn - starts_crispr
    only_cris     = starts_crispr - starts_dyn
    intersect     = starts_dyn & starts_crispr

    # 2. Bar chart of true exclusives + intersection
    plt.figure(figsize=(8, 6))
    plt.bar(
        ['Dynamic only', 'Crisprme only', 'Both'],
        [len(only_dyn), len(only_cris), len(intersect)],
        color="skyblue"
    )
    plt.title("Unique Start Values")
    plt.ylabel("Unique Value Count")
    plt.savefig(f"{Output_folder}/unique_start_values.png")
    plt.close()

    # 3. 2‑set Venn diagram
    plt.figure(figsize=(6, 6))
    venn2(
        subsets=(len(only_dyn), len(only_cris), len(intersect)),
        set_labels=('Dynamic', 'CRISPRme')
    )
    plt.title("Overlap of Start Positions")
    plt.savefig(f"{Output_folder}/start_positions_venn.png")
    plt.close()

    files_labels  = {
    'both_matched.csv':   'BOTH',
    'crisprme_only.csv':  'CRISPRME',
    'dynamic_only.csv':   'DYNAMIC'
    }
    combined_fname = 'combined.csv'

    # --- LOAD, TAG, COLLECT ---
    dfs = []
    for fname, label in files_labels.items():
        path = os.path.join(Output_folder, fname)
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing expected file: {path}")
        df = pd.read_csv(path, dtype=str)          # read everything as string to avoid dtype mismatches
        df['Source'] = label
        dfs.append(df)

    # --- CONCAT & SAVE ---
    combined = pd.concat(dfs, ignore_index=True, sort=False)
    out_path = os.path.join(Output_folder, combined_fname)
    combined.to_csv(out_path, index=False)
    
    print(f"Saved combined file to {out_path}")
    files_labels  = {
        'both_matched.csv':   'BOTH',
        'crisprme_only.csv':  'CRISPRME',
        'dynamic_only.csv':   'DYNAMIC'
    }
    combined_fname = 'combined.csv'

    # --- LOAD, TAG, COLLECT ---
    dfs = []
    for fname, label in files_labels.items():
        path = os.path.join(Output_folder, fname)
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing expected file: {path}")
        df = pd.read_csv(path, dtype=str)          # read everything as string to avoid dtype mismatches
        df['Source'] = label
        dfs.append(df)

    # --- CONCAT & SAVE ---
    combined = pd.concat(dfs, ignore_index=True, sort=False)
    out_path = os.path.join(Output_folder, combined_fname)
    combined.to_csv(out_path, index=False)
    print(f"Saved combined file to {out_path}")


if __name__ == "__main__":
    main()

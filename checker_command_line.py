#!/usr/bin/env python3
import numpy as np
import csv
import argparse
from Bio import SeqIO

# ----------------------------
# (1) IUPAC dict and matcher — unchanged
iupac_dict = {
    "A": {"A"}, "C": {"C"}, "G": {"G"},
    "T": {"T", "U"}, "U": {"T", "U"},
    "R": {"A", "G"}, "Y": {"C", "T"},
    "S": {"G", "C"}, "W": {"A", "T"},
    "K": {"G", "T"}, "M": {"A", "C"},
    "B": {"C", "G", "T"}, "D": {"A", "G", "T"},
    "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"}
}

def iupac_match(a, b):
    return bool(iupac_dict.get(a.upper(), set()) & iupac_dict.get(b.upper(), set()))

# ----------------------------
# (2) Needleman–Wunsch for all alignments — unchanged
def needleman_wunsch_all_alignments(
    RNA, DNA,
    match=0, mismatch=-1, gap=-1,
    max_DNA_gaps=1, max_RNA_gaps=1, max_mismatches=4
):
    m, n = len(RNA), len(DNA)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    traceback = [[[] for _ in range(n + 1)] for _ in range(m + 1)]
    for i in range(1, m + 1):
        dp[i][0] = dp[i - 1][0] + gap
        traceback[i][0].append(('up', i - 1, 0))
    for j in range(1, n + 1):
        dp[0][j] = dp[0][j - 1] + gap
        traceback[0][j].append(('left', 0, j - 1))

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if iupac_match(RNA[i - 1], DNA[j - 1]):
                diag_score = dp[i - 1][j - 1] + match
            else:
                diag_score = dp[i - 1][j - 1] + mismatch
            up_score = dp[i - 1][j] + gap
            # no gap‐penalty on last row for left moves:
            left_score = dp[i][j - 1] if i == m else dp[i][j - 1] + gap

            max_score = max(diag_score, up_score, left_score)
            dp[i][j] = max_score
            if diag_score == max_score:
                traceback[i][j].append(('diag', i - 1, j - 1))
            if up_score == max_score:
                traceback[i][j].append(('up', i - 1, j))
            if left_score == max_score:
                traceback[i][j].append(('left', i, j - 1))

    alignments = []
    def backtrack(i, j, aligned1, aligned2):
        if i == 0 and j == 0:
            RNA_GAPS = aligned1.strip("-").count("-")
            DNA_GAPS = aligned2.count("-")
            Missmatches = 0 - dp[m][n] - RNA_GAPS - DNA_GAPS
            if (aligned1[-1] != "-" and
                RNA_GAPS <= max_RNA_gaps and
                DNA_GAPS <= max_DNA_gaps and
                Missmatches <= max_mismatches):
                alignments.append((
                    aligned1[::-1],
                    aligned2[::-1],
                    Missmatches,
                    RNA_GAPS,
                    DNA_GAPS
                ))
            return
        for direction, prev_i, prev_j in traceback[i][j]:
            if direction == 'diag':
                backtrack(prev_i, prev_j,
                          aligned1 + RNA[i - 1],
                          aligned2 + DNA[j - 1])
            elif direction == 'up':
                backtrack(prev_i, prev_j,
                          aligned1 + RNA[i - 1],
                          aligned2 + '-')
            else:  # left
                backtrack(prev_i, prev_j,
                          aligned1 + '-',
                          aligned2 + DNA[j - 1])

    minscore = -max_mismatches - max_DNA_gaps - max_RNA_gaps
    if dp[m][n] >= minscore:
        backtrack(m, n, '', '')
        if alignments:
            return alignments

# ----------------------------
# (3) clear_csv — unchanged
def clear_csv(file_path):
    with open(file_path, 'w'):
        pass

# ----------------------------
# (4) reverse_complement — unchanged
def reverse_complement(seq):
    complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'R': 'Y', 'Y': 'R',
        'S': 'S', 'W': 'W',
        'K': 'M', 'M': 'K',
        'B': 'V', 'V': 'B',
        'D': 'H', 'H': 'D',
        'N': 'N'
    }
    return "".join(complement.get(b, b) for b in reversed(seq.upper()))

# ----------------------------
# (5) scan — only changed to:
#    • accept max_mismatches and output_file
#    • forward those into the alignment call
def scan(
    RNA, DNA, Strand,
    max_dna_gaps=1, max_rna_gaps=1, max_mismatches=4,
    chr_num=22, output_file='results.csv'
):
    with open(output_file, 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for starting_pos in range(0, len(DNA) - len(RNA) + max_dna_gaps + 1):
            if starting_pos % 10000 == 0:
                print(f"Position {starting_pos}")
            slice_end = min(starting_pos + len(RNA) + max_rna_gaps, len(DNA))
            DNA_Slice = DNA[starting_pos:slice_end]
            if "N" in DNA_Slice:
                continue
            aligns = needleman_wunsch_all_alignments(
                RNA, DNA_Slice,
                max_DNA_gaps=max_dna_gaps,
                max_RNA_gaps=max_rna_gaps,
                max_mismatches=max_mismatches
            )
            if not aligns:
                continue
            for aln in aligns:
                r_aln, d_aln, mm, rg, dg = aln
                actual_len = len(r_aln.rstrip('-'))
                if Strand == "+":
                    start = starting_pos
                    end   = starting_pos + actual_len
                else:
                    start = len(DNA) - starting_pos - actual_len + dg
                    end   = len(DNA) - starting_pos + dg
                writer.writerow([
                    f"chr{chr_num}",
                    r_aln[:actual_len],
                    d_aln[:actual_len],
                    Strand,
                    start,
                    end,
                    mm, rg, dg
                ])

# ----------------------------
# (6) CLI entrypoint
if __name__ == '__main__':
    p = argparse.ArgumentParser(
        description="Scan a FASTA for matches of an RNA sequence (w/ IUPAC ambiguity)."
    )
    p.add_argument('--fasta',        required=True,
                   help="Path to enriched FASTA file")
    p.add_argument('--rna',          required=True,
                   help="RNA sequence (IUPAC allowed)")
    p.add_argument('--max-mismatches', type=int, default=4,
                   help="Maximum mismatches allowed")
    p.add_argument('--max-dna-gaps',   type=int, default=1,
                   help="Maximum DNA gaps")
    p.add_argument('--max-rna-gaps',   type=int, default=1,
                   help="Maximum RNA gaps")
    p.add_argument('--chr', dest='chr_num', type=int, required=True,
                   help="Chromosome number (for output prefix)")
    p.add_argument('--output',       required=True,
                   help="CSV output file path")
    args = p.parse_args()

    # read sequences
    record = SeqIO.read(args.fasta, "fasta")
    DNA    = str(record.seq)
    RNAseq = args.rna

    # prepare output
    clear_csv(args.output)
    with open(args.output, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            "CHR", "RNA", "DNA", "Strand",
            "Start", "END",
            "mismatches", "gaps in RNA", "gaps in DNA"
        ])

    # scan forward and reverse strands
    scan(
        RNAseq, DNA, "+",
        args.max_dna_gaps,
        args.max_rna_gaps,
        args.max_mismatches,
        args.chr_num,
        args.output
    )
    scan(
        RNAseq, reverse_complement(DNA), "-",
        args.max_dna_gaps,
        args.max_rna_gaps,
        args.max_mismatches,
        args.chr_num,
        args.output
    )

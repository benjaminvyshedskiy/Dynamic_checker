import csv
from Bio import SeqIO

# IUPAC matching rules
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
    return bool(iupac_dict[a.upper()] & iupac_dict[b.upper()])

def reverse_complement(seq):
    comp = {
        'A':'T','T':'A','C':'G','G':'C',
        'R':'Y','Y':'R','S':'S','W':'W',
        'K':'M','M':'K','B':'V','V':'B',
        'D':'H','H':'D','N':'N'
    }
    return "".join(comp.get(b, b) for b in reversed(seq.upper()))

def needleman_wunsch_all_alignments(
    RNA, DNA,
    max_mismatches=4,
    max_DNA_gaps=1,
    max_RNA_gaps=1
):
    m, n = len(RNA), len(DNA)

    # dp_score and dp_cells as before
    dp_score = [[float('inf')]*(n+1) for _ in range(m+1)]
    dp_cells = [[{} for _ in range(n+1)] for _ in range(m+1)]

    # Base case
    dp_score[0][0] = 0
    dp_cells[0][0] = {(0,0,0): []}

    # First column (DNA gaps)
    for i in range(1, m+1):
        cnt = (0, i, 0)
        if i <= max_DNA_gaps:
            dp_score[i][0] = i
            dp_cells[i][0] = {cnt: [('up', i-1, 0)]}

    # First row (RNA gaps)
    for j in range(1, n+1):
        cnt = (0, 0, j)
        if j <= max_RNA_gaps:
            dp_score[0][j] = j
            dp_cells[0][j] = {cnt: [('left', 0, j-1)]}

    # Fill DP with threshold pruning
    for i in range(1, m+1):
        for j in range(1, n+1):
            candidates = []
            # Diagonal
            for (mm, dg, rg) in dp_cells[i-1][j-1]:
                mis = 0 if iupac_match(RNA[i-1], DNA[j-1]) else 1
                new = (mm+mis, dg, rg)
                if new[0] <= max_mismatches:
                    candidates.append((dp_score[i-1][j-1] + mis,
                                       new, 'diag', i-1, j-1))
            # Up
            for (mm, dg, rg) in dp_cells[i-1][j]:
                new = (mm, dg+1, rg)
                if new[1] <= max_DNA_gaps:
                    candidates.append((dp_score[i-1][j] + 1,
                                       new, 'up', i-1, j))
            # Left
            for (mm, dg, rg) in dp_cells[i][j-1]:
                if i == m:
                    new, cost = (mm, dg, rg), dp_score[i][j-1]
                else:
                    new, cost = (mm, dg, rg+1), dp_score[i][j-1] + 1
                if new[2] <= max_RNA_gaps:
                    candidates.append((cost, new, 'left', i, j-1))

            if not candidates:
                continue

            best = min(c[0] for c in candidates)
            dp_score[i][j] = best
            cell = {}
            for cost, cnt, move, pi, pj in candidates:
                if cost == best:
                    cell.setdefault(cnt, []).append((move, pi, pj))
            dp_cells[i][j] = cell

    # Backtrack and recompute counts from the actual alignment strings
    result = []
    seen = set()

    def backtrack(i, j, cnt, a1, a2):
        if i == 0 and j == 0:
            aln1, aln2 = a1[::-1], a2[::-1]
            # recompute to avoid any propagation issues
            mismatches = sum(
                0 if iupac_match(x, y) else 1
                for x, y in zip(aln1, aln2)
                if x != '-' and y != '-'
            )
            dna_gaps = aln2.count('-')
            rna_gaps = aln1.rstrip('-').count('-')
            key = (aln1, aln2, mismatches, dna_gaps, rna_gaps)
            if key not in seen:
                seen.add(key)
                result.append((aln1, aln2, mismatches, rna_gaps, dna_gaps))
            return

        for move, pi, pj in dp_cells[i][j][cnt]:
            if move == 'diag':
                mis = 0 if iupac_match(RNA[i-1], DNA[j-1]) else 1
                prev = (cnt[0] - mis, cnt[1], cnt[2])
                backtrack(pi, pj, prev,
                          a1 + RNA[i-1], a2 + DNA[j-1])
            elif move == 'up':
                prev = (cnt[0], cnt[1] - 1, cnt[2])
                backtrack(pi, pj, prev,
                          a1 + RNA[i-1], a2 + '-')
            else:  # left
                if i == m:
                    prev = cnt
                else:
                    prev = (cnt[0], cnt[1], cnt[2] - 1)
                backtrack(pi, pj, prev,
                          a1 + '-', a2 + DNA[j-1])

    # start backtracking for each final count
    for final_cnt, ptrs in dp_cells[m][n].items():
        backtrack(m, n, final_cnt, "", "")

    return result or None

def clear_csv(path):
    with open(path, 'w'):
        pass

def scan(
    RNA, DNA, Strand,output_file,
    CHR=22,
    max_mismatches=4,
    max_DNA_gaps=1,
    max_RNA_gaps=1
):
    
    with open(output_file,'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        if Strand == "+": 
            writer.writerow([
            "CHR","RNA","DNA","Strand",
            "Start","END","mismatches",
            "gaps_in_RNA","gaps_in_DNA"
            ])
        L = len(RNA)
        for start in range(0, len(DNA) - L + max_DNA_gaps + 1):
            if start%10000==0:print(start)
            window = DNA[start : start + L + max_RNA_gaps]
            if "N" in window:
                continue
            hits = needleman_wunsch_all_alignments(
                RNA, window,
                max_mismatches, max_DNA_gaps, max_RNA_gaps
            )
            if not hits:
                continue
            for aRNA, aDNA, mm, rg, dg in hits:
                tl = len(aRNA.rstrip('-'))
                if Strand == "+":
                    writer.writerow([
                        f"chr{CHR}", aRNA[:tl], aDNA[:tl],
                        Strand, start, start + tl,
                        mm, rg, dg
                    ])
                else:
                    writer.writerow([
                        f"chr{CHR}", aRNA[:tl], aDNA[:tl],
                        Strand,
                        len(DNA) - start - tl + dg,
                        len(DNA) - start + dg,
                        mm, rg, dg
                    ])

if __name__ == "__main__":
    #Path to file
    record = SeqIO.read("chr22.enriched.fa", "fasta")
    DNA = str(record.seq)
    #Write guide
    Guide = "CTAACAGTTGCTTTTATCACNGG"
    output_file = "result.csv"
    clear_csv(output_file)
    scan(Guide, DNA, Strand="+",output_file=output_file, CHR=22)
    scan(Guide, reverse_complement(DNA), Strand="-", output_file=output_file,CHR=22)

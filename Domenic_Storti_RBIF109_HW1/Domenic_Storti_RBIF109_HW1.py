#!/usr/bin/env python3
# RBIF109 - Week 1 Homework
# Author: Domenic Storti

import random


# ---------- Tools (format and debugging) ----------

def _clean(seq):
    """Make uppercase and strip spaces and newlines so that comparisons are in same format."""
    return "".join(str(seq).upper().split())


def reverse_complement(dna):
    """Return reverse complement of a DNA string using dictionary and loop.
    """
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    dna = _clean(dna)
    rev = ""
    for i in range(len(dna) - 1, -1, -1):
        base = dna[i]
        rev += comp.get(base, 'N')  # Use 'N' for unknown
    return rev


# ---------- Q1: Edit distance ----------

def edit_distance_equal_length(seq_a, seq_b):
    """Computes the number of positions where the two sequences differ.
    Raise error if lengths differ.
    Utilizes for-loop over indices and keeps mismatch count.
    """
    a = _clean(seq_a)
    b = _clean(seq_b)
    if len(a) != len(b):
        raise ValueError("Sequences must be the same length")
    d = 0
    for i in range(len(a)):
        if a[i] != b[i]:
            d += 1
    return d


# ---------- Q2: in silico PCR ----------

def in_silico_pcr(template_dna, primer_fwd, primer_rev):
    """ Finds the forward primer directly on the template, then finds the reverse complement of the reverse
    primer on the same template, and checks that the forward site occurs before the reverse-complement site.
    If both are present in the correct orientation, it returns the substring that includes both primer
    sequences.
    Otherwise, returns None.
    """
    template = _clean(template_dna)
    fwd = _clean(primer_fwd)
    rev = _clean(primer_rev)

    fwd_index = template.find(fwd)
    if fwd_index == -1:
        return None

    rev_rc = reverse_complement(rev)
    rev_index = template.find(rev_rc)
    if rev_index == -1:
        return None

    if rev_index <= fwd_index:
        return None

    end_index = rev_index + len(rev_rc)
    return template[fwd_index:end_index]


# ---------- Q3: Most frequent k-mer in list of sequences ----------

def most_frequent_kmer_in_list(seq_list, k):
    """ Counts all length‑k substrings across an entire list of sequences.
    Returns the set of most frequent k‑mers, plus the whole counts dictionary for reference.
    Uses a dictionary to tally counts and simple loops to iterate k‑length across each string.
    """
    counts = {}
    for seq in seq_list:
        s = _clean(seq)
        if k > len(s):
            continue
        for i in range(0, len(s) - k + 1):
            kmer = s[i:i + k]
            if kmer not in counts:
                counts[kmer] = 1
            else:
                counts[kmer] += 1

    if not counts:
        return [], counts

    max_count = 0
    for v in counts.values():
        if v > max_count:
            max_count = v

    winners = []
    for kmer, v in counts.items():
        if v == max_count:
            winners.append((kmer, v))

    winners.sort(key=lambda x: x[0])
    return (winners, counts)


# ---------- Q4: DNA -> Protein ----------

GENETIC_CODE = {
    # Standard genetic code
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*",
    "TGC": "C", "TGT": "C", "TGA": "*", "TGG": "W"
}

STOP_CODONS = {"TAA", "TAG", "TGA"}
START = "ATG"


def translate_first_orf(dna):
    """Scans for the first 'ATG' (start) in frame 0 on the + strand, then translates codons in steps of
    3 until a stop codon ('TAA','TAG','TGA') or end of string.
    Returns the amino-acid string. If no complete ORF is found, returns ''.
    Uses a provided dictionary for the standard genetic code.
    """
    s = _clean(dna)
    start = s.find(START)
    if start == -1:
        return ""

    protein = ""
    i = start
    while i + 2 < len(s):
        codon = s[i:i + 3]
        aa = GENETIC_CODE.get(codon, "X")  # X = Unknown codon
        if codon in STOP_CODONS:
            return protein
        protein += aa
        i += 3
    return ""


def simulate_random_dna(n):
    """Generates a random DNA string of length n using A/T/G/C only"""
    alphabet = ['A', 'T', 'G', 'C']
    out = []
    for _ in range(n):
        out.append(random.choice(alphabet))
    return "".join(out)


# ---------- Demo / debugging ----------

def demo():
    print("----------")
    print("RBIF109 Week 1 Homework - Domenic Storti")
    print("----------")

    # Q1 demo
    seq_a = "ATGCC"
    seq_b = "ACGCT"
    d = edit_distance_equal_length(seq_a, seq_b)
    print(f"[Q1] Edit distance for {seq_a} vs {seq_b} -> {d}")

    # Q2 demo
    template = "AAACTGATGCGTACGTTAGCATGACCTAGGCTAACGTTT"
    fwd = "ATGCG"
    rev = "AACT"
    amplicon = in_silico_pcr(template, fwd, rev)
    print(f"[Q2] Amplicon -> {amplicon if amplicon else 'No product'}")

    # Q3 demo
    seqs = ["ATGCATGC", "CATGCGGG", "ATGCCCATG"]
    winners, counts = most_frequent_kmer_in_list(seqs, k=3)
    print(f"[Q3] Most frequent 3-mers across {len(seqs)} sequences -> {winners}")

    # Q4 demo
    translations = []
    for i in range(10):
        dna = simulate_random_dna(200)
        prot = translate_first_orf(dna)
        translations.append((i + 1, dna, prot, len(prot)))
    found = sum(1 for (_, _, p, _) in translations if p != "")
    print(f"[Q4] Found {found} complete ORFs out of 10 random sequences (frame 0, + strand only).")
    return translations


if __name__ == "__main__":
    translations = demo()
    with open("Domenic_Storti_RBIF109_HW1_outputs.txt", "w") as fout:
        fout.write("RBIF109 - Week 1 Homework Outputs\n")
        fout.write("Student: Domenic Storti\n\n")

        # Q1
        fout.write("[Q1] Edit distance example with Seq_a = ATGCC, Seq_b = ACGCT -> 2\n\n")

        # Q2
        template = "AAACTGATGCGTACGTTAGCATGACCTAGGCTAACGTTT"
        fwd = "ATGCG"
        rev = "AACT"
        amplicon = in_silico_pcr(template, fwd, rev)
        fout.write(f"[Q2] in silico PCR (template={template}, fwd={fwd}, rev={rev})\n")
        fout.write(f"     Amplicon: {amplicon if amplicon else 'No product'}\n\n")

        # Q3
        seqs = ["ATGCATGC", "CATGCGGG", "ATGCCCATG"]
        winners, counts = most_frequent_kmer_in_list(seqs, k=3)
        fout.write(f"[Q3] Most frequent 3-mers across {len(seqs)} sequences\n")
        fout.write(f"     Winners: {winners}\n")
        fout.write(f"     Total unique 3-mers counted: {len(counts)}\n\n")

        # Q4
        found = sum(1 for (_, _, p, _) in translations if p != "")
        fout.write(f"[Q4] Translating first in-frame ORF (ATG..stop) on + strand, frame 0\n")
        fout.write(f"     Complete ORFs found: {found} / 10\n")
        for (idx, dna, prot, plen) in translations[:3]:  # keep it short in outputs file
            fout.write(f"     Example {idx}: protein length = {plen} aa; protein head = '{prot[:15]}'\n")

    print("----------")
    print("Outputs written to: Domenic_Storti_RBIF109_HW1_outputs.txt")
    print("----------")

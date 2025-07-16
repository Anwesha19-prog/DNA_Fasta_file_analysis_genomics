# # -*- coding: utf-8 -*-
# """
# Created on Wed Jul 16 09:27:03 2025

# @author: DELL
# """

from Bio import SeqIO

def find_orfs(seq, frame=0): #finds origin of open reading frames
    seq = seq[frame:]         #An ORF is a sequence of DNA that starts with a start codon and ends with a stop codon, potentially coding for a protein. 
    stop_codons = {'TAA', 'TAG', 'TGA'}
    orfs = []
    start_pos = None
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if start_pos is None:
            if codon == 'ATG':  # start codon
                start_pos = i
        else:
            if codon in stop_codons:
                orf_len = i + 3 - start_pos
                orfs.append((start_pos + frame + 1, orf_len))  # +1 for 1-based position
                start_pos = None
    # If ORF runs to end of seq without stop codon, consider it as well
    if start_pos is not None:
        orfs.append((start_pos + frame + 1, len(seq) - start_pos))
    return orfs

def count_repeats(sequences, repeat_length):
    repeats = {}
    for seq in sequences:
        s = str(seq)
        for i in range(len(s) - repeat_length + 1):
            kmer = s[i:i+repeat_length]
            repeats[kmer] = repeats.get(kmer, 0) + 1 #It retrieves the current count of the k-mer using .get(), adds 1 to it, and assigns the result back to
    return repeats

def analyze_fasta(fasta_file):
    record_count = 0
    longest_seq_len = 0
    shortest_seq_len = float('inf')
    longest_orf_rf2 = 0
    longest_orf_rf3_start = 0
    longest_orf_rf3_len = 0
    longest_orf_any_frame = 0

    seqs = []
    seq_dict = {}  # store seqs by ID for specific analysis
    repeats_6 = {}
    repeats_12 = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        record_count += 1
        seq = str(record.seq).upper()
        seq_len = len(seq)
        seqs.append(seq)
        seq_dict[record.id] = seq

        # Update longest and shortest seq lengths
        longest_seq_len = max(longest_seq_len, seq_len)
        shortest_seq_len = min(shortest_seq_len, seq_len)

        # Longest ORF in reading frame 2 (frame index 1)
        orfs_rf2 = find_orfs(seq, frame=1)
        if orfs_rf2:
            max_orf_rf2 = max(orfs_rf2, key=lambda x: x[1])[1]
            longest_orf_rf2 = max(longest_orf_rf2, max_orf_rf2)

        # Longest ORF in reading frame 3 and its start pos (frame index 2)
        orfs_rf3 = find_orfs(seq, frame=2)
        if orfs_rf3:
            max_orf_rf3 = max(orfs_rf3, key=lambda x: x[1])
            if max_orf_rf3[1] > longest_orf_rf3_len:
                longest_orf_rf3_len = max_orf_rf3[1]
                longest_orf_rf3_start = max_orf_rf3[0]

        # Longest ORF in any forward reading frame
        for frame in range(3):
            orfs = find_orfs(seq, frame=frame)
            if orfs:
                max_orf = max(orfs, key=lambda x: x[1])[1]
                longest_orf_any_frame = max(longest_orf_any_frame, max_orf)

    # Longest forward ORF in specific sequence id
    target_id = "gi|142022655|gb|EQ086233.1|16"
    longest_orf_target_seq = 0
    if target_id in seq_dict:
        target_seq = seq_dict[target_id]
        for frame in range(3):
            orfs = find_orfs(target_seq, frame=frame)
            if orfs:
                max_orf = max(orfs, key=lambda x: x[1])[1]
                longest_orf_target_seq = max(longest_orf_target_seq, max_orf)
    else:
        longest_orf_target_seq = None  # ID not found

    # Most frequent repeat of length 6 and its count
    repeats_6 = count_repeats(seqs, 6)
    max_count_6 = max(repeats_6.values())
    most_freq_6 = [k for k,v in repeats_6.items() if v == max_count_6]

    # Repeats of length 12, how many occur max times
    repeats_12 = count_repeats(seqs, 12)
    max_count_12 = max(repeats_12.values())
    max_12_count = sum(1 for v in repeats_12.values() if v == max_count_12)

    # Count occurrences of specified length-7 repeats
    repeats_7 = count_repeats(seqs, 7)
    candidates_7 = ["AATGGCA", "CATCGCC", "CGCGCCG", "TGCGCGC"]
    max_candidate_count = 0
    max_candidate = None
    for c in candidates_7:
        count_c = repeats_7.get(c, 0)
        if count_c > max_candidate_count:
            max_candidate_count = count_c
            max_candidate = c

    return {
        "record_count": record_count,
        "longest_seq_len": longest_seq_len,
        "shortest_seq_len": shortest_seq_len,
        "longest_orf_rf2": longest_orf_rf2,
        "longest_orf_rf3_start": longest_orf_rf3_start,
        "longest_orf_rf3_len": longest_orf_rf3_len,
        "longest_orf_any_frame": longest_orf_any_frame,
        "longest_orf_target_seq": longest_orf_target_seq,
        "most_freq_6_count": max_count_6,
        "max_12_count": max_12_count,
        "max_candidate_7": max_candidate,
        "max_candidate_7_count": max_candidate_count,
    }

# Usage example:
fasta_file = r"E:\Coursera_Certs\dna2.fasta"  # replace with your file path
results = analyze_fasta(fasta_file)

print(f"Number of records: {results['record_count']}")
print(f"Length of the longest sequence: {results['longest_seq_len']}")
print(f"Length of the shortest sequence: {results['shortest_seq_len']}")
print(f"Length of the longest ORF in reading frame 2: {results['longest_orf_rf2']}")
print(f"Starting position of the longest ORF in reading frame 3: {results['longest_orf_rf3_start']}")
print(f"Length of the longest ORF in reading frame 3: {results['longest_orf_rf3_len']}")
print(f"Length of the longest ORF in any forward reading frame: {results['longest_orf_any_frame']}")
print(f"Longest forward ORF length in { 'gi|142022655|gb|EQ086233.1|16' }: {results['longest_orf_target_seq']}")
print(f"Most frequent repeat of length 6 occurs {results['most_freq_6_count']} times")
print(f"Number of different 12-mers occurring max times ({results['max_12_count']}): {results['max_12_count']}")
print(f"Repeat of length 7 with max occurrences: {results['max_candidate_7']} (occurs {results['max_candidate_7_count']} times)")


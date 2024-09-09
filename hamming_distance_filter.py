#!/usr/bin/env python

from Bio import SeqIO
import matplotlib.pyplot as plt

def hamming_distance(seq1, seq2, count_gaps_Ns = True):
    # Check if the sequences are of equal length
    if len(seq1) != len(seq2):
        raise ValueError("Input sequences must have the same length")

    # Initialize the Hamming distance
    hamming_distance = 0

    # Iterate through the characters in the sequences
    for char1, char2 in zip(seq1, seq2):
        # Check if both characters are not gaps and are different
        if count_gaps_Ns:
            if char1 != char2:
                hamming_distance += 1
        else:
            if char1 != '-' and char2 != '-' and char1 != char2:
                if char1 != 'N' and char2 != 'N':
                    hamming_distance += 1

    return hamming_distance

# def hamming_distance(seq1, seq2):
#     return sum(1 for a, b in zip(seq1, seq2) if a != b)

def filter_sequences_by_threshold(records, reference_seq, threshold, count_gaps):
    filtered_records = []
    for record in records:
        distance = hamming_distance(reference_seq, str(record.seq), count_gaps_Ns = count_gaps)
        if distance <= threshold:
            filtered_records.append(record)
    return filtered_records

def plot_histogram(hamming_distances):
    plt.hist(hamming_distances, bins=20, edgecolor='black')
    plt.xlabel('Hamming Distance')
    plt.ylabel('Frequency')
    plt.title('Hamming Distance Histogram')
    plt.show()

def write_fasta_file(output_filename, records):
    with open(output_filename, 'w') as output_file:
        SeqIO.write(records, output_file, 'fasta')

def main(input_filename, reference_seq_name, threshold, output_filename, count_gaps):
    records = list(SeqIO.parse(input_filename, 'fasta'))
    
    reference_record = next(record for record in records if record.id == reference_seq_name)
    reference_seq = str(reference_record.seq)

    filtered_records = filter_sequences_by_threshold(records, reference_seq, threshold, count_gaps)
    write_fasta_file(output_filename, filtered_records)
    
    print(f'Total sequences: {len(records)}')
    print(f'Sequences below threshold ({threshold}): {len(filtered_records)}')

if __name__ == '__main__':
    #input_filename = 'remove_N_perct_thresh_2_gisaid_GA_2020-12-01_to_2021-03-31_downloaded230501 assembled to GA-EHC-705E.cleaned.merged.fasta'  # Provide your input FASTA file name
    #input_filename = 'gisaid_GA_2020-12-01_to_2021-03-31_downloaded230501 assembled to GA-EHC-705E.cleaned.merged 231127.fasta'
    input_filename = 'gisaid_GA_2020-12-01_to_2021-03-31_downloaded230501 assembled to GA-EHC-705E.cleaned.merged_231127.fasta'
    reference_seq_name = 'hCoV-19/USA/GA-EHC-705E/2021|EPI_ISL_6913932|2021-02-01'  # Provide the name of the reference sequence in the FASTA file
    threshold = 10  # Define your Hamming distance threshold
    output_filename = 'hamming_ignoreN_705E_below_'+ str(threshold)+input_filename  # Provide desired output file name
    count_gaps = False

    main(input_filename, reference_seq_name, threshold, output_filename, count_gaps)

import os.path
from collections import Counter
from Bio import SeqIO


def statistical_sequence_count(seq_record):
    sequence_count = sequence_counter[seq_record.seq]
    return sequence_count


def compute_sequence_length_difference(seq_record, target_length):
    sequence_length = len(seq_record.seq)
    difference = abs(sequence_length - target_length)
    return difference


def sort_by_sequence_count(input_fasta_file, output_fasta_file):
    records = list(SeqIO.parse(input_fasta_file, "fasta"))
    global sequence_counter
    sequence_counter = Counter(record.seq for record in records)
    sorted_records = sorted(records, key=lambda x: statistical_sequence_count(x), reverse=True)
    with open(output_fasta_file, "w") as output_handle:
        for record in sorted_records:
            sequence_without_newline = str(record.seq).replace("\n", "")
            output_handle.write(f">{record.id}\n{sequence_without_newline}\n")


def sort_by_length_difference(input_fasta_file, output_fasta_file, target_length):
    records = list(SeqIO.parse(input_fasta_file, "fasta"))
    sorted_records = sorted(records, key=lambda x: compute_sequence_length_difference(x, target_length))
    with open(output_fasta_file, "w") as output_handle:
        for record in sorted_records:
            sequence_without_newline = str(record.seq).replace("\n", "")
            output_handle.write(f">{record.id}\n{sequence_without_newline}\n")


def sequencing_preprocessing(input_file, target_length):
    file_path, file_extension = os.path.splitext(input_file)
    sorted_file = file_path + '_sorted.fasta'
    sort_by_sequence_count(input_file, sorted_file)
    sort_by_length_difference(sorted_file, sorted_file, target_length)
    return sorted_file


def fasta_to_list(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        for line in file:
            line = line.rstrip()
            sequence = ''
            if line.startswith('>'):
                continue
            else:
                sequence += line
            sequences.append(sequence)
    return sequences


def read_file(original_file, sequencing_file, correct_length):
    original_sequences = fasta_to_list(original_file)
    sequencing_sequences = fasta_to_list(sequencing_file)
    min_length = correct_length - 10
    max_length = correct_length + 10
    reads = [seq for seq in sequencing_sequences if len(seq) >= min_length and len(seq) <= max_length]
    return original_sequences, reads


def write_fasta(final_clusters, reads, output_file):
    with open(output_file, 'w') as file:
        for i, cluster in enumerate(final_clusters):
            file.write(f'>Cluster{i + 1}\n')
            for seq in cluster:
                file.write(f'{reads[seq]}\n')


if __name__ == '__main__':
    input_file = "../datasets/test_data/test_cluster_50_with_labels.fasta"
    target_length = 110
    sorted_file = sequencing_preprocessing(input_file, target_length)
    print("Sort success!")
    print("path: ", sorted_file)

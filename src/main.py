from sequencing_preprocessing import *
from lsh_sketch import *
from fuzzy_clustering import *
from cluster_merging_refinement import *
from clusters_msa import *
from compute_performances import *
from tqdm import tqdm
import time


def main():
    encoding_file = '../datasets/test_data/test_origin_50.fasta'
    sequencing_file = '../datasets/test_data/test_cluster_50_with_labels.fasta'
    k = 14
    gap = 2
    drift = 2
    sketch_size = 5
    star_position = 0
    good_bad = 3
    correct_length = 110
    update_match_threshold = 3 / sketch_size
    span = (2 * drift + k) * sketch_size + gap * (sketch_size - 1)
    end_position = star_position + span
    output_file = '../results/clustering_result.txt'
    result_file = '../results/result.txt'
    print('****************************************Starting the main programme!***************************************')

    s1 = time.time()
    sequencing_sorted_file = sequencing_preprocessing(sequencing_file, correct_length)
    s2 = time.time()
    original_sequences, reads = read_file(encoding_file, sequencing_sorted_file, correct_length)
    print("Data read successfully!")

    s3 = time.time()
    reads_lsh_sketches = []
    for seq in tqdm(reads, desc='Generate hash sketches', ncols=100):
        read_lsh_index = segment_lsh_index(seq[star_position:end_position], k, gap, drift)
        reads_lsh_sketches.append(read_lsh_index)
    s4 = time.time()
    diff_list = generate_difference_list(k)
    print("Sketches build successfully!")

    s5 = time.time()
    clusters_dict = clustering_variable(reads_lsh_sketches, diff_list, sketch_size, drift)
    s6 = time.time()
    print("Sequences clustering successfully!")

    good_clusters = [[c] + list(set(clusters_dict[c])) for c in clusters_dict if len(clusters_dict[c]) >= good_bad]
    bad_clusters = [[c] + list(set(clusters_dict[c])) for c in clusters_dict if len(clusters_dict[c]) < good_bad]
    s7 = time.time()
    for cluster in tqdm(good_clusters, desc='Update representation', ncols=100):
        cluster = update_representation(cluster, reads_lsh_sketches, diff_list, update_match_threshold)
    s8 = time.time()
    merge_clusters = cluster_merging(good_clusters, reads_lsh_sketches, diff_list, sketch_size, drift)
    s9 = time.time()
    # final_clusters = good_clusters + bad_clusters  # not merging and refinement
    # final_clusters = merge_clusters + bad_clusters  # only merging
    # final_clusters = cluster_refinement(good_clusters, bad_clusters, reads_lsh_sketches, diff_list, sketch_size, drift)  # only refinement
    final_clusters = cluster_refinement(merge_clusters, bad_clusters, reads_lsh_sketches, diff_list, sketch_size,
                                        drift)
    s10 = time.time()
    write_fasta(final_clusters, reads, output_file)
    print("Cluster update representative and merging successfully!")

    s11 = time.time()
    msa_results = clusters_msa(final_clusters, reads, 15)
    s12 = time.time()
    candidates = generate_candidates(msa_results)
    print("MSA and candidates generate successfully!")

    true_labels = getting_true_labels(sequencing_sorted_file, correct_length)  # true labels
    pred_labels = getting_cluster_labels(final_clusters, len(reads))  # cluster labels
    nmi = compute_NMI(true_labels, pred_labels)
    ami = compute_AMI(true_labels, pred_labels)
    ari = compute_ARI(true_labels, pred_labels)
    fmi = compute_FMI(true_labels, pred_labels)
    homo = compute_HOMO(true_labels, pred_labels)
    comp = compute_COMP(true_labels, pred_labels)
    # v_measure = compute_V_measure(true_labels, pred_labels)
    # purity = compute_purity(true_labels, pred_labels)
    # accuracy = compute_accuracy(true_labels, pred_labels)
    # intra_distances = compute_intra_cluster_levenshtein(final_clusters, reads)
    # intra_distances = sum(intra_distances) / len(intra_distances)
    # inter_distances = compute_inter_cluster_levenshtein(final_clusters, reads)
    # inter_distances = sum(inter_distances) / len(inter_distances)
    fraction, redundancy = fraction_recovered(candidates, original_sequences)
    re_rate = reconstruction_rate(candidates, original_sequences)

    with open(result_file, 'w') as rfile:
        rfile.write(f'The number of original sequences is: {len(original_sequences)}\n')
        rfile.write(f'The number of sequencing sequences is: {len(reads)}\n')
        rfile.write(f'The number of Good Clusters is: {len(good_clusters)}\n')
        rfile.write(f'The number of Bad Clusters is: {len(bad_clusters)}\n')
        rfile.write(f'The number of merged Clusters is: {len(merge_clusters)}\n')
        rfile.write(f'The number of Final Clusters is: {len(final_clusters)}\n')
        rfile.write(f'Clustering time: {round((s6 - s5) + (s8 - s7) + (s9 - s8) + (s10 - s9), 2)} s \n')
        rfile.write(f'MSA time: {round(s12 - s11, 2)} s \n')
        rfile.write(f'NMI: {round(nmi, 3)}\n')
        rfile.write(f'AMI: {round(ami, 3)}\n')
        rfile.write(f'ARI: {round(ari, 3)}\n')
        rfile.write(f'FMI: {round(fmi, 3)}\n')
        rfile.write(f'HOMO: {round(homo, 3)}\n')
        rfile.write(f'COMP: {round(comp, 3)}\n')
        rfile.write(f'Sequence recovery rate: {round(fraction, 4)}\n')
        rfile.write(f'Sequence reconstruction rate: {round(re_rate, 4)}\n')

    print('Sucessfully!')


if __name__ == '__main__':
    main()

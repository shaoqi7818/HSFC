import time
from sequencing_preprocessing import *
from lsh_sketch import *
from fuzzy_clustering import *
from cluster_merging_refinement import *
from clusters_msa import *
from compute_performances import *
from tqdm import tqdm


def clustering_reconstruction(encoding_file, sequencing_file):
    k = 14
    gap = 5
    drift = 2
    sketch_size = 6
    star_position = 34
    good_bad = 3
    correct_length = 200
    update_match_threshold = 3 / sketch_size
    span = (2 * drift + k) * sketch_size + gap * (sketch_size - 1)
    end_position = star_position + span
    print('1.Sequencing preprocessing AND Read sequencing files')
    sequencing_sorted_file = sequencing_preprocessing(sequencing_file, correct_length)
    original_sequences, reads = read_file(encoding_file, sequencing_sorted_file, correct_length)
    print('2.Generate hash signatures for each sequence')
    reads_lsh_index = [segment_lsh_index(seq[star_position:end_position], k, gap, drift) for seq in reads]
    diff_list = generate_difference_list(k)
    reads_lsh_sketches = list_to_tuple(reads_lsh_index)
    print('3.Greedy clustering')
    clusters_dict = clustering_variable(reads_lsh_sketches, diff_list, sketch_size, drift)
    print('4.Update clusters representative and Cluster merging')
    good_clusters = [[c] + list(set(clusters_dict[c])) for c in clusters_dict if len(clusters_dict[c]) >= good_bad]
    bad_clusters = [[c] + list(set(clusters_dict[c])) for c in clusters_dict if len(clusters_dict[c]) < good_bad]
    for cluster in tqdm(good_clusters, desc='Update representation', ncols=100):
        cluster = update_representation(cluster, reads_lsh_sketches, diff_list, update_match_threshold)
    merge_clusters = cluster_merging(good_clusters, reads_lsh_sketches, diff_list, sketch_size, drift)
    final_clusters = cluster_refinement(merge_clusters, bad_clusters, reads_lsh_sketches, diff_list, sketch_size, drift)
    print('5.Multiple Sequence Alignment AND Majority voting for Candidates')
    msa_results = clusters_msa(final_clusters, reads, 15)
    candidates = generate_candidates(msa_results)
    return original_sequences, candidates


if __name__ == '__main__':
    encoding_file = '../datasets/simulate_data/origin.fasta'
    simulate_error_path = '../datasets/simulate_data'
    error_rates = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
    iterations = 10
    result_recovery_fraction = {}
    result_reconstruction_rate = {}
    average_result_recovery_fraction = {}
    average_result_reconstruction_rate = {}
    total_time = []
    for error_rate in tqdm(error_rates, desc='Calculate performance', ncols=100):
        recovery_fractions = []
        reconstruction_rates = []
        for it in range(iterations):
            print(f'----------------开始执行————错误率：{error_rate}，第{it}次迭代----------------')
            sequencing_file = f'{simulate_error_path}/error_{error_rate}/simulate_error{error_rate}_iteration{it + 1}.fasta'

            t1 = time.time()
            original_sequences, candidates = clustering_reconstruction(encoding_file, sequencing_file)
            t2 = time.time()
            total_time.append(round(t2 - t1, 2))

            recovery, redundancy = fraction_recovered(candidates, original_sequences)
            reconstruction = reconstruction_rate(candidates, original_sequences)
            recovery_fractions.append(round(recovery, 4))
            reconstruction_rates.append(round(reconstruction, 4))
        average_recovery_fraction = sum(recovery_fractions) / iterations
        average_reconstruction_rate = sum(reconstruction_rates) / iterations
        average_result_recovery_fraction[error_rate] = round(average_recovery_fraction, 4)
        average_result_reconstruction_rate[error_rate] = round(average_reconstruction_rate, 4)
        result_recovery_fraction[error_rate] = recovery_fractions
        result_reconstruction_rate[error_rate] = reconstruction_rates

    print('total time:', sum(total_time), total_time)

    print('---------------------------------------------')
    for key, value in result_recovery_fraction.items():
        print(f'{key}错误率下的序列回收分数为：{value}')
    print('---------------------------------------------')
    for key, value in result_reconstruction_rate.items():
        print(f'{key}错误率下的序列重建率为：{value}')

    print('---------------------------------------------')
    for key, value in average_result_recovery_fraction.items():
        print(f'{key}错误率下的平均序列回收分数为：{value}')
    print('---------------------------------------------')
    for key, value in average_result_reconstruction_rate.items():
        print(f'{key}错误率下的平均序列重建率为：{value}')

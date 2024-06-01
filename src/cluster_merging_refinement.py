from tqdm import tqdm
from collections import Counter


def fuzzy_match(seq_lsh_index_1, seq_lsh_index_2, difference_list):
    match_list = []
    for subseq_index_1, subseq_index_2 in zip(seq_lsh_index_1, seq_lsh_index_2):
        sub_difference = [abs(index_1 - index_2) for index_1 in subseq_index_1 for index_2 in subseq_index_2]
        similarity = set(sub_difference) & set(difference_list)
        if similarity:
            match_list.append(1)
        else:
            match_list.append(0)
    return match_list


def calculate_jaccard_index(match_list):
    jaccard_index = match_list.count(1) / len(match_list)
    return jaccard_index


def most_frequent_element(lst):
    filtered_list = [x for x in lst if x != -1]
    counter = Counter(filtered_list)
    most_common_element = counter.most_common(1)
    return most_common_element[0][0] if most_common_element else None


def update_representation(cluster, reads_lsh_sketches, diff_list, update_match_threshold):
    match_num = {}
    for c in cluster:
        match_num[c] = 0

    for i in range(0, len(cluster)):
        for j in range(i + 1, len(cluster)):
            match_list = fuzzy_match(reads_lsh_sketches[cluster[i]], reads_lsh_sketches[cluster[j]], diff_list)
            if calculate_jaccard_index(match_list) >= update_match_threshold:
                match_num[cluster[i]] += 1
                match_num[cluster[j]] += 1

    max_key = max(match_num, key=match_num.get)
    max_index = cluster.index(max_key)
    cluster[0], cluster[max_index] = cluster[max_index], cluster[0]
    return cluster


def cluster_merging(good_clusters, reads_lsh_sketches, diff_list, sketch_size, drift):
    lsh_dict_list = [{} for _ in range(sketch_size)]
    merge_clusters = {}
    for index, good in enumerate(tqdm(good_clusters, desc='Cluster merging', ncols=100)):
        sketch_diff = []
        for signs in reads_lsh_sketches[good[0]]:
            signs_diff = []
            for sig in signs:
                sig_diff = set([(sig + diff) for diff in diff_list] + [(sig - diff) for diff in diff_list])
                signs_diff += sig_diff
            sketch_diff.append(signs_diff)

        cluster_flag = []
        for i, signs_diff in enumerate(sketch_diff):
            for _ in signs_diff:
                if _ in lsh_dict_list[i]:
                    cluster_flag.append(lsh_dict_list[i][_])
                    break
            else:
                cluster_flag.append(-1)
                continue

        max_flag = most_frequent_element(cluster_flag)

        if max_flag:
            merge_clusters[max_flag] = merge_clusters[max_flag] + good
        else:
            merge_clusters[index] = good
            for i in range(sketch_size):
                lsh_dict_list[i][reads_lsh_sketches[good[0]][i][drift]] = index

    merged_clusters = list(merge_clusters.values())
    return merged_clusters


def cluster_refinement(good_clusters, bad_clusters, reads_lsh_sketches, diff_list, sketch_size, drift):
    lsh_dict_list = [{} for _ in range(sketch_size)]
    for ind, good in enumerate(good_clusters):
        for i, signs in enumerate(reads_lsh_sketches[good[0]]):
            lsh_dict_list[i][signs[drift]] = ind

    remain_clusters = []
    for index, bad in enumerate(tqdm(bad_clusters, desc='Cluster refinement', ncols=100)):
        sketch_diff = []
        for signs in reads_lsh_sketches[bad[0]]:
            signs_diff = []
            for sig in signs:
                sig_diff = set([(sig + diff) for diff in diff_list] + [(sig - diff) for diff in diff_list])
                signs_diff += sig_diff
            sketch_diff.append(signs_diff)

        cluster_flag = []
        for i, signs_diff in enumerate(sketch_diff):
            for _ in signs_diff:
                if _ in lsh_dict_list[i]:
                    cluster_flag.append(lsh_dict_list[i][_])
                    break
            else:
                cluster_flag.append(-1)
                continue

        max_flag = most_frequent_element(cluster_flag)

        if max_flag:
            good_clusters[max_flag] += bad
        else:
            remain_clusters.append(bad)

    refine_clusters = good_clusters + remain_clusters
    return refine_clusters

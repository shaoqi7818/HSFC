from collections import Counter
from tqdm import tqdm


def most_frequent_element(lst):
    filtered_list = [x for x in lst if x != -1]
    counter = Counter(filtered_list)
    most_common_element = counter.most_common(1)
    return most_common_element[0][0] if most_common_element else None


def clustering_variable(reads_lsh_sketches, diff_list, sketch_size, drift):
    lsh_dict_list = [{} for _ in range(sketch_size)]
    clusters_dict = {}
    update_flag = {}
    for index, read_sketch in enumerate(tqdm(reads_lsh_sketches, desc='Clustering', ncols=100)):
        sketch_diff = []
        for signs in read_sketch:
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

        if max_flag is None:
            clusters_dict[index] = []
            update_flag[index] = 0
            for i in range(sketch_size):
                lsh_dict_list[i][read_sketch[i][drift]] = index
        else:
            clusters_dict[max_flag] += [index]

    return clusters_dict

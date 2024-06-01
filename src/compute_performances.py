import re
import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn import metrics
import Levenshtein
from tqdm import tqdm


def levenshtein_distance(seq1, seq2):
    return Levenshtein.distance(seq1, seq2)


def compute_intra_cluster_levenshtein(clusters, reads):
    average_distances = []
    for cluster in clusters:
        distance = []
        for i in range(len(cluster)):
            for j in range(i + 1, len(cluster)):
                distance.append(levenshtein_distance(reads[cluster[i]], reads[cluster[j]]))
        if len(distance) != 0:
            average_distance = sum(distance) / len(distance)
        else:
            average_distance = 0
        average_distances.append(average_distance)
    return average_distances


def compute_inter_cluster_levenshtein(clusters, reads):
    average_distances = []
    for i in range(len(clusters)):
        distance = 0
        for j in range(len(clusters)):
            if i != j:
                distance += levenshtein_distance(reads[clusters[i][0]], reads[clusters[j][0]])
        average_distance = distance / len(clusters) - 1
        average_distances.append(average_distance)
    return average_distances


def compute_SI(x, labels):
    si = metrics.silhouette_score(x, labels)
    return si


def compute_CHI(x, labels):
    chi = metrics.calinski_harabasz_score(x, labels)
    return chi


def compute_NCC(x, labels):
    m = x.shape[0]
    n = x.shape[1]
    Y = np.zeros((m, m))
    for r in range(m):
        for s in range(m):
            if labels[r] == labels[s]:
                Y[r, s] = 1
    drs = np.zeros((m, m))
    for r in range(m):
        for s in range(m):
            for att in range(n):
                if x[r, att] != x[s, att]:
                    drs[r, s] += 1
    ncc = 0.0
    for r in range(m):
        for s in range(m):
            if r != s:
                ncc += (n - 2 * drs[r, s]) * Y[r, s] + drs[r, s]
    return ncc


def getting_cluster_labels(clusters, length):
    pred_dict = {}
    for i, clu in enumerate(clusters):
        for cl in clu:
            pred_dict[cl] = i
    pred_labels = []
    i = 0
    j = -1
    while i <= length - 1:
        if i in pred_dict:
            pred_labels.append(pred_dict[i])
            i += 1
        else:
            pred_labels.append(j)
            j -= 1
            i += 1
    return pred_labels


def getting_true_labels(input_file, correct_length):
    true_labels = []
    min_length = correct_length - 10
    max_length = correct_length + 10
    with open(input_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                cluster_label = int(re.search(r'Cluster(\d+)', line).group(1))
                true_labels.append(cluster_label)
            else:
                if min_length <= len(line) <= max_length:
                    continue
                else:
                    true_labels.pop()
                    continue
    return true_labels


def compute_NMI(true_labels, pred_labels):
    nmi = metrics.normalized_mutual_info_score(true_labels, pred_labels)
    return nmi


def compute_AMI(true_labels, pred_labels):
    ami = metrics.adjusted_mutual_info_score(true_labels, pred_labels)
    return ami


def compute_ARI(true_labels, pred_labels):
    ari = metrics.adjusted_rand_score(true_labels, pred_labels)
    return ari


def compute_FMI(true_labels, pred_labels):
    fmi = metrics.fowlkes_mallows_score(true_labels, pred_labels)
    return fmi


def compute_HOMO(true_labels, pred_labels):
    homo = metrics.homogeneity_score(true_labels, pred_labels)
    return homo


def compute_COMP(true_labels, pred_labels):
    comp = metrics.completeness_score(true_labels, pred_labels)
    return comp


def compute_V_measure(true_labels, pred_labels):
    v_measure = metrics.v_measure_score(true_labels, pred_labels)
    return v_measure


def compute_purity(true_labels, pred_labels):
    clusters = np.unique(pred_labels)
    true_labels = np.reshape(true_labels, (-1, 1))
    pred_labels = np.reshape(pred_labels, (-1, 1))
    count = []

    for c in clusters:
        idx = np.where(pred_labels == c)[0]
        labels_tmp = true_labels[idx, :].reshape(-1)
        count.append(np.bincount(labels_tmp).max())

    return np.sum(count) / true_labels.shape[0]


def compute_accuracy(true_labels, pred_labels):
    unique_true_labels = np.unique(true_labels)
    unique_predicted_labels = np.unique(pred_labels)

    contingency_matrix = np.zeros((len(unique_true_labels), len(unique_predicted_labels)))

    for i, true_label in enumerate(unique_true_labels):
        true_mask = (true_labels == true_label)
        for j, predicted_label in enumerate(unique_predicted_labels):
            predicted_mask = (pred_labels == predicted_label)
            contingency_matrix[i, j] = np.sum(np.logical_and(true_mask, predicted_mask))

    row_ind, col_ind = linear_sum_assignment(-contingency_matrix)

    accuracy = contingency_matrix[row_ind, col_ind].sum() / len(true_labels)

    return accuracy


def fraction_recovered(candidates, orig_seqs):
    count_dict = {}
    redundancy = 0

    for seq in orig_seqs:
        count_dict[seq] = 0
    for cand in candidates:
        if cand in count_dict:
            count_dict[cand] += 1

    fraction = sum([count_dict[seq] > 0 for seq in count_dict]) / len(count_dict)
    if fraction > 0:
        redundancy = sum([count_dict[seq] for seq in count_dict]) / len(count_dict) / fraction

    return fraction, redundancy


def reconstruction_rate(candidates, orig_seqs):
    reconstruct_dict = {}
    for seq in orig_seqs:
        reconstruct_dict[seq] = 0
    for cand in tqdm(candidates, desc='Reconstruction rate', ncols=100):
        if cand in reconstruct_dict:
            reconstruct_dict[cand] += 1
        else:
            min_distance = float('inf')
            target_seq = None
            for seq in orig_seqs:
                distance = levenshtein_distance(cand, seq)
                if distance < min_distance:
                    min_distance = distance
                    target_seq = seq
            re_rate = 1 - min_distance / len(target_seq)
            reconstruct_dict[target_seq] += re_rate

    average_re_rate = sum([reconstruct_dict[seq] for seq in reconstruct_dict]) / len(candidates)

    return average_re_rate


if __name__ == "__main__":
    print(1)

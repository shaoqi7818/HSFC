def list_to_tuple(lst):
    return tuple(tuple(tuple(inner) for inner in middle) for middle in lst)


def generate_drifts_list(x):
    list = []
    i = 0
    k = -x - 1
    while True:
        k = k + 1
        list.append(k)
        i += 1
        if k == x:
            break
    return list


def segment_lsh_index(sequence, k, gap, drift):
    seq = []
    seq_index = []
    for i in range(drift, len(sequence) - k + 1, k + gap + 2 * drift):
        subseq = []
        for d in generate_drifts_list(drift):
            subseq_ = sequence[i + d: i + d + k]
            subseq.append(subseq_)
        seq.append(subseq)
        subseq_index = []
        for subseq_unit in subseq:
            index = 0
            for j, base in enumerate(subseq_unit):
                if base == 'A':
                    x = 0
                elif base == 'C':
                    x = 1
                elif base == 'G':
                    x = 2
                elif base == 'T':
                    x = 3
                else:
                    index = -1
                    break
                index += x * (4 ** j)
            subseq_index.append(index)
        seq_index.append(subseq_index)
    return seq_index


def generate_difference_list(k):
    difference_list = [0]
    for i in range(k):
        result = [x * 4 ** i for x in [1, 2, 3]]
        difference_list.extend(result)
    return difference_list


if __name__ == '__main__':
    print(1)

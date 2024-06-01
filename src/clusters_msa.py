import operator
import subprocess
from random import shuffle
from Bio import AlignIO
from tqdm import tqdm


def muscle(cluster):
    input_alignment = "../output/msa/clm.fasta"
    output_alignment = "../output/msa/clmout.fasta"

    file = open(input_alignment, "w")
    for i, c in enumerate(cluster):
        file.write(">S%d\n" % i)
        file.write(c)
        file.write("\n")
    file.close()

    muscle_exe = "muscle5"
    muscle_command = [muscle_exe, "-align", input_alignment, "-output", output_alignment]
    # muscle_command = [muscle_exe, "-super5", input_alignment, "-output", output_alignment]
    try:
        subprocess.run(muscle_command, check=True)
        print("MUSCLE alignment completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running MUSCLE: {e}")

    msa_sequences = AlignIO.read(output_alignment, "fasta")
    aligned_cluster = []
    for i in msa_sequences:
        aligned_cluster += [i.seq]
    return aligned_cluster


def clusters_msa(clusters, reads, maxsize=15):
    msa_results = []
    for i, cluster_indexes in enumerate(tqdm(clusters, ncols=100)):
        cluster = [reads[index] for index in cluster_indexes]
        if len(cluster) < 5:
            continue
        elif len(cluster) <= maxsize:
            aligned_cluster = muscle(cluster)
            msa_results.append(aligned_cluster)
        else:
            for j in range(5):
                shuffle(cluster)
                aligned_cluster = muscle(cluster[:maxsize])
                msa_results.append(aligned_cluster)
        if i % 1000 == 0:
            print("%", round(i * 100 / len(clusters), 2), "of the clusters are aligned.")
    return msa_results


def sequence_voting(msa, weight=0.4):
    res = ""
    for i in range(len(msa[0])):
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, '-': 0, 'N': 0}
        for j in range(len(msa)):
            counts[msa[j][i]] += 1
        counts['-'] *= weight
        mv = max(counts.items(), key=operator.itemgetter(1))[0]
        if mv != '-':
            res += mv
    return res


def generate_candidates(msa_results):
    candidates = []
    for msa in msa_results:
        candidates.append(sequence_voting(msa, weight=0.5))
    return candidates


if __name__ == "__main__":
    pass

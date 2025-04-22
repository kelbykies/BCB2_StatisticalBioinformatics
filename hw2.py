
import numpy as nm
from Bio import SeqIO

# 1.b
def gen_next_base(seq, t_prob):
    """

    :param seq: sequence that is being generated
    :param t_prob: Is a row in the transition matrix; a 1x4 vector
    :return: seq
    """
    y = nm.random.uniform(0, 1)
    if y < t_prob[0]:
        seq = seq + "A"
    elif y < (t_prob[0] + t_prob[1]):
        seq = seq + "T"
    elif y < (t_prob[0] + t_prob[1] + t_prob[2]):
        seq = seq + "G"
    else:
        seq = seq + "C"

    return seq


# Transition Matrix
p = [[0.4, 0.2, 0.2, 0.2],
     [0.2, 0.4, 0.2, 0.2],
     [0.2, 0.2, 0.4, 0.2],
     [0.2, 0.2, 0.2, 0.4]]

alpha = [.95, 0.05, 0, 0]

generated_seq = ""

for i in range(0, 10000):
    # initial state
    if i == 0:
        # Generate number from [0,1], x
        x = nm.random.uniform(0, 1)
        if x < alpha[0]:
            generated_seq = generated_seq + "A"
        elif x < (alpha[0] + alpha[1]):
            generated_seq = generated_seq + "T"
        elif x < (alpha[0] + alpha[1] + alpha[2]):
            generated_seq = generated_seq + "G"
        else:
            generated_seq = generated_seq + "C"
    else:
        # Find the previous seq base
        prev_base = generated_seq[i-1]
        y = nm.random.uniform(0, 1)
        if prev_base == "A":
            # Check the first row
            generated_seq = gen_next_base(generated_seq, p[0])
        elif prev_base == "T":
            generated_seq = gen_next_base(generated_seq, p[1])
        elif prev_base == "G":
            generated_seq = gen_next_base(generated_seq, p[2])
        else:
            generated_seq = gen_next_base(generated_seq, p[3])
        if i % 70 == 69:
            generated_seq = generated_seq + "\n"


# Print to a FASTA
with open("Background_Seq.fa", "w+") as outfile:
    outfile.write(">Background_seq\n" + generated_seq)
outfile.close()


##########################################################


# 2a: Derivation can be seen on writeup
# MLE of p_hat(i) = number of i's / total number of nucleotides
# Where i can be A,C,G or T for each position of the motif.

import numpy as nm
from Bio import SeqIO


def count_nuc(motif, count_matrix):
    """
    Count the number of A's, C's, G's and T's that appear for each position in the motif
    :param motif: Nucleotide sequence
    :param count_matrix: Stores tallies of the amount of A's, C's, G's and T's
    :return: count_matrix
    """

    for i in range(0,len(motif)):
        nuc = motif[i]

        if nuc == "A":
            count_matrix[i][0] += 1

        elif nuc == "C":
            count_matrix[i][1] += 1
        elif nuc == "G":
            count_matrix[i][2] += 1
        else:
            count_matrix[i][3] += 1
    return count_matrix


# number of reads in Fasta File
total_reads = 0
motif_length = 9
# 2D matrix that will hold mles for paramters(pA, pC, pG, pT)
mle_matrix = nm.zeros((motif_length, 4))

# Loop through the Fasta file seq. by seq.
for record in SeqIO.parse("henri_motif.fa", "fasta"):
    total_reads += 1
    mle_matrix = count_nuc(record.seq, mle_matrix)

# Estimate parameters
mle_matrix = mle_matrix/total_reads

print(mle_matrix)



# 2b










import numpy as np

def needleman_wunsch(seq1, seq2, match_reward=2, mismatch_penalty=-1, gap_penalty=-1):
    # Create a matrix to store the scores
    n = len(seq1)
    m = len(seq2)
    score_matrix = np.zeros((n + 1, m + 1))

    # Initialize the scoring matrix with gap penalties
    for i in range(1, n + 1):
        score_matrix[i][0] = i * gap_penalty
    for j in range(1, m + 1):
        score_matrix[0][j] = j * gap_penalty

    # Fill the scoring matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = score_matrix[i - 1][j - 1] + (match_reward if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty
            score_matrix[i][j] = max(match, delete, insert)

    # Traceback to create the alignment
    align1, align2 = '', ''
    i, j = n, m

    while i > 0 and j > 0:
        score_current = score_matrix[i][j]
        score_diagonal = score_matrix[i - 1][j - 1]
        score_up = score_matrix[i - 1][j]
        score_left = score_matrix[i][j - 1]

        if score_current == score_diagonal + (match_reward if seq1[i - 1] == seq2[j - 1] else mismatch_penalty):
            align1 += seq1[i - 1]
            align2 += seq2[j - 1]
            i -= 1
            j -= 1
        elif score_current == score_up + gap_penalty:
            align1 += seq1[i - 1]
            align2 += '-'
            i -= 1
        elif score_current == score_left + gap_penalty:
            align1 += '-'
            align2 += seq2[j - 1]
            j -= 1

    # Add remaining characters/gaps if needed
    while i > 0:
        align1 += seq1[i - 1]
        align2 += '-'
        i -= 1
    while j > 0:
        align1 += '-'
        align2 += seq2[j - 1]
        j -= 1

    # Reverse the alignments since we traced back
    align1 = align1[::-1]
    align2 = align2[::-1]

    # Get final alignment score
    final_score = score_matrix[n][m]

    return score_matrix, align1, align2, final_score

# Test example
seq1 = "GATTACA"
seq2 = "GCATGCU"

# Run Needleman-Wunsch algorithm
score_matrix, alignment1, alignment2, final_score = needleman_wunsch(seq1, seq2)

print("Score Matrix:")
print(score_matrix)

print("\nAlignment:")
print(alignment1)
print(alignment2)

print("\nFinal Alignment Score:", final_score)

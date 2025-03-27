import numpy as np
import matplotlib.pyplot as plt

def generate_dot_plot(seq1, seq2):
    """Generates a dot plot for two sequences."""
    n = len(seq1)
    m = len(seq2)
    dot_matrix = np.zeros((n, m))

    # Fill the dot matrix
    for i in range(n):
        for j in range(m):
            if seq1[i] == seq2[j]:
                dot_matrix[i][j] = 1  # Mark a dot where sequences match

    return dot_matrix

def plot_dot_plot(dot_matrix, seq1, seq2):
    """Plots the dot plot using Matplotlib."""
    plt.figure(figsize=(10, 10))
    plt.imshow(dot_matrix, cmap='Greys', interpolation='nearest')
    plt.xticks(ticks=np.arange(len(seq2)), labels=list(seq2))
    plt.yticks(ticks=np.arange(len(seq1)), labels=list(seq1))
    plt.xlabel('Sequence 2')
    plt.ylabel('Sequence 1')
    plt.title('Dot Plot')
    plt.show()

def scoring_matrix(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """Calculates the score of an alignment between two sequences."""
    n = len(seq1)
    m = len(seq2)
    score_matrix = np.zeros((n+1, m+1))

    # Initialize the score matrix with gap penalties
    for i in range(1, n+1):
        score_matrix[i][0] = i * gap
    for j in range(1, m+1):
        score_matrix[0][j] = j * gap

    # Fill the score matrix
    for i in range(1, n+1):
        for j in range(1, m+1):
            if seq1[i-1] == seq2[j-1]:
                score = match  # Match score
            else:
                score = mismatch  # Mismatch penalty

            score_matrix[i][j] = max(
                score_matrix[i-1][j-1] + score,  # Diagonal (match/mismatch)
                score_matrix[i-1][j] + gap,      # Gap in seq2 (insertion)
                score_matrix[i][j-1] + gap       # Gap in seq1 (deletion)
            )

    return score_matrix

def print_alignment_score(score_matrix):
    """Prints the final alignment score."""
    print("Alignment Score:")
    print(f"Score = {score_matrix[-1][-1]}")

# Bigger example usage
seq1 = "ACTGACCTGAC"
seq2 = "AGTCCGACACTG"

# Generate the dot plot
dot_matrix = generate_dot_plot(seq1, seq2)

# Plot the dot plot
plot_dot_plot(dot_matrix, seq1, seq2)

# Generate the scoring matrix
score_matrix = scoring_matrix(seq1, seq2)

# Print the score
print_alignment_score(score_matrix)


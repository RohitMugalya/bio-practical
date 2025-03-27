def match_mismatch_score(seq1, seq2, match_score=1, mismatch_penalty=-1):
    """
    Function to compute the match/mismatch score between two sequences.

    Parameters:
    - seq1: First sequence (string)
    - seq2: Second sequence (string)
    - match_score: Score for a match (default = 1)
    - mismatch_penalty: Penalty for a mismatch (default = -1)

    Returns:
    - total_score: Total score based on match/mismatch
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length")

    total_score = 0

    # Compare each base pair in the aligned sequences
    for base1, base2 in zip(seq1, seq2):
        if base1 == base2:
            total_score += match_score  # Match score
        else:
            total_score += mismatch_penalty  # Mismatch penalty

    return total_score

# Example usage
seq1 = "AGCTTAGC"
seq2 = "AGCTCAGC"

# Calculate the match/mismatch score with default scoring scheme
score = match_mismatch_score(seq1, seq2)

print(f"Match/Mismatch score between '{seq1}' and '{seq2}': {score}")

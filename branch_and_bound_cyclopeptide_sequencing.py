# Dictionary of amino acid masses
amino_acid_masses = {
    'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99,
    'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
    'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131,
    'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186
}

def linear_spectrum(peptide):
    """Calculate the linear spectrum (prefix masses) of a peptide."""
    prefix_mass = [0]  # starting mass is 0
    for i in range(len(peptide)):
        prefix_mass.append(prefix_mass[i] + amino_acid_masses[peptide[i]])

    spectrum = [0]
    for i in range(len(prefix_mass)):
        for j in range(i+1, len(prefix_mass)):
            spectrum.append(prefix_mass[j] - prefix_mass[i])

    return sorted(spectrum)

def cyclic_spectrum(peptide):
    """Calculate the cyclic spectrum of a cyclic peptide."""
    n = len(peptide)
    prefix_mass = [0] * (n + 1)
    for i in range(n):
        prefix_mass[i+1] = prefix_mass[i] + amino_acid_masses[peptide[i]]

    peptide_mass = prefix_mass[n]
    spectrum = [0]

    # Calculate subpeptide masses
    for i in range(n):
        for j in range(i + 1, n + 1):
            spectrum.append(prefix_mass[j] - prefix_mass[i])
            if i > 0 and j < n:
                spectrum.append(peptide_mass - (prefix_mass[j] - prefix_mass[i]))

    return sorted(spectrum)

def is_consistent(peptide, spectrum):
    """Check if the linear spectrum of the peptide is consistent with the given spectrum."""
    peptide_spectrum = linear_spectrum(peptide)
    peptide_spectrum_counts = {x: peptide_spectrum.count(x) for x in peptide_spectrum}
    spectrum_counts = {x: spectrum.count(x) for x in spectrum}

    for mass in peptide_spectrum_counts:
        if peptide_spectrum_counts[mass] > spectrum_counts.get(mass, 0):
            return False
    return True

def branch_and_bound_cyclopeptide_sequencing(spectrum):
    """Branch and Bound method for cyclopeptide sequencing."""
    candidate_peptides = [[]]  # Start with an empty peptide
    final_peptides = []

    parent_mass = max(spectrum)

    while candidate_peptides:
        # Expand peptides by adding one more amino acid to each candidate
        candidate_peptides = [peptide + [aa] for peptide in candidate_peptides for aa in amino_acid_masses]

        # Copy list to avoid modification while iterating
        for peptide in candidate_peptides[:]:
            peptide_mass = sum(amino_acid_masses[aa] for aa in peptide)

            if peptide_mass == parent_mass:
                # If the mass matches, check the cyclic spectrum
                if cyclic_spectrum(peptide) == spectrum:
                    final_peptides.append(peptide)
                candidate_peptides.remove(peptide)
            elif peptide_mass > parent_mass or not is_consistent(peptide, spectrum):
                # Prune the peptide if its mass exceeds or it's inconsistent with the spectrum
                candidate_peptides.remove(peptide)

    return final_peptides

def convert_to_string(peptide):
    """Convert the peptide list back to a string of amino acids."""
    return ''.join(peptide)

# Example usage
spectrum = [0, 113, 128, 186, 241, 299, 314, 427]
result = branch_and_bound_cyclopeptide_sequencing(spectrum)

# Print the found cyclic peptides
print("Cyclic peptides found:")
for peptide in result:
    print(convert_to_string(peptide))

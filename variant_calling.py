def variant_calling(reference, sample):
    variants = []
    ref_len = len(reference)
    sample_len = len(sample)
    i = 0
    j = 0

    while i < ref_len and j < sample_len:
        if reference[i] == sample[j]:
            # Bases match, move to the next position
            i += 1
            j += 1
        else:
            # A variant has been detected
            if reference[i] != sample[j]:
                # SNP (Single Nucleotide Polymorphism)
                variants.append({
                    "type": "SNP",
                    "position": i + 1,
                    "ref_base": reference[i],
                    "sample_base": sample[j]
                })
                i += 1
                j += 1
            elif i < ref_len and j < sample_len and reference[i] != sample[j]:
                # Insertions/Deletions (indels)
                if i + 1 < ref_len and reference[i + 1] == sample[j]:
                    # Deletion in the sample
                    variants.append({
                        "type": "Deletion",
                        "position": i + 1,
                        "deleted_base": reference[i]
                    })
                    i += 1
                elif j + 1 < sample_len and reference[i] == sample[j + 1]:
                    # Insertion in the sample
                    variants.append({
                        "type": "Insertion",
                        "position": i + 1,
                        "inserted_base": sample[j]
                    })
                    j += 1

    # Check for remaining insertions or deletions
    if i < ref_len:
        for k in range(i, ref_len):
            variants.append({
                "type": "Deletion",
                "position": k + 1,
                "deleted_base": reference[k]
            })

    if j < sample_len:
        for k in range(j, sample_len):
            variants.append({
                "type": "Insertion",
                "position": ref_len + k + 1,
                "inserted_base": sample[k]
            })

    return variants


# Example usage:
reference_genome = "ACGTACGTACGT"
sample_genome = "ACGTTGGTACG"

# Call variants
variants = variant_calling(reference_genome, sample_genome)

# Print identified variants
print("Identified variants:")
for variant in variants:
    if variant["type"] == "SNP":
        print(f"SNP at position {variant['position']}: {variant['ref_base']} -> {variant['sample_base']}")
    elif variant["type"] == "Insertion":
        print(f"Insertion at position {variant['position']}: {variant['inserted_base']} inserted")
    elif variant["type"] == "Deletion":
        print(f"Deletion at position {variant['position']}: {variant['deleted_base']} deleted")

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def blast_sequence(seq):
    result_handle = NCBIWWW.qblast("blastn", "nt", seq)
    blast_record = NCBIXML.read(result_handle)

    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            print("*Alignment*")
            print("Sequence:", alignment.title)
            print("Length:", alignment.length)
            print("Score:", hsp.score)
            print("E-value:", hsp.expect)

query_seq = "ATGCGTACGTAGCTAGCTGACTGATCGTAGCTAGCTGACGTAGCTAGCATCGTACG"
blast_sequence(query_seq)

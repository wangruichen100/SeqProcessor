import typer

from Bio.Seq import Seq


def fasta_read(file_path):
    sequences = {}
    with open(file_path, 'r') as file:
        sequence_name = None
        sequence = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence_name is not None:
                    sequences[sequence_name] = sequence
                sequence_name = line[1:]
                sequence = ''
            else:
                sequence += line
        # Adding the last sequence
        if sequence_name is not None:
            sequences[sequence_name] = sequence

    return sequences

def fasta_read2(file_path):
    names = []
    sequences = []
    current_name = None
    current_sequence = ""
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_name:
                    names.append(current_name)
                    sequences.append(current_sequence)
                current_name = line[1:]
                current_sequence = ""
            else:
                current_sequence += line
        if current_name:
            names.append(current_name)
            sequences.append(current_sequence)
    return names, sequences


def translate_dna_to_protein(dna_sequence, codon_table):
    # Create a Seq object from the DNA sequence
    dna_seq = Seq(dna_sequence)
    
    # Translate the DNA sequence to protein sequence using the specified codon table
    protein_seq = dna_seq.translate(table=codon_table)
    
    return protein_seq

def translate_sequence(nucleotide_sequence):
    codon_table = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
    }
    
    nucleotide_sequence = nucleotide_sequence.upper()
    protein_sequence = ""
    for i in range(0, len(nucleotide_sequence), 3):
        codon = nucleotide_sequence[i:i+3]
        if codon in codon_table:
            protein_sequence += codon_table[codon]
        else:
            protein_sequence += "?"
    return protein_sequence
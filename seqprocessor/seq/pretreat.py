import typer

from Bio import SeqIO
from enum import Enum

from seqprocessor.utils import fasta_read, fasta_read2



pretreat_app = typer.Typer(help="Sequence pretreatment")

# Modify sequences with duplicate IDs.
def handle_duplicate_ids(names, sequences, delete_format):
    """
    Handle duplicate IDs in sequences.

    Args:
        names (list): List of sequence IDs.
        sequences (list): List of sequences.
        delete_format (DeleteDuplicateFormat): Format for handling duplicates.

    Returns:
        tuple: Tuple containing new list of sequence IDs and sequences.
    """
    new_names = []
    new_sequences = []
    id_counts = {}
    for name, sequence in zip(names, sequences):
        if name in id_counts:
            if delete_format == DeleteDuplicateFormat.f1:
                # If duplicate IDs should be deleted, add a suffix to the name
                id_counts[name] += 1
                new_name = f"{name}_{id_counts[name]}"
                new_names.append(new_name)
                new_sequences.append(sequence)
        else:
            id_counts[name] = 0
            new_names.append(name)
            new_sequences.append(sequence)
    return new_names, new_sequences

class DeleteDuplicateFormat(str, Enum):
    f1 = "delete"  # "Delete duplicate sequences."
    f2 = "rename"  # "Rename duplicate sequence names."

@pretreat_app.command(name="id_duplicate", 
                      help="Handle sequences with duplicate IDs.")
def id_duplicate(
    file_path: str = typer.Option(..., "--input", "-i", help="A FASTA file"),
    out_path: str = typer.Option(..., "--out", "-o", help="Output file path"),
    delete_format: DeleteDuplicateFormat = typer.Option(DeleteDuplicateFormat.f1, "--format", "-f",  help="Delete format"),
):
    """
    Handle sequences with duplicate IDs.

    Args:
        file_path (str): Path to the input FASTA file.
        out_path (str): Path to the output file.
        delete_format (DeleteDuplicateFormat): Format for handling duplicates.
    """
    names, sequences = fasta_read2(file_path)
    names, sequences = handle_duplicate_ids(names, sequences, delete_format)
    
    with open(out_path, "w") as outfile:
        for name, sequence in zip(names, sequences):
            outfile.write(f">{name}\n{sequence}\n")


# Handling characters outside of ATCG in the sequences.
class ProcessFormat(str, Enum):
    f1 = "f1"  # "Delete characters outside of ATCG."
    f2 = "f2"  # "Replace characters outside of ATCG with N."
    f3 = "f3"  # "Replace characters outside of ATCG with -."

@pretreat_app.command(name="seq_standardize", 
                      help="Process characters outside of ATCG in the FASTA file sequences")
def seq_standardize(
    file_path: str = typer.Option(..., "--input", "-i", help="A FASTA file"),
    out_path: str = typer.Option(..., "--out", "-o", help="Output file path"),
    standardize_format: ProcessFormat = typer.Option(ProcessFormat.f1, "--format", "-f",  help="Standardize format"),
):
    """
    Process characters outside of ATCG in the FASTA file sequences.

    Args:
        file_path (str): Path to the input FASTA file.
        out_path (str): Path to the output file.
        standardize_format (ProcessFormat): Format for standardizing sequences.

    Returns:
        None
    """
    # Read sequences from the input FASTA file
    sequences = fasta_read(file_path)
    records = []

    # Process each sequence
    for name, seq in sequences.items():
        # Convert sequence to uppercase
        seq = seq.upper()
        if standardize_format == ProcessFormat.f1:
            # Delete characters outside of ATCG
            seq = ''.join(c for c in seq if c in 'ATCGatcg')
        elif standardize_format == ProcessFormat.f2:
            # Replace characters outside of ATCG with N
            seq = ''.join(c if c in 'ATCGatcg' else 'N' for c in seq)
        elif standardize_format == ProcessFormat.f3:
            # Replace characters outside of ATCG with -
            seq = ''.join(c if c in 'ATCGatcg' else '-' for c in seq)

        # Append the processed sequence to the records list
        records.append({"Name": name, "Sequence": seq})

    # Write the processed sequences to the output file
    with open(out_path, 'w') as outfile:
        for record in records:
            outfile.write(f'>{record["Name"]}\n{record["Sequence"]}\n')


# Changing the case of sequences.
class CaseFormat(str, Enum):
    f1 = "upper" #"Convert to uppercase."
    f2 = "lower" #"Convert to lowercase."

@pretreat_app.command(name="seq_case", 
                  help="Convert to uppercase or lowercase.")
def seq_case(
    file_path: str = typer.Option(..., "--input", "-i", help="A FASTA file"),
    out_path: str = typer.Option(..., "--out", "-o", help="Output file path"),
    case_format: CaseFormat = typer.Option(CaseFormat.f1, "--format", "-f",  help="Uppercase or lowercase"),
):
    sequences = fasta_read(file_path)
    records = []
    for name, seq in sequences.items():
        if case_format == CaseFormat.f1:
            # Convert to uppercase
            seq = seq.upper()
        elif case_format == CaseFormat.f2:
            # Convert to lowercase
            seq = seq.lower()

        records.append({"Name": name, "Sequence": seq})

    with open(out_path, 'w') as outfile:
        for record in records:
            outfile.write(f'>{record["Name"]}\n{record["Sequence"]}\n')
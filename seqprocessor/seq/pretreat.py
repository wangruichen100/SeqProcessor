import typer
import pandas as pd

from Bio import SeqIO
from enum import Enum
from collections import Counter

from seqprocessor.utils import fasta_read, fasta_read2
from seqprocessor.options import OutputFormat


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
    f1 = "upper"  # "Convert to uppercase."
    f2 = "lower"  # "Convert to lowercase."

@pretreat_app.command(name="seq_case", 
                      help="Convert sequences to uppercase or lowercase.")
def seq_case(
    file_path: str = typer.Option(..., "--input", "-i", help="Path to the input FASTA file."),
    out_path: str = typer.Option(..., "--out", "-o", help="Path to the output file."),
    case_format: CaseFormat = typer.Option(CaseFormat.f1, "--format", "-f",  help="Uppercase or lowercase"),
):
    """
    Convert sequences to uppercase or lowercase.

    Args:
        file_path (str): Path to the input FASTA file.
        out_path (str): Path to the output file.
        case_format (CaseFormat): Format for case conversion.

    Returns:
        None
    """
    # Read sequences from the input FASTA file
    sequences = fasta_read(file_path)
    records = []

    # Process each sequence
    for name, seq in sequences.items():
        if case_format == CaseFormat.f1:
            # Convert to uppercase
            seq = seq.upper()
        elif case_format == CaseFormat.f2:
            # Convert to lowercase
            seq = seq.lower()

        # Append the processed sequence to the records list
        records.append({"Name": name, "Sequence": seq})

    # Write the processed sequences to the output file
    with open(out_path, 'w') as outfile:
        for record in records:
            outfile.write(f'>{record["Name"]}\n{record["Sequence"]}\n')

# Changing the name of sequences.
@pretreat_app.command(name="name_change", 
                      help="Replace sequence names according to the corresponding table.")
def name_change(
    file_path: str = typer.Option(..., "--input", "-i", help="Path to the input FASTA file."),
    info_path: str = typer.Option(..., "--info", help="Path to the input Info file."),
    out_path: str = typer.Option(..., "--out", "-o", help="Path to the output file."),
    info_format: OutputFormat = typer.Option(OutputFormat.csv, "--format", "-f",  help="Info file format"),
):
    """
    Replace sequence names according to the corresponding table.

    Args:
        file_path (str): Path to the input FASTA file.
        info_path (str): Path to input Info file.
        out_path (str): Path to the output file.
        info_format (info_format): Format for Info file format.

    Returns:
        None
    """
    names, sequences = fasta_read2(file_path)
    
    if info_format == OutputFormat.csv:
        df_info = pd.read_csv(info_path)
    elif info_format == OutputFormat.table:
        df_info = pd.read_csv(info_path, sep="\t")
    elif info_format == OutputFormat.excel:
        df_info = pd.read_excel(info_path)
    
    old_name = df_info.iloc[:, 0].to_list()
    new_name = df_info.iloc[:, 1].to_list()
    name_dict = dict(zip(old_name,new_name))

    with open(out_path, 'w') as outfile:
        for name, sequence in zip(names, sequences):
            name_ = name_dict[name]
            outfile.write(f'>{name_}\n{sequence}\n')

# Filtering low-quality sequences.
@pretreat_app.command(name="quality_control", 
                      help="Filtering low-quality sequences.")
def quality_control(
    file_path: str = typer.Option(..., "--input", "-i", help="Path to the input FASTA file."),
    out_path: str = typer.Option(..., "--out", "-o", help="Path to the output file."),
    qc_percentage: float = typer.Option(..., "--qc", "-q", help="Qc percentage, 0.0 to 1.0."),
):

    names, sequences = fasta_read2(file_path)

    with open(out_path, 'w') as outfile:
        for name, sequence in zip(names, sequences):
            sequence = sequence.upper()
            seq_length = len(sequence)
            seq_count = Counter(sequence)
            atcg_count = seq_count["A"]+seq_count["T"]+seq_count["C"]+seq_count["G"]
            # Checking if percentage of A, T, C, and G bases is greater than or equal to qc_percentage
            if atcg_count/seq_length >= qc_percentage:
                outfile.write(f'>{name}\n{sequence}\n')

if __name__ == "__main__":
    pretreat_app()
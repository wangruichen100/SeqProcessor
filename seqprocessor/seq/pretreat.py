import typer

from enum import Enum

from seqprocessor.utils import fasta_read



pretreat_app = typer.Typer(help="Sequence pretreatment")

class ProcessFormat(str, Enum):
    f1 = "f1" #"Delete characters outside of ATCG."
    f2 = "f2" #"Replace characters outside of ATCG with N."
    f3 = "f3" #"Replace characters outside of ATCG with -."

@pretreat_app.command(name="fasta_standardize", 
                  help="Process characters outside of ATCG in the FASTA file sequences")
def fasta_standardize(
    file_path: str = typer.Option(..., "--input", "-i", help="A FASTA file"),
    out_path: str = typer.Option(..., "--out", "-o", help="Output file path"),
    standardize_format: ProcessFormat = typer.Option(ProcessFormat.f1, "--format", "-f",  help="Standardize format"),
):
    sequences = fasta_read(file_path)
    records = []

    for name, seq in sequences.items():
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

        records.append({"Name": name, "Sequence": seq})

    with open(out_path, 'w') as outfile:
        for record in records:
            outfile.write(f'>{record["Name"]}\n{record["Sequence"]}\n')


class CaseFormat(str, Enum):
    f1 = "upper" #"Convert to uppercase."
    f2 = "lower" #"Convert to lowercase."

@pretreat_app.command(name="fasta_case", 
                  help="Convert to uppercase or lowercase.")
def fasta_standardize(
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
import typer
import pandas as pd

from Bio import SeqIO

from seqprocessor.utils import fasta_read
from seqprocessor.options import OutputFormat


info_app = typer.Typer(help="Get information")

def process_genebank_record(record):
    feature = record.features[0] if record.features else None
    qualifiers = {
        "id": record.name,
        "length": len(record.seq),
        "organism": feature.qualifiers.get("organism", ["no"])[0],
        "strain": feature.qualifiers.get("strain", ["no"])[0],
        "isolate": feature.qualifiers.get("isolate", ["no"])[0],
        "country": feature.qualifiers.get("country", ["no"])[0],
        "host": feature.qualifiers.get("host", ["no"])[0],
        "date": feature.qualifiers.get("collection_date", ["no"])[0],
        "note": feature.qualifiers.get("note", ["no"])[0],
        "taxonomy": record.annotations.get("taxonomy", "no"),
    }
    return qualifiers

@info_app.command(name="genbank_profile", help="Get information from GenBank")
def genbank_profile(
    file_path: str = typer.Option(..., "--input", "-i", help="A GenBank file"),
    out_path: str = typer.Option(..., "--out", "-o", help="Output file path"),
    out_format: OutputFormat = typer.Option(OutputFormat.csv, "--format", "-f", help="Output file format"),
):
    try:
        records = [process_genebank_record(record) for record in SeqIO.parse(file_path, "genbank")]
    except Exception as e:
        typer.echo(f"Error processing GenBank file: {e}")
        raise typer.Abort()

    df = pd.DataFrame(records)

    try:
        if out_format == OutputFormat.csv:
            df.to_csv(out_path, index=False)
        elif out_format == OutputFormat.table:
            df.to_csv(out_path, sep="\t", index=False)
        elif out_format == OutputFormat.excel:
            df.to_excel(out_path, index=False)
    except Exception as e:
        typer.echo(f"Error writing output file: {e}")
        raise typer.Abort()

@info_app.command(name="fasta_profile", help="Get name, sequence, and length from FASTA")
def fasta_profile(
    file_path: str = typer.Option(..., "--input", "-i", help="A FASTA file"),
    out_path: str = typer.Option(..., "--out", "-o", help="Output file path"),
    out_format: OutputFormat = typer.Option(OutputFormat.csv, "--format", "-f", help="Output file format"),
):
    try:
        sequences = fasta_read(file_path)
        records = [{"Name": name, "Sequence": seq, "Length": len(seq)} for name, seq in sequences.items()]
    except Exception as e:
        typer.echo(f"Error processing FASTA file: {e}")
        raise typer.Abort()

    df = pd.DataFrame(records)

    try:
        if out_format == OutputFormat.csv:
            df.to_csv(out_path, index=False)
        elif out_format == OutputFormat.table:
            df.to_csv(out_path, sep="\t", index=False)
        elif out_format == OutputFormat.excel:
            df.to_excel(out_path, index=False)
    except Exception as e:
        typer.echo(f"Error writing output file: {e}")
        raise typer.Abort()
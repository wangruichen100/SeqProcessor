import typer

from seqprocessor.utils import fasta_read2, translate_sequence

translate_app = typer.Typer(help="Translate FASTA file to protein sequences")

@translate_app.command(name="nt_aa", help="Translate FASTA file to protein sequences.")
def nt_aa(
    file_path: str = typer.Option(..., "--input", "-i", help="A FASTA file"),
    out_path: str = typer.Option(..., "--out", "-o", help="Output file path"),
):
    # Open the output file for writing
    with open(out_path, "w") as output_file:
        # Iterate over each sequence in the input FASTA file
        names, sequences = fasta_read2(file_path)
        for name, sequence in zip(names,sequences):
            # Translate the DNA sequence to protein sequence
            protein_seq = translate_sequence(str(sequence))
            # Write the protein sequence to the output file
            output_file.write(f">{name}\n{protein_seq}\n")

    
    

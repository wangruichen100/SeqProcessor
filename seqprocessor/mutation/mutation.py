import typer

from Bio import SeqIO
from collections import Counter
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import pandas as pd
import os
from tqdm import tqdm

from seqprocessor.utils import fasta_read, fasta_read2
from seqprocessor.options import OutputFormat

mutation_app = typer.Typer(help="Mutation analysis")

@mutation_app.command(name="pair_mutation", help="Compare pair mutations by ref sequence.")
def pair_mutation(
    file_path: str = typer.Option(..., "--input", "-i", help="A FASTA file"),
    out_format: OutputFormat = typer.Option(OutputFormat.csv, "--format", "-f", help="Output file format"),
):
    # Read sequences from the FASTA file
    names, sequences = fasta_read2(file_path)
    ref_sequence = sequences[0]  # Reference sequence is the first sequence in the file

    all_differences = []
    # Compare each sequence to the reference sequence
    for name, sequence in zip(names[1:], sequences[1:]):
        differences = find_differences(ref_sequence, sequence, name)
        all_differences.extend(differences)
    
    # Create a DataFrame from the differences
    df = pd.DataFrame(all_differences, columns=["Site", "Ref seq", "Aligned seq", "Aligned seq name"])
    # Pivot the DataFrame to wide format
    wide_df = df.pivot(index='Site', columns='Aligned seq name', values='Aligned seq')
    
    # Get unique sites from the DataFrame and sort them
    ref_sites = sorted(set(df["Site"]))
    # Get reference letters for each site
    ref_letters = [ref_sequence[i - 1] for i in ref_sites]
    # Insert reference letters as the first column in the wide DataFrame
    wide_df.insert(0, 'Ref seq', ref_letters)

    # Fill NaN values in wide DataFrame with reference letters
    for index, row in wide_df.iterrows():
        if pd.isna(row['Ref seq']):  # Skip if the value in column 'Ref seq' is NaN
            continue
        for col in wide_df.columns[1:]:  # Iterate over other columns except 'Ref seq'
            if pd.isna(row[col]):  # Replace NaN values in the current column with the corresponding value from column 'Ref seq'
                wide_df.at[index, col] = row['Ref seq']
    
    # Save the results to the specified output format
    out_path = file_path.replace(".*","")
    if out_format == OutputFormat.csv:
        df.to_csv(f"{out_path}_l.csv", index=False)
        wide_df.to_csv(f"{out_path}_w.csv")
    elif out_format == OutputFormat.table:
        df.to_csv(f"{out_path}_l.tab", sep="\t", index=False)
        wide_df.to_csv(f"{out_path}_w.tab", sep="\t")
    elif out_format == OutputFormat.excel:
        df.to_excel(f"{out_path}_l.xlsx", index=False)
        wide_df.to_excel(f"{out_path}_w.xlsx")

def find_differences(ref_seq, seq, name):
    differences = []
    for i, (letter1, letter2) in enumerate(zip(ref_seq, seq)):
        if letter1 != letter2:
            differences.append([i+1, letter1, letter2, name])
    return differences



@mutation_app.command(name="group_aamutation", help="Compare amino acid mutations in different classifications.")
def group_aamutation(
    file_path: str = typer.Option(..., "--input", "-i", help="A FASTA file"),
    info_path: str = typer.Option(..., "--info", "-info", help="Information file path"),
    out_path: str = typer.Option(..., "--out", "-o", help="Output file path"),
    id_column: str = typer.Option("id", "--id", "-id", help="ID column"),
    type_column: str = typer.Option("type", "--type", "-type", help="Classification column"),
    picture_type: str = typer.Option("bar", "--ptype", "-pt", help="Picture format"),
):
    seq_record = fasta_read(file_path)
    df = pd.read_excel(info_path)

    os.makedirs(f"{out_path}/{type_column}", exist_ok=True)

    for _, group_df in df.groupby(type_column):
        group_name = group_df.iloc[0][type_column]
        group_fasta_path = f"{out_path}/{type_column}/{group_name}.fas"

        with open(group_fasta_path, "w") as f:
            for seq_id in group_df[id_column]:
                f.write(f">{seq_id}\n{seq_record[seq_id]}\n")

        create_plots_aa(group_fasta_path, group_name, out_path, picture_type)

def create_plots_aa(input_fasta, group_name, out_path, picture_type):
    pdf_path = f"{out_path}/{group_name}.pdf"

    if os.path.exists(pdf_path):
        os.remove(pdf_path)

    seq_data = [list(seq.seq) for seq in SeqIO.parse(input_fasta, "fasta")]
    df = pd.DataFrame(seq_data)

    freq_df = []
    for col in df.columns:
        stat = Counter(df[col])
        freq = [stat.get(aa, 0) / len(df) for aa in 'GAVLIFWYDNKEQMSTCPHR']
        freq_df.append(freq)

    freq_df = pd.DataFrame(freq_df, columns=list('GAVLIFWYDNKEQMSTCPHR'))

    with PdfPages(pdf_path) as pdf:
        for i, row in tqdm(freq_df.iterrows(), desc="Plotting"):
            if picture_type == "bar":
                row.plot(kind='bar', stacked=True, colormap='tab20', yticks=[0.5, 1])
            elif picture_type == "line":
                row.plot(kind='line', colormap='tab20', yticks=[0.5, 1])

            plt.ylim((0, 1.1))
            plt.legend(prop={'size': 8}, bbox_to_anchor=(1.12, 1), loc='upper right', fontsize=5)
            plt.title(f'Site {i + 1}')
            plt.xlabel(None)
            pdf.savefig()
            plt.clf()

@mutation_app.command(name="group_ntmutation", help="Compare nucleotide mutations in different classifications.")
def group_ntmutation(
    file_path: str = typer.Option(..., "--input", "-i", help="A FASTA file"),
    info_path: str = typer.Option(..., "--info", "-info", help="Information file path"),
    out_path: str = typer.Option(..., "--out", "-o", help="Output file path"),
    id_column: str = typer.Option("id", "--id", "-id", help="ID column"),
    type_column: str = typer.Option("type", "--type", "-type", help="Classification column"),
    picture_type: str = typer.Option("bar", "--ptype", "-pt", help="Picture format"),
):
    seq_record = fasta_read(file_path)
    df = pd.read_excel(info_path)

    os.makedirs(f"{out_path}/{type_column}", exist_ok=True)

    for _, group_df in df.groupby(type_column):
        group_name = group_df.iloc[0][type_column]
        group_fasta_path = f"{out_path}/{type_column}/{group_name}.fas"

        with open(group_fasta_path, "w") as f:
            for seq_id in group_df[id_column]:
                f.write(f">{seq_id}\n{seq_record[seq_id]}\n")

        create_plots_nt(group_fasta_path, group_name, out_path, picture_type)

def create_plots_nt(input_fasta, group_name, out_path, picture_type):
    pdf_path = f"{out_path}/{group_name}.pdf"

    if os.path.exists(pdf_path):
        os.remove(pdf_path)

    seq_data = [list(seq.seq) for seq in SeqIO.parse(input_fasta, "fasta")]
    df = pd.DataFrame(seq_data)

    freq_df = []
    for col in df.columns:
        stat = Counter(df[col])
        freq = [stat.get(aa, 0) / len(df) for aa in 'ATCG']
        freq_df.append(freq)

    freq_df = pd.DataFrame(freq_df, columns=list('ATCG'))

    with PdfPages(pdf_path) as pdf:
        for i, row in tqdm(freq_df.iterrows(), desc="Plotting"):
            if picture_type == "bar":
                row.plot(kind='bar', stacked=True, colormap='Set2', yticks=[0.5, 1])
            elif picture_type == "line":
                row.plot(kind='line', colormap='Set2', yticks=[0.5, 1])

            plt.ylim((0, 1.1))
            plt.legend(prop={'size': 8}, bbox_to_anchor=(1.12, 1), loc='upper right', fontsize=5)
            plt.title(f'Site {i + 1}')
            plt.xlabel(None)
            pdf.savefig()
            plt.clf()

if __name__ == "__main__":
    mutation_app()

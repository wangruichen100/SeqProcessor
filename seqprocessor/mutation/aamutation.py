import typer

from Bio import SeqIO
from collections import Counter
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import pandas as pd
import os
from tqdm import tqdm

from seqprocessor.utils import fasta_read

aamutation_app = typer.Typer(help="Amino acid mutation analysis")

@aamutation_app.command(name="aa_mutation", help="Compare amino acid mutations in different classifications.")
def aa_mutation(
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

        create_plots(group_fasta_path, group_name, out_path, picture_type)

def create_plots(input_fasta, group_name, out_path, picture_type):
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

if __name__ == "__main__":
    aamutation_app()

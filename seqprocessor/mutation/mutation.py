import typer

from Bio import SeqIO
from collections import Counter
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

from seqprocessor.utils import fasta_read, fasta_read2, tab21_black
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
    seq_type = sorted(set(df[type_column].to_list()))
    os.makedirs(f"./{out_path}/{type_column}", exist_ok=True)

    for seq_type_ in seq_type:
        df_type = df[df[type_column]==seq_type_]
        # print(seq_type_)
        with open(f"./{out_path}/{type_column}/{seq_type_}.fas", "w") as f:
            for name_ in df_type[id_column].to_list():          
                f.write(f">{name_}\n{seq_record[name_]}\n")

    pdf = f"{out_path}/{type_column}/{type_column}.pdf"
    file_path = f"{out_path}/{type_column}"
    file_name = seq_type

    site_all = []
    aa_all = []
    stat_df_all = []
    freq_df_all = []
    group_all = []

    for name_ in file_name:
        print(name_)
        site, aa, stat_df, freq_df, group_ = aatj("{0}/{1}.fas".format(file_path, name_), name_)

        site_all += site
        aa_all += aa
        stat_df_all += stat_df
        freq_df_all += freq_df
        group_all += group_

    df = pd.DataFrame()
    df["site"] = site_all
    df["aa"] = aa_all
    df["stat"] = stat_df_all
    df["freq"] = freq_df_all
    df["group"] = group_all
    df = df.fillna(0)
    df.to_csv(f"{out_path}/{type_column}/{type_column}.csv",index = False)

    with PdfPages(pdf) as pdf:
        for i_ in range(df["site"].min(), df["site"].max() + 1):
            print(i_ + 1)
            df_site = df[df["site"] == i_]

            pivot_df = df_site.pivot(index='group', columns='aa', values='freq')
            if picture_type == "bar":
                pivot_df.plot(kind='bar', stacked=True, colormap=tab21_black(), yticks=[0.5, 1])
            elif picture_type == "line":
                pivot_df.plot(kind='line', colormap=tab21_black(), yticks=[0.5, 1])
            plt.ylim((0, 1.1))
            plt.legend(prop={'size': 8}, bbox_to_anchor=(1.12, 1),loc='upper right', fontsize=5)
            plt.title('Site' + str(i_+1))
            plt.xlabel(None)
            # plt.show()
            pdf.savefig()
            plt.clf()
            plt.close('all')

def aatj(input_fasta, group):
    site = []
    aa = []
    stat_df = []
    freq_df = []

    seq_data = [list(i.seq) for i in SeqIO.parse(input_fasta, "fasta")]
    df = pd.DataFrame(seq_data)

    for column in df.columns:
        aa_ = ['G', 'A', 'V', 'L', 'I', 'F', 'W', 'Y', 'D', 'N', 'E', 'K', 'Q', 'M', 'S', 'T', 'C', 'P', 'H', 'R', "-"]
        aa += aa_
        site += [column] * len(aa_)

        stat = Counter(df[column])
        aa_dict = {aa: stat.get(aa, 0) for aa in aa_}
        stat_df.extend(aa_dict.values())

        freq = np.array(list(aa_dict.values())) / sum(aa_dict.values())
        freq_df.extend(freq)

    group_ = [group] * len(site)

    return site, aa, stat_df, freq_df, group_


@mutation_app.command(name="group_ntmutation", help="Compare nucleotide in different classifications.")
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
    seq_type = sorted(set(df[type_column].to_list()))
    os.makedirs(f"./{out_path}/{type_column}", exist_ok=True)

    for seq_type_ in seq_type:
        df_type = df[df[type_column]==seq_type_]
        # print(seq_type_)
        with open(f"./{out_path}/{type_column}/{seq_type_}.fas", "w") as f:
            for name_ in df_type[id_column].to_list():          
                f.write(f">{name_}\n{seq_record[name_].upper()}\n")

    pdf = f"{out_path}/{type_column}/{type_column}.pdf"
    file_path = f"{out_path}/{type_column}"
    file_name = seq_type

    site_all = []
    nt_all = []
    stat_df_all = []
    freq_df_all = []
    group_all = []

    for name_ in file_name:
        print(name_)
        site, nt, stat_df, freq_df, group_ = nttj("{0}/{1}.fas".format(file_path, name_), name_)

        site_all += site
        nt_all += nt
        stat_df_all += stat_df
        freq_df_all += freq_df
        group_all += group_

    df = pd.DataFrame()
    df["site"] = site_all
    df["nt"] = nt_all
    df["stat"] = stat_df_all
    df["freq"] = freq_df_all
    df["group"] = group_all
    df = df.fillna(0)
    df.to_csv(f"{out_path}/{type_column}/{type_column}.csv",index = False)

    with PdfPages(pdf) as pdf:
        for i_ in range(df["site"].min(), df["site"].max() + 1):
            print(i_ + 1)
            df_site = df[df["site"] == i_]

            pivot_df = df_site.pivot(index='group', columns='nt', values='freq')
            if picture_type == "bar":
                pivot_df.plot(kind='bar', stacked=True, colormap="Set1", yticks=[0.5, 1])
            elif picture_type == "line":
                pivot_df.plot(kind='line', colormap="Set1", yticks=[0.5, 1])
            plt.ylim((0, 1.1))
            plt.legend(prop={'size': 8}, bbox_to_anchor=(1.12, 1),loc='upper right', fontsize=5)
            plt.title('Site' + str(i_+1))
            plt.xlabel(None)
            # plt.show()
            pdf.savefig()
            plt.clf()
            plt.close('all')

def nttj(input_fasta, group):
    site = []
    nt = []
    stat_df = []
    freq_df = []

    seq_data = [list(i.seq) for i in SeqIO.parse(input_fasta, "fasta")]
    df = pd.DataFrame(seq_data)

    for column in df.columns:
        nt_ = ['A',"T","C","G","-"]
        nt += nt_
        site += [column] * len(nt_)

        stat = Counter(df[column])
        nt_dict = {nt: stat.get(nt, 0) for nt in nt_}
        stat_df.extend(nt_dict.values())

        freq = np.array(list(nt_dict.values())) / sum(nt_dict.values())
        freq_df.extend(freq)

    group_ = [group] * len(site)

    return site, nt, stat_df, freq_df, group_


if __name__ == "__main__":
    mutation_app()

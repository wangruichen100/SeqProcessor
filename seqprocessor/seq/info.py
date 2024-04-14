import typer
import pandas as pd

from Bio import SeqIO

from seqprocessor.utils import fasta_read
from seqprocessor.options import OutputFormat


info_app = typer.Typer(help="Get information")

@info_app.command(name="fasta_profile", help="Get name, seq and length")
def fasta_profile(
    file_path: str = typer.Option(..., "--input", "-i", help="A fasta file"),
    out_path: str = typer.Option(..., "--out", "-o", help="Out file"),
    out_format: OutputFormat = typer.Option(OutputFormat.csv, "--format", "-f", help="Out file format"),
):
    sequences = fasta_read(file_path)
    seq_name = []
    seq_seq = []
    seq_length = []
    for num, (name, sequence) in enumerate(sequences.items()):
        print(num, name)
        length = len(sequence)
        seq_name.append(name)
        seq_seq.append(sequence)
        seq_length.append(length)

    df = pd.DataFrame()
    df["Name"] = seq_name
    df["Sequence"] = seq_seq
    df["Length"] = seq_length

    if out_format == OutputFormat.csv:
        df.to_csv(out_path)

    elif out_format == OutputFormat.table:
        df.to_csv(out_path, sep="\t")

    elif out_format == OutputFormat.excel:
        df.to_excel(out_path,index=False)
    

@info_app.command(name="genebank_profile", help="Get information from genebank")
def genebank_profile(
    file_path: str = typer.Option(..., "--input", "-i", help="A genbank file"),
    out_path: str = typer.Option(..., "--out", "-o", help="Out file"),
    out_format: OutputFormat = typer.Option(OutputFormat.csv, "--format", "-f", help="Out file format"),
):
    sequences = SeqIO.parse(file_path,"genbank")
    name = []
    lenght = []
    organism = []
    strain = []
    isolate = []
    country = []
    host = []
    coll_date = []
    note = []
    taxonomy = []

    for num,i in enumerate(sequences):
        print(num+1, i.name)
        name.append(i.name)
        lenght.append(len(i.seq))

        if "strain" in i.features[0].qualifiers.keys():
            strain.append(i.features[0].qualifiers["strain"][0])
        if not ("strain" in i.features[0].qualifiers.keys()):
            strain.append("no")

        if "organism" in i.features[0].qualifiers.keys():
            organism.append(i.features[0].qualifiers["organism"][0])
        if not ("organism" in i.features[0].qualifiers.keys()):
            organism.append("no")

        if "isolate" in i.features[0].qualifiers.keys():
            isolate.append(i.features[0].qualifiers["isolate"][0])
        if not ("isolate" in i.features[0].qualifiers.keys()):
            isolate.append("no")

        if "country" in i.features[0].qualifiers.keys():
            country.append(i.features[0].qualifiers["country"][0])
        if not ("country" in i.features[0].qualifiers.keys()):
            country.append("no")

        if "host" in i.features[0].qualifiers.keys():
            host.append(i.features[0].qualifiers["host"][0])
        if not ("host" in i.features[0].qualifiers.keys()):
            host.append("no")

        if "collection_date" in i.features[0].qualifiers.keys():
            coll_date.append(i.features[0].qualifiers["collection_date"][0])
        if not ("collection_date" in i.features[0].qualifiers.keys()):
            coll_date.append("no")
            
        if "note" in i.features[0].qualifiers.keys():
            note.append(i.features[0].qualifiers["note"][0])
        if not ("note" in i.features[0].qualifiers.keys()):
            note.append("no")
            
        if "taxonomy" in i.annotations.keys():
            taxonomy.append(i.annotations["taxonomy"])
        if not ("taxonomy" in i.annotations.keys()):
            taxonomy.append("no")
            
    df = pd.DataFrame()
    df["id"] = name
    df["lenght"] = lenght
    df["organism"] = organism
    df["strain"] = strain
    df["isolate"] = isolate
    df["country"] = country
    df["host"] = host
    df["date"] = coll_date
    df["note"] = note
    df["taxonomy"] = taxonomy

    if out_format == OutputFormat.csv:
        df.to_csv(out_path)

    elif out_format == OutputFormat.table:
        df.to_csv(out_path, sep="\t")

    elif out_format == OutputFormat.excel:
        df.to_excel(out_path,index=False)

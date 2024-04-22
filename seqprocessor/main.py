import typer

from seqprocessor.seq.info import info_app
from seqprocessor.seq.pretreat import pretreat_app
from seqprocessor.seq.align import align_app
from seqprocessor.seq.phylo import phylo_app

app = typer.Typer()

app.add_typer(info_app, name="info_app", help="Get information")
app.add_typer(pretreat_app, name="pretreat_app", help="Sequence pretreatment")
app.add_typer(align_app, name="align_app", help="Multiple sequence alignment")
app.add_typer(phylo_app, name="phylo_app", help="Reconstructing Phylogenetic Tree")

if __name__ == "__main__":
    app()


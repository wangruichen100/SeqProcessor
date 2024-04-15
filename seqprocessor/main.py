import typer

from seqprocessor.seq.info import info_app
from seqprocessor.seq.pretreat import pretreat_app

app = typer.Typer()

app.add_typer(info_app, name="info_app", help="Get information")
app.add_typer(pretreat_app, name="pretreat_app", help="Sequence pretreatment")

if __name__ == "__main__":
    app()


import typer

from seqprocessor.seq.info import info_app

app = typer.Typer()

app.add_typer(info_app, name="info_app", help="get information")

if __name__ == "__main__":
    app()


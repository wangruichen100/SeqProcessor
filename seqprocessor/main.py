import typer

from seqprocessor.seq.info import info_app

app = typer.Typer()

app.add_typer(info_app, name="info_app", help="Get information")

if __name__ == "__main__":
    app()


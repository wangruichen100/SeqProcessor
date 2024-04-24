import subprocess
import typer

align_app = typer.Typer(help="Multiple sequence alignment")

# Command to perform multiple sequence alignment with MAFFT
# Check if MAFFT is installed
def is_mafft_installed():
    try:
        # Execute conda list command and capture the output
        output = subprocess.check_output(["conda", "list"], universal_newlines=True)
        # Check if "mafft" is in the output
        return "mafft" in output
    except subprocess.CalledProcessError:
        # If the command fails, Conda may not be installed or configured correctly
        return False

# Install MAFFT using Conda
def install_mafft_with_conda():
    try:
        # Try installing MAFFT 7.525 using the correct Conda installation command
        subprocess.run(["conda", "install", "bioconda::mafft", "-y"], check=True)
        return True
    except subprocess.CalledProcessError:
        return False

@align_app.command(name="align_mafft",
                   help="Multiple sequence alignment by MAFFT.")
def align_mafft(
    file_path: str = typer.Option(..., "--input", "-i", help="Path to the input FASTA file."),
    out_path: str = typer.Option(..., "--out", "-o", help="Path to the output file."),
    
):
    print("Check if MAFFT is installed")
    if not is_mafft_installed():
        print("MAFFT is not installed. Attempting to install with Conda...")
        if not install_mafft_with_conda():
            print("Failed to install MAFFT with Conda.")
            return
        print("MAFFT installed successfully with Conda.")

    # Execute the MAFFT command
    try:
        subprocess.run(["mafft", file_path], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print("Error running MAFFT:", e)
        return

    # Write the result to the output file
    try:
        with open(out_path, "w") as f:
            subprocess.run(["mafft", file_path], check=True, stdout=f)
        print("MAFFT run successful. Results saved in", out_path)
    except subprocess.CalledProcessError as e:
        print("Error writing output file:", e)

if __name__ == "__main__":
    align_app()

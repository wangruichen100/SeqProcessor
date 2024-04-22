import subprocess
import typer

phylo_app = typer.Typer(help="Reconstructing Phylogenetic Tree")

# Check if IQ-TREE is installed
def is_iqtree_installed():
    try:
        # Execute conda list command and capture the output
        output = subprocess.check_output(["conda", "list"], universal_newlines=True)
        # Check if "IQ-TREE" is in the output
        return "iqtree" in output
    except subprocess.CalledProcessError:
        # If the command fails, Conda may not be installed or configured correctly
        return False

# Install IQ-TREE using Conda
def install_iqtree_with_conda():
    try:
        # Try installing IQ-TREE using the correct Conda installation command
        subprocess.run(["conda", "install", "bioconda::iqtree", "-y"], check=True)
        return True
    except subprocess.CalledProcessError:
        return False

# Command to perform phylogenetic tree reconstruction with IQ-TREE
@phylo_app.command(name="tree_iqtree",
                   help="Build phylogenetic tree using IQ-TREE.")
def tree_iqtree(
    file_path: str = typer.Option(..., "--input", "-i", help="Path to the input FASTA file."),
    model: str = typer.Option("MFP", "--model", "-m", help="Model to use."),
    threads: str = typer.Option("AUTO", "--threads", "-T", help="Threads to use."),
):
    print("Check if IQ-TREE is installed")
    if not is_iqtree_installed():
        print("IQ-TREE is not installed. Attempting to install with Conda...")
        if not install_iqtree_with_conda():
            print("Failed to install IQ-TREE with Conda.")
            return
        print("IQ-TREE installed successfully with Conda.")

    # Execute IQ-TREE command to build the tree
    # "iqtree -s input.fas -m MFP -bb 1000 -T AUTO"
    try:
        subprocess.run(["iqtree", "-s", file_path, "-m", model,"-bb", "1000", "-T", threads], check=True)
        print("Tree construction completed.")
    except subprocess.CalledProcessError as e:
        print("Error building tree with IQ-TREE:", e)
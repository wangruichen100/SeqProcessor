import typer


def fasta_read(file_path):
    sequences = {}
    with open(file_path, 'r') as file:
        sequence_name = None
        sequence = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence_name is not None:
                    sequences[sequence_name] = sequence
                sequence_name = line[1:]
                sequence = ''
            else:
                sequence += line
        # Adding the last sequence
        if sequence_name is not None:
            sequences[sequence_name] = sequence

    return sequences

def fasta_read2(file_path):
    names = []
    sequences = []
    current_name = None
    current_sequence = ""
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_name:
                    names.append(current_name)
                    sequences.append(current_sequence)
                current_name = line[1:]
                current_sequence = ""
            else:
                current_sequence += line
        if current_name:
            names.append(current_name)
            sequences.append(current_sequence)
    return names, sequences


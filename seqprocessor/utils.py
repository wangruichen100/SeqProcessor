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


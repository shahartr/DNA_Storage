
def run():
    with open('final_20_length.txt', 'r') as f:
        primers = [line.strip() for line in f]
    #create fasta file from the 22k primers
    with open('final_20_length.fasta', 'w') as f:
        for i, primer in enumerate(primers):
            f.write(f">primer_{i}\n{primer}\n")


if __name__ == '__main__':
    run()
# Create a FASTA file from the "final 20-length.txt" primers



#load primer sequences from a file (each in line) and export fasta file
if __name__ == '__main__':
    primers=[]
    with open('final_20_length25.txt', 'r') as f:
        for line in f.readlines():
            primers.append(line.strip())

    with open('primers.fasta', 'w') as f:
        for i, primer in enumerate(primers):
            f.write('>primer'+str(i)+'\n')
            f.write(primer+'\n')

    print('Primers saved to primers.fasta')



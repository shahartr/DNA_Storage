from Bio.SeqUtils import MeltingTemp as mt
from Bio.Blast import NCBIWWW, NCBIXML
from primer3 import calcHairpin, calcHomodimer


def run():
    count = 0
    with open('../../../../../data/results.txt', 'r') as infile:
        results = infile.readlines()
        print(len(results))
    for result in results:
        line = result.strip().split('\t')
        identity = float(line[2])
        alignment_length = int(line[3])
        e_value = float(line[10])
        if identity > 85 and alignment_length > 15 and e_value < 0.001:
            count += 1

    print(f"Number of primers left: {count}")


def check_tm(primer):
    tm = mt.Tm_NN(primer)
    if tm >= 55 and tm <= 60:
        return True
    else:
        return False


def remoteBlast():
    with open('final_20_length000002.txt', 'r') as file:
        primers = file.readlines()
    print(len(primers))
    final_primers = []
    count = 0
    for primer in primers:
        result_handle = NCBIWWW.qblast("blastn", "nt", primer)
        with open("../../../../../data/results.xml", "w") as out_handle:
            out_handle.write(result_handle.read())
        result_handle.close()
        with open("../../../../../data/results.xml") as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        if hsp.expect < 0.001:#hsp is high-scoring segment pair
                            if check_tm(primer):
                                final_primers.append(primer)
                            else:
                                count += 1


if __name__ == '__main__':
    #run()
    remoteBlast()

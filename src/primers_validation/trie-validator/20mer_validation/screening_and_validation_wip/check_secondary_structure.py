from nupack import *
Min_GC = 0.45
Max_GC = 0.55
Min_TM = 55
Max_TM = 65



def secondary_structure_found(primer):
    # Check for hairpin and homodimer formation using NUPACK secondary structure
    # prediction
    results = mfe([primer], material='dna', temp=55.0)
    if results:
        energy = results[0].energy
        if energy < -3.0:
            return True
    return False


def melting_temperature(primer):
    # Calculate the melting temperature of the primer
    complex_types = [Complex([primer], 1)]
    complex_analyses = mfe(complex_types, material='dna', temp=37.0)

    for analysis in complex_analyses:
        return analysis.melting_temperature

def in_mt_range(primer):
    tm = melting_temperature(primer)
    return Min_TM <= tm <= Max_TM

def run():
    with open('../../../../../data/valid_primers_so_far.txt', 'r') as f:
        primers = [line.strip() for line in f]

    final_primers = []
    count_primers_disqualify_by_structure = 0
    count_primers_disqualify_by_tm = 0
    for primer in primers:
        if not secondary_structure_found(primer):
            if not in_mt_range(primer):
                final_primers.append(primer)
            else:
                count_primers_disqualify_by_tm += 1
        else:
            count_primers_disqualify_by_structure += 1

    print(f"Disqualified by secondary structure: {count_primers_disqualify_by_structure} primers")
    print(f"Disqualified by melting temperature: {count_primers_disqualify_by_tm} primers")
    print(f"Number of valid primers: {len(final_primers)}")

if __name__ == '__main__':
    run()
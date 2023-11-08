from Bio import SeqIO
import os

def read_fasta(filename):
    sequences = {}
    for record in SeqIO.parse(filename, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def smith_waterman(seq1, seq2, gap, match, mismatch):
    m, n = len(seq1), len(seq2)
    matrix = [[0] * (n + 1) for _ in range(m + 1)]

    max_score = 0
    max_i, max_j = 0, 0

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i - 1] == seq2[j - 1]:
                score = matrix[i - 1][j - 1] + match
            else:
                score = max(
                    matrix[i - 1][j] - gap,
                    matrix[i][j - 1] - gap,
                    matrix[i - 1][j - 1] + mismatch,
                    0
                )

            matrix[i][j] = score
            if score > max_score:
                max_score = score
                max_i, max_j = i, j

    doneseq1, doneseq2 = [], []
    i, j = max_i, max_j

    while i > 0 and j > 0 and matrix[i][j] > 0:
        if seq1[i - 1] == seq2[j - 1]:
            doneseq1.insert(0, seq1[i - 1])
            doneseq2.insert(0, seq2[j - 1])
            i -= 1
            j -= 1
        elif matrix[i][j] == matrix[i - 1][j] - gap:
            doneseq1.insert(0, seq1[i - 1])
            doneseq2.insert(0, '-')
            i -= 1
        elif matrix[i][j] == matrix[i][j - 1] - gap:
            doneseq1.insert(0, '-')
            doneseq2.insert(0, seq2[j - 1])
            j -= 1
        else:
            doneseq1.insert(0, seq1[i - 1])
            doneseq2.insert(0, seq2[j - 1])
            i -= 1
            j -= 1

    return ''.join(doneseq1), ''.join(doneseq2), max_score

def main():
    while True:
        fasta_filename = input("Podaj nazwę pliku FASTA (z rozszerzeniem fa lub fasta): ")

        if not os.path.exists(fasta_filename):
            print(f'Plik {fasta_filename} nie istnieje.')
            choice = input("Czy chcesz wprowadzić nazwę pliku ponownie? (0 -> tak, 1 -> nie): ")
            if choice == '1':
                break
        else:
            gap= get_numeric_input("Podaj wartość Gap: ")
            mismatch = get_numeric_input("Podaj wartość Mismatch: ")
            match = get_numeric_input("Podaj wartość Match: ")

            sequences = read_fasta(fasta_filename)

            if len(sequences) != 2:
                print('ERROR -> plik FASTA powinien zawierać dwie sekwencje.')
            else:
                seq1 = sequences[list(sequences.keys())[0]]
                seq2 = sequences[list(sequences.keys())[1]]

                aligned_seq1, aligned_seq2, max_score = smith_waterman(seq1, seq2, gap, match, mismatch)

                with open('alignment_result.txt', 'w') as output_file:
                    output_file.write(f'Sekwencja 1: {aligned_seq1}\n')
                    output_file.write(f'Sekwencja 2: {aligned_seq2}\n\n')
                    output_file.write(f'Max Score: {max_score}\n')

                print('Wynik dopasowania zapisany w alignment_result.txt')
                break
def get_numeric_input(a):
    while True:
        user_input = input(a)
        try:
            numeric_value = float(user_input)
            return numeric_value
        except ValueError:
            print("Wprowadź prawidłową wartość liczbową (ułamki po kropce).")

if __name__ == '__main__':
    main()

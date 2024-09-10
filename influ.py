from Bio import SeqIO
import matplotlib.pyplot as plt
import os

# Wczytanie sekwencji z pliku FASTA
file_path = "./data/InfluA.fasta"  # ścieżka do pliku
fasta = list(SeqIO.parse(file_path, format= "fasta"))

# Tworzenie słownika
seqs = {}
for entry in fasta:
    seqs[entry.id] = entry.seq
#9700..11110 pos bialko NA od ref NC_026433

#wyznaczenie pozycji startu i konca genu ref po przyrownaniu
#loop seq - gap
def gapped_pos(seq, pos):
    non_gap = 0
    gaps = 0
    for nt in seq:
        if nt != '-':
            non_gap += 1
        else:
            gaps += 1
        if non_gap == pos:
            return pos + gaps


start = gapped_pos(seqs['Influenza-A-H1N1_31'], 9700)
#11595
end = gapped_pos(seqs['Influenza-A-H1N1_31'], 11110)
#13839 15431
'''
##### Dla aminokwasowej sekwencji
# Wyznaczenie fragmentów genu HA (bez przerw "-")[21584-1:25414]
'''
NA = {}
for seqrecord in fasta:
    NA_seq = seqrecord.seq[start-1:end] # Wybieramy odpowiedni fragment genu NA
    cleaned_seq = NA_seq.replace('-', '')  # Usuwanie przerw "-"
    NA[seqrecord.id] = cleaned_seq


# Przetłumaczenie sekwencji nukleotydowych na aminokwasy i zapisanie do pliku FASTA
with open('NA_translated.fasta', 'w') as f:  # Nowy plik wynikowy
    for NA_id, NA_seq in NA.items():
        translated_seq = NA_seq.translate()  # Tłumaczenie na aminokwasy
        f.write('>' + NA_id + '\n')  # Zapisanie nazwy sekwencji
        f.write(str(translated_seq) + '\n')  # Zapisanie przetłumaczonej sekwencji aminokwasowej

print("Plik 'NA_translated.fasta' został utworzony.")


# Wczytanie sekwencji z pliku FASTA
file_path = "./data/influ-NA-ali.fasta"  # ścieżka do pliku
NA_aa = list(SeqIO.parse(file_path, format="fasta"))

# Tworzenie słownika sekwencji bezpośrednio
seqs = {}
for entry in NA_aa:
    seqs[entry.id] = entry.seq
def get_aa_mutations(initial, variant):
    out = []
    seqs = list(zip(initial, variant))
    for pos, aa in enumerate(seqs):
        if aa[0] != aa[1]:
            out.append(aa[0].upper() + str(pos) + aa[1].upper())
    return out


for item in seqs:
    print(item + ' ' + str(len(get_aa_mutations(seqs['Influenza-A-H1N1_31'], seqs[item]))))

# Ścieżka do folderu, gdzie mają zostać zapisane pliki
output_folder = "./result/plots"


# Funkcja do wyodrębnienia sekwencji na podstawie nazwy wariantu
def get_variant_sequences(seqs, variant_keyword):
    return [seq_id for seq_id in seqs if variant_keyword in seq_id]


# Wyodrębnienie sekwencji dla każdego wariantu
a_seqs = get_variant_sequences(seqs, "Influenza-A")
b_seqs = get_variant_sequences(seqs, "Influenza-B")


# Funkcja rysująca wykresy dla danej grupy sekwencji i zapisująca do plików
def plot_and_save_mutation_chart(group, group_name):
    plt.figure(figsize=(11.7, 8.3))  # Ustawienie rozmiaru wykresu w formacie A4

    for y, item in enumerate(group):
        plt.plot((0, len(seqs['Influenza-A-H1N1_31'])), (y, y), color='lightgrey')  # linia dla każdej sekwencji
        plt.text(-600, y + .5, item, va='center', ha='left')  # wyświetlenie nazwy sekwencji

        # Dodanie mutacji do wykresu
        for yy, mutation in enumerate(get_aa_mutations(seqs['Influenza-A-H1N1_31'], seqs[item])):
            pos = int(mutation[1:-1])  # pozycja mutacji
            aa_change = mutation[-1]  # zmiana aminokwasu

            # Oznaczamy delecje jako Δ, jeśli pozycja nie zawiera aminokwasu w wariancie
            if aa_change == '-':
                aa_change = 'Δ'

            # Rysowanie mutacji w różnych pozycjach (aby teksty się nie pokrywały)
            if yy % 3 == 0:
                plt.text(pos, y - .2, aa_change, va='center', ha='center')
            elif yy % 2 == 0:
                plt.text(pos, y, aa_change, va='center', ha='center')
            else:
                plt.text(pos, y + .2, aa_change, va='center', ha='center')

    # Ustawienia zakresu osi
    plt.xlim(-600, len(seqs['Influenza-A-H1N1_31']) + 100)  # margines na lewo dla nazw
    plt.ylim(-.75, len(group) - .25)

    # Tytuł wykresu
    plt.title(f"Mutacje dla {group_name}")

    # Zapis wykresu do pliku w formacie JPEG i TIFF
    jpeg_filename = os.path.join(output_folder, f"{group_name}_mutations-NA.jpeg")
    tiff_filename = os.path.join(output_folder, f"{group_name}_mutations-NA.tiff")

    plt.savefig(jpeg_filename, format='jpeg')
    plt.savefig(tiff_filename, format='tiff')

    plt.close()  # Zamknięcie wykresu po zapisaniu


# Generowanie wykresów i zapis do plików
plot_and_save_mutation_chart(a_seqs, "grypy typu A")
plot_and_save_mutation_chart(b_seqs, "grypy typu B")


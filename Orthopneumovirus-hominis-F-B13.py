from Bio import SeqIO
import matplotlib.pyplot as plt
import os

# Wczytanie sekwencji z pliku FASTA
file_path = "./data/ortho-align.fasta"  # ścieżka do pliku
fasta = list(SeqIO.parse(file_path, format="fasta"))


seqs = {}
for entry in fasta:
    seqs[entry.id] = entry.seq
    

# Funkcja porównująca mutacje do sekwencji referencyjnej i wypisująca różnice oraz liczbę mutacji
# Ustawienie referencyjnej sekwencji
reference_seq = seqs['Orthopneumovirus-hominis-B_13']
def get_mutations(reference, variant):
    # Zamiana sekwencji na listę krotek (pary: nt z referencji, nt z wariantu)
    seqs = list(zip(reference, variant))
    mutations = []
    for pos, nt in enumerate(seqs):  # Iteracja po pozycjach i nukleotydach
        if nt[0] != nt[1]:  # Jeśli nukleotyd w referencji i wariancie są różne
            mutation = f"{nt[0].upper()}{pos + 1}{nt[1].upper()}"  # Zapis różnicy
            mutations.append(mutation)

    return mutations, len(mutations)


# Porównanie wszystkich sekwencji z referencją i wypisanie liczby mutacji
for entry in fasta[1:]:
    print(f"Mutacje w sekwencji: {entry.id}")
    mutations, mutation_count = get_mutations(reference_seq, entry.seq)
    if mutations:
        print(f"Liczba mutacji: {mutation_count}")
        print("Mutacje: ", mutations)
    else:
        print("Brak mutacji")


###### Funkcja obliczająca procentową zgodność z referencją
def calculate_identity(reference, variant):
    matches = sum(nt[0] == nt[1] for nt in zip(reference, variant))  # Zliczamy zgodne pozycje
    return matches / len(reference) * 100  # Obliczamy procent zgodności


# Obliczenie procentowej zgodności wszystkich sekwencji z referencją
identities = []
for entry in fasta[1:]:
    identity = calculate_identity(reference_seq, entry.seq)
    identities.append((entry.id, identity))

# Sortowanie sekwencji na podstawie zgodności
sorted_identities = sorted(identities, key=lambda x: x[1], reverse=True)

# Wypisanie 5 najbardziej i 5 najmniej podobnych sekwencji
print("5 najbardziej podobnych sekwencji:")
for seq_id, identity in sorted_identities[:6]:
    print(f"{seq_id}: {identity:.2f}%")

print("\n5 najmniej podobnych sekwencji:")
for seq_id, identity in sorted_identities[-5:]:
    print(f"{seq_id}: {identity:.2f}%")


#### Białko F od ref NC_038235, region 5661..7385
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

# Sprawdzenie pozycji startu i końca genu F w referencyjnej sekwencji Orthopneumovirus-hominis-B_13
print(gapped_pos(seqs['Orthopneumovirus-hominis-B_13'], 5661))  # 5810
print(gapped_pos(seqs['Orthopneumovirus-hominis-B_13'], 7385))  # 7535


# Wyznaczanie fragmentów genu F (bez przerw "-")
pF = {}
for seq_id, seq in seqs.items():
    pF_seq = seq[5850-1:7574]  # Fragment genu F
    cleaned_seq = pF_seq.replace('-', '')  # Usunięcie przerw "-"
    pF[seq_id] = cleaned_seq  # Zapisanie oczyszczonej sekwencji

# Przetłumaczenie sekwencji nukleotydowych na aminokwasy i zapisanie do pliku FASTA
with open('result/pF_translated.fasta', 'w') as f:  # Nowy plik wynikowy
    for pF_id, pF_seq in pF.items():
        translated_seq = pF_seq.translate()  # Tłumaczenie na aminokwasy
        f.write('>' + pF_id + '\n')  # Zapisanie nazwy sekwencji
        f.write(str(translated_seq) + '\n')  # Zapisanie przetłumaczonej sekwencji aminokwasowej
    f.close()
print("Plik 'pF_translated.fasta' został utworzony.")


file_path2 = "./data/orthoF-B-ali.fasta"  # ścieżka do pliku
pF_aa = list(SeqIO.parse(file_path2, format= "fasta"))
pF_aa

seqs = {}
for entry in pF_aa:
    seqs[entry.id] = entry.seq

def get_aa_mutations(initial, variant):
    out = []
    seqs = list(zip(initial, variant))
    for pos, aa in enumerate(seqs):
        if aa[0] != aa[1]:
            out.append(aa[0].upper() + str(pos) + aa[1].upper())
    return out
    

for item in seqs:
    print(item + ' '+str(len(get_aa_mutations(seqs['Orthopneumovirus-hominis-B_13'], seqs[item]))))

# Ścieżka do folderu, gdzie mają zostać zapisane pliki
output_folder = "./result/plots"

# Funkcja do wyodrębnienia sekwencji na podstawie nazwy wariantu
def get_variant_sequences(seqs, variant_keyword):
    return [seq_id for seq_id in seqs if variant_keyword in seq_id]

# Wyodrębnienie sekwencji dla każdego wariantu
a_seqs = get_variant_sequences(seqs, "Orthopneumovirus-hominis-A")
b_seqs = get_variant_sequences(seqs, "Orthopneumovirus-hominis-B")


# Funkcja rysująca wykresy dla danej grupy sekwencji i zapisująca do plików
def plot_and_save_mutation_chart(group, group_name):
    plt.figure(figsize=(11.7, 8.3))  # Ustawienie rozmiaru wykresu w formacie A4
    
    for y, item in enumerate(group):
        plt.plot((0, len(seqs['Orthopneumovirus-hominis-B_13'])), (y, y), color='lightgrey')  # linia dla każdej sekwencji
        plt.text(-300, y + .5, item, va='center', ha='left')  # wyświetlenie nazwy sekwencji
        
        # Dodanie mutacji do wykresu
        for yy, mutation in enumerate(get_aa_mutations(seqs['Orthopneumovirus-hominis-B_13'], seqs[item])):
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
    plt.xlim(-300, len(seqs['Orthopneumovirus-hominis-B_13']) + 100)  # margines na lewo dla nazw
    plt.ylim(-.75, len(group) - .25)
    
    # Tytuł wykresu
    plt.title(f"Mutacje dla typu {group_name}")

    # Zapis wykresu do pliku w formacie JPEG i TIFF
    jpeg_filename = os.path.join(output_folder, f"{group_name}_mutations.jpeg")
    tiff_filename = os.path.join(output_folder, f"{group_name}_mutations.tiff")
    
    plt.savefig(jpeg_filename, format='jpeg')
    plt.savefig(tiff_filename, format='tiff')
    
    plt.close()  # Zamknięcie wykresu po zapisaniu

# Generowanie wykresów i zapis do plików
plot_and_save_mutation_chart(a_seqs, "Orthopneumovirus-hominis-A")
plot_and_save_mutation_chart(b_seqs, "Orthopneumovirus-hominis-B")






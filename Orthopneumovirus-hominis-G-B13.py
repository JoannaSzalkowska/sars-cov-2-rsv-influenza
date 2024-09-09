
from Bio import SeqIO
import matplotlib.pyplot as plt
import os

# Wczytanie sekwencji z pliku FASTA
file_path = "./data/ortho-align.fasta"  # ścieżka do pliku
fasta = list(SeqIO.parse(file_path, format= "fasta"))


reference_seq = fasta[30].seq

# Funkcja porównująca mutacje do sekwencji referencyjnej i wypisująca różnice oraz liczbę mutacji
def get_mutations(reference, variant):
    # Zamiana sekwencji na listę krotek (pary: nt z referencji, nt z wariantu)
    seqs = list(zip(reference, variant))  
    mutations = []
    for pos, nt in enumerate(seqs):# Iteracja po pozycjach i nukleotydach
        if nt[0] != nt[1]:# Jeśli nukleotyd w referencji i wariancie są różne
            mutation = f"{nt[0].upper()}{pos + 1}{nt[1].upper()}"# Zapis różnicy
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

##### Funkcja porównująca dwie sekwencje i zwracająca tylko procentową zgodność
def compare_sequences(seq1, seq2):
    matches = 0
    total_positions = min(len(seq1), len(seq2))  # porównujemy tylko do najkrótszej sekwencji
    for nt1, nt2 in zip(seq1, seq2):
        if nt1 == nt2:
            matches += 1
    identity = (matches / total_positions) * 100
    return identity

# Lista do przechowywania wyników porównań
similarities = []

# Porównanie każdej sekwencji z każdą i zapisanie wyników
for i, seq1 in enumerate(fasta):
    for j, seq2 in enumerate(fasta[i + 1:], start=i + 1):
        identity = compare_sequences(seq1.seq, seq2.seq)
        similarities.append((seq1.id, seq2.id, identity))

# Sortowanie wyników na podstawie procentowej zgodności (od największej do najmniejszej)
similarities_sorted = sorted(similarities, key=lambda x: x[2], reverse=True)

# Wyświetlenie 10 najbardziej podobnych par sekwencji
print("10 najbardziej podobnych par sekwencji:")
for pair in similarities_sorted[:10]:
    print(f"{pair[0]} vs {pair[1]}: {pair[2]:.2f}%")

# Wyświetlenie 10 najmniej podobnych par sekwencji
print("\n10 najmniej podobnych par sekwencji:")
for pair in similarities_sorted[-10:]:
    print(f"{pair[0]} vs {pair[1]}: {pair[2]:.2f}%")
    
#### sekwencja aminokwasowa
seqs = {}
for entry in fasta:
    seqs[entry.id] = entry

##### białko G od ref NC_038235 lub cds 4688..5584

# Funkcja wyznaczająca pozycję startu i końca genu z uwzględnieniem przerw
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

# Wyznaczanie pozycji startu i końca genu referencyjnego dla sekwencji B
start =gapped_pos(seqs['Orthopneumovirus-hominis-B_13'], 4688)#4732
end = gapped_pos(seqs['Orthopneumovirus-hominis-B_13'], 5584)#5765

# Wyznaczanie fragmentów genu G (bez przerw "-")
pG= {}
for seqrecord in fasta:
    pG_seq = seqrecord.seq[start-1:end]  # Wybieramy odpowiedni fragment genu G
    cleaned_seq = pG_seq.replace('-', '')  # Usuwanie przerw "-"   
    pG[seqrecord.id] = cleaned_seq

# Przetłumaczenie sekwencji nukleotydowych na aminokwasy i zapisanie do pliku FASTA
with open('result/G_translated.fasta', 'w') as f:  # Nowy plik wynikowy
    for pG_id, pG_seq in pG.items():
        translated_seq = pG_seq.translate()  # Tłumaczenie na aminokwasy
        f.write('>' + pG_id + '\n')  # Zapisanie nazwy sekwencji
        f.write(str(translated_seq) + '\n')  # Zapisanie przetłumaczonej sekwencji aminokwasowej
    f.close()
print("Plik 'pG_translated.fasta' został utworzony.")

       
file_path2 = "./data/ortho-prot.fasta"
pG_aa = list(SeqIO.parse(file_path2, format= "fasta"))

# Funkcja do wyznaczania mutacji
seqs = {}
for entry in pG_aa:
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
    
# Funkcja do wyodrębnienia sekwencji na podstawie nazwy wariantu
def get_variant_sequences(seqs, variant_keyword):
    return [seq_id for seq_id in seqs if variant_keyword in seq_id]

# Wyodrębnienie sekwencji dla każdego wariantu
a_seqs = get_variant_sequences(seqs, "Orthopneumovirus-hominis-A")
b_seqs = get_variant_sequences(seqs, "Orthopneumovirus-hominis-B")
    
# Ścieżka do folderu, gdzie mają zostać zapisane pliki
output_folder = "./result/plots"

# Funkcja do rysowania wykresów mutacji i zapisywania plików
def plot_and_save_mutation_chart(group, group_name):
    plt.figure(figsize=(11.7, 8.3))  # Ustawienie rozmiaru wykresu w formacie A4
    
    for y, item in enumerate(group):
        plt.plot((0, len(seqs['Orthopneumovirus-hominis-B_13'])), (y, y), color='lightgrey')
        plt.text(-280, y + .5, item, va='center', ha='left')
        
        # Dodanie mutacji do wykresu
        for yy, mutation in enumerate(get_aa_mutations(seqs['Orthopneumovirus-hominis-B_13'], seqs[item])):
            pos = int(mutation[1:-1])
            aa_change = mutation[-1]
            if aa_change == '-': # delecje
                aa_change = 'Δ'
                
            if yy % 3 == 0:
                plt.text(pos, y - .2, aa_change, va='center', ha='center')
            elif yy % 2 == 0:
                plt.text(pos, y, aa_change, va='center', ha='center')
            else:
                plt.text(pos, y + .2, aa_change, va='center', ha='center')

    plt.xlim(-280, len(seqs['Orthopneumovirus-hominis-B_13']) + 100)
    plt.ylim(-.75, len(group) - .25)
    
    plt.title(f"Mutacje dla typu {group_name}")
    
    jpeg_filename = os.path.join(output_folder, f"{group_name}_mutations333.jpeg")
    tiff_filename = os.path.join(output_folder, f"{group_name}_mutations333.tiff")
    
    plt.savefig(jpeg_filename, format='jpeg')
    plt.savefig(tiff_filename, format='tiff')
    
    plt.close()
    
# Generowanie wykresów
plot_and_save_mutation_chart(get_variant_sequences(seqs, "Orthopneumovirus-hominis-A"), "Orthopneumovirus-hominis-A")
plot_and_save_mutation_chart(get_variant_sequences(seqs, "Orthopneumovirus-hominis-B"), "Orthopneumovirus-hominis-B")

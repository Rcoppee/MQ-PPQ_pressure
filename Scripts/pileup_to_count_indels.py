import pandas as pd
import re

def parse_pileup_with_symbols(pileup_file):
    # Initialisation de la liste pour stocker les résultats
    data = []

    # Ouverture du fichier pileup
    with open(pileup_file, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            chromosome = fields[0]
            position = int(fields[1])
            ref_base = fields[2].upper()  # Base de référence (A, T, C, G)
            reads_count = int(fields[3])
            reads_bases = fields[4]

            # Compter les occurrences des symboles + et -
            plus_count = reads_bases.count('+')
            minus_count = reads_bases.count('-')

            # Identifier le symbole dominant
            symbol_counts = {'+': plus_count, '-': minus_count}
            if plus_count == 0 and minus_count == 0:
                dominant_symbol = 'no'  # Si aucun symbole, le dominant est vide
                dominant_proportion = 0.0
            else:
                dominant_symbol = max(symbol_counts, key=symbol_counts.get)
                dominant_count = symbol_counts[dominant_symbol]
                dominant_proportion = (dominant_count / reads_count) * 100

            # Ajouter les informations à la liste
            data.append([chromosome, position, reads_count, plus_count, minus_count, dominant_symbol, dominant_proportion])

    # Conversion en DataFrame Pandas pour un tableau structuré
    df = pd.DataFrame(data, columns=['Chromosome', 'Position', 'Reads Count', 'Plus Count', 'Minus Count', 'Dominant Symbol', 'Dominant Proportion (%)'])
    return df

# Exemple d'utilisation
pileup_file = 'F78-all.pileup'
df = parse_pileup_with_symbols(pileup_file)

# Sauvegarder le tableau en CSV
df.to_csv('indels_F78.csv', index=False)

# Afficher les premières lignes du tableau
print(df.head())


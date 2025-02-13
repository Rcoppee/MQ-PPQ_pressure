import pandas as pd
import re

def parse_pileup(pileup_file):
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

            # Remplacer les . et , par la base de référence
            reads_bases = reads_bases.replace('.', ref_base).replace(',', ref_base)

            # Supprimer les marqueurs d'indels
            #reads_bases = re.sub(r'[\^][^\^]*[\$]', '', reads_bases)  # Remove start (^) and end ($) of reads
            #reads_bases = re.sub(r'[\+\-][0-9]+[ACGTNacgtn]+', '', reads_bases)  # Remove indels

            #print(reads_bases)



            # Compter les occurrences de chaque base (A, T, C, G)
            a_count = reads_bases.upper().count('A')
            t_count = reads_bases.upper().count('T')
            c_count = reads_bases.upper().count('C')
            g_count = reads_bases.upper().count('G')

            # Identifier le nucléotide dominant
            base_counts = {'A': a_count, 'T': t_count, 'C': c_count, 'G': g_count}
            dominant_base = max(base_counts, key=base_counts.get)
            dominant_count = base_counts[dominant_base]

            # Calculer la proportion du nucléotide dominant
            if reads_count > 0:
                dominant_proportion = (dominant_count / reads_count) * 100
            else:
                dominant_proportion = 0.0

            # Ajouter les informations à la liste
            data.append([chromosome, position, reads_count, a_count, t_count, c_count, g_count, dominant_base, dominant_proportion])


            #if(position==998):
            #    break;

    # Conversion en DataFrame Pandas pour un tableau structuré
    df = pd.DataFrame(data, columns=['Chromosome', 'Position', 'Reads Count', 'A Count', 'T Count', 'C Count', 'G Count', 'Dominant Base', 'Dominant Proportion (%)'])
    return df

# Exemple d'utilisation
pileup_file = 'F78-all.pileup'
df = parse_pileup(pileup_file)

# Sauvegarder le tableau en CSV
df.to_csv('pileup_all-78.csv', index=False)

# Afficher les premières lignes du tableau
print(df.head())

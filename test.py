from Bio import Entrez, SeqIO


# data import
Entrez.email = "your_email@example.com"

handle = Entrez.efetch(
    db="nucleotide",
    id="NM_000527.5",
    rettype="fasta",
    retmode="text"
)

record = SeqIO.read(handle, "fasta")
full_sequence = str(record.seq).lower()

# print(full_sequence[:100])

# Extract coding sequence
cds = full_sequence[86:2669]
print("CDS length:", len(cds))

# codon table
codon_table = {
    'ttt':'F','ttc':'F','tta':'L','ttg':'L',
    'tct':'S','tcc':'S','tca':'S','tcg':'S',
    'tat':'Y','tac':'Y','taa':'*','tag':'*',
    'tgt':'C','tgc':'C','tga':'*','tgg':'W',

    'ctt':'L','ctc':'L','cta':'L','ctg':'L',
    'cct':'P','ccc':'P','cca':'P','ccg':'P',
    'cat':'H','cac':'H','caa':'Q','cag':'Q',
    'cgt':'R','cgc':'R','cga':'R','cgg':'R',

    'att':'I','atc':'I','ata':'I','atg':'M',
    'act':'T','acc':'T','aca':'T','acg':'T',
    'aat':'N','aac':'N','aaa':'K','aag':'K',
    'agt':'S','agc':'S','aga':'R','agg':'R',

    'gtt':'V','gtc':'V','gta':'V','gtg':'V',
    'gct':'A','gcc':'A','gca':'A','gcg':'A',
    'gat':'D','gac':'D','gaa':'E','gag':'E',
    'ggt':'G','ggc':'G','gga':'G','ggg':'G'
}

# Translation function
def translate(dna):
    protein = ""
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i+3]
        protein += codon_table.get(codon, "?")
    return protein

# rs5926
# Original base
mutation_pos = 1920 - 1  # 0-based index
print("Original Nucleobases: ", cds[mutation_pos])

# Oringal codon
codon_start = mutation_pos - (mutation_pos % 3)
original_codon = cds[codon_start:codon_start+3]
print("Original codon:", original_codon)

print("\nThe SNP rs5926 is a mutation in which the last nucleotide of the codon AAC is changed to either G or T.")

# Create mutation
# convert them to list to modify
mut_cds_G = list(cds)
mut_cds_T = list(cds)

# Apply mutations
mut_cds_G[mutation_pos] = "g"
mut_cds_T[mutation_pos] = "t"

# Convert back to string
mut_cds_G = "".join(mut_cds_G)
mut_cds_T = "".join(mut_cds_T)

# translate
protein_normal = translate(cds)
protein_G = translate(mut_cds_G)
protein_T = translate(mut_cds_T)

# example region around mutation
print("example region around mutation")
print("Normal: ", cds[mutation_pos-10:mutation_pos+10])
print("C>G   : ", mut_cds_G[mutation_pos-10:mutation_pos+10])
print("C>T   : ", mut_cds_T[mutation_pos-10:mutation_pos+10])


# identify affected amino acid
aa_pos = mutation_pos // 3

print("\nAmino acid position:", aa_pos + 1)

print("Normal:", protein_normal[aa_pos-10:aa_pos+10])
print("C>G   :", protein_G[aa_pos-10:aa_pos+10])
print("C>T   :", protein_T[aa_pos-10:aa_pos+10])

# Classify mutation
def classify(normal, mutated):
    if normal == mutated:
        return "No change"
    elif mutated == "*":
        return "Nonsense"
    else:
        return "Missense"

print("\nMutation classification:")
print("C>G :", classify(protein_normal[aa_pos], protein_G[aa_pos]))
print("C>T :", classify(protein_normal[aa_pos], protein_T[aa_pos]))
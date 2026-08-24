import os
print(os.getcwd())
with open("sample.fasta", "r") as file:
    lines = file.readlines()

dna = ""

for line in lines:
    if not line.startswith(">"):
        dna += line.strip()

dna = dna.upper()

a = dna.count("A")
t = dna.count("T")
g = dna.count("G")
c = dna.count("C")

length = len(dna)

gc_content = ((g + c) / length) * 100

complement = dna.translate(str.maketrans("ATGC", "TACG"))
reverse_complement = complement[::-1]
rna = dna.replace("T", "U")
start_position = rna.find("AUG")
stop_codon = ["UAA", "UAG", "UGA"]
stop_position = -1
for codon in stop_codon:
    position = rna.find(codon)
    if position != -1:
        stop_position = position
        break
orf = ""
if start_position != -1 and stop_position != -1:
    orf =rna[start_position:stop_position+3] 
codon_table = {
    "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
    "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",

    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",

    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",

    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"
}
orf_protien = ""
for i in range(0, len(orf) - 2, 3):
    codon = orf[i:i+3]
    amino_acid = codon_table.get(codon,
 "?")
    orf_protien += amino_acid
protien = "" 

for i in range(0, len(rna) - 2, 3):
    codon = rna[i:i+3]
    amino_acid = codon_table.get(codon,
"?")
    protien += amino_acid

print("\n=============================="
)   
print("           DNA ANALYSIS")
print("==============================")
print("Sequence:" , dna)
print("Length:" , length)
print("A:", a)
print("T:", t)
print("G:", g)
print("C:", c)
print("GC Content:", round(gc_content,
2), "%")
print("Reverse Complement:",
reverse_complement)
print("RNA:", rna)
print("protien:", protien)
print("start codon position:",
start_position)
print("stop position codon:",
stop_position)
print("ORF:", orf)
print("ORF Sequence:", orf)
print("Protien Sequence:", orf_protien)
print("==============================")
print("      ANALYSIS COMPLETE"       )
print("==============================")
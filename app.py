from flask import Flask, render_template, request

app = Flask(__name__)


codon_table = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",

    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",

    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",

    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G"
}


def analyze_dna(dna):

    dna = dna.upper()
    dna = dna.replace(" ", "").replace("\n", "").replace("\r", "")

    valid = set("ATGC")
    if any(base not in valid for base in 
    dna):
         return None, "Invalid DNA sequence. Use only A, T, G, and C."

    if not dna:
        return None, "No DNA sequence found."

    if not set(dna).issubset(valid):
        return None, "Invalid DNA sequence."

    length = len(dna)

    a = dna.count("A")
    t = dna.count("T")
    g = dna.count("G")
    c = dna.count("C")

    gc = ((g + c) / length) * 100

    rna = dna.replace("T", "U")

    complement_table = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }

    complement = "".join(
        complement_table[x] for x in dna
    )

    reverse_complement = complement[::-1] 

    start = dna.find("ATG")

    stop = -1

    if start != -1:

        stop_codons = ["TAA", "TAG", "TGA"]

        for i in range(start + 3, length - 2, 3):

            codon = dna[i:i+3]

            if codon in stop_codons:
                stop = i
                break


    orf = ""

    if start != -1 and stop != -1:

        orf = dna[start:stop+3]


    protein = ""

    if orf:

        for i in range(0, len(orf)-2, 3):

            codon = orf[i:i+3]

            protein += codon_table.get(codon, "?")


    result = {

        "sequence": dna,
        "length": length,

        "A": a,
        "T": t,
        "G": g,
        "C": c,

        "GC": round(gc, 2),
        "AT": round(((a + t) / length) * 100, 2),

        "RNA": rna,

        "complement": complement,

        "reverse_complement": reverse_complement,

        "start_position": start if start != -1 else "Not found",

        "stop_position": stop if stop != -1 else "Not found",

        "orf": orf if orf else "Not found",

        "protein": protein if protein else "Not found"

    }


    return result, None



@app.route("/", methods=["GET", "POST"])
def home():

    result = None
    error = None
    filename = None


    if request.method == "POST":

        uploaded = request.files.get("fasta_file")


        if uploaded and uploaded.filename:

            filename = uploaded.filename

            content = uploaded.read().decode("utf-8")

            lines = content.splitlines()

            dna = "".join(
                line.strip()
                for line in lines
                if not line.startswith(">")
            )


        else:

            dna = request.form.get("dna", "")


        result, error = analyze_dna(dna)



    return render_template(
        "index.html",
        result=result,
        error=error,
        filename=filename
    )



if __name__ == "__main__":

    app.run(host="0.0.0.0", port=5000,
    debug=True) 
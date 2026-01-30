# Algorithm/main.py
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from builder.plasmid_builder import assemble_plasmid

def main():
    if len(sys.argv) != 3:
        print("Usage: python main.py <input_fasta> <design_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    design_file = sys.argv[2]

    plasmid_sequence = assemble_plasmid(fasta_file, design_file)

    record = SeqRecord(
        Seq(plasmid_sequence),
        id="Synthetic_Plasmid",
        description="Plasmid assembled using ORI prediction and design rules"
    )

    SeqIO.write(record, "Output.fa", "fasta")
    print("Output plasmid saved as Output.fa")

if __name__ == "__main__":
    main()

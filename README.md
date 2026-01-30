# Plasmid Construction Pipeline

This project is a computational pipeline for designing and constructing plasmids from a given input FASTA file and a design file. The pipeline identifies the origin of replication (ORI), incorporates user-defined genes and restriction sites, and outputs a final plasmid sequence.

## Project Structure

```
.
├── Algorithm
│   ├── builder
│   │   ├── plasmid_builder.py
│   │   ├── restriction_sites.json
│   │   └── restriction_sites.py
│   ├── main.py
│   └── ori_finder.py
├── Antibiotic_Marker
│   ├── Ampicillin.fa
│   ├── Chloramphenicol.fa
│   └── Kanamycin.fa
├── Input
│   ├── Design_pUC19.txt
│   ├── markers.tab
│   └── pUC19.fa
├── LICENSE
└── README.md
```

## How it Works

The pipeline consists of the following main components:

1.  **ORI Finder (`ori_finder.py`)**: This module is responsible for identifying the origin of replication (ORI) in the input genome. It uses a two-step process:
    *   **GC Skew Analysis**: It first calculates the cumulative GC skew of the sequence to find the approximate location of the ORI.
    *   **MEME-like Motif Discovery**: It then uses a MEME-like algorithm to search for conserved motifs within the region identified by the GC skew analysis. This helps to refine the ORI prediction.

2.  **Plasmid Builder (`plasmid_builder.py`)**: This module assembles the final plasmid sequence. It performs the following steps:
    *   It calls the ORI finder to get the ORI sequence.
    *   It parses a user-provided design file to determine which genes (e.g., antibiotic resistance markers) and restriction sites to include in the plasmid.
    *   It loads the sequences of the specified genes from the `Antibiotic_Marker` directory.
    *   It removes any internal restriction sites from the plasmid sequence to ensure that the final plasmid can be correctly assembled.
    *   It adds the desired restriction sites to form a multiple cloning site (MCS).

3.  **Main Script (`main.py`)**: This is the entry point of the pipeline. It takes the input FASTA file and the design file as command-line arguments, orchestrates the plasmid construction process, and saves the final plasmid sequence to a file named `Output.fa`.

## How to Run

1.  **Prerequisites**: Make sure you have Python and the Biopython library installed:
    ```bash
    pip install biopython
    ```

2.  **Run the pipeline**:
    ```bash
    python Algorithm/main.py <input_fasta> <design_file>
    ```
    For example:
    ```bash
    python Algorithm/main.py Input/pUC19.fa Input/Design_pUC19.txt
    ```

3.  **Output**: The pipeline will generate a file named `Output.fa` containing the sequence of the newly constructed plasmid.

## Input Files

*   **Input FASTA (`<input_fasta>`)**: A FASTA file containing the DNA sequence of the plasmid or genome from which the ORI will be extracted.
*   **Design File (`<design_file>`)**: A text file that specifies the components to be included in the plasmid. Each line in the file should be in the format `component_type,component_name`. For example:
    ```
    BamHI_site,BamHI
    AmpR_gene,Ampicillin
    ```

## Included Data

*   **`Antibiotic_Marker` directory**: Contains FASTA files for common antibiotic resistance genes (Ampicillin, Chloramphenicol, Kanamycin).
*   **`Input` directory**:
    *   `pUC19.fa`: The sequence of the pUC19 cloning vector.
    *   `Design_pUC19.txt`: An example design file for constructing a plasmid based on pUC19.
    *   `markers.tab`: A table of common molecular biology markers.
*   **`Algorithm/builder/restriction_sites.json`**: A JSON file containing a list of common restriction enzymes and their recognition sequences.
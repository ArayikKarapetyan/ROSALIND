import requests
from Bio import SeqIO
from io import StringIO
import re

def fetch_protein_sequence_biopython(protein_id):
    """Fetch protein sequence using BioPython with ID cleaning."""
    # STRIP THE ID: 'P01866_GCB_MOUSE' -> 'P01866'
    accession = protein_id.split('_')[0]
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    
    try:
        # Increased timeout and updated URL to the current REST API
        response = requests.get(url, timeout=15)
        if response.status_code == 200:
            fasta_data = StringIO(response.text)
            record = SeqIO.read(fasta_data, "fasta")
            return str(record.seq)
        return None
    except Exception as e:
        # Hidden print to avoid messing up Rosalind output format
        return None

def find_n_glycosylation_motif_biopython(sequence):
    """Find motif using regex with lookahead for overlapping matches."""
    positions = []
    # N{P}[ST]{P} -> N, (not P), (S or T), (not P)
    pattern = re.compile(r'(?=(N[^P][ST][^P]))')
    
    for match in pattern.finditer(sequence):
        positions.append(match.start() + 1)
    
    return positions

# Load your dataset
try:
    with open("Dataset.txt", "r") as file:
        protein_ids = [line.strip() for line in file if line.strip()] 
except FileNotFoundError:
    print("File not found.")
    protein_ids = []

# Process and Print in Rosalind Format
for protein_id in protein_ids:
    seq = fetch_protein_sequence_biopython(protein_id)
    if seq:
        matches = find_n_glycosylation_motif_biopython(seq)
        if matches:
            # Print original ID from dataset and the found positions
            print(protein_id)
            print(" ".join(map(str, matches)))
import re
import requests

def fetch_protein_sequence(protein_id):
    """
    Fetches the protein sequence from UniProt.
    Extracts the Accession ID to ensure the URL is valid even for IDs 
    containing underscores (e.g., P01133_EGF_HUMAN).
    """
    # Use only the accession part (before the underscore) for the URL
    accession = protein_id.split('_')[0]
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    
    try:
        response = requests.get(url, timeout=15)
        if response.status_code == 200:
            lines = response.text.strip().split('\n')
            # Join all lines except the header
            return "".join(lines[1:])
        else:
            print(f"Failed to fetch {protein_id} (Status: {response.status_code})")
            return None
    except Exception as e:
        print(f"Error fetching {protein_id}: {e}")
        return None

def find_glycosylation_motif(sequence):
    """
    Finds the N-glycosylation motif N{P}[ST]{P} using a lookahead 
    assertion to capture overlapping occurrences.
    """
    # Pattern: N followed by (not P), then (S or T), then (not P)
    # (?=...) is a lookahead that finds the match without 'consuming' characters
    pattern = r'(?=(N[^P][ST][^P]))'
    
    # We find all matches and return their 1-based start positions
    return [match.start() + 1 for match in re.finditer(pattern, sequence)]

def main():
    # Read the dataset from your provided file name
    # Based on 
    try:
        with open("Dataset.txt", "r") as file:
            protein_ids = [line.strip() for line in file if line.strip()]
    except FileNotFoundError:
        print("Dataset file not found.")
        return

    for protein_id in protein_ids:
        sequence = fetch_protein_sequence(protein_id)
        
        if sequence:
            positions = find_glycosylation_motif(sequence)
            
            # Rosalind requires: Only print ID and positions if the motif is found
            if positions:
                print(protein_id)
                print(*(positions))

if __name__ == "__main__":
    main()
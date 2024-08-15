import pandas as pd
import json



def process_gene_names(row):
    """Iterates through a tsv formatted list of genes from uniprot to create a list of unique genes.

    Args:
        row (pd.Series): The data of the row as a series

    Returns:
        list: a list containing all unique gene names
    """    
    all_genes = []
    for col in ['Gene Names (primary)', 'Gene Names']:
        if isinstance(row[col], str):
            all_genes.extend(row[col].split())  # Split into individual genes
    return list(set(all_genes))  # Return unique gene names

if __name__ == "__main__":
    # Extracted from Uniprot - genes for Homo sapiens
    df = pd.read_csv("./uniprot_human.tsv", sep="\t") 
    # Apply the processing function
    df['combined_genes'] = df.apply(process_gene_names, axis=1)

    # Create the dictionary
    gene_entry_dict = {}
    for index, row in df.iterrows():
        entry = row['Entry']
        for gene in row['combined_genes']:
            gene_entry_dict[gene] = entry 

    # Save as JSON file where key-value pair is
    with open('./gene_entry.json', 'w') as outfile:
        json.dump(gene_entry_dict, outfile, indent=4)
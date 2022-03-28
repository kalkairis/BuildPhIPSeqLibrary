import logging
import os
import urllib.request

import pandas as pd

from Figures.figures_config import allergens_raw_data, allergens_dir, infectious_dir, infectious_raw_data


def get_uniprot_fasta(sequence_id):
    url = f"http://www.uniprot.org/uniprot/?query={sequence_id}&format=fasta"
    with urllib.request.urlopen(url) as r:
        fasta = r.read().decode('utf-8').strip()
    return ''.join(fasta.split('\n')[1:])


def get_sequence_and_convert_to_input_format(dir_path, file_path):
    df = pd.read_csv(file_path, header=[0, 1], low_memory=False)
    uniprot_to_fasta_df = pd.DataFrame({'Parent_Protein_IRI': df[('Epitope', 'Parent Protein IRI')].unique()}).dropna()
    # Leave only uniprot references
    uniprot_to_fasta_df = uniprot_to_fasta_df[
        uniprot_to_fasta_df.Parent_Protein_IRI.str.startswith('http://www.uniprot.org/uniprot/')]
    uniprot_to_fasta_df['sequence_ID'] = uniprot_to_fasta_df['Parent_Protein_IRI'].apply(lambda url: url.split('/')[-1])
    uniprot_to_fasta_df['AA_sequence'] = uniprot_to_fasta_df.sequence_ID.apply(get_uniprot_fasta)
    input_path = os.path.join(dir_path, 'Input')
    os.makedirs(input_path, exist_ok=True)
    uniprot_to_fasta_df.to_csv(os.path.join(input_path, 'input_df.csv'))


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    logging.info("Converting allergens")
    get_sequence_and_convert_to_input_format(allergens_dir, allergens_raw_data)
    logging.info("Converting infectious diseases")
    get_sequence_and_convert_to_input_format(infectious_dir, infectious_raw_data)

"""
This code maps the origin as well as the possible mappings of each oligo.
Note that this is NOT part of the main library building code.
"""
import time


def map_single_oligo_to_sequences_list(oligo, sequences_dict):
    ret = []
    if '*' in oligo:
        oligo = oligo[:oligo.index('*')]
    for sequence_id, sequence in sequences_dict.items():
        for position in re.finditer(oligo, sequence):
            ret.append((sequence_id, position.span()[0]))
    return ret

# def map_sequences(oligos_aa_sequences, id_to_sequences):
#     print("Merge oligos sequences, this takes time", time.ctime())
#     #TODO: read existing oligos. base_df = read_oligo_sequences_to_file()
#
#     # Stage 1: add oligos already in origins into origins list and remove from new sequences
#     existing_oligos = oligos_aa_sequences.index.intersection(base_df['oligo_aa_sequence'].values)
#     for existing_oligo in existing_oligos:
#         base_df.loc[base_df['oligo_aa_sequence'].eq(existing_oligo).idxmax(), 'origins'] += \
#             oligos_aa_sequences.loc[existing_oligo]['origins']
#     oligos_aa_sequences.drop(index=existing_oligos, inplace=True)
#
#     # Stage 2: map newly added sequences into existing oligos list
#     base_df['mapped'] += base_df['oligo_aa_sequence'].apply(
#         lambda oligo: map_single_oligo_to_sequences_list(oligo, id_to_sequences))
#
#     # Stage 3: map newly added oligos to all sequences (from sequences IDs file)
#     all_ids_to_sequences = read_sequence_ids_file()['AA_sequence'].to_dict()
#     oligos_aa_sequences.reset_index(inplace=True)
#     oligos_aa_sequences['mapped'] = oligos_aa_sequences['oligo_aa_sequence'].apply(
#         lambda oligo: map_single_oligo_to_sequences_list(oligo, all_ids_to_sequences))
#     running_index = 0 if len(base_df) == 0 else base_df.index.str.split('_').str[1].astype(int).max() + 1
#     oligos_aa_sequences['oligo_id'] = list(
#         map(lambda idx: f'oligo_{idx}', range(running_index, running_index + len(oligos_aa_sequences))))
#     oligos_aa_sequences.set_index('oligo_id', inplace=True)
#
#     # Stage 4: concatenate dfs
#     ret = pd.concat([base_df, oligos_aa_sequences], axis=0, ignore_index=False)
#     ret.to_csv(OLIGO_SEQUENCES_FILE)
#     return ret, oligos_aa_sequences
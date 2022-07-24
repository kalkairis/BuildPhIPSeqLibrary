import matplotlib.patches as mpatches
import pandas
import os

from statannotations.Annotator import Annotator
import seaborn as sns
import numpy
import matplotlib.pyplot as plt
import matplotlib

from BuildPhIPSeqLibrary.config import AMINO_INFO
from BuildPhIPSeqLibrary.read_pipeline_files import read_barcoded_nucleotide_files


def is_in(origins, check_in):
    for seq_ID, p in origins:
        if seq_ID in check_in:
            return True
    return False

def create_subfigure_2_a(ax, fig, outer_spec):
    codon_usage = AMINO_INFO.set_index('codon')[['amino_acid', 'corrected_relative_frequency']]

    path = os.path.join(os.path.dirname(__file__), 'data')
    file_template = os.path.join(path, "%s_barcoded_nuc_file.csv")
    for f in ['infectious', 'allergens']:
        nuc_seqs = read_barcoded_nucleotide_files(file_template % f)['nuc_sequence']
        codon_usage[f] = 0
        lseq = len(nuc_seqs.iloc[0])
        for i in range(0, lseq, 3):
            vs = nuc_seqs.str[i:i + 3].value_counts()
            codon_usage.loc[vs.index, f] += vs

    stops = codon_usage[codon_usage.amino_acid == '*']
    codon_usage = codon_usage[codon_usage.amino_acid != '*']

    for AA, df in codon_usage.groupby('amino_acid'):
        for f in ['infectious', 'allergens']:
            codon_usage.loc[df.index, f + "_prop"] = df[f] / df[f].sum()
    # for f in ['infectious', 'allergens']:
    #     worse_codon = ((codon_usage.corrected_relative_frequency - codon_usage[f + "_prop"]) /
    #                    codon_usage.corrected_relative_frequency).abs().sort_values().index[-1]
    #     print("worse:")
    #     print(codon_usage[codon_usage.amino_acid == codon_usage.loc[worse_codon, 'amino_acid']])

    codon_usage = pandas.concat([codon_usage, stops[['amino_acid']]]).fillna(0)

    letter_order = ['T', 'C', 'A', 'G']
    palette = {'Theoretical': 'limegreen', 'Infectious': 'orangered', 'Allergens': 'gold'}
    plt.rcParams['font.size'] = 16
    # fig = plt.figure(figsize=(22, 21))
    handles = []
    for k, v in palette.items():
        handles.append(mpatches.Patch(color=v, label=k))
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.02),
               handles=handles, ncol=3, fontsize=22)
    ax.axis('off')
    letter_box_annot_opts = {'fontsize': 24, 'xy': (0.5, 0.5), 'xycoords': 'axes fraction', 'va': 'center', 'ha': 'center'}
    spec = outer_spec.subgridspec(5, 6, width_ratios=[1, 5, 5, 5, 5, 1], height_ratios=[1, 5, 5, 5, 5], wspace=0, hspace=0)
    shared_axis = None
    for col, col_letter in enumerate(letter_order):
        for row, row_letter in enumerate(letter_order):
            ax = fig.add_subplot(spec[row + 1, col + 1])
            ax.set_facecolor('azure')
            ax.set_xticks([])
            ax.set_yticks([])
            inner_spec = spec[row + 1, col + 1].subgridspec(1, 2, width_ratios=[5, 1], wspace=0)

            ax = fig.add_subplot(inner_spec[0], sharex=shared_axis, sharey=shared_axis, xlim=(0, 1))
            if shared_axis is None:
                shared_axis = ax
            if row < 3:
                ax.get_xaxis().set_visible(False)
            ax.yaxis.set_ticks([])

            sub_df = codon_usage[codon_usage.index.str.startswith(row_letter + col_letter)].set_index('amino_acid',
                                                                                                      append=True)
            sub_df = sub_df[['corrected_relative_frequency', 'infectious_prop', 'allergens_prop']].rename(
                columns={'corrected_relative_frequency': 'Theoretical', 'infectious_prop': 'Infectious',
                         'allergens_prop': 'Allergens'}).stack().reset_index()
            sub_df.columns = ['codon', 'amino_acid', 'src', 'prop']
            sns.barplot(data=sub_df, x='prop', y='codon', hue='src', ax=ax,
                        order=list(map(lambda l: row_letter + col_letter + l, letter_order)),
                        palette=palette)
            ax.set_facecolor('azure')
            box_pairs = []
            for c in sub_df.amino_acid.unique():
                box_pairs.append((tuple(sub_df[sub_df.amino_acid.eq(c)].iloc[0][['codon', 'src']].values.tolist()),
                                  tuple(sub_df[sub_df.amino_acid.eq(c)].iloc[-1][['codon', 'src']].values.tolist())))
            annot = Annotator(ax, box_pairs, data=sub_df, x='prop', y='codon', hue='src',
                              order=list(map(lambda l: row_letter + col_letter + l, letter_order)), orient='h')
            annot.configure(test=None, loc='outside', fontsize=30).set_custom_annotations(sub_df.amino_acid.unique())
            annot.annotate()
            ax.legend_.remove()
            ax.yaxis.set_ticks([])
            ax.set_ylabel('')
            ax.set_xlabel('')
            ax.spines['right'].set_visible(False)

            if row == 3:
                ax.set_xticks([tick for tick in ax.get_xticks() if tick < 1.1])

    for row, letter in enumerate(letter_order):
        ax = fig.add_subplot(spec[row + 1, 0], facecolor='steelblue')
        ax.annotate(letter, **letter_box_annot_opts)
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])

        ax = fig.add_subplot(spec[row + 1, -1], facecolor='skyblue')
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        inner_spec = spec[row + 1, -1].subgridspec(4, 1, hspace=0)
        for i, letter_2 in enumerate(letter_order):
            ax = fig.add_subplot(inner_spec[i, 0], facecolor='skyblue')
            ax.annotate(letter_2, **letter_box_annot_opts)
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            ax.axis('off')

    for col, letter in enumerate(letter_order):
        ax = fig.add_subplot(spec[0, col + 1], facecolor='steelblue')
        ax.annotate(letter, **letter_box_annot_opts)
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])

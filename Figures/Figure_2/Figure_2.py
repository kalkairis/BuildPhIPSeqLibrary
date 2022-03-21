import matplotlib.pyplot as plt
import pandas as pd
import os
import seaborn as sns
import numpy as np
import string

from Figures.Figure_2.Figure_2a import create_subfigure_2_a


def subfigure_c(ax):
    # Simulations
    df = pd.read_csv(os.path.join(os.path.dirname(__file__), 'data', 'simulations.csv'), index_col=0)
    df = df[df.barcode_length_in_aa < 5]
    df['log_num_oligos'] = np.log10(df['num_oligos'])
    sns.lineplot(data=df, x='barcode_length_in_aa', y='log_num_oligos', hue='num_barcodes', err_style='bars', ci='sd',
                 palette={3: sns.color_palette('Set2')[1], 5: sns.color_palette('Set2')[2], 1: sns.color_palette('Set2')[0]}, ax=ax)
    plt.xticks(range(1, 5), ["%d aa" % i for i in range(1, 5)])
    plt.xlabel("Number AAs per barcode part")
    plt.ylabel("Log average number of random AAs\nuntil first no-coding")
    plt.yticks(range(1, 7), ["10", "100", "1000", "$10^{4}$", "$10^{5}$", "$10^{6}$"])
    plt.legend(
        ["1 barcode part (no correction)", "3 barcode parts, correct 1 error", "5 barcode parts, correct 2 errors"])
    # plt.tight_layout()
    # plt.savefig("Fig2c.png")
    # plt.close()


def subfigure_b(ax):
    # Simulations
    df = pd.read_csv(os.path.join(os.path.dirname(__file__), 'data', 'sim_divide.csv'))
    sns.barplot(data=df, x='len_in_AAs', y='0', hue='num_parts',
                     palette={3: sns.color_palette('Set2')[1], 5: sns.color_palette('Set2')[2], 7: sns.color_palette('Set2')[4]}, ax=ax)
    plt.yscale('log')
    plt.xlabel("Total number AAs for all barcode parts")
    plt.ylabel("Log average number of radon AAs\nuntil first no-coding")
    ax.legend_.set_title('Number of barcode parts')
    ax.legend_.texts[0].set_text("3 parts, correct 1 error")
    ax.legend_.texts[1].set_text("5 parts, correct 2 error")
    ax.legend_.texts[2].set_text("7 parts, correct 3 error")
    # plt.tight_layout()
    # plt.savefig("Fig2b.png")
    # plt.close()


def create_data_b():
    dfs = {}
    for i in range(5, 11):
        dfs[i] = pd.read_csv("/home/sigall/PycharmProjects/BuildPhIPSeqLibrary/Figures/Figure_4/all_opts_%d.csv" % i,
                             index_col=0)
    df = pd.concat(dfs)
    df.index.names = ['len_in_AAs', 'num_parts']
    df.columns.names = ['sim']
    df.stack().to_csv(os.path.join(os.path.dirname(__file__), 'data', 'sim_divide.csv'))


if __name__ == "__main__":
    fig = plt.figure(figsize=(37, 21))
    spec = fig.add_gridspec(1, 2, width_ratios=[27, 11], wspace=0.2)

    # sub figure a
    ax = fig.add_subplot(spec[0,0])
    text_size = 30
    ax.text(-0.01, 1.05, string.ascii_lowercase[0], transform=ax.transAxes, size=text_size, weight='bold')
    create_subfigure_2_a(ax, fig, spec[0, 0])

    # create_data_b()
    # Create subfigures b and c
    right_spec = spec[0, 1].subgridspec(2, 1)
    ax = fig.add_subplot(right_spec[0, 0])
    ax.text(-0.2, 1.1, string.ascii_lowercase[1], transform=ax.transAxes, size=text_size, weight='bold')
    subfigure_b(ax)
    ax = fig.add_subplot(right_spec[1, 0])
    ax.text(-0.2, 1.1, string.ascii_lowercase[2], transform=ax.transAxes, size=text_size, weight='bold')
    subfigure_c(ax)
    plt.savefig("Figure_2.png")

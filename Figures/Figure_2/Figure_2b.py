import matplotlib.pyplot as plt
import pandas as pd
import os
import seaborn as sns
import numpy as np


def subfigure_c():
    # Simulations
    df = pd.read_csv(os.path.join(os.path.dirname(__file__), 'data', 'simulations.csv'), index_col=0)
    df = df[df.barcode_length_in_aa < 5]
    df['log_num_oligos'] = np.log10(df['num_oligos'])
    sns.lineplot(data=df, x='barcode_length_in_aa', y='log_num_oligos', hue='num_barcodes', err_style='bars', ci='sd',
                 palette={3: sns.color_palette('Set2')[1], 5: sns.color_palette('Set2')[2], 1: sns.color_palette('Set2')[0]})
    plt.xticks(range(1, 5), ["%d aa" % i for i in range(1, 5)])
    plt.xlabel("Num AAs per barcode part")
    plt.ylabel("Log average number of random AAs\nuntil first no-coding")
    plt.yticks(range(1, 7), ["10", "100", "1000", "$10^{4}$", "$10^{5}$", "$10^{6}$"])
    plt.legend(
        ["1 barcode part (no correction)", "3 barcode parts, correct 1 error", "5 barcode parts, correct 2 errors"])
    plt.tight_layout()
    plt.savefig("Fig2c.png")
    plt.close()


def subfigure_b():
    # Simulations
    df = pd.read_csv(os.path.join(os.path.dirname(__file__), 'data', 'sim_divide.csv'))
    ax = sns.barplot(data=df, x='len_in_AAs', y='0', hue='num_parts',
                     palette={3: sns.color_palette('Set2')[1], 5: sns.color_palette('Set2')[2], 7: sns.color_palette('Set2')[4]})
    plt.yscale('log')
    plt.xlabel("Total num AAs for all barcode parts")
    plt.ylabel("Log average number of radon AAs\nuntil first no-coding")
    ax.legend_.set_title('Number of barcode parts')
    ax.legend_.texts[0].set_text("3 parts, correct 1 error")
    ax.legend_.texts[1].set_text("5 parts, correct 2 error")
    ax.legend_.texts[2].set_text("7 parts, correct 3 error")
    plt.tight_layout()
    plt.savefig("Fig2b.png")
    plt.close()


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

    # פfig = plt.figure(figsize=(22, 32))

    # create_data_b()
    subfigure_b()
    subfigure_c()

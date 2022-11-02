import matplotlib.pyplot as plt

from BuildPhIPSeqLibrary.config import AMINO_INFO

if __name__ == "__main__":
    amino_acid_df = AMINO_INFO[['amino_acid', 'codon']].copy()
    for i in range(3):
        amino_acid_df[f'codon_{i}'] = amino_acid_df['codon'].str.get(i)
    codon_values = {'T': 0, 'C': 1, 'A': 2, 'G': 3}

    amino_acid_df.sort_values(by=['codon_0', 'codon_1', 'codon_2'], key=lambda c: c.map(codon_values), inplace=True)

    fig, ax = plt.subplots()
    size = 0.3
    colors = ['#315B96', '#A7384C', '#22448D', '#BF82BB']
    outer = amino_acid_df.groupby(['codon_0', 'codon_1', 'codon_2']).count().reset_index().sort_values(
        by=['codon_0', 'codon_1', 'codon_2'], key=lambda c: c.map(codon_values))
    ax.pie(outer.amino_acid.values.flatten(),
           radius=1,
           # labeldistance=1 - (size/2),
           startangle=90,
           # labels=outer['codon_2'],
           colors=colors,
           wedgeprops=dict(width=0.2, edgecolor='w'))

    center = amino_acid_df.groupby(['codon_0', 'codon_1']).count().reset_index().sort_values(by=['codon_0', 'codon_1'],
                                                                                             key=lambda c: c.map(
                                                                                                 codon_values))
    ax.pie(center.amino_acid.values.flatten(),
           radius=0.8,
           labels=center['codon_1'],
           startangle=90,
           colors=colors,
           labeldistance=0.8,
           textprops={'size': 14, 'va': 'center', 'ha': 'center', 'color': 'white'},
           wedgeprops=dict(width=0.5, edgecolor='w'))

    inner = amino_acid_df.groupby('codon_0').count().reset_index().sort_values(by=['codon_0'],
                                                                               key=lambda c: c.map(codon_values))
    ax.pie(inner.amino_acid.values.flatten(),
           radius=0.5,
           labels=inner['codon_0'],
           startangle=90,
           labeldistance=0.5,
           textprops={'size': 24, 'va': 'center', 'ha': 'center', 'color': 'white'},
           colors=colors,
           wedgeprops=dict(width=0.5, edgecolor='w'))

    adj_check = (amino_acid_df.amino_acid != amino_acid_df.amino_acid.shift()).cumsum()
    amino_acid_wedges = amino_acid_df.groupby(['amino_acid', adj_check], as_index=False, sort=False).count()
    ax.pie(

        amino_acid_wedges.codon.values.flatten(),
        radius=1.4,
        labels=amino_acid_wedges.amino_acid,
        colors=['lightgrey'],
        startangle=90,
        labeldistance=0.85,
        textprops={'size': 11, 'va': 'center', 'ha': 'center', 'color': 'black'},
        wedgeprops=dict(width=0.4, edgecolor='grey'))
    plt.savefig('figure_1_amino_acid_circle.png')

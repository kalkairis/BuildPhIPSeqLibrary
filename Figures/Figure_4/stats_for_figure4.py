import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time
from BuildPhIPSeqLibrary.config import AMINO_ACIDS

hamming_distance = 1
df = pd.DataFrame({'length_of_barcode': np.arange(3, 100, 3)})
df['num_of_oligos_random'] = 4**df.length_of_barcode

from_web_table = """UUU F 0.57 19.7 (   101)
UCU S 0.11  5.7 (    29)
UAU Y 0.53 16.8 (    86)
UGU C 0.42  5.9 (    30)
UUC F 0.43 15.0 (    77)
UCC S 0.11  5.5 (    28)
UAC Y 0.47 14.6 (    75)
UGC C 0.58  8.0 (    41)
UUA L 0.15 15.2 (    78)
UCA S 0.15  7.8 (    40)
UAA * 0.64  1.8 (     9)
UGA * 0.36  1.0 (     5)
UUG L 0.12 11.9 (    61)
UCG S 0.16  8.0 (    41)
UAG * 0.00  0.0 (     0)
UGG W 1.00 10.7 (    55)
CUU L 0.12 11.9 (    61)
CCU P 0.17  8.4 (    43)
CAU H 0.55 15.8 (    81)
CGU R 0.36 21.1 (   108)
CUC L 0.10 10.5 (    54)
CCC P 0.13  6.4 (    33)
CAC H 0.45 13.1 (    67)
CGC R 0.44 26.0 (   133)
CUA L 0.05  5.3 (    27)
CCA P 0.14  6.6 (    34)
CAA Q 0.30 12.1 (    62)
CGA R 0.07  4.3 (    22)
CUG L 0.46 46.9 (   240)
CCG P 0.55 26.7 (   137)
CAG Q 0.70 27.7 (   142)
CGG R 0.07  4.1 (    21)
AUU I 0.58 30.5 (   156)
ACU T 0.16  8.0 (    41)
AAU N 0.47 21.9 (   112)
AGU S 0.14  7.2 (    37)
AUC I 0.35 18.2 (    93)
ACC T 0.47 22.8 (   117)
AAC N 0.53 24.4 (   125)
AGC S 0.33 16.6 (    85)
AUA I 0.07  3.7 (    19)
ACA T 0.13  6.4 (    33)
AAA K 0.73 33.2 (   170)
AGA R 0.02  1.4 (     7)
AUG M 1.00 24.8 (   127)
ACG T 0.24 11.5 (    59)
AAG K 0.27 12.1 (    62)
AGG R 0.03  1.6 (     8)
GUU V 0.25 16.8 (    86)
GCU A 0.11 10.7 (    55)
GAU D 0.65 37.9 (   194)
GGU G 0.29 21.3 (   109)
GUC V 0.18 11.7 (    60)
GCC A 0.31 31.6 (   162)
GAC D 0.35 20.5 (   105)
GGC G 0.46 33.4 (   171)
GUA V 0.17 11.5 (    59)
GCA A 0.21 21.1 (   108)
GAA E 0.70 43.7 (   224)
GGA G 0.13  9.2 (    47)
GUG V 0.40 26.4 (   135)
GCG A 0.38 38.5 (   197)
GAG E 0.30 18.4 (    94)
GGG G 0.12  8.6 (    44)"""
from_web_table = pd.DataFrame({'whole_str': from_web_table.split('\n')})
from_web_table['codon'] = from_web_table.whole_str.str.split().str[0]
from_web_table['amino_acid'] = from_web_table.whole_str.str.split().str[1]
from_web_table['fraction'] = from_web_table.whole_str.str.split().str[2].astype(float)
from_web_table['frequency_per_thousand'] = from_web_table.whole_str.str.split().str[3].astype(float)

frequency_per_aa = from_web_table.groupby('amino_acid').frequency_per_thousand.sum().drop(index='*').to_frame()
frequency_per_aa['frequency_per_thousand'] /= frequency_per_aa['frequency_per_thousand'].sum()

num_options_per_aa = from_web_table.groupby('amino_acid').codon.count()

MAX_TRY = 10**7
stats = {}

num_bars = 5
for barcode_aa_length in range(1, 7):
    rounds = []
    if barcode_aa_length < 5:
        ROUNDS = 1000
    else:
        ROUNDS = 100
    for round_i in range(ROUNDS):
        if (round_i % 5) == 0:
            print(f"For {barcode_aa_length}: in round {round_i}", time.ctime())
        res = {}
        for n in range(num_bars):
            res[n] = {}
        flag_MAX = True
        for i in range(MAX_TRY):
            for n in range(num_bars):
                sequence = np.random.choice(frequency_per_aa.index, size=barcode_aa_length,
                                            p=frequency_per_aa['frequency_per_thousand'])
                num_options = res[n].get(''.join(sequence), np.prod(list(map(num_options_per_aa.get, sequence))))
                if num_options == 0:
                    flag_MAX = False
                    break
                res[n][''.join(sequence)] = num_options - 1
            if not flag_MAX:
                break
        if flag_MAX:
            rounds.append(np.nan)
        else:
            rounds.append(i)
    rounds = pd.Series(rounds)
    rounds.to_csv(f"round_for_{barcode_aa_length}_{num_bars}.csv")
    print("For (%d * %d): %d of %d at max (%d)" % (barcode_aa_length, num_bars, rounds.isna().sum(), ROUNDS, MAX_TRY))
    if rounds.isna().sum() < ROUNDS:
        print("mean of non - max % g(median % d)" % (rounds.mean(), rounds.median()))
    stats[barcode_aa_length] = [rounds.isna().sum(), ROUNDS, rounds.mean(), rounds.median()]
print(pd.DataFrame(stats))

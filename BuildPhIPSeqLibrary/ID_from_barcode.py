import argparse
import pandas
import os
import numpy

from BuildPhIPSeqLibrary.config import BARCODE_NUC_LENGTHS, MAPPED_OLIGO_SEQUENCES_FILE, BARCODE_IN_5_PRIME_END

from BuildPhIPSeqLibrary.read_pipeline_files import read_barcoded_nucleotide_files, read_sequence_ids_file
from BuildPhIPSeqLibrary.construct_nucleotide_sequences import get_barcode_from_nuc_seq


def get_maps(ID):
    if os.path.exists(MAPPED_OLIGO_SEQUENCES_FILE):
        return pandas.read_csv(MAPPED_OLIGO_SEQUENCES_FILE, index_col=0).loc[ID]
    return None


def find_ID(inp, ham, used_barcode, allow_indel, df_barcodes, print_flag):
    # allow no indels - as barcode was planned
    start_location = 0
    IDs = {}
    for i, barcode_length in enumerate(used_barcode):
        part = inp[start_location:start_location + barcode_length]
        if not BARCODE_IN_5_PRIME_END:
            part = part[::-1]
        if part in df_barcodes[f'barcode_{i}'].values:
            IDs[(0, i)] = df_barcodes[df_barcodes[f'barcode_{i}'] == part].index[0]
        start_location += barcode_length
    vIDs = pandas.Series(IDs, dtype=str).value_counts()
    if len(vIDs) > 0:
        errs = len(used_barcode) - vIDs.iloc[0]
        if errs <= ham:
            if print_flag:
                print("barcode identified with %d errors" % errs)
            return vIDs.index[0], 0, errs
    if ham == 0:
        if print_flag:
            print("Can't ID barcode")
        return None, 0, 0

    if not allow_indel:
        return None, 0, 0

    # allow single indel - because usually possible
    for shft in [-1, 1]:
        start_location = shft
        for i, barcode_length in enumerate(used_barcode):
            if (start_location >= 0) and ((start_location+barcode_length) < len(inp)):
                part = inp[start_location:start_location + barcode_length]
                if not BARCODE_IN_5_PRIME_END:
                    part = part[::-1]
                if part in df_barcodes[f'barcode_{i}'].values:
                    IDs[(shft, i)] = df_barcodes[df_barcodes[f'barcode_{i}'] == part].index[0]
            start_location += barcode_length
    df_IDs = pandas.Series(IDs, dtype=str)
    found = []
    for ID in df_IDs.unique():
        errs = len(used_barcode) - df_IDs[df_IDs == ID].index.get_level_values(1).nunique()
        if errs <= ham:
            found.append(ID)
    if len(found) == 1:
        errs = len(used_barcode) - df_IDs[df_IDs == found[0]].index.get_level_values(1).nunique()
        if print_flag:
            print("barcode identified with 1 indel and %d errors" % max(0, (errs-1)))
        return found[0], 1, max(0, (errs-1))
    if (len(found) > 1) or (ham == 1):
        if print_flag:
            print("Can't ID barcode with 1 indel")
        return None, 0, 0

    # allow 2 indels - should we?
    for shft in [-2, 2]:
        start_location = shft
        for i, barcode_length in enumerate(used_barcode):
            if (start_location >= 0) and ((start_location+barcode_length) < len(inp)):
                part = inp[start_location:start_location + barcode_length]
                if not BARCODE_IN_5_PRIME_END:
                    part = part[::-1]
                if part in df_barcodes[f'barcode_{i}'].values:
                    IDs[(shft, i)] = df_barcodes[df_barcodes[f'barcode_{i}'] == part].index[0]
            start_location += barcode_length
    df_IDs = pandas.Series(IDs)
    found = []
    for ID in df_IDs.unique():
        errs = len(used_barcode) - df_IDs[df_IDs == ID].index.get_level_values(1).nunique()
        if errs <= ham:
            found.append(ID)
    if len(found) == 1:
        errs = len(used_barcode) - df_IDs[df_IDs == found[0]].index.get_level_values(1).nunique()
        if print_flag:
            print("barcode identified with 2 indel and %d errors" % max(0, (errs-2)))
        return found[0], 2, max(0, (errs-2))

    if print_flag:
        print("Can't ID barcode with up to 2 indels")
    return None, 0, 0


def out_sources(ID, print_flag):
    maps = get_maps(ID)
    if maps is None:
        if print_flag:
            print("No mapping file exists. Can't locate origins of IDed %s" % ID)
        return 0, 0
    origins = eval(maps.origins)
    maps = eval(maps.mapped)
    seqs = read_sequence_ids_file()
    if print_flag:
        print("\nThe origins of IDed oligo %s:" % ID)
    for o in origins:
        if print_flag:
            print("Position %d of %s (file %s)" % (o[1], seqs.loc[o[0]].sequence_ID, seqs.loc[o[0]].input_file))
        maps.pop(maps.index(o))
    if len(maps) > 0:
        if print_flag:
            print("And also maps to:")
        for o in maps:
            if print_flag:
                print("Position %d of %s (file %s)" % (o[1], seqs.loc[o[0]].sequence_ID, seqs.loc[o[0]].input_file))
    return origins, maps


def find_and_output(inp, allow_indel, df_barcodes=None, print_flag=False):
    if df_barcodes is None:
        df_barcodes = read_barcoded_nucleotide_files()

    len_barcode = sum(BARCODE_NUC_LENGTHS)
    if len(inp) < len_barcode:
        print("Only partial barcode was given")
    found = False
    for ham in range(int((len(BARCODE_NUC_LENGTHS) - 1) / 2), -1, -1):
        used_barcode = BARCODE_NUC_LENGTHS[:2 * ham + 1]
        if sum(used_barcode) <= len(inp):
            if print_flag:
                print("Checking for barcode of length %d (of %d given), at hamming distance <= %d" % (sum(used_barcode),
                                                                                                      len(inp), ham))
            found = True
            break
    assert found, "Can't ID with less then 1 part of the barcode. Should be %s." % (str(BARCODE_NUC_LENGTHS))
    for i, _ in enumerate(used_barcode):
        assert ('barcode_%d' % i) in df_barcodes.columns, "Given barcode file had no barcode_%d column" % i

    ID, indels, errs = find_ID(inp, ham, used_barcode, allow_indel, df_barcodes, print_flag)
    if ID is None:
        return ID, indels, errs, [], []

    origins, maps = out_sources(ID, print_flag)
    return ID, indels, errs, origins, maps


def run_test(num_rounds):
    df_barcodes = read_barcoded_nucleotide_files()

    cnt = [0, 0]
    for i in range(num_rounds):
        oli = numpy.random.choice(df_barcodes.index)
        if BARCODE_IN_5_PRIME_END:
            input = df_barcodes.loc[oli].nuc_sequence[:sum(BARCODE_NUC_LENGTHS)]
        else:
            input = df_barcodes.loc[oli].nuc_sequence[-sum(BARCODE_NUC_LENGTHS):][::-1]
        ID, indels, errs, _, _ = find_and_output(input, False, df_barcodes)
        if ID == oli:
            cnt[0] += 1
            if (indels == 0) and (errs == 0):
                cnt[1] += 1

        pos = numpy.random.randint(0, sum(BARCODE_NUC_LENGTHS))
        ID, indels, errs, _, _ = find_and_output(input[:pos] + "N" + input[pos+1:], False, df_barcodes)
        if ID == oli:
            cnt[0] += 1
            if (indels == 0) and (errs == 1):
                cnt[1] += 1

        pos2 = numpy.random.randint(0, sum(BARCODE_NUC_LENGTHS))
        tmp_input = input[:pos] + "N" + input[pos + 1:]
        tmp_input = tmp_input[:pos2] + "N" + tmp_input[pos2 + 1:]
        ID, indels, errs, _, _ = find_and_output(tmp_input, False, df_barcodes)
        if ID == oli:
            cnt[0] += 1
            if (indels == 0) and (errs <= 2):
                cnt[1] += 1

        ID, indels, errs, _, _ = find_and_output(input[:pos] + input[pos + 1:] + "N", True, df_barcodes)
        if ID == oli:
            cnt[0] += 1
            if ((indels <= 1) and ((indels + errs) <= 1)) or ((indels == 0) and (errs <=2)):
                cnt[1] += 1

        tmp_input = input[:pos] + "N" + input[pos + 1:]
        tmp_input = tmp_input[:pos2] + tmp_input[pos2 + 1:] + "N"
        ID, indels, errs, _, _ = find_and_output(tmp_input, True, df_barcodes)
        if ID == oli:
            cnt[0] += 1
            if (indels <= 1) and ((indels + errs) <= 2):
                cnt[1] += 1

        tmp_input = input[:pos] + input[pos + 1:] + "N"
        tmp_input = tmp_input[:pos2] + tmp_input[pos2 + 1:] + "N"
        ID, indels, errs, _, _ = find_and_output(tmp_input, True, df_barcodes)
        if ID == oli:
            cnt[0] += 1
            if (indels <= 2) and ((indels + errs) <= 2):
                cnt[1] += 1

    print(cnt, "of %d" % (num_rounds*6))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ID the source of a barcode.")
    parser.add_argument('barcode', type=str, help='the barcode to identify')
    parser.add_argument('--allow_indel', dest='allow_indel', default=False, action='store_true')
    args = parser.parse_args()

    # run_test(10)
    ID, indels, errs, origins, other_maps = find_and_output(args.barcode, args.allow_indel, print_flag=True)
    print()

import pandas as pd
import numpy as np
from operator import attrgetter
import dask.dataframe as dd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--pgen_pdata_file', help='a file containing counts of public sequences in positive-labeled '
                                              'repertoires and all repertoires, generation probability of sequences '
                                              'an indication whether the sequence is ground truth signal or not. '
                                              'These information are expected to be supplied through the following '
                                              'fields respectively: "pos_count", "public_count", "pgen", "ref_seq"')
parser.add_argument('--outfile', help='the output file with a mapping between discerete intervals of generation '
                                      'probability and population incidence')
parser.add_argument('--total_sample_size', help='Total sample size to be used in computation of population incidence '
                                                'of sequences', type=int)
parser.add_argument('--poslabel_sample_size', help='Sample size of just the positive-labaled repertoires to be used '
                                                   'in computation of population incidence of signal sequences', type=int)
parser.add_argument('--ground_truth_true', help='an optional flag when enabled generates the mapping for signal '
                                                'sequences. If not supplied generates for all the remaining public '
                                                'sequences', action="store_true",)
args = parser.parse_args()


def _read_and_validate(pgen_pdata_file):
    df = pd.read_csv(pgen_pdata_file, index_col=0, header=0, sep='\t')
    obligatory_columns = ['pos_count', 'public_count', 'pgen', 'ref_seq']
    assert all(col_name in df.columns for col_name in obligatory_columns), "Expected required fields are missing"
    return df


def _compute_count_proportions(df, total_sample_size, poslabel_sample_size, ground_truth_true=False):
    ddf = dd.from_pandas(df, npartitions=10)
    df_relevant = ddf[ddf['ref_seq'] == ground_truth_true]
    if not ground_truth_true:
        df_relevant = df_relevant.map_partitions(
            lambda dat: dat.assign(count_prop=dat.public_count / total_sample_size))
    else:
        df_relevant = df_relevant.map_partitions(
            lambda dat: dat.assign(count_prop=dat.pos_count / poslabel_sample_size))
    df_relevant = df_relevant[['pgen', 'count_prop']].compute()
    return df_relevant


def _bin_sequences(df):
    df['count_prop_bins'] = pd.cut(df['count_prop'],
                                   bins=np.unique(np.nanpercentile(df['count_prop'], np.array(
                                       [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99, 100]))),
                                   include_lowest=True)
    df['pgen_bins'] = pd.cut(df['pgen'],
                             bins=np.array(
                                 [0, 1e-20, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-09, 1e-08, 1e-02]),
                             include_lowest=True)
    return df


def _calculate_bin_counts(df):
    pgen_grouped = df.groupby('pgen_bins', observed=True)
    bin_counts = pgen_grouped['count_prop_bins'].value_counts().rename('values')
    bin_counts = bin_counts.to_frame().reset_index()
    return bin_counts


def _calculate_probabilities(bin_counts):
    bin_counts['prob'] = bin_counts['values'] / bin_counts.groupby('pgen_bins')['values'].transform('sum')
    bin_counts['pgen_left'] = bin_counts['pgen_bins'].map(attrgetter('left'))
    bin_counts['pgen_right'] = bin_counts['pgen_bins'].map(attrgetter('right'))
    bin_counts['sample_size_prop_left'] = bin_counts['count_prop_bins'].map(attrgetter('left'))
    bin_counts['sample_size_prop_right'] = bin_counts['count_prop_bins'].map(attrgetter('right'))
    return bin_counts


def _write_bin_counts_to_file(bin_counts, outfile):
    bin_counts_write = bin_counts[
        ['pgen_left', 'pgen_right', 'sample_size_prop_left', 'sample_size_prop_right', 'prob']]
    bin_counts_write.to_csv(outfile, sep='\t', index=None)


def generate_pgen_count_map(pgen_pdata_file, outfile, total_sample_size, poslabel_sample_size, ground_truth_true=False):
    df = _read_and_validate(pgen_pdata_file)
    df_relevant = _compute_count_proportions(df, total_sample_size, poslabel_sample_size,
                                             ground_truth_true=ground_truth_true)
    df_relevant = _bin_sequences(df_relevant)
    bin_counts = _calculate_bin_counts(df_relevant)
    bin_counts = _calculate_probabilities(bin_counts)
    _write_bin_counts_to_file(bin_counts, outfile)


def execute():
    generate_pgen_count_map(args.pgen_pdata_file, args.outfile, args.total_sample_size, args.poslabel_sample_size,
                            args.ground_truth_true)

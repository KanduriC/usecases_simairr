import os
import pandas as pd
import argparse
import numpy as np
import glob
from multiprocessing import Pool
from scipy.stats import binom

parser = argparse.ArgumentParser()
parser.add_argument('--files_path',
                    help='the path where pgen files and corresponding sequence count files are found',
                    required=True)
parser.add_argument('--n_threads',
                    help='number of threads to be used for parallel processing (usually as many threads as pgen files',
                    type=int, required=True)
parser.add_argument('--metadata_file',
                    help='the path to metadata file containing subject_id column and a label field of bool dtype',
                    required=True)
parser.add_argument('--label_field',
                    help='column in metadata file that corresponds to label field of bool dtype', required=True)
parser.add_argument('--avg_rep_size',
                    help='average repertoire size used when computing probability based on binomial distribution',
                    type=int, required=True)
args = parser.parse_args()


def compute_pval(pos_count, pgen, total_pos_count, avg_rep_size):
    p_occur_once = binom.sf(0, avg_rep_size, pgen)
    p_occur_by_chance = binom.sf(pos_count, total_pos_count, p_occur_once)
    return p_occur_by_chance


def _merge_pgen_seqcounts(pgen_file_path, seqcounts_file_path):
    pgen_file = pd.read_csv(pgen_file_path, header=None, sep='\t', error_bad_lines=False)
    counts_file = pd.read_csv(seqcounts_file_path, header=None, sep='\t', error_bad_lines=False)
    pgen_file.columns = ['aa_seq_pgen', 'pgen']
    counts_file.columns = ['aa_seq_counts', 'v_gene', 'j_gene', 'pos_count', 'neg_count']
    pgen_file['row_index_pgen'] = np.arange(len(pgen_file))
    counts_file['row_index_counts'] = np.arange(len(counts_file))
    merged_df = pd.merge(counts_file, pgen_file, how="left", left_on=['aa_seq_counts', 'row_index_counts'],
                         right_on=['aa_seq_pgen', 'row_index_pgen'])
    merged_df = merged_df[['aa_seq_counts', 'v_gene', 'j_gene', 'pos_count', 'neg_count', 'pgen']]
    return merged_df


class ComputePval:
    def __init__(self, files_path, n_threads, avg_rep_size, metadata_file, label_field):
        self.files_path = files_path
        self.n_threads = n_threads
        self.avg_rep_size = avg_rep_size
        self.metadata_file = metadata_file
        self.label_field = label_field

    def concat_compute_pval(self, seqcounts_file_path):
        pgen_file_path = os.path.join(os.path.dirname(seqcounts_file_path),
                                      "pgen_" + os.path.basename(seqcounts_file_path))
        merged_df = _merge_pgen_seqcounts(pgen_file_path, seqcounts_file_path)
        merged_df['pval_pos'] = merged_df.apply(
            lambda x: compute_pval(x['pos_count'], x['pgen'], len(self.pos_label_examples),
                                   self.avg_rep_size), axis=1)
        merged_df['pval_neg'] = merged_df.apply(
            lambda x: compute_pval(x['neg_count'], x['pgen'], len(self.neg_label_examples),
                                   self.avg_rep_size), axis=1)
        pgen_pdata_file_path = os.path.join(os.path.dirname(pgen_file_path),
                                            'pdata_' + os.path.basename(pgen_file_path))
        merged_df.to_csv(pgen_pdata_file_path, sep='\t', index=None)

    def multi_compute_pval(self):
        self._set_labels()
        count_files = glob.glob(self.files_path + "/public_seq_*", recursive=True)
        pool = Pool(self.n_threads)
        pool.map(self.concat_compute_pval, count_files)

    def _set_labels(self):
        metadata = pd.read_csv(self.metadata_file, header=0, index_col=None, error_bad_lines=False)
        self.pos_label_examples = metadata[metadata[self.label_field]]["subject_id"].to_list()
        self.neg_label_examples = metadata[~metadata[self.label_field]]["subject_id"].to_list()


def execute():
    pval_test = ComputePval(files_path=args.files_path, n_threads=args.n_threads, avg_rep_size=args.avg_rep_size,
                            metadata_file=args.metadata_file, label_field=args.label_field)
    pval_test.multi_compute_pval()

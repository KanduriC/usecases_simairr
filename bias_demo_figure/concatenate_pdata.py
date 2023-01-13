import os
import glob
import argparse
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests

parser = argparse.ArgumentParser()
parser.add_argument('--dir_path', help='path to pdata_pgen_* files')
parser.add_argument('--reference_seq_file', help='path to file containing true reference sequences')
args = parser.parse_args()


def concatenate_pdata_files(dir_path, reference_seq_file):
    ref_seq = pd.read_csv(reference_seq_file, header=0, index_col=None, error_bad_lines=False, sep='\t')
    found_files = glob.glob(dir_path + "/pdata_pgen_*", recursive=True)
    li = []
    for i, filename in enumerate(found_files):
        print("processing file number:", filename)
        fn = pd.read_csv(filename, header=0, sep='\t', error_bad_lines=False)
        li.append(fn)
    outfile_df = pd.concat(li, axis=0)
    outfile_df['qval_pos'] = multipletests(outfile_df['pval_pos'], method="fdr_bh")[1]
    outfile_df['qval_neg'] = multipletests(outfile_df['pval_neg'], method="fdr_bh")[1]
    outfile_df = pd.merge(outfile_df, ref_seq, how="left", indicator='ref_seq')
    outfile_df['ref_seq'] = np.where(outfile_df.ref_seq == "both", True, False)
    outfile_df['public_count'] = outfile_df['neg_count'] + outfile_df['pos_count']
    public_df = outfile_df[outfile_df['public_count'] > 1]
    private_df = outfile_df[outfile_df['public_count'] == 1]
    public_df.to_csv(os.path.join(dir_path, "concatenated_pdata_pgen_public.tsv"), sep='\t')
    private_df.to_csv(os.path.join(dir_path, "concatenated_pdata_pgen_private.tsv"), sep='\t')


def execute():
    concatenate_pdata_files(args.dir_path, args.reference_seq_file)

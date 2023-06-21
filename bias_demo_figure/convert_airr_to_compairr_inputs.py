import argparse
import pandas as pd
import numpy as np
import re

parser = argparse.ArgumentParser()
parser.add_argument('--metadata_file', help='metadata file containing a filename column with path to files')
parser.add_argument('--repertoire_concat_file', help='Output file with concatenated repertoires')
parser.add_argument('--seq_concat_file', help='Output file with concatenated sequences of repertoires')
args = parser.parse_args()


def concatenate_airr_files(metadata_file, repertoire_concat_file, seq_concat_file):
    meta = pd.read_csv(metadata_file, header=0)
    outfile_df = _read_and_concatenate_repertoires(meta)
    outfile_df['duplicate_count'] = 1
    if 'v_genes' in outfile_df.columns:
        relevant_cols = ['repertoire_id', 'sequence_id', 'duplicate_count', 'v_genes', 'j_genes', 'junction', 'junction_aa']
        modified_colnames = ['repertoire_id', 'sequence_id', 'duplicate_count', 'v_call', 'j_call', 'junction', 'junction_aa']
    else:
        relevant_cols = ['repertoire_id', 'sequence_id', 'duplicate_count', 'v_call', 'j_call', 'junction', 'junction_aa']
        modified_colnames = ['repertoire_id', 'sequence_id', 'duplicate_count', 'v_call', 'j_call', 'junction',
                             'junction_aa']
    outfile_df_rep = outfile_df[relevant_cols].copy()
    outfile_df_rep.columns = modified_colnames
    outfile_df_rep.to_csv(repertoire_concat_file, index=None, sep="\t")
    if 'v_genes' in outfile_df.columns:
        outfile_df = outfile_df.drop_duplicates(subset=['junction_aa', 'v_genes', 'j_genes'])
    else:
        outfile_df = outfile_df.drop_duplicates(subset=['junction_aa', 'v_call', 'j_call'])
    outfile_df['repertoire_id'] = 0
    outfile_df['sequence_id'] = range(1, 1 + len(outfile_df))
    outfile_df_seq = outfile_df[relevant_cols].copy()
    outfile_df_seq.columns = modified_colnames
    outfile_df_seq.to_csv(seq_concat_file, index=None, sep="\t")


def _read_and_concatenate_repertoires(meta):
    li = []
    for i, filename in enumerate(meta['filename']):
        print("processing file number:", i)
        fn = pd.read_csv(filename, header=0, sep='\t')
        fn.insert(loc=0, column='repertoire_id', value=meta['subject_id'][i])
        if 'v_genes' in fn.columns:
            fn[['v_subgroups', "v_genes"]] = fn[['v_subgroups', "v_genes"]].fillna('TRBJ2-2P')
            fn['v_genes'] = np.vectorize(_replace_trb_gene)(fn["v_genes"], fn["v_subgroups"])
        li.append(fn)
    outfile_df = pd.concat(li, axis=0, ignore_index=True)
    return outfile_df


def _replace_trb_gene(v_gene_col, j_gene_col):
    return re.sub(r"^TRBJ.*|''", j_gene_col, v_gene_col)


def execute():
    concatenate_airr_files(args.metadata_file, args.repertoire_concat_file, args.seq_concat_file)
import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--metadata_file', help='metadata file containing a filename column with path to files')
parser.add_argument('--repertoires_path', help='Path containing repertoire files')
args = parser.parse_args()


def generate_compairr_input_from_olga_format(metadata_file, repertoires_path):
    repertoire_concat_file = os.path.join(os.path.dirname(metadata_file), "concatenated_repertoire_file.tsv")
    seq_concat_file = os.path.join(os.path.dirname(metadata_file), "concatenated_sequence_file.tsv")
    meta = pd.read_csv(metadata_file, header=0)
    li = _read_all_repertoires(meta, repertoires_path)
    outfile_df = pd.concat(li, axis=0, ignore_index=True)
    outfile_df['sequence_id'] = range(1, 1 + len(outfile_df))
    outfile_df = outfile_df[
        ['repertoire_id', 'sequence_id', 'duplicate_count', 'v_call', 'j_call', 'junction', 'junction_aa']]
    outfile_df.to_csv(repertoire_concat_file, index=None, sep="\t")
    outfile_df = outfile_df.drop_duplicates(subset=['junction_aa', 'v_call', 'j_call'])
    outfile_df['repertoire_id'] = 0
    outfile_df['sequence_id'] = range(1, 1 + len(outfile_df))
    outfile_df.to_csv(seq_concat_file, index=None, sep="\t")


def _read_all_repertoires(meta, repertoires_path):
    li = []
    for i, filename in enumerate(meta['filename']):
        print("processing file number:", i)
        fn = pd.read_csv(os.path.join(repertoires_path, filename), header=None, index_col=None, sep='\t')
        fn.columns = ['junction', 'junction_aa', 'v_call', 'j_call']
        fn['repertoire_id'] = os.path.splitext(filename)[0]
        fn['duplicate_count'] = 1
        li.append(fn)
    return li


def execute():
    generate_compairr_input_from_olga_format(args.metadata_file, args.repertoires_path)

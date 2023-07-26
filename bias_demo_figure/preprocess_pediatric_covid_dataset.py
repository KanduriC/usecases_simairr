import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--metadata_file', help='metadata file containing a filename column with path to files')
parser.add_argument('--repertoire_concat_file', help='Output file with concatenated repertoires')
parser.add_argument('--seq_concat_file', help='Output file with concatenated sequences of repertoires')
args = parser.parse_args()


def generate_compairr_input(metadata_file, repertoire_concat_file, seq_concat_file):
    meta = pd.read_csv(metadata_file, header=0)
    outfile_df = _read_and_concatenate_repertoires(meta)
    outfile_df.to_csv(repertoire_concat_file, index=None, sep="\t")
    outfile_df = outfile_df.drop_duplicates(subset=['junction_aa', 'v_call', 'j_call'])
    outfile_df['repertoire_id'] = 0
    outfile_df['sequence_id'] = range(1, 1 + len(outfile_df))
    outfile_df.to_csv(seq_concat_file, index=None, sep="\t")


def _read_and_concatenate_repertoires(meta):
    li = []
    for i, filename in enumerate(meta['filename']):
        print("processing file number:", i)
        fn = pd.read_csv(filename, header=0, sep='\t')
        fn['sequence_id'] = range(0, len(fn))
        fn['duplicate_count'] = 1
        fn.insert(loc=0, column='repertoire_id', value=meta['subject_id'][i])
        fn = drop_missing_amino_acids(fn, amino_acid_col="aminoAcid")
        fn = clean_and_fill(fn, "vGeneName", "vFamilyName")
        fn = clean_and_fill(fn, "jGeneName", "jFamilyName")
        relevant_colnames = ["repertoire_id", "sequence_id", "duplicate_count", "vGeneName", "jGeneName", "nucleotide",
                             "aminoAcid"]
        modified_colnames = ['repertoire_id', 'sequence_id', 'duplicate_count', 'v_call', 'j_call', 'junction',
                             'junction_aa']
        fn = fn[relevant_colnames]
        fn.columns = modified_colnames
        fn = drop_missing_amino_acids(fn, amino_acid_col="v_call")
        fn = drop_missing_amino_acids(fn, amino_acid_col="j_call")
        fn = fn[fn['v_call'] != 'TCRBV21-01']
        fn = fn[fn['v_call'] != 'TCRBV07-05']
        li.append(fn)
    outfile_df = pd.concat(li, axis=0, ignore_index=True)
    return outfile_df


def drop_missing_amino_acids(df, amino_acid_col):
    df.dropna(subset=[amino_acid_col], inplace=True)
    return df


def clean_and_fill(df, column1, column2):
    df[column1] = df[column1].str.split('/').str[0]
    df[column1] = df[column1].fillna(df[column2])
    return df


def execute():
    generate_compairr_input(args.metadata_file, args.repertoire_concat_file, args.seq_concat_file)

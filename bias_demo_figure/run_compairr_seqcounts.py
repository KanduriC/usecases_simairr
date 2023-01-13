import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--concatenated_rep_file',
                    help='concatenated repertoire file in AIRR format as desired by compairr')
parser.add_argument('--concatenated_seq_file', help='concatenated sequence file in AIRR format as desired by compairr')
parser.add_argument('--compairr_path', help='path to compairr binary file')
parser.add_argument('--batch_size', help='number of rows (sequences) to process in chunks', type=int, required=True)
parser.add_argument('--ignore_genes', help='binary (0 or 1) indicating whether to ignore or not ignore genes', type=int,
                    required=True)
args = parser.parse_args()


def run_compairr_seqcount(concatenated_rep_file, concatenated_seq_file, compairr_path, batch_size, ignore_genes=0):
    """
    cuts the concatenated sequence file into chunks and runs compairr on each chunk against the concatenated repertoire
    file

    """
    seqfile_path = os.path.dirname(concatenated_seq_file)
    chunkfiles_path = os.path.join(seqfile_path, "seqfile_chunks")
    if not os.path.exists(chunkfiles_path):
        os.makedirs(chunkfiles_path)
    output_files_path = os.path.join(seqfile_path, "compairr_outputfiles")
    if not os.path.exists(output_files_path):
        os.makedirs(output_files_path)
    for i, chunk in enumerate(pd.read_csv(concatenated_seq_file, chunksize=batch_size, header=0, sep='\t')):
        chunk_seq_file = os.path.join(chunkfiles_path, 'chunk{}.tsv'.format(i))
        chunk.to_csv(chunk_seq_file, index=False, sep='\t')
        compairr_output_file = os.path.join(output_files_path, 'compairr_design_mat_{}.txt'.format(i))
        compairr_log_file = os.path.join(output_files_path, 'compairr_logfile_{}.txt'.format(i))
        if ignore_genes == 0:
            os.system(
                compairr_path + ' -x ' + chunk_seq_file + ' ' + concatenated_rep_file + ' -f -g -o ' + compairr_output_file + ' -l ' + compairr_log_file + ' -u -t 8')
        else:
            os.system(
                compairr_path + ' -x ' + chunk_seq_file + ' ' + concatenated_rep_file + ' -f -o ' + compairr_output_file + ' -l ' + compairr_log_file + ' -u -t 8')


def execute():
    run_compairr_seqcount(args.concatenated_rep_file, args.concatenated_seq_file, args.compairr_path, args.batch_size,
                          args.ignore_genes)
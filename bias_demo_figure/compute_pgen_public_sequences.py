import glob
import os
import argparse
import pandas as pd
from multiprocessing import Pool
from scipy.stats import binom

parser = argparse.ArgumentParser()
parser.add_argument('--design_mats_path',
                    help='the path to sequence presence matrix files written by compairr',
                    required=True)
parser.add_argument('--seq_concat_file', help='Output file with concatenated sequences of repertoires', required=True)
parser.add_argument('--legal_vgenes_file', help='a text file with legal V gene names', required=True)
parser.add_argument('--legal_jgenes_file', help='a text file with legal J gene names', required=True)
parser.add_argument('--n_threads',
                    help='number of threads to be used for parallel processing (usually as many thread as design matrix'
                         'files', type=int, required=True)
parser.add_argument('--metadata_file',
                    help='the path to metadata file containing subject_id column and a label field of bool dtype',
                    required=True)
parser.add_argument('--label_field',
                    help='column in metadata file that corresponds to label field of bool dtype', required=True)
parser.add_argument('--avg_rep_size',
                    help='average repertoire size used when computing probability based on binomial distribution',
                    type=int, required=True)
parser.add_argument('--public_counts_in',
                    help='if 0, public counts based on pos_label, if 1, public counts based on neg_labels. If 2, '
                         'based on both',
                    type=int, required=True)
parser.add_argument('--compute_pgen',
                    help='default is True to compute pgen (assumes training data)', required=True)
parser.add_argument('--vj_gene_conditioning',
                    help='if 1, sequences are conditioned on V and J genes before computing pgen, else no conditioning',
                    type=int, required=True)
args = parser.parse_args()


def _get_legal_genes(legal_genes_file):
    genes = pd.read_csv(legal_genes_file, header=None, sep='\t')
    genes.columns = ['gene']
    gene_list = genes['gene'].to_list()
    return gene_list


def _compute_pval(pos_count, pgen, total_pos_count, avg_rep_size):
    p_occur_once = binom.sf(0, avg_rep_size, pgen)
    p_occur_by_chance = binom.sf(pos_count, total_pos_count, p_occur_once)
    return p_occur_by_chance


def _get_output_file_path(compairr_design_mat):
    return os.path.join(os.path.dirname(compairr_design_mat),
                        'public_seq_' + os.path.basename(compairr_design_mat))


def _read_design_matrix(compairr_design_mat):
    return pd.read_csv(compairr_design_mat, header=0, index_col=0, sep='\t', error_bad_lines=False)


def _compute_counts(design_mat, label_examples):
    return design_mat[label_examples].astype(bool).sum(axis=1)


def _combine_counts(counts_df_pos, counts_df_neg):
    return pd.concat([counts_df_pos, counts_df_neg], axis=1)


def _write_to_file(genefilt_df, out_file_path):
    genefilt_df.to_csv(out_file_path, sep='\t', index=None, header=None)


def _get_pgen_file_path(out_file_path):
    return os.path.join(os.path.dirname(out_file_path), 'pgen_' + os.path.basename(out_file_path))


def _run_olga_tool(command):
    exit_code = os.system(command)
    if exit_code != 0:
        raise RuntimeError(f"Running olga tool failed:{command}.")


def _read_pgen_file(pgen_file_path):
    pgen_file = pd.read_csv(pgen_file_path, header=None, index_col=0, sep='\t')
    pgen_file.columns = ['pgen']
    return pgen_file


class ComputePgenPublic:
    def __init__(self, design_mats_path, seq_concat_file_df, legal_vgenes_file, legal_jgenes_file, n_threads,
                 metadata_file, label_field, avg_rep_size,
                 public_counts_in=1, compute_pgen=True, vj_gene_conditioning=1):
        self.design_mats_path = design_mats_path
        self.seq_concat_file_df = seq_concat_file_df
        self.legal_vgenes_file = legal_vgenes_file
        self.legal_jgenes_file = legal_jgenes_file
        self.n_threads = n_threads
        self.metadata_file = metadata_file
        self.label_field = label_field
        self.pos_label = self.label_field + '_true'
        self.neg_label = self.label_field + '_false'
        self.avg_rep_size = avg_rep_size
        self.public_counts_in = public_counts_in
        self.compute_pgen = compute_pgen
        self.vj_gene_conditioning = vj_gene_conditioning

    def compute_pgen_public_sequences(self, compairr_design_mat):
        vgene_list, jgene_list = self._read_legal_genes()
        out_file_path = _get_output_file_path(compairr_design_mat)
        design_mat = _read_design_matrix(compairr_design_mat)
        metadata = self._read_metadata()
        pos_label_examples = self._get_pos_label_examples(metadata)
        neg_label_examples = self._get_neg_label_examples(metadata)
        counts_df_pos = _compute_counts(design_mat, pos_label_examples)
        counts_df_neg = _compute_counts(design_mat, neg_label_examples)
        counts_df = _combine_counts(counts_df_pos, counts_df_neg)
        counts_df = self._add_labels(counts_df)
        counts_df = self._filter_counts(counts_df)
        merged_counts_df = self._merge_counts_with_sequences(counts_df)
        genefilt_df = self._filter_by_gene_list(merged_counts_df, vgene_list, jgene_list)
        _write_to_file(genefilt_df, out_file_path)
        if self.compute_pgen:
            pgen_file_path = _get_pgen_file_path(out_file_path)
            command = self._get_olga_command(pgen_file_path, out_file_path)
            _run_olga_tool(command)
            pgen_file = _read_pgen_file(pgen_file_path)
            self._compute_and_combine_pdata(genefilt_df, pgen_file, pgen_file_path, pos_label_examples, neg_label_examples)

    def multi_compute_pgen_public_sequences(self):
        found_files = glob.glob(self.design_mats_path + "/*_design_mat_*", recursive=True)
        pool = Pool(self.n_threads)
        pool.map(self.compute_pgen_public_sequences, found_files)

    def _read_legal_genes(self):
        return tuple(_get_legal_genes(gene_file) for gene_file in [self.legal_vgenes_file, self.legal_jgenes_file])

    def _read_metadata(self):
        return pd.read_csv(self.metadata_file, header=0, index_col=None, error_bad_lines=False)

    def _get_pos_label_examples(self, metadata):
        return metadata[metadata[self.label_field]]["subject_id"].to_list()

    def _get_neg_label_examples(self, metadata):
        return metadata[~metadata[self.label_field]]["subject_id"].to_list()

    def _add_labels(self, counts_df):
        counts_df.columns = [self.pos_label, self.neg_label]
        return counts_df

    def _filter_counts(self, counts_df):
        if self.public_counts_in == 0:
            counts_df = counts_df.loc[counts_df[self.pos_label] > 1]
        elif self.public_counts_in == 1:
            counts_df = counts_df.loc[counts_df[self.neg_label] > 1]
        elif self.public_counts_in == 2:
            counts_df = counts_df[counts_df[[self.pos_label, self.neg_label]].sum(axis=1) > 1]
        else:
            counts_df = counts_df
        return counts_df

    def _merge_counts_with_sequences(self, counts_df):
        merged_counts_df = counts_df.join(self.seq_concat_file_df)
        merged_counts_df = merged_counts_df[['junction_aa', 'v_call', 'j_call', self.pos_label, self.neg_label]]
        return merged_counts_df

    def _filter_by_gene_list(self, merged_counts_df, vgene_list, jgene_list):
        if self.vj_gene_conditioning == 1:
            vfilt_df = merged_counts_df[merged_counts_df['v_call'].isin(vgene_list)]
            return vfilt_df[vfilt_df['j_call'].isin(jgene_list)]
        return merged_counts_df

    def _get_olga_command(self, pgen_file_path, out_file_path):
        if self.vj_gene_conditioning == 1:
            return f'olga-compute_pgen --humanTRB -i {out_file_path} -o {pgen_file_path} --seq_type_out aaseq --v_in 1 --j_in 2'
        return f'olga-compute_pgen --humanTRB -i {out_file_path} -o {pgen_file_path}'

    def _compute_and_combine_pdata(self, genefilt_df, pgen_file, pgen_file_path, pos_label_examples, neg_label_examples):
        pgen_pdata_file_path = os.path.join(os.path.dirname(pgen_file_path),
                                            'pdata_' + os.path.basename(pgen_file_path))
        genefilt_df = genefilt_df.set_index('junction_aa')
        pgen_pdata = pd.concat([genefilt_df, pgen_file], axis=1)
        pgen_pdata['pval_pos'] = pgen_pdata.apply(
            lambda x: _compute_pval(x[self.pos_label], x['pgen'], len(pos_label_examples),
                                    self.avg_rep_size), axis=1)
        pgen_pdata['pval_neg'] = pgen_pdata.apply(
            lambda x: _compute_pval(x[self.neg_label], x['pgen'], len(neg_label_examples),
                                    self.avg_rep_size), axis=1)
        pgen_pdata.to_csv(pgen_pdata_file_path, sep='\t')


def execute():
    concat_seqfile_df = pd.read_csv(args.seq_concat_file, header=0, sep='\t', index_col='sequence_id')
    compute_pgen_wf = ComputePgenPublic(design_mats_path=args.design_mats_path, seq_concat_file_df=concat_seqfile_df,
                                        legal_vgenes_file=args.legal_vgenes_file,
                                        legal_jgenes_file=args.legal_jgenes_file,
                                        n_threads=args.n_threads, metadata_file=args.metadata_file,
                                        label_field=args.label_field, avg_rep_size=args.avg_rep_size,
                                        public_counts_in=args.public_counts_in, compute_pgen=args.compute_pgen,
                                        vj_gene_conditioning=args.vj_gene_conditioning)
    compute_pgen_wf.multi_compute_pgen_public_sequences()

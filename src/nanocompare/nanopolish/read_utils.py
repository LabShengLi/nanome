import pandas as pd
import numpy as np


def read_nanopolish_meth_freq_output(infn):
    """
    Sample input file:

    head K562.methylation_frequency.tsv

    chromosome	start	end	num_motifs_in_group	called_sites	called_sites_methylated	methylated_frequency	group_sequence
    chr1	10524	10524	1	1	1	1.000	TGCTCCGCCTT
    chr1	10541	10541	1	1	1	1.000	ACCACCGAAAT
    chr1	10562	10562	1	1	1	1.000	split-group
    chr1	10570	10570	1	1	1	1.000	split-group
    chr1	10576	10576	1	1	1	1.000	split-group
    chr1	10578	10578	1	1	1	1.000	split-group
    chr1	10588	10588	1	1	1	1.000	split-group
    chr1	10630	10630	1	1	1	1.000	split-group
    chr1	10632	10632	1	1	1	1.000	split-group

    :param infn:
    :return:
    """
    freq_dtypes = {'chromosome'      : np.object,
            'start'                  : np.int64,
            'end'                    : np.int64,
            'num_motifs_in_group'    : np.int64,
            'called_sites'           : np.int64,
            'called_sites_methylated': np.int64,
            'methylated_frequency'   : np.float64,
            'group_sequence'         : np.object}
    df = pd.read_csv(infn, sep="\t", header=0, dtype=freq_dtypes)
    df = df.set_index(['chromosome', 'start', 'end'])
    return df


def read_nanopolish_meth_call_output(infn):
    """
    Nanopolish header are as follows:
    chromosome	strand	start	end	read_name	log_lik_ratio	log_lik_methylated	log_lik_unmethylated	num_calling_strands	num_motifs	sequence
    chr20	+	5000104	5000104	d763e7f0-c28f-413d-aa80-345eeb054f62	-2.81	-87.00	-84.18	1	1	TAATACGTTCA
    chr20	+	5000224	5000224	d763e7f0-c28f-413d-aa80-345eeb054f62	-1.56	-98.39	-96.83	1	1	GGCAACGTCAC
    chr20	+	5000236	5000236	d763e7f0-c28f-413d-aa80-345eeb054f62	-2.78	-115.09	-112.31	1	1	TGCAACGGTGA

    :param infn:
    :return:
    """
    freq_dtypes = {'chromosome'   : np.object,
            'strand'              : np.object,
            'start'               : np.int64,
            'end'                 : np.int64,
            'read_name'           : np.object,
            'log_lik_ratio'       : np.float64,
            'log_lik_methylated'  : np.float64,
            'log_lik_unmethylated': np.float64,
            'num_calling_strands' : np.int64,
            'num_motifs'          : np.int64,
            'sequence'            : np.object}
    df = pd.read_csv(infn, sep="\t", header=0, dtype=freq_dtypes)
    return df


def test_nanop_freq():
    infn = "/projects/liuya/workspace/nanopolish_test/methylation_example/methylation_frequency.tsv"

    infn = '/projects/liuya/workspace/tcgajax/nanocompare/nanopolish/K562.methylation_frequency.tsv'
    df = read_nanopolish_meth_freq_output(infn)

    ind1 = ('chr1', 10524, 10524)
    if ind1 in df.index:
        row_ind1 = df.loc[ind1, :]
        print(row_ind1)

    print(df)
    pass


def test_nanop_call():
    infn = "/projects/liuya/workspace/nanopolish_test/methylation_example/methylation_calls.tsv"
    df = read_nanopolish_meth_call_output(infn)
    print(df)


if __name__ == '__main__':
    test_nanop_call()
    pass

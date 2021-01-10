"""
Load all important data for plotting purpose, includes:

1. Performance data
2. Correlation data
3. Running time and resource usage data

etc.
"""

import pandas as pd

from global_config import *
from nanocompare.collect_data import collect_wide_format_newly_exp
from nanocompare.nanocompare_global_settings import locations_category, map_from_tool_to_abbr, locations_category2, locations_singleton2


def load_running_time_and_mem_usage():
    """
    Load running time and mem usage DF
    :return:
    """
    infn = os.path.join(pkl_base_dir, 'nanocompare', 'nanocompare_running_time.pkl')
    return pd.read_pickle(infn)
    pass


def load_all_perf_data():
    """
    Load the existing get all results, find all Dataset and tagname folder to get all tsv, before refine's data

            2020-01-11 15:47:11,451 - [nanocompare_global_settings.py:92] - INFO:                                                prefix  ...                Dataset
        0   K562_WGBS_joined/K562_WGBS_joined.Nanopolish_c...  ...       K562_WGBS_joined
        1   K562_WGBS_joined/K562_WGBS_joined.Nanopolish_c...  ...       K562_WGBS_joined
        2   K562_WGBS_joined/K562_WGBS_joined.Nanopolish_c...  ...       K562_WGBS_joined
        3   K562_WGBS_joined/K562_WGBS_joined.Nanopolish_c...  ...       K562_WGBS_joined
        6   K562_WGBS_joined/K562_WGBS_joined.Nanopolish_c...  ...       K562_WGBS_joined
        ..                                                ...  ...                    ...
        6   HL60_AML_oxBsseq_cut5/HL60_AML_oxBsseq_cut5.De...  ...  HL60_AML_oxBsseq_cut5
        8   HL60_AML_oxBsseq_cut5/HL60_AML_oxBsseq_cut5.De...  ...  HL60_AML_oxBsseq_cut5
        14  HL60_AML_oxBsseq_cut5/HL60_AML_oxBsseq_cut5.De...  ...  HL60_AML_oxBsseq_cut5
        20  HL60_AML_oxBsseq_cut5/HL60_AML_oxBsseq_cut5.De...  ...  HL60_AML_oxBsseq_cut5
        21  HL60_AML_oxBsseq_cut5/HL60_AML_oxBsseq_cut5.De...  ...  HL60_AML_oxBsseq_cut5
        [261 rows x 20 columns]

        df3.info()
        <class 'pandas.core.frame.DataFrame'>
        Int64Index: 261 entries, 0 to 21
        Data columns (total 20 columns):
        prefix               261 non-null object
        coord                261 non-null object
        accuracy             261 non-null float64
        roc_auc              261 non-null float64
        precision_5C         261 non-null float64
        recall_5C            261 non-null float64
        F1_5C                261 non-null float64
        Csites               261 non-null int64
        precision_5mC        261 non-null float64
        recall_5mC           261 non-null float64
        F1_5mC               261 non-null float64
        mCsites              261 non-null int64
        referenceCpGs        261 non-null int64
        corrMix              212 non-null float64
        Corr_mixedSupport    261 non-null int64
        corrAll              236 non-null float64
        Corr_allSupport      261 non-null int64
        Tool                 261 non-null object
        Location             261 non-null object
        Dataset              261 non-null object
        dtypes: float64(10), int64(5), object(5)
        memory usage: 42.8+ KB

        df3.Dataset.unique()
        Out[5]:
        array(['K562_WGBS_joined', 'K562_RRBS_joined', 'APL_BSseq_cut10',
               'APL_Bsseq', 'APL_oxBSseq_cut5', 'HL60_RRBS_joined',
               'HL60_AML_Bsseq_cut5', 'HL60_AML_oxBsseq_cut5'], dtype=object)

    :return:
    """
    outfn = os.path.join(pkl_base_dir, 'nanocompare', "Nanocompare_performance_alldata.pkl")
    df = pd.read_pickle(outfn)
    return df


def load_all_perf_data_for_dataset(dsname):
    """
    Load a specific dataset's performance original data, such as NA19240_RRBS_joined, APL_BSseq_cut10
    :param dsname:
    :return:
    """
    dfall = load_all_perf_data()
    df1 = dfall[dfall['Dataset'] == dsname]
    df1['Tool'] = df1['Tool'].apply(map_from_tool_to_abbr)
    return df1


def load_all_perf_data_for_dataset_list(dslist=['K562_WGBS_joined', 'APL_BSseq_cut10', 'HL60_AML_Bsseq_cut5', 'NA19240_RRBS_joined']):
    """
    Load a specific dataset's performance original data, such as NA19240_RRBS_joined, APL_BSseq_cut10
    :param dsname:
    :return:
    """
    dfall = load_all_perf_data()
    df1 = dfall[dfall['Dataset'].isin(dslist)]
    df1['Tool'] = df1['Tool'].apply(map_from_tool_to_abbr)
    return df1


def get_one_dsname_perf_data(dsname='K562_WGBS_joined'):
    """
    Get one data sets performance results as the simple case study
    :param dsname:
    :return:
    """
    df = load_refined_data()
    listOfSelected_cat = ["GW", "cpgIslandExt", "promoters_500bp", "exonFeature", "intergenic"]

    measure = 'F1_5mC'
    df1 = df[(df['Dataset'] == dsname) & (df['Measurement'] == measure)]
    df1 = df1[df['Location'].isin(listOfSelected_cat)]
    return df1


# Deprecated
def load_refined_data():
    """
    refined all data for plots

            df1
        Out[3]:
                           Dataset         Location  ... Measurement Performance
        0         K562_WGBS_joined               GW  ...       F1_5C    0.973656
        1         K562_WGBS_joined       singletons  ...       F1_5C    0.972089
        2         K562_WGBS_joined    nonsingletons  ...       F1_5C    0.978164
        3         K562_WGBS_joined     cpgIslandExt  ...       F1_5C    0.981434
        4         K562_WGBS_joined      exonFeature  ...       F1_5C    0.982214
        ..                     ...              ...  ...         ...         ...
        517  HL60_AML_oxBsseq_cut5      exonFeature  ...      F1_5mC    0.573355
        518  HL60_AML_oxBsseq_cut5       intergenic  ...      F1_5mC    0.739625
        519  HL60_AML_oxBsseq_cut5  promoters_500bp  ...      F1_5mC    0.426156
        520  HL60_AML_oxBsseq_cut5       discordant  ...      F1_5mC    0.642202
        521  HL60_AML_oxBsseq_cut5       concordant  ...      F1_5mC    0.723352
        [522 rows x 5 columns]
        df1.Dataset.unique()
        Out[4]:
        array(['K562_WGBS_joined', 'K562_RRBS_joined', 'APL_BSseq_cut10',
               'APL_Bsseq', 'APL_oxBSseq_cut5', 'HL60_RRBS_joined',
               'HL60_AML_Bsseq_cut5', 'HL60_AML_oxBsseq_cut5'], dtype=object)
        df1.Location.unique()
        Out[5]:
        array(['GW', 'singletons', 'nonsingletons', 'cpgIslandExt', 'exonFeature',
               'intergenic', 'promoters_500bp', 'discordant', 'concordant'],
              dtype=object)
        df1.Measurement.unique()
        Out[6]: array(['F1_5C', 'F1_5mC'], dtype=object)

    :return:
    """
    outfn = os.path.join(pkl_base_dir, 'nanocompare', "Nanocompare_performance_refined.pkl")
    df = pd.read_pickle(outfn)
    return df


def load_refined_perf_data_by_measurement(meas_list=['F1_5C', 'F1_5mC']):
    """
    load some performance index of all data, melt performance measurements for subsequent analysis and plotting

    Sample

    df = load_perf_data_refine_by_measurement(['F1_5C', 'F1_5mC', 'accuracy','roc_auc'])
    logger.info(f"df={df}")

    :param meas_list:
    :return:
    """
    alld = load_all_perf_data()

    sel_col = ['Dataset', 'Location', 'Tool'] + meas_list
    refine_alld = alld[sel_col]

    refine_alld = pd.melt(refine_alld, id_vars=['Dataset', 'Location', 'Tool'], var_name='Measurement', value_name='Performance')

    return refine_alld


def load_refined_perf_data_by_measurement2(meas_list=['F1_5C', 'F1_5mC']):
    """
    load some performance index of all data, melt performance measurements for subsequent analysis and plotting

    Sample

    df = load_perf_data_refine_by_measurement(['F1_5C', 'F1_5mC', 'accuracy','roc_auc'])
    logger.info(f"df={df}")

    :param meas_list:
    :return:
    """
    alld = collect_wide_format_newly_exp()

    sel_col = ['Dataset', 'Location', 'Tool'] + meas_list
    refine_alld = alld[sel_col]

    refine_alld = pd.melt(refine_alld, id_vars=['Dataset', 'Location', 'Tool'], var_name='Measurement', value_name='Performance')

    return refine_alld


def get_performance_from_datasets_and_locations(ds_list=['K562_WGBS_joined', 'APL_BSseq_cut10', 'HL60_AML_Bsseq_cut5'], location_list=locations_category, meas_list=['F1_5C', 'F1_5mC']):
    """
    Get a list of datasets results with refined, (melt to Performance and Measurement)
    :param ds_list: a list of datasets
    :return:
    """
    df = load_refined_perf_data_by_measurement(meas_list=meas_list)

    dfsel = df[df['Dataset'].isin(ds_list)]
    dfsel = dfsel[dfsel['Location'].isin(location_list)]

    dfsel['Tool'] = dfsel['Tool'].apply(map_from_tool_to_abbr)

    return dfsel


def load_box_plot_all_data():
    """
    data from report dump folder

    df2
    :return:
    """
    infn = os.path.join(pkl_base_dir, 'nanocompare', 'box_plots_all_data.joinedReports.pkl')
    df = pd.read_pickle(infn)

    df['BasecallTool'] = df['BasecallTool'].apply(map_from_tool_to_abbr)

    return df


def load_sing_nonsing_count_df():
    infn = os.path.join(pkl_base_dir, "nanocompare", "singleton_nonsingleton_count_with_locations.csv")
    df = pd.read_csv(infn, header=0, index_col=0)
    logger.debug(f"df={df}")
    return df


def save_pivoted_sing_nonsing_count_df():
    df = load_sing_nonsing_count_df()

    df = df.drop('ds_source', axis=1)

    df1 = df.pivot(index='dsname', columns='location')

    df1 = df1.iloc[:, 0:3]
    df1.columns = df1.columns.droplevel()

    df1 = df1.rename(columns={'absolute': 'Singletons', 'concordant': 'Concordant', 'discordant': 'Discordant'})
    df1['Non-Singletons'] = df1['Concordant'] + df1['Discordant']

    df1 = df1[['Singletons', 'Non-Singletons', 'Concordant', 'Discordant']]

    fout = os.path.join(pic_base_dir, 'dataset.singleton.nonsingleton.counts.xlsx')
    df1.to_excel(fout)


def load_long_format_perf_data_from_newly_exp_by_measure(meas_list=['F1_5C', 'F1_5mC']):
    """
    load some performance index of all data, melt performance measurements for subsequent analysis and plotting into long format

    Sample

    df = load_perf_data_refine_by_measurement(['F1_5C', 'F1_5mC', 'accuracy','roc_auc'])
    logger.info(f"df={df}")

    Dataframe as long format:
    'Dataset', 'Location', 'Tool', 'Measurement', 'Performance'

    :param meas_list:
    :return:
    """
    wide_df = collect_wide_format_newly_exp()  # wide format

    sel_col = ['Dataset', 'Location', 'Tool'] + meas_list
    wide_df = wide_df[sel_col]

    long_df = pd.melt(wide_df, id_vars=['Dataset', 'Location', 'Tool'], var_name='Measurement', value_name='Performance')

    return long_df


def get_long_format_perf_within_measures_and_locations(location_list=locations_category2, meas_list=['F1_5C', 'F1_5mC']):
    """
    Get a list of datasets results with refined, (melt to Performance and Measurement)
    :param ds_list: a list of datasets
    :return:
    """
    df = load_long_format_perf_data_from_newly_exp_by_measure(meas_list=meas_list)

    dfsel = df[df['Location'].isin(location_list)]

    # dfsel['Tool'] = dfsel['Tool'].apply(map_from_tool_to_abbr)

    return dfsel


def load_corr_data_tsv_fns():
    """
    Load the full file and path names of correlation on four tools tsv results
    :return:
    """
    # nanocompare_results_dir = "/projects/li-lab/yang/NanoCompare-paper/results/corr-plot-data"

    base_dir = os.path.join(data_base_dir, 'corr-plot-data')

    files = ["Methylation_correlation_plotting_data.APL_oxBS_cut10.tsv", "Methylation_correlation_plotting_data.APL_WGBS_cut10.tsv", "Methylation_correlation_plotting_data.HL60_RRBS_rep_ENCFF000MDA_Bismark.tsv", "Methylation_correlation_plotting_data.HL60_RRBS_rep_ENCFF000MDF_Bismark.tsv",
            "Methylation_correlation_plotting_data.K562_WGBS_rep_ENCFF721JMB.tsv",
            "Methylation_correlation_plotting_data.K562_WGBS_rep_ENCFF867JRG.tsv"]

    files = ["Methylation_correlation_plotting_data.APL_WGBS_cut10.tsv", "Methylation_correlation_plotting_data.HL60_RRBS_rep_ENCFF000MDF_Bismark.tsv",
            "Methylation_correlation_plotting_data.K562_WGBS_rep_ENCFF867JRG.tsv"]

    filesFull = [os.path.join(base_dir, fn) for fn in files]

    na19240files = ['Methylation_correlation_plotting_data.NA19240_RRBS_Origin_time_01-17-23-27-10-489023.tsv',
            'Methylation_correlation_plotting_data.NA19240_RRBS_ENCFF000LZS_time_01-18-01-55-06-136184.tsv',
            'Methylation_correlation_plotting_data.NA19240_RRBS_ENCFF000LZT_time_01-17-23-27-21-284827.tsv',
            ]

    na19240files = ['Methylation_correlation_plotting_data.NA19240_RRBS_Origin_time_01-17-23-27-10-489023.tsv'
            ]
    na19240List = [os.path.join(base_dir, fn) for fn in na19240files]

    ret_list = ['/projects/li-lab/yang/results/2020-12-21/K562_WGBS_Joined/Meth_corr_plot_data-K562_WGBS_Joined-bsCov5-minCov4-time-2020-12-21-15-57.tsv']

    ret_list = ['/projects/li-lab/yang/results/2021-01-08/K562_WGBS_Joined_NewRuns/Meth_corr_plot_data-K562_WGBS_Joined_NewRuns-bsCov10-minCov4-baseCount0.tsv']

    # return filesFull

    return ret_list


if __name__ == '__main__':

    pass


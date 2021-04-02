"""
Plots all Nanocompare paper out
"""
import argparse
import glob
import pickle
from collections import defaultdict

import matplotlib.patches as mpatches
import numpy as np
import seaborn as sns
# from ggplot import ggplot, aes, geom_line, geom_abline
from scipy.stats import pearsonr
from sklearn import metrics
from sklearn.metrics import precision_recall_curve

from nanocompare.collect_data import collect_singleton_vs_nonsingleton_df, save_wide_format_performance_results
from nanocompare.global_settings import *
from nanocompare.load_data import *


def single_ds_5mc_5c_performance(dsname="HL60_AML_Bsseq_cut5"):
    df_sel = get_one_dsname_perf_data(dsname)

    style = 'Location'
    filled_markers = ('o', 'P', '>', 's', 'X')
    filled_markers = ['o', 'P', '>', 's', 'X']

    ax = sns.scatterplot(x="Location", y="Performance", hue='Tool', data=df_sel, style=style, markers=filled_markers, s=100, alpha=0.85, linewidths=None)

    # ax = sns.stripplot(x="Location", y="Performance", hue='Tool', size=8, data=df_sel)

    # sns.catplot(x="time", y="pulse", hue="kind", data=exercise)

    # ax = sns.catplot(x="Location", y="Performance", hue='Tool', data=df_sel)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    outfn = os.path.join(pic_base_dir, f"singleds_nanocompare_location_f1_5c_5mc_performance_time_{current_time_str()}.png")
    ax.figure.savefig(outfn, format='png', dpi=600)
    logger.info(f"save to {outfn}")
    plt.show()


def set_style(font_scale=1.2):
    # sns.set()
    # sns.set(font_scale=1.15)
    # sns.set_style("white")

    # This sets reasonable defaults for font size for
    # a figure that will go in a paper
    sns.set_context("paper")
    sns.set(font_scale=font_scale)

    # Set the font to be serif, rather than sans
    # sns.set(font='serif')

    # Make the background white, and specify the
    # specific font family
    # sns.set_style("white", {
    #         "font.family": "serif",
    #         "font.serif" : ["Times", "Palatino", "serif"]
    #         })
    sns.set_style("white")


def cut_abbr_dsname(dsname):
    """
    Cut the abbr from dsname with other useful info
    :param dsname:
    :return:
    """
    return dsname[0:dsname.find('_')]


def set_labels_fig_3a(grid, ds_list, location_list, meas_list=['F1_5mC', 'F1_5C']):
    """
    Set labels of figure 3A, FacetGrid
    :param grid:
    :return:
    """

    ncol = len(ds_list)
    [plt.setp(ax.texts, text="") for ax in grid.axes.flat]
    # grid.set_titles(row_template='{row_name}', col_template='{col_name}')
    # grid.set_xticklabels(rotation=45)

    grid.set_titles(row_template='', col_template='')
    # grid.set_xticklabels([])
    # grid.set(xticks=[])

    # Iterate through each axis
    for axi, ax in enumerate(grid.axes.flat):

        # only set first col figure's Y axis label
        if axi % ncol == 0:
            ax.set_ylabel(meas_list[axi // ncol])
        else:
            ax.set_ylabel("")

        if axi < ncol:
            ax.set_title(ds_list[axi])

        ax.set_xticklabels([])

    grid.fig.tight_layout()

    legend_data = grid._legend_data

    # modify legend label to another one
    # legend_data["Tool"] = legend_data.pop("Tool")

    # keep this label order to print dict of legend_data
    # label_order = ['Tool'] + tools_abbr + ['Location'] + location_list

    grid.add_legend(label_order=legend_data)
    # plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5))
    pass


def scatter_facet_grid_multiple_ds_5mc_5c_performance(df, location_list=locations_singleton, meas_list=['F1_5mC', 'F1_5C']):
    """
    Generate three datasets 5mc and 5C results in a FacetGrid
    :return:
    """

    # Set the style for plotting
    set_style(font_scale=1.5)

    # style = 'Location'

    # figure specific settings

    if len(location_list) == 4:
        filled_markers = ('o', 'P', '>', 's')
    else:
        filled_markers = ('o', 'P', '>', 's', 'X', '^')

    marker_size = 120

    perf_order = meas_list

    tools_showing_order = ['DeepSignal', 'Tombo', 'Nanopolish', 'DeepMod', 'Megalodon']

    grid = sns.FacetGrid(df, row='Measurement', col='Dataset', row_order=perf_order, margin_titles=True)

    # grid.map_dataframe(sns.catplot, x="Location", y="Performance", hue='Tool', hue_order=tools_abbr)

    grid.map_dataframe(sns.scatterplot, x="Location", y="Performance", hue='Tool', style="Location", markers=filled_markers, hue_order=tools_showing_order, s=marker_size, x_jitter=3, y_jitter=5)

    ds_list = df['Dataset'].unique()
    # Set plot labels
    set_labels_fig_3a(grid, ds_list=ds_list, location_list=location_list, meas_list=meas_list)

    # Save plot and show if possible
    outfn = os.path.join(pic_base_dir, f"scatter_facet_grid_multiple_ds_location_{location_list[0]}_{meas_list[0]}_performance.png")

    # grid.savefig(outfn, format='png', bbox_inches='tight', dpi=600)

    plt.savefig(outfn, format='png', bbox_inches='tight', dpi=600)

    plt.show()
    plt.close()
    logger.info(f"save to {outfn}")


def corr_grid_plot_for_fig5a(infn):
    """
    Plot the grid of corr COE, distribution and scatter plot based on input files
    :return:
    """

    ## Load data into df
    logger.debug(infn)
    df = pd.read_csv(infn, sep=',')

    sel_col = []
    rename_dict = {}
    for col in df.columns:
        if str(col).endswith('_freq'):
            sel_col.append(col)

            new_col_name = get_tool_name(str(col).replace('_freq', ''))
            rename_dict.update({str(col): new_col_name})
    df = df[sel_col]
    df = df.rename(columns=rename_dict)
    num_col = len(df.columns)
    df = df.iloc[:, list(range(1, num_col)) + [0]]
    logger.debug(f'Load data = {len(df)}')

    # Rearrange columns by correlation from high to low
    corrdf = df.corr()
    # logger.info(corrdf)

    coe_series = corrdf.iloc[:, -1]
    # logger.info(coe_series[:-1].sort_values(ascending=False))

    orderedColumns = coe_series[:-1].sort_values(ascending=False).index.tolist()

    df = df.loc[:, ['BGTruth'] + orderedColumns]

    ## Plot correlation grid figure
    basefn = os.path.basename(infn)
    outfileName = "{}.jpg".format(basefn.replace(".csv", ""))
    outfn = os.path.join(args.o, outfileName)

    plt.clf()

    fig, ax = plt.subplots()

    fig.set_size_inches(10, 10)

    left, width = .25, .5
    bottom, height = .25, .5
    right = left + width
    top = bottom + height

    position = 1
    for y in range(1, num_col + 1):
        for x in range(1, num_col + 1):
            if x == y:
                # Diagonal, distribution lines:
                plt.subplot(num_col, num_col, position)
                #             df[fields[x-1]].hist(bins=100)
                params = {'legend.fontsize': 16, 'legend.handlelength': 0, 'legend.handletextpad': 0, 'legend.fancybox': True}
                plt.rcParams.update(params)
                ax = sns.kdeplot(df.iloc[:, x - 1], shade=True, color="black", legend=True)
                leg = ax.legend(labels=[str(df.columns[x - 1])], fontsize=14)
                for item in leg.legendHandles:
                    item.set_visible(False)

                # ax.annotate("abcd", xy=(0.2, 0.8), xycoords='axes fraction',
                #             fontsize=14)
                #
                # ax.xaxis.set_major_locator(plt.MaxNLocator(3))
                # ax.yaxis.set_major_locator(plt.MaxNLocator(3))
                ax.set_yticklabels([])
                ax.set_xticklabels([])
                ax.set_ylabel('')
                ax.set_xlabel('')

            elif x > y:
                # upper triangle, COE number:
                ax2 = plt.subplot(num_col, num_col, position)

                corrValue = pearsonr(df.iloc[:, x - 1], df.iloc[:, y - 1])
                corrValueStr = "{0:2.2f}".format(corrValue[0])

                ax2.text(0.5 * (left + right), 0.5 * (bottom + top), corrValueStr,
                         horizontalalignment='center',
                         verticalalignment='center',
                         fontsize=25, color='black',
                         transform=ax2.transAxes)

                ax2.set_yticklabels([])
                ax2.set_xticklabels([])

            elif x < y:
                # lower triangle, scatter plot, using hexbin, x-column label, y-row label:
                ax3 = plt.subplot(num_col, num_col, position)

                # plt.hexbin(df.iloc[:, x - 1], df.iloc[:, y - 1], gridsize=(gridRes, gridRes), cmap=agressiveHot)  # plt.cm.gray_r )
                if len(df) < 100:
                    gridRes = 15  # for HL60
                else:
                    gridRes = 40

                mincnt = None
                plt.hexbin(df.iloc[:, x - 1], df.iloc[:, y - 1], gridsize=(gridRes, gridRes), cmap='Blues', bins='log', mincnt=mincnt)

                # ax3.set_yticklabels([])
                # ax3.set_xticklabels([])
                plt.xticks([0, 1], fontsize=8)
                plt.yticks([0, 1], fontsize=8)

            position += 1

    fig.savefig(outfn, dpi=600, bbox_inches='tight')
    plt.show()
    plt.close()
    logger.info(f"save to {outfn}")


def smooth_scatter_cor_plot():
    # files = ["Methylation_correlation_plotting_data.APL_oxBS_cut10.tsv", "Methylation_correlation_plotting_data.APL_WGBS_cut10.tsv", "Methylation_correlation_plotting_data.HL60_RRBS_rep_ENCFF000MDA_Bismark.tsv", "Methylation_correlation_plotting_data.HL60_RRBS_rep_ENCFF000MDF_Bismark.tsv",
    #         "Methylation_correlation_plotting_data.K562_WGBS_rep_ENCFF721JMB.tsv",
    #         "Methylation_correlation_plotting_data.K562_WGBS_rep_ENCFF867JRG.tsv"]
    # cor_tsv_fields = ["DeepSignal_freq", "Tombo_freq", "Nanopolish_freq", "DeepMod_freq", "DeepMod_clust_freq", "BSseq"]

    # fields = ["DeepSignal_freq", "Tombo_freq", "Nanopolish_freq", "DeepMod_freq",  "BSseq"]

    optionalSuffix = current_time_str()

    for infn in load_corr_data_tsv_fns():
        df = pd.read_csv(infn, sep='\t')
        df = df.rename(columns=dict_cor_tsv_to_abbr())

        basename = os.path.basename(infn)
        outfileName = "{}_time_{}.png".format(basename.replace(".tsv", ""), optionalSuffix)
        outfn = os.path.join(pic_base_dir, outfileName)
        position = 1

        logger.debug(f"df={df}")
        df.info()

        xdata = df['BSseq']
        ydata = df['DeepSignal']

        ydata = df['Tombo']

        # smooth_scatter_call_r(x=xdata, y=ydata, x_label='BGTruth', y_label=f"DeepSignal", outdir=pic_base_dir, is_show=True)

        scatter_plot_x_y_smoothlike(xdata=xdata, ydata=ydata)

        break

        plt.clf()

        fig, ax = plt.subplots()
        fig.set_size_inches(10, 10)

        left, width = .25, .5
        bottom, height = .25, .5
        right = left + width
        top = bottom + height

        gridRes = 30

        for y in range(1, len(cor_tsv_fields) + 1):
            for x in range(1, len(cor_tsv_fields) + 1):
                if x == y:
                    # Diagonal:
                    plt.subplot(len(cor_tsv_fields), len(cor_tsv_fields), position)
                    #             df[fields[x-1]].hist(bins=100)
                    ax = sns.kdeplot(df[cor_tsv_fields_abbr[x - 1]], shade=True, color="black")

                    ax.set_yticklabels([])
                    ax.set_xticklabels([])

                elif x > y:
                    # upper triangle:
                    ax2 = plt.subplot(len(cor_tsv_fields), len(cor_tsv_fields), position)

                    corrValue = pearsonr(df[cor_tsv_fields_abbr[x - 1]], df[cor_tsv_fields_abbr[y - 1]])
                    corrValueStr = "{0:2.2f}".format(corrValue[0])

                    ax2.text(0.5 * (left + right), 0.5 * (bottom + top), corrValueStr,
                             horizontalalignment='center',
                             verticalalignment='center',
                             fontsize=25, color='black',
                             transform=ax2.transAxes)

                    ax2.set_yticklabels([])
                    ax2.set_xticklabels([])

                elif x < y:
                    # lower triangle:
                    ax3 = plt.subplot(len(cor_tsv_fields), len(cor_tsv_fields), position)

                    plt.hexbin(df[cor_tsv_fields_abbr[x - 1]], df[cor_tsv_fields_abbr[y - 1]], gridsize=(gridRes, gridRes), cmap=agressiveHot)  # plt.cm.gray_r )

                    ax3.set_yticklabels([])
                    ax3.set_xticklabels([])

                position += 1

        fig.savefig(outfn, dpi=600, bbox_inches='tight')
        plt.show()
        logger.info(f"save to {outfn}")
        # break


def get_simple_title_from_ds(str):
    return str[0:str.find('_')]
    pass


def cor_box_plot():
    """
    Box plot the correlation of all or mixed CpGs
    :return:
    """

    folders = ["AML_Bsseq_cut10", "AML_Bsseq_cut5", "AML_oxBSseq_cut10", "AML_oxBSseq_cut5", "APL_BSseq_cut10", "APL_BSseq_cut5", "APL_oxBSseq_cut10", "APL_oxBSseq_cut5", "HL60_AML_Bsseq_cut5", "HL60_AML_oxBsseq_cut5", "HL60_RRBS_rep_ENCFF000MDA_Bismark", "HL60_RRBS_rep_ENCFF000MDF_Bismark", "HL60_RRBS_rep_ENCFF001TNE", "HL60_RRBS_rep_ENCFF001TNF",
            "K562_RRBS_rep_ENCFF001TOL", "K562_RRBS_rep_ENCFF001TOM", "K562_WGBS_rep_ENCFF721JMB", "K562_WGBS_rep_ENCFF721JMB_chr20", "K562_WGBS_rep_ENCFF721JMB_chr20_except", "K562_WGBS_rep_ENCFF721JMB_chr20_except_2", "K562_WGBS_rep_ENCFF867JRG"]
    folders = [
            "APL_BSseq_cut10"]  # , "APL_BSseq_cut5", "APL_oxBSseq_cut10", "APL_oxBSseq_cut5", "HL60_AML_Bsseq_cut5", "HL60_AML_oxBsseq_cut5", "HL60_RRBS_rep_ENCFF000MDA_Bismark", "HL60_RRBS_rep_ENCFF000MDF_Bismark", "HL60_RRBS_rep_ENCFF001TNE", "HL60_RRBS_rep_ENCFF001TNF", "K562_RRBS_rep_ENCFF001TOL", "K562_RRBS_rep_ENCFF001TOM", "K562_WGBS_rep_ENCFF721JMB", "K562_WGBS_rep_ENCFF867JRG"]

    # 'K562_WGBS_joined', 'APL_BSseq_cut10', 'HL60_AML_Bsseq_cut5', 'NA19240_RRBS_joined'
    folders = ["NA19240_RRBS_joined", 'K562_WGBS_joined', 'APL_BSseq_cut10', 'HL60_AML_Bsseq_cut5']

    stat = "corrAll"

    for folder in folders:
        print(folder)
        df = load_all_perf_data_for_dataset(folder)

        # indir = os.path.join(nanocompare_basedir, 'reports', folder)
        # data = pp.load_data("/home/rosikw/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/reports/{}".format(folder))
        # df = pp.load_data(indir)

        logger.debug(f"df={df}")

        # data.to_csv("{}.combined.tsv".format(folder), sep='\t', index=False)

        #######

        plt.clf()

        # df = pd.read_csv("{}.combined.tsv".format(folder), delimiter="\t")
        fig = plt.gcf()
        fig.set_size_inches(4, 4)
        plt.rcParams.update({'font.size': 14})
        # sns.boxplot(x="method", y="AUC", hue="BinaryLabels", data=df, palette="Set1")
        # ax = sns.boxplot(x="BasecallTool", y=stat, data=df, palette="Set1", order=["Tombo_calls", "Nanopolish_calls", "DeepSignal_calls", "DeepMod_calls"])

        ax = sns.boxplot(x="Tool", y=stat, data=df, order=tools_abbr)

        # ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
        ax.set_xticklabels([], rotation=45)

        plt.ylim(0, 1)

        title_str = get_simple_title_from_ds(folder)
        plt.title(title_str)

        outfn = os.path.join(pic_base_dir, f"{folder}.{stat}.time.{current_time_str()}.BoxPlot.png")
        plt.savefig(outfn, bbox_inches='tight', dpi=600, format="png")
        logger.info(f"save to {outfn}")
        plt.show()
        plt.close()

        # break

    stat = "corrMix"

    for folder in folders:
        print(folder)
        df = load_all_perf_data_for_dataset(folder)

        plt.clf()

        fig = plt.gcf()
        fig.set_size_inches(4, 4)
        plt.rcParams.update({'font.size': 14})

        ax = sns.boxplot(x="Tool", y=stat, data=df, order=tools_abbr)
        # ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
        ax.set_xticklabels([], rotation=45)

        title_str = get_simple_title_from_ds(folder)

        plt.ylim(0, 1)
        plt.title(title_str)

        current_palette = sns.color_palette()
        pal_plot = [mpatches.Patch(color=current_palette[i], label=tools_abbr[i]) for i in range(4)]
        plt.legend(handles=pal_plot, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title='Tool')

        outfn = os.path.join(pic_base_dir, f"{folder}.{stat}.time.{current_time_str()}.BoxPlot.png")
        plt.savefig(outfn, bbox_inches='tight', dpi=600, format="png")
        logger.info(f"save to {outfn}")
        plt.show()
        plt.close()
        # break

    pass


def plot_runnnig_time_and_mem_usage():
    df = load_running_time_and_mem_usage()

    pal_plot = get_pal_4()

    set_style(font_scale=1.5)
    sns.set_style('ticks')
    # fig, ax = plt.subplots()
    fig = plt.gcf()
    fig.set_size_inches(15, 4)
    # plt.rcParams.update({'font.size': 18})

    # Subplot 1
    plt.subplot(1, 2, 1)

    ax1 = sns.barplot(x="mem", y="tool", data=df, order=tools_abbr)

    ax1.set(xlabel='RAM memory used (GB)', ylabel='')
    sns.despine()
    plt.tight_layout()

    # Subplot 2
    plt.subplot(1, 2, 2)
    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25, wspace=0.35)
    ax2 = sns.barplot(x="time", y="tool", data=df, order=tools_abbr)

    plt.locator_params(axis='x', nbins=4)
    plt.ylabel("")
    plt.xlabel("Total CPU time consumed (seconds)")
    sns.despine()

    plt.legend(handles=pal_plot, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title='Tool')

    plt.tight_layout()

    outfn = os.path.join(pic_base_dir, f"nanocompare_resources_usage_{current_time_str()}.png")
    plt.savefig(outfn, dpi=600, bbox_inches='tight', format="png")
    logger.info(f"save to {outfn}")

    plt.show()
    plt.close()


def gen_figure_3a_4a():
    """
    Generate Figure 3-a and 4-a in paper
    Scatter plot on FacetGrid
    :return:
    """
    measure_list = [['F1_5mC', 'F1_5C']]

    # measure_list = [['F1_5mC', 'F1_5C'],
    #         ['Accuracy', 'ROC_AUC'],
    #         ['Precision_5mC', 'Recall_5mC'],
    #         ['Precision_5C', 'Recall_5C']]

    for meas_list in measure_list:  # meas_list = ['precision_5mC', 'precision_5C']
        df = get_long_format_perf_within_measures_and_locations(location_list=locations_singleton, meas_list=meas_list)
        scatter_facet_grid_multiple_ds_5mc_5c_performance(df, location_list=locations_singleton, meas_list=meas_list)

        df = get_long_format_perf_within_measures_and_locations(location_list=locations_category, meas_list=meas_list)
        scatter_facet_grid_multiple_ds_5mc_5c_performance(df, location_list=locations_category, meas_list=meas_list)


def gen_figure_3b_4b():
    """
    Generate figure 3b and 4b in paper
    Box plot of all results
    :return:
    """

    metrics = ["F1_5mC", "F1_5C", "Accuracy", "ROC_AUC", 'Precision_5mC', 'Precision_5C', 'Recall_5mC', 'Recall_5C']

    # metrics = ["F1_5mC", "F1_5C"]
    # metrics = ["Accuracy", "ROC_AUC"]

    box_plots_two_locations(metrics=metrics)


def gen_figure_5a(infn):
    """
    Generate figure 5-a in paper
    Correlation plots
    :return:
    """
    corr_grid_plot_for_fig5a(infn)
    pass


def gen_figure_5b():
    """
    Correlation plots
    :return:
    """
    cor_box_plot()
    pass


def gen_figure_5c():
    """
    Bar plots of running time and mem usage
    :return:
    """
    plot_runnnig_time_and_mem_usage()


def pie_plot_for_ds(dsname="NA19240"):
    df = load_sing_nonsing_count_df()
    #
    df1 = df[df['dsname'] == dsname]
    #
    logger.debug(f"df1={df1}")
    #
    nsing = df1[df1['location'].isin(['absolute'])]['count'].sum()
    ncond = df1[df1['location'] == 'concordant']['count'].sum()
    ndisc = df1[df1['location'] == 'discordant']['count'].sum()

    # df = load_count_ds_original()
    #
    # nsing = df.loc[dsname, ('singletons', 'total.sites')]
    # ncond = df.loc[dsname, ('concordant', 'total.sites')]
    # ndisc = df.loc[dsname, ('discordant', 'total.sites')]

    piedata = [nsing, ncond, ndisc]
    pielabel = ['Singletons', "Concordant", "Discordant"]

    # plot_pie(piedata, pielabel)

    plot_pie_chart(piedata, pielabel, dsname=dsname)

    pass


def pie_plot_all():
    dslist = ['NA19240', 'K562', "HL60", 'APL']

    for dsname in dslist:
        pie_plot_for_ds(dsname=dsname)


# def plot_pie(data, labels):
#     # explode = (0, 0, 0, 0.1, 0, 0)
#     plt.pie(data, labels=labels, autopct='%1.1f%%', shadow=False, startangle=150)
#     plt.title("Pie chart")
#     plt.show()


def plot_pie_chart(data, labels, dsname="NA19240"):
    set_style(font_scale=1.0)
    fig, ax = plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))

    percent = np.array(data) / np.sum(data) * 100

    logger.debug(f"percent={percent}")

    labels_legend = [f"{labels[k]}: {percent[k]:.1f}%" for k in range(len(labels))]

    # wedges, texts = ax.pie(data)
    # wedges, texts = ax.pie(data, wedgeprops=dict(width=0.5), startangle=-40)

    explode = (0, 0.1, 0.1)  # only "explode" the 2nd slice (i.e. 'Hogs')

    wedges, texts = ax.pie(data, explode=explode, labels=data, startangle=90)

    if False:
        bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
        kw = dict(arrowprops=dict(arrowstyle="-"),
                  bbox=bbox_props, zorder=0, va="center")

        for i, p in enumerate(wedges):
            ang = (p.theta2 - p.theta1) / 2. + p.theta1
            y = np.sin(np.deg2rad(ang))
            x = np.cos(np.deg2rad(ang))
            horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
            connectionstyle = "angle,angleA=0,angleB={}".format(ang)
            kw["arrowprops"].update({"connectionstyle": connectionstyle})
            ax.annotate(f"{percent[i]:.1f}%", xy=(x, y), xytext=(1.35 * np.sign(x), 1.4 * y),
                        horizontalalignment=horizontalalignment, **kw)

    # ax.set_title("Matplotlib bakery: A donut")
    ax.legend(wedges, labels_legend,
              title="",
              loc="center left",
              bbox_to_anchor=(1, 0, 0.5, 1))

    ax.set_title(f"{dsname}")

    # plt.legend(wedges, labels_legend, loc="lower left")

    ax.figure.tight_layout()
    outfn = os.path.join(pic_base_dir, f"pie_plot_{dsname}_time_{current_time_str()}.png")
    ax.figure.savefig(outfn, format='png', dpi=600)

    logger.info(f"save to {outfn}")

    plt.show()


def parse_arguments():
    parser = argparse.ArgumentParser(description='Plot out in Nano-compare paper.')
    parser.add_argument("cmd", help="name of command, lung or lesion")
    parser.add_argument('-i', nargs='+', help='list of input files', default=[])
    parser.add_argument('-o', type=str, help="output dir", default=pic_base_dir)
    parser.add_argument('--tagname', type=str, help="tagname of files", default=None)

    args = parser.parse_args()
    return args


def plot_performance_curves(ret, outdir, tagname="tagname"):
    figure_size = (4, 4)

    title_font_size = 18
    label_font_size = 16

    plt.clf()
    plt.figure(figsize=figure_size)
    for toolname, toolcolor in zip(ToolNameList, ToolsColorList):
        ytrue = ret[f'{toolname}_true']
        ypred = ret[f'{toolname}_pred']
        yscore = ret[f'{toolname}_score']

        # logger.debug(f'ytrue={len(ytrue)}, ypred={len(ypred)}, yscore={len(yscore)}')
        fpr, tpr, threshold = metrics.roc_curve(ytrue, yscore)
        roc_auc = metrics.auc(fpr, tpr)
        plt.plot(fpr, tpr, toolcolor, label=f'{toolname}={roc_auc:.2f}')
    plt.legend(loc='lower right')
    plt.plot([0, 1], [0, 1], 'r--')
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.ylabel('True Positive Rate', fontsize=label_font_size)
    plt.xlabel('False Positive Rate', fontsize=label_font_size)
    plt.title('ROC Curves', fontsize=title_font_size)
    outfn = os.path.join(outdir, f'{tagname}.roc.curves.jpg')
    plt.savefig(outfn, format='png', bbox_inches='tight', dpi=600)

    plt.show()
    plt.close()

    ## PR curves
    plt.clf()
    plt.figure(figsize=figure_size)
    for toolname, toolcolor in zip(ToolNameList, ToolsColorList):
        ytrue = ret[f'{toolname}_true']
        yscore = ret[f'{toolname}_score']

        precision, recall, _ = precision_recall_curve(ytrue, yscore)
        average_precision = metrics.average_precision_score(ytrue, yscore)
        plt.plot(recall, precision, toolcolor, label=f'{toolname}={average_precision:.2f}')
    plt.legend(loc='lower left')
    plt.xlim([0, 1])
    plt.xlabel('Recall', fontsize=label_font_size)
    plt.ylabel('Precision', fontsize=label_font_size)
    plt.ylim([0.0, 1.05])
    plt.title('Precision-Recall Curves', fontsize=title_font_size)

    outfn = os.path.join(outdir, f'{tagname}.pr.curves.jpg')
    plt.savefig(outfn, format='png', bbox_inches='tight', dpi=600)

    plt.show()
    plt.close()

    pass


if __name__ == '__main__':
    set_log_debug_level()

    args = parse_arguments()
    logger.debug(args)

    if args.cmd == 'fig5a':
        ## find ${resultDir} -name 'Meth_corr_plot_data*.csv' \
        ##      -type f -exec python plot_figure.py fig5a -i {} \;
        ## python plot_figure.py fig5a -i /projects/li-lab/Nanopore_compare/result/MethCorr-NA19240_RRBS/Meth_corr_plot_data_joined-NA19240_RRBS-bsCov5-minToolCov3-baseFormat1.csv
        for fn in args.i:
            gen_figure_5a(fn)
    elif args.cmd == 'fig5c':
        gen_figure_5c()
    elif args.cmd == 'fig5b':
        gen_figure_5b()
    elif args.cmd == 'fig3a':
        ## python plot_figure.py fig3a
        gen_figure_3a_4a()
    elif args.cmd == 'fig3b':
        gen_figure_3b_4b()
    elif args.cmd == 'all':
        gen_figure_3a_4a()
        gen_figure_3b_4b()
        gen_figure_5a()
    elif args.cmd == 'export-data':
        # using global setting variable runPrefixDict to locate all data
        ## python plot_figure.py export-data -i /projects/li-lab/Nanopore_compare/result/MethPerf-APL_RRBS_CPG /projects/li-lab/Nanopore_compare/result/MethPerf-HL60_RRBS_CPG /projects/li-lab/Nanopore_compare/result/MethPerf-K562_WGBS_CPG --tagname CPG

        ## python plot_figure.py export-data -i /projects/li-lab/yang/results/2021-02-22/Meth
        if len(args.i) == 0:
            run_prefix = runPrefixDict
        else:
            run_prefix = defaultdict()
            for dirname in args.i:
                run_prefix[os.path.basename(dirname)] = dirname
        save_wide_format_performance_results(run_prefix, args.o, args.tagname)

        pattern1 = "*.summary.singleton.nonsingleton.cov1.csv"
        df = collect_singleton_vs_nonsingleton_df(run_prefix, pattern=pattern1)
        outfn = os.path.join(args.o, f'dataset.singleton.vs.non-singleton{f"-{args.tagname}" if args.tagname else ""}.cov1.csv')
        df.to_csv(outfn)

        pattern2 = "*.summary.singleton.nonsingleton.cov5.csv"
        df = collect_singleton_vs_nonsingleton_df(run_prefix, pattern=pattern2)
        outfn = os.path.join(args.o, f'dataset.singleton.vs.non-singleton{f"-{args.tagname}" if args.tagname else ""}.cov5.csv')
        df.to_csv(outfn)

        logger.info(f'save stats of singleton and non-singleton to {outfn}')
    elif args.cmd == 'export-curve-data':
        ## python plot_figure.py export-curve-data -i /projects/li-lab/Nanopore_compare/result/MethPerf-HL60_RRBS /projects/li-lab/Nanopore_compare/result/MethPerf-K562_WGBS

        ## python plot_figure.py export-curve-data -i /projects/li-lab/Nanopore_compare/result/MethPerf-APL_RRBS_CPG /projects/li-lab/Nanopore_compare/result/MethPerf-HL60_RRBS_CPG /projects/li-lab/Nanopore_compare/result/MethPerf-K562_WGBS_CPG --tagname CPG
        outdir = os.path.join(args.o, f'plot-curve-data{f"-{args.tagname}" if args.tagname else ""}')
        os.makedirs(outdir, exist_ok=True)

        if len(args.i) == 0:
            args.i = list(runPrefixDict.values())
            logger.info(f'No input, we use default args.i={args.i}')

        for bdir in args.i:
            runPrefix = os.path.basename(bdir)
            cnt = 0
            for coordinate_bed_name in coordDictToStandardName.keys():
                curve_data = defaultdict(list)
                for toolname in ToolNameList:
                    pattern_str = os.path.join(bdir, 'performance?results', 'curve_data', f'*.{toolname}.*{coordinate_bed_name}.curve_data.pkl')
                    fnlist = glob.glob(pattern_str)
                    if len(fnlist) != 1:
                        raise Exception(f'Can not locate curve_data for Tool={toolname}, at Coord={coordinate_bed_name}, with pattern={pattern_str}, find results={fnlist}. Please check if MethPerf results folder={bdir} specified is correct.')
                    cnt += 1
                    with open(fnlist[0], 'rb') as f:
                        ret = pickle.load(f)
                        # logger.debug(ret)
                        curve_data[f'{toolname}_true'] = ret['yTrue']  # Ground truth label
                        curve_data[f'{toolname}_pred'] = ret['yPred']  # Prediction label
                        curve_data[f'{toolname}_score'] = ret['yScore']  # Prediction score, used for roc curve plotting
                outfn = os.path.join(outdir, f'{runPrefix}.plot.curve.data.ytrue.ypred.yscore.{coordDictToStandardName[coordinate_bed_name]}.pkl')
                with open(outfn, 'wb') as f:
                    pickle.dump(curve_data, f)
                # logger.info(f'save to {outfn}')

                outfn = os.path.join(outdir, f'{runPrefix}.plot.curve.data.ytrue.ypred.yscore.{coordDictToStandardName[coordinate_bed_name]}.dat')
                with open(outfn, "w") as f:
                    for toolname in ToolNameList:
                        f.write(f"{toolname}_true:")
                        outstr = ','.join([str(value) for value in curve_data[f'{toolname}_true']])
                        f.write(outstr)
                        f.write("\n")

                        f.write(f"{toolname}_pred:")
                        outstr = ','.join([str(value) for value in curve_data[f'{toolname}_pred']])
                        f.write(outstr)
                        f.write("\n")

                        f.write(f"{toolname}_score:")
                        outstr = ','.join([f'{value:.4f}' for value in curve_data[f'{toolname}_score']])
                        f.write(outstr)
                        f.write("\n")
                # logger.info(f'save to {outfn}')
            logger.info(f'For runPrefix={runPrefix} at dir={bdir}, total files={cnt}')
    elif args.cmd == 'plot-curve-data':
        ## find /projects/li-lab/Nanopore_compare/result/plot-curve-data -name '*.pkl' -exec python plot_figure.py plot-curve-data -i {} \;
        outdir = os.path.join(args.o, f'curves-figures{f"-{args.tagname}" if args.tagname else ""}')
        os.makedirs(outdir, exist_ok=True)
        for fn in args.i:
            logger.debug(f'Plot data from fn={fn}')
            with open(fn, 'rb') as infn:
                basename = os.path.basename(fn)
                bn = os.path.splitext(basename)[0]
                ret = pickle.load(infn)
                # logger.debug(ret.keys())
                plot_performance_curves(ret, outdir, tagname=bn)
    elif args.cmd == 'guppy-qos':
        ## python plot_figure.py guppy-qos -i /fastscratch/liuya/nanocompare/HL60-Runs/HL60-N50-basecall /fastscratch/liuya/nanocompare/K562-Runs/K562-N50-basecall /projects/li-lab/Nanopore_compare/nanopore_fast5/NA19240-N300-basecall /fastscratch/liuya/nanocompare/APL-Runs/APL-N50-basecall
        for basedir in args.i:
            fpatstr = os.path.join(basedir, '*', 'sequencing_summary.txt')
            flist = glob.glob(fpatstr)
            # logger.debug(flist)
            dflist = []
            for fn in flist:
                df = pd.read_csv(fn, sep='\t')
                dflist.append(df)
            outdf = pd.concat(dflist)
            logger.info(outdf)

            bfn = os.path.basename(basedir)
            outfn = os.path.join(pic_base_dir, f'{bfn}.sequencing_summary.txt')
            outdf.to_csv(outfn, index=False, sep='\t')



    else:
        raise Exception(f'Command={args.cmd} is not support')

    logger.info('DONE')

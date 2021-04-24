"""
Plots all Nanocompare paper out
"""
import argparse
import pickle
from collections import defaultdict

import glob
import matplotlib.pyplot as plt
#### Additional libraries:
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr
from sklearn import metrics
from sklearn.metrics import precision_recall_curve

from nanocompare.eval_common import correlation_report_on_regions, logger, set_log_debug_level, pic_base_dir
from nanocompare.global_settings import *


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


def grid_plot_correlation_matrix_for_fig5a(infn):
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

    # Rename to BS-seq name based on WR's suggestions
    df = df.rename(columns={'BGTruth': 'BS-seq'})

    ## Plot correlation grid figure
    basefn = os.path.basename(infn)
    outfileName = "{}.jpg".format(basefn.replace(".csv", ""))
    outfn = os.path.join(outdir, outfileName)

    plt.clf()

    fig, ax = plt.subplots()

    fig.set_size_inches(10, 10)

    left, width = .25, .5
    bottom, height = .25, .5
    right = left + width
    top = bottom + height

    position = 1
    for yrow in range(1, num_col + 1):  # row
        for xcol in range(1, num_col + 1):  # column, each tool
            logger.debug("Please wait, processing column xcol {}, yrow {}.".format(xcol, yrow))  # adding this so that one would be able to estimate how much time one has for a coffe
            if xcol == yrow:
                # Diagonal, distribution lines:
                plt.subplot(num_col, num_col, position)
                params = {'legend.fontsize'   : 50,
                        'legend.handlelength' : 0,
                        'legend.handletextpad': 0,
                        'legend.fancybox'     : False,
                        'legend.loc'          : 'upper right',
                        'legend.framealpha'   : 0,
                        'legend.borderaxespad': 0}
                plt.rcParams.update(params)
                ax = sns.kdeplot(df.iloc[:, xcol - 1], shade=False, color="black", legend=True)
                leg = ax.legend(labels=[str(df.columns[xcol - 1])], fontsize=10)
                for item in leg.legendHandles:
                    item.set_visible(False)

                # ax.annotate("abcd", xy=(0.2, 0.8), xycoords='axes fraction',
                #             fontsize=14)
                ax.set_yticklabels([])
                ax.set_xticklabels([])
                ax.set_ylabel('')
                ax.set_xlabel('')
                ax.tick_params(axis=u'both', which=u'both', length=0)

            elif xcol > yrow:
                # upper triangle, COE number:
                ax2 = plt.subplot(num_col, num_col, position)

                corrValue = pearsonr(df.iloc[:, xcol - 1], df.iloc[:, yrow - 1])
                corrValueStr = "{0:.3f}".format(corrValue[0])

                #                 print(xcol, yrow)
                if yrow == 1:
                    ax2.text(0.5 * (left + right), 0.5 * (bottom + top), corrValueStr,
                             horizontalalignment='center',
                             verticalalignment='center',
                             fontsize=19, color='#7B2323', weight='bold',
                             transform=ax2.transAxes)
                else:
                    ax2.text(0.5 * (left + right), 0.5 * (bottom + top), corrValueStr,
                             horizontalalignment='center',
                             verticalalignment='center',
                             fontsize=19, color='black',
                             transform=ax2.transAxes)

                ax2.set_yticklabels([])
                ax2.set_xticklabels([])
                ax2.tick_params(axis=u'both', which=u'both', length=0)

            elif xcol < yrow:
                # lower triangle, scatter plot, using hexbin, x-column label, y-row label:
                ax3 = plt.subplot(num_col, num_col, position)

                ### Kernel density plot:
                BGres = 300j  # 150j # 75j # 300j ## for debugging change this parameter to "50j" or less

                m1_sign = np.array(list(df.iloc[:, xcol - 1]))
                m2_sign = np.array(list(df.iloc[:, yrow - 1]))
                xmin = ymin = -0.05
                xmax = ymax = 1.05

                X, Y = np.mgrid[xmin:xmax:BGres, ymin:ymax:BGres]
                positions = np.vstack([X.ravel(), Y.ravel()])
                values = np.vstack([m1_sign, m2_sign])
                kernel = stats.gaussian_kde(values)
                Z = np.reshape(kernel(positions).T, X.shape)

                ax3.imshow(np.rot90(Z), cmap=plt.cm.Spectral_r, extent=[xmin, xmax, ymin, ymax], aspect="auto")
                # ax3.imshow(np.rot90(Z), cmap=plt.cm.gist_stern, extent=[xmin, xmax, ymin, ymax], aspect="auto") # An example how to switch from Spectral_r to gist_stern colormap
                # ax3.plot(m1_sign, m2_sign, marker=',', color="white", linestyle='None', markersize=0.25, alpha=0.25) # uncomment this line if you want to generate version with points marker very lightly on the figure

                ### Scatter plot (previous implementation for the plot - i am keeping this for now but should be cleaned if we decide to go with KDE version)
                # if len(df) < 100:
                #                     gridRes = 15  # for HL60
                #                 else:
                #                     gridRes = 50
                #
                #                 mincnt = 1
                #
                #                 if len(df) < 100:  # for HL60 plot
                #                     plt.hexbin(df.iloc[:, xcol - 1], df.iloc[:, yrow - 1], gridsize=(gridRes, gridRes), cmap='Blues', bins='log', mincnt=None)
                #                     pass
                #                 else:
                #                     plt.hexbin(df.iloc[:, xcol - 1], df.iloc[:, yrow - 1], gridsize=(gridRes, gridRes), cmap='Blues', bins='log', mincnt=mincnt)

                if yrow == num_col:  # last row scatter plot shows x ticks
                    plt.xticks([0, 1], fontsize=10)
                    ax3.set_xticklabels(["0%", "100%"], fontsize=10)
                    for label, alnType in zip(ax3.get_xticklabels(), ['left', 'right']):
                        label.set_horizontalalignment(alnType)
                else:
                    plt.xticks([], fontsize=10)

                if xcol == 1:  # first column scatter plot shows y ticks
                    plt.yticks([0, 1], fontsize=10)
                    ax3.set_yticklabels(["0%", "100%"], fontsize=10)
                    for label, alnType in zip(ax3.get_yticklabels(), ['bottom', 'top']):
                        label.set_verticalalignment(alnType)
                #
                else:
                    plt.yticks([], fontsize=10)

            position += 1

    fig.savefig(outfn, dpi=300, bbox_inches='tight')
    fig.savefig(outfn.replace("jpg", ".pdf"), dpi=300, bbox_inches='tight')  # Generate also PDF version
    plt.show()
    plt.close()
    logger.info(f"save to {outfn}")


def gen_figure_5a(infn):
    """
    Generate figure 5-a in paper
    Correlation plots
    :return:
    """
    grid_plot_correlation_matrix_for_fig5a(infn)
    pass


def plot_ROC_PR_curves(ret, outdir, tagname="tagname"):
    """
    Plot ROC AUC and PR curves
    :param ret:
    :param outdir:
    :param tagname:
    :return:
    """
    figure_size = (4, 4)

    # title_font_size = 18
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
    # plt.title('ROC Curves', fontsize=title_font_size)
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
    # plt.title('Precision-Recall Curves', fontsize=title_font_size)

    outfn = os.path.join(outdir, f'{tagname}.pr.curves.jpg')
    plt.savefig(outfn, format='png', bbox_inches='tight', dpi=600)

    plt.show()
    plt.close()


def collect_singleton_vs_nonsingleton_df(runPrefix, pattern="*.summary.singleton.nonsingleton.cov1.csv"):
    dflist = []
    # logger.debug(runPrefix)
    for run1 in runPrefix:
        filepat = os.path.join(runPrefix[run1], pattern)
        # logger.debug(run1)
        flist = glob.glob(filepat)
        if len(flist) != 1:
            raise Exception(f"Too much/No summary of singleton vs non-singleton for {runPrefix} in folder {runPrefix[runPrefix]} with pattern={pattern}, len={len(flist)}")
        # logger.debug(f'Get file:{flist[0]}')
        df = pd.read_csv(flist[0], index_col=0)
        dflist.append(df)
    retdf = pd.concat(dflist)
    retdf.index.name = 'Dataset'
    return retdf


def collect_performance_report_as_df(runPrefix):
    """
    create report from list of runPrefix, return specified columns
    :return:
    """
    dflist = []
    for runKey in runPrefix:
        logger.debug(f'runPrefix={runKey}')
        pattern = os.path.join(runPrefix[runKey], "performance?results", "*.performance.report.csv")
        files = glob.glob(pattern)
        for infn in files:
            df = pd.read_csv(infn, index_col=0, sep=",")
            dflist.append(df)
            logger.debug(f'Collect data from {runKey}:{os.path.basename(infn)} = {len(df)}')
    if len(dflist) == 0:
        raise Exception(f"Can not find report at {runKey} in folder {runPrefix[runKey]} with pattern={pattern}")
    combdf = pd.concat(dflist, ignore_index=True)
    logger.info(f'We collected total {len(dflist)} files with {len(combdf)} records')
    return combdf


def load_wide_format_performance_results(runPrefix, sel_locations=locations_category + locations_singleton):
    """
    Collect the currently new performance of exp results for paper
    :return:
    """
    df = collect_performance_report_as_df(runPrefix)
    retdf = df[df['Location'].isin(sel_locations)]

    logger.debug(f"collect_newly_exp_data, wide-format seldf={len(retdf)} using locations={sel_locations}")

    return retdf


def save_wide_format_performance_results(runPrefix, outdir, tagname):
    """
    Save all performance report results into a csv
    :return:
    """
    df = load_wide_format_performance_results(runPrefix)
    outfn = os.path.join(outdir, f'performance-results{f"-{tagname}" if tagname else ""}.csv')
    df.to_csv(outfn)
    logger.info(f'save to {outfn}')


def parse_arguments():
    parser = argparse.ArgumentParser(description='Plot and export data for Nanocompare paper.')
    parser.add_argument("cmd", help="name of command, fig5a, export-data, etc.")
    parser.add_argument('-i', nargs='+', help='list of input files', default=[])
    parser.add_argument('-o', type=str, help="output dir", default=None)  # TODO: check all correct when change this to None
    parser.add_argument('--tagname', type=str, help="tagname of files", default=None)
    parser.add_argument('--beddir', type=str, help="bed file dir used for finding Concordant and Discordant", default=None)
    parser.add_argument('--cutoff', type=int, help="the cutoff used on bed file", default=1)

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    set_log_debug_level()

    args = parse_arguments()
    logger.debug(args)

    if args.cmd == 'fig5a':
        ## find ${resultDir} -name 'Meth_corr_plot_data*.csv' \
        ##      -type f -exec python plot_figure.py fig5a -i {} \;
        ## python plot_figure.py fig5a -i /projects/li-lab/Nanopore_compare/result/MethCorr-NA19240_RRBS/Meth_corr_plot_data_joined-NA19240_RRBS-bsCov5-minToolCov3-baseFormat1.csv

        if not args.o:
            outdir = pic_base_dir
        for fn in args.i:
            gen_figure_5a(fn)
    elif args.cmd == 'export-corr-data':
        # python plot_figure.py export-corr-data -i /projects/li-lab/yang/results/2021-04-02-methcorr/MethCorr-HL60_RRBS_2Reps --beddir /projects/li-lab/yang/results/2021-04-07/MethPerf-cut5
        dflist = []
        for indir in args.i:
            fnlist = glob.glob(os.path.join(indir, 'Meth_corr_plot_data_joined-*.csv'))
            if len(fnlist) != 1:
                raise Exception(f'Found more fnlist={fnlist}')
            logger.info(f'Find file: {fnlist[0]}')

            basefn = os.path.basename(fnlist[0])
            tagname = basefn.replace('Meth_corr_plot_data_joined-', '')
            dsname = tagname[:tagname.find('_')]
            df = correlation_report_on_regions(fnlist[0], beddir=args.beddir, dsname=dsname, outdir=args.o)
            dflist.append(df)
        outdf = pd.concat(dflist)
        outfn = os.path.join(args.o, 'All.corrdata.coe.pvalue.each.regions.xlsx')
        outdf.to_excel(outfn)
        logger.info(f'save to {outfn}')
    elif args.cmd == 'export-data':
        # using global setting variable runPrefixDict to locate all data
        ## python plot_figure.py export-data -i /projects/li-lab/Nanopore_compare/result/MethPerf-APL_RRBS_CPG /projects/li-lab/Nanopore_compare/result/MethPerf-HL60_RRBS_CPG /projects/li-lab/Nanopore_compare/result/MethPerf-K562_WGBS_CPG --tagname CPG

        ## python plot_figure.py export-data -i /projects/li-lab/yang/results/2021-02-22/Meth
        if len(args.i) == 0:
            raise Exception(f'No run dir specified')
        else:
            run_prefix = defaultdict()
            for dirname in args.i:
                run_prefix[os.path.basename(dirname)] = dirname
        logger.info(run_prefix)
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
            raise Exception("No run dir specified")

        for bdir in args.i:
            runPrefix = os.path.basename(bdir)
            cnt = 0
            for coordinate_bed_name in location_filename_to_abbvname.keys():
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
                outfn = os.path.join(outdir, f'{runPrefix}.plot.curve.data.ytrue.ypred.yscore.{location_filename_to_abbvname[coordinate_bed_name]}.pkl')
                with open(outfn, 'wb') as f:
                    pickle.dump(curve_data, f)
                # logger.info(f'save to {outfn}')

                outfn = os.path.join(outdir, f'{runPrefix}.plot.curve.data.ytrue.ypred.yscore.{location_filename_to_abbvname[coordinate_bed_name]}.dat')
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
                plot_ROC_PR_curves(ret, outdir, tagname=bn)
    elif args.cmd == 'guppy-qos':
        ## collect basecall output for summary results, used for qos
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
    elif args.cmd == 'bed-to-bedGraph':
        infn = args.i[0]
        cutoff = args.cutoff
        df = pd.read_csv(infn, sep='\t', header=None)
        logger.info(len(df))

        df = df[df.iloc[:, 7] >= cutoff]
        logger.info(f'After cutoff={cutoff}, len(df)={len(df)}')
        if not args.o:
            outdir = pic_base_dir
            outfn = os.path.join(outdir, f'{os.path.basename(infn).replace(".cov1.bed", "")}.cov{cutoff}.bedGraph')
        else:
            outfn = args.o
        df = df.iloc[:, [0, 1, 2, 6]]
        df.to_csv(outfn, sep='\t', header=False, index=None)
        logger.info(f'save to {outfn}')

        from pybedtools import BedTool

        bed1 = BedTool(outfn)
        bed1 = bed1.sort()
        outfn2 = outfn.replace('.bedGraph', '.sorted.bedGraph')
        bed1.saveas(outfn2)
        logger.info(f'after sort bed, save to {outfn2}')

        pass
    else:
        raise Exception(f'Command={args.cmd} is not support')

    logger.info('Plot data script DONE')

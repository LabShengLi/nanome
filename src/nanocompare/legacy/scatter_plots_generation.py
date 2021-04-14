# import performance_plots as pp
# %matplotlib inline
#
import pandas as pd
import seaborn as sns
import glob
from lilab.tcga.global_tcga import *
from nanocompare.legacy.collect_data import collect_data_selected_locations
from nanocompare.global_settings import nanocompare_basedir, locations_category
from nanocompare.legacy.load_data import load_all_perf_data, load_refined_data


def load_data_general(path, extension="tsv", prefix="APL_BSseq_cut10/APL_Bsseq_cut10."):
    # all required files should be placed in one directory
    # print(path)
    pattern = path + "/*." + extension
    # print(pattern)
    files = glob.glob(pattern)

    dfs = []

    for f in files:
        df = pd.read_csv(f, index_col=0, sep="\t")
        dfs.append(df)

    cobmined_data = pd.concat(dfs)

    tool_names = cobmined_data["prefix"].str.split(".", expand=True)

    #     tool_names = cobmined_data["prefix"].str.split(".",n=1,expand=True).str.split("_",n=1,expand=True)#,n=1,expand=True)
    #     print(tool_names)
    cobmined_data["Tool"] = tool_names[1]

    cobmined_data = cobmined_data.replace(to_replace="False", value="x.x.GW")
    #     cobmined_data=cobmined_data.replace(to_replace=prefix,value="x.")
    #     print(cobmined_data)
    cobmined_data = cobmined_data.replace(to_replace="hg38_singletons.bed", value="x.x.singletons")
    cobmined_data = cobmined_data.replace(to_replace="hg38_nonsingletons.bed", value="x.x.nonsingletons")
    cobmined_data = cobmined_data.replace(to_replace="ONT.hg38.promoterFeature.flank_100.bed", value="x.x.promoters_100bp")
    cobmined_data = cobmined_data.replace(to_replace="ONT.hg38.promoterFeature.flank_1000.bed", value="x.x.promoters_1000bp")
    cobmined_data = cobmined_data.replace(to_replace="ONT.hg38.promoterFeature.flank_200.bed", value="x.x.promoters_200bp")
    cobmined_data = cobmined_data.replace(to_replace="ONT.hg38.promoterFeature.flank_2000.bed", value="x.x.promoters_2000bp")
    cobmined_data = cobmined_data.replace(to_replace="ONT.hg38.promoterFeature.flank_500.bed", value="x.x.promoters_500bp")
    cobmined_data = cobmined_data.replace(to_replace="ONT.hg38.promoterFeature.flank_750.bed", value="x.x.promoters_750bp")

    #     if "hg38_singletons.mixed.bed

    #     location = cobmined_data["coord"].str.split(".",expand=True)
    location = cobmined_data["coord"].str.replace(prefix, "x.")
    location = location.str.split(".", expand=True)

    #     tool = cobmined_data["prefix"].str.split(".",n=4,expand=True)

    #     print(location)
    cobmined_data["Location"] = location[2]
    #     cobmined_data["BasecallTool"] = tool[1]

    ### Because I have toooooo many labels, I am removing some categories:
    cobmined_data = cobmined_data[cobmined_data["Location"] != "promoters_100bp"]
    cobmined_data = cobmined_data[cobmined_data["Location"] != "promoters_200bp"]
    cobmined_data = cobmined_data[cobmined_data["Location"] != "promoters_750bp"]
    cobmined_data = cobmined_data[cobmined_data["Location"] != "promoters_1000bp"]
    cobmined_data = cobmined_data[cobmined_data["Location"] != "promoters_2000bp"]

    return cobmined_data


def performance_scatter_general1(df, x, y, hue=None, style=None, outfile=None):
    # https://matplotlib.org/3.1.1/api/markers_api.html
    filled_markers = ('o', 'X', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X', 'o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X')
    fig, ax = plt.subplots()
    fig.set_size_inches(11.7, 11.27)
    ax = sns.scatterplot(x=x, y=y, data=df, hue=hue, style=style, markers=filled_markers)
    # ax = fig.gca()
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    if outfile != None:
        fig.savefig(outfile, dpi=300, bbox_inches='tight')


def performance_scatter_general(df, x, y, hue=None, style=None, outfile=None, filled_markers=None):  # ('o', '<', '>', 'X', 'D', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X', 'o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X')):
    # https://matplotlib.org/3.1.1/api/markers_api.html
    #     filled_markers = ('o', '<', '>', 'X', 'D', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X', 'o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X')
    fig, ax = plt.subplots()
    fig.set_size_inches(11.7, 11.27)
    ax = sns.scatterplot(x=x, y=y, data=df, hue=hue, style=style, markers=filled_markers, s=500)
    # ax = fig.gca()
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax.set_xlabel(x, fontsize=40)
    ax.set_ylabel(y, fontsize=40)
    ax.tick_params(labelsize=30)

    ax.set_ylim(0, 1)
    ax.set_xlim(0, 1)

    if outfile != None:
        fig.savefig(outfile, dpi=300, bbox_inches='tight')


def scatter_plot_5mc_5c():
    for prefix in ["K562_WGBS_joined", "APL_BSseq_cut10", "APL_Bsseq", "K562_RRBS_joined", "HL60_RRBS_joined", "APL_oxBSseq_cut5"]:  # , "AML_Bsseq"]:
        data = collect_data_selected_locations(f"{nanocompare_basedir}/reports/{prefix}", sel_locations=["GW", "cpgIslandExt", "promoters_500bp", "exonFeature", "intergenic"])
        # print(data)#.head()

        # pd.set_option('display.max_rows', None)

        # performance_scatter(data, "recall_5C", "recall_5mC", hue="Tool", style="Location")
        outfile = os.path.join(pic_base_dir, f"{prefix}.selectedCategories.pdf")

        performance_scatter_general(data, "F1_5C", "F1_5mC", hue="Tool", style="Location", outfile=outfile, filled_markers=('o', 'P', '>', 's', 'X'))

        data = collect_data_selected_locations(f"{nanocompare_basedir}/reports/{prefix}", sel_locations=["singletons", "nonsingletons", "discordant", "concordant"])
        # print(data)#.head()
        # print(data.info())
        # print(data)
        # performance_scatter(data, "recall_5C", "recall_5mC", hue="Tool", style="Location")
        outfile = os.path.join(pic_base_dir, f"{prefix}.selectedCategories.nonsingletons_subtypes.pdf")
        performance_scatter_general(data, "F1_5C", "F1_5mC", hue="Tool", style="Location", outfile=outfile, filled_markers=('s', 'D', 'v', '^'))

        print(prefix, "DONE!")


def refine_and_save_data():
    alld = load_all_perf_data()

    refine_alld = alld[['Dataset', 'Location', 'Tool', 'F1_5C', 'F1_5mC']]

    refine_alld = pd.melt(refine_alld, id_vars=['Dataset', 'Location', 'Tool'], var_name='Measurement', value_name='Performance')

    outfn = os.path.join(pic_base_dir, "Nanocompare_performance_refined.xlsx")
    refine_alld.to_excel(outfn)

    outfn = os.path.join(pic_base_dir, "Nanocompare_performance_refined.pkl")
    refine_alld.to_pickle(outfn)

    return refine_alld


def simple_step1_1ds_plot():
    df = load_refined_data()

    dsname = 'K562_WGBS_joined'
    measure = 'F1_5mC'
    df1 = df[(df['Dataset'] == dsname) & (df['Measurement'] == measure)]
    df1 = df1[df['Location'].isin(locations_category)]
    logger.info(f"df1={df1}")

    sns_plot = sns.catplot(x="Location", y="Performance", hue='Tool', data=df1)

    fig = sns_plot.get_figure()
    outfn = os.path.join(pic_base_dir, f"{dsname}_{measure}_catplot.png")
    fig.savefig(outfn, dpi=300, format='png')


if __name__ == '__main__':
    set_log_debug_level()
    # set_log_info_level()

    # add_dataset_report_to_all()

    # alld = save_alldata()

    # df = refine_and_save_data()

    # df = load_refined_data()
    # simple_step1_1ds_plot()

    # refine_and_save_data()

    # scatter_plot_5mc_5c()

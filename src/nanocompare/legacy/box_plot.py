import importlib
import nanocompare.legacy.performance_plots as pp
import seaborn as sns
from nanocompare.nanocompare_global_settings import nanocompare_basedir
from lilab.tcga.global_tcga import *
from lilab.tcga.utils import current_time_str

# %matplotlib inline

importlib.reload(pp)

folder = "reports_dump"
prefix = current_time_str()
metric = "roc_auc"
metric = "accuracy"
metrics = ["accuracy", "roc_auc", "F1_5C", "F1_5mC"]

data = pp.load_data(f"{nanocompare_basedir}/reports/{folder}")

regions = ["GW", "singletons", "nonsingletons", "discordant", "concordant", "cpgIslandExt", "cpgShoresExt", "cpgShelvesExt", "promoterFeature", "exonFeature", "intronFeature", "geneFeature", "intergenic"]
regions = ["GW", "singletons", "nonsingletons", "discordant", "concordant", "cpgIslandExt", "promoterFeature", "exonFeature", "intronFeature", "intergenic"]

for metric in metrics:
    plt.clf()
    fig = plt.gcf()
    fig.set_size_inches(20, 5)
    plt.rcParams.update({'font.size': 15})

    for i, region in enumerate(regions):
        order = i + 1
        df_filtered = data[data.Location == region]
        df_filtered = df_filtered[df_filtered.accuracy > 0]
        df_filtered.head()

        plt.subplot(1, len(regions), order)

        # # sns.boxplot(x="method", y="AUC", hue="BinaryLabels", data=df, palette="Set1")
        # ax = sns.violinplot(x="BasecallTool", y="accuracy", data=df_filtered, palette="Set1", order=["Tombo_calls", "Nanopolish_calls", "DeepSignal_calls", "DeepMod_calls"], width=0.3)
        ax = sns.boxplot(x="BasecallTool", y=metric, data=df_filtered, palette="Set1", order=["Tombo_calls", "Nanopolish_calls", "DeepSignal_calls", "DeepMod_calls"], width=0.85)
        # ax = sns.boxplot(x="BasecallTool", y="accuracy", data=df_filtered, palette="Set1", order=["Tombo_calls", "Nanopolish_calls", "DeepSignal_calls", "DeepMod_calls"])
        ax = sns.swarmplot(x="BasecallTool", y=metric, data=df_filtered, order=["Tombo_calls", "Nanopolish_calls", "DeepSignal_calls", "DeepMod_calls"], color="yellow")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        plt.ylim(0, 1)
        if order > 1:
            ax.set_yticklabels([])
            plt.ylabel("")
        ax.set_xticklabels(["Tombo", "Nanopolish", "DeepSignal", "DeepMod"])
        if region == "cpgIslandExt":
            plt.title("CGI")
        elif region == "promoterFeature":
            plt.title("Promoters")
        elif region == "intronFeature":
            plt.title("Introns")
        elif region == "exonFeature":
            plt.title("Exons")
        elif region == "singletons":
            plt.title("singletons")
        elif region == "nonsingletons":
            plt.title("NS")
        elif region == "discordant":
            plt.title("dist-NS")
        elif region == "concordant":
            plt.title("con-NS")
        elif region == "intergenic":
            plt.title("Intergenic")
        else:
            plt.title(region)

    outfn = os.path.join(pic_base_dir, "{}.{}.{}_BoxPlot.png".format(prefix, folder, metric))
    plt.savefig(outfn, bbox_inches='tight', dpi=300)
    logger.info(f"save to {outfn}")

outfn = os.path.join(pic_base_dir, "{}.joinedReports.tsv".format(prefix))
data.to_csv(outfn, sep='\t', index=False)
logger.info(f"save to {outfn}")

# for metric in metrics:
#     plt.clf()
#     fig = plt.gcf()
#     fig.set_size_inches(20, 5)
#     plt.rcParams.update({'font.size': 15})
#
#     for i, region in enumerate(regions):
#         order = i + 1
#         df_filtered = data[data.Location == region]
#         df_filtered = df_filtered[df_filtered.accuracy > 0]
#         df_filtered.head()
#
#         plt.subplot(1, len(regions), order)
#
#         # # sns.boxplot(x="method", y="AUC", hue="BinaryLabels", data=df, palette="Set1")
#         # ax = sns.violinplot(x="BasecallTool", y="accuracy", data=df_filtered, palette="Set1", order=["Tombo_calls", "Nanopolish_calls", "DeepSignal_calls", "DeepMod_calls"], width=0.3)
#         ax = sns.barplot(x="BasecallTool", y=metric, data=df_filtered, palette="Set1", order=["Tombo_calls", "Nanopolish_calls", "DeepSignal_calls", "DeepMod_calls"])  # , width=0.85)
#         # ax = sns.boxplot(x="BasecallTool", y="accuracy", data=df_filtered, palette="Set1", order=["Tombo_calls", "Nanopolish_calls", "DeepSignal_calls", "DeepMod_calls"])
#         #         ax = sns.swarmplot(x="BasecallTool", y=metric, data=df_filtered, order=["Tombo_calls", "Nanopolish_calls", "DeepSignal_calls", "DeepMod_calls"], color="yellow")
#         ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
#         plt.ylim(0, 1)
#         if order > 1:
#             ax.set_yticklabels([])
#             plt.ylabel("")
#         ax.set_xticklabels(["Tombo", "Nanopolish", "DeepSignal", "DeepMod"])
#         if region == "cpgIslandExt":
#             plt.title("CGI")
#         elif region == "promoterFeature":
#             plt.title("Promoters")
#         elif region == "intronFeature":
#             plt.title("Introns")
#         elif region == "exonFeature":
#             plt.title("Exons")
#         elif region == "singletons":
#             plt.title("singletons")
#         elif region == "nonsingletons":
#             plt.title("NS")
#         elif region == "discordant":
#             plt.title("dist-NS")
#         elif region == "concordant":
#             plt.title("con-NS")
#         elif region == "intergenic":
#             plt.title("Intergenic")
#         else:
#             plt.title(region)
#
#     plt.savefig("{}.{}.{}_BarPlot.pdf".format(prefix, folder, metric), bbox_inches='tight', dpi=300)
#     print("{}.{}.{}_BarPlot.pdf PLOTTED".format(prefix, folder, metric))
#
# data.to_csv("{}.joinedReports.BarPlot.tsv".format(prefix), sep='\t', index=False)

# df_filtered = data[data.Location == "GW"]
# df_filtered.head()

# fig = plt.gcf()
# fig.set_size_inches(10,10)
# plt.rcParams.update({'font.size': 18})
# # # sns.boxplot(x="method", y="AUC", hue="BinaryLabels", data=df, palette="Set1")
# # ax = sns.violinplot(x="BasecallTool", y="accuracy", data=df_filtered, palette="Set1", order=["Tombo_calls", "Nanopolish_calls", "DeepSignal_calls", "DeepMod_calls"], width=0.3)
# ax = sns.boxplot(x="BasecallTool", y="accuracy", data=df_filtered, palette="Set1", order=["Tombo_calls", "Nanopolish_calls", "DeepSignal_calls", "DeepMod_calls"], width=0.3)
# # ax = sns.boxplot(x="BasecallTool", y="accuracy", data=df_filtered, palette="Set1", order=["Tombo_calls", "Nanopolish_calls", "DeepSignal_calls", "DeepMod_calls"])
# ax = sns.swarmplot(x="BasecallTool", y="accuracy", data=df_filtered, order=["Tombo_calls", "Nanopolish_calls", "DeepSignal_calls", "DeepMod_calls"], color="yellow")
# ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
# plt.title(folder)

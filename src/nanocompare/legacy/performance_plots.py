import pandas as pd
import seaborn as sns
import glob
import os


def load_data(path, extension="tsv"):
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

    tool_names = cobmined_data["prefix"].str.split("_", n=1, expand=True)
    cobmined_data["Tool"] = tool_names[0]

    cobmined_data = cobmined_data.replace(to_replace="False", value="x.x.GW")

    cobmined_data = cobmined_data.replace(to_replace="hg38_singletons.bed", value="x.x.singletons")
    cobmined_data = cobmined_data.replace(to_replace="hg38_nonsingletons.bed", value="x.x.nonsingletons")

    # 	if "hg38_singletons.mixed.bed

    location = cobmined_data["coord"].str.split(".", n=4, expand=True)

    tool = cobmined_data["prefix"].str.split(".", n=4, expand=True)

    # 	print(location)
    cobmined_data["Location"] = location[2]
    cobmined_data["BasecallTool"] = tool[1]
    return cobmined_data


def performance_scatter(df, x, y, hue=None, style=None, outfile=None):
    ax = sns.scatterplot(x=x, y=y, data=df, hue=hue, style=style)
    # ax = fig.gca()
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

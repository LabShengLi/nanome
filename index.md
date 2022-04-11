## Tutorial of DNA methylation calling for ONT data
In this totorial, you will learn how to do methylation calling on Oxford Nanopore sequencing data by latest tools. Please create a fresh new folder to run following commands.

### 1. Software installation
#### 1.1 Install Guppy
Instal the latest version of Guppy:
```
wget https://americas.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_6.1.1_linux64.tar.gz
tar -xzf ont-guppy-cpu_6.1.1_linux64.tar.gz && \
    rm -f ont-guppy-cpu_6.1.1_linux64.tar.gz

ont-guppy-cpu/bin/guppy_basecaller  -v
```

#### 1.2 Install Conda
If you do not have conda, please follow this link(https://docs.conda.io/en/latest/miniconda.html) to install conda.

#### 1.3 Install Megalodon
[Megalodon](https://github.com/nanoporetech/megalodon) is a popular and latest ONT developed methylation-calling tool. It can be installed in conda enviroment.
```
conda create --name megalodon python=3.9
conda activate megalodon

pip install megalodon
megalodon -v
```




You can use the [editor on GitHub](https://github.com/TheJacksonLaboratory/nanome/edit/gh-pages/index.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [Basic writing and formatting syntax](https://docs.github.com/en/github/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/TheJacksonLaboratory/nanome/settings/pages). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.

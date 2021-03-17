bootstrap: docker
from: continuumio/miniconda3

%files
    /projects/li-lab/yang/workspace/nano-compare/environment.yml

%post
    /home/liuya/anaconda3/bin/conda env create -f /projects/li-lab/yang/workspace/nano-compare/environment.yml

%runscript
    exec /home/liuya/anaconda3/envs/$(head -n 1 /projects/li-lab/yang/workspace/nano-compare/environment.yml | cut -f 2 -d ' ')/bin/"$@"

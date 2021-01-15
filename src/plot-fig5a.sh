#!/bin/bash

set -x

resultDir=/projects/li-lab/yang/results/2021-01-15/

find ${resultDir} -name 'Meth_corr_plot_data*.csv' \
		-type f -exec python plot_figure.py fig5a -i {} \;



#!/bin/bash

awk ' BEGIN { OFS="\t" } {print $1,$2,$3,$6,$5 }' $1 |sed '1d'>${1}.splice_plot_ready.canonical.bed


#!/bin/bash

awk 'BEGIN {OFS="\t"} { if ($4 ~ /circ/) print $1,$2,$3,$6,$7,$0 }' $1 > ${1}.build_db.bed

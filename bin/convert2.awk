awk 'BEGIN { FS = "\t" }; {print $1,$2,$3,$6,$4}' $1>${1}.circles.bed

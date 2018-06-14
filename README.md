# circleVis #
Visualize coverage, canonical, and backsplice junctions.
## Requirements ##
circleVis works with Python 3.0 and newer.
## Dependencies ##
* Matplotlib
* Numpy
## Installation ##
To install circleVis:

```
NOT IMPLEMENTED YET: pip install circleVis
```

Or:

```
git clone https://github.com/flemingtonlab/circleVis.git
```

## Example ##
To run the example dataset:

```
git clone https://github.com/flemingtonlab/circleVis.git 

cd circleVis/example 

wget ftp://ftp.ensembl.org/pub/release-92/gtf/homo_sapiens/Homo_sapiens.GRCh38.92.gtf.gz .

gunzip Homo_sapiens.GRCh38.92.gtf.gz 

python ../bin/circbuild --gtf Homo_sapiens.GRCh38.92.gtf --wigneg chr19a_negative_strand.wig --wigpos chr19a_positive_strand.wig --splicejunction chr19_canonical_junctions.bed --circlejunction chr19_backsplice_junctions.bed --output example

python ../bin/circplot -db example.db -g PIP5K1C
```

## Authors ##
Created by Nathan Ungerleider and Erik Flemington

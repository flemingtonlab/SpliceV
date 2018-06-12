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
pip install circleVis
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

wget ftp://ftp.ensembl.org/pub/release-92/gtf/homo_sapiens . 

python build_db.py --gtf Homo_sapiens.GRCh38.92.gtf --wigneg chr19a_negative_strand.wig --wigpos chr19_positive_strand --splicejunction chr19_canonical_junctions.bed --circlejunction chr19_backsplice_junctions.bed --output example

python plot_transcript.py -db example.db -g PIP5K1C
```

## Authors ##
Created by Nathan Ungerleider and Erik Flemington

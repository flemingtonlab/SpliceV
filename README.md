# SpliceV #
Visualize coverage, canonical, and backsplice junctions.

![Example plot](https://github.com/flemingtonlab/SpliceV/blob/master/etc/example.png)

## Documentation ##
See https://splicev.readthedocs.io/en/master/

## Example pipeline ##
See https://github.com/flemingtonlab/SpliceV/blob/master/docs/example.pdf

This will generate figure 1B and 1C from our manuscript (DOI pending)

## Requirements ##
SpliceV works with Python 2.7 and 3.0+.
## Dependencies ##
* Matplotlib
* Numpy
* pysam
## Installation ##
To install SpliceV:

```
pip install SpliceV
```

Or:

```
git clone https://github.com/flemingtonlab/SpliceV.git
```

## Example ##
To run the example dataset:

```
git clone https://github.com/flemingtonlab/SpliceV.git 

cd SpliceV/example 

python ../bin/SpliceV -b example.vta1.bam -gtf vta1.gtf -sj example.canonical.bed -bsj example.circles.bed -g VTA1 -f 3 -is 3

```

The sample data comes from Akata cells (a B Cell Lymphoma line) that were treated with the exonuclease RNase R prior to sequencing. RNase R exclusively digests RNA with a free end, helping enrich circularized RNA abundance in the sample. 

These example commands will generate the following plot:
![User example plot](https://github.com/flemingtonlab/SpliceV/blob/master/etc/vta1.png)

This plot reveals a prominant circle from exon 2 through exon 4 (evidenced by the back-splice junction reads which are plotted as curves below the axis. This is also demonstrated by the higher level of coverage in exon 2-4, shown by the relative intensity of color).

![User example plot explained](https://github.com/flemingtonlab/SpliceV/blob/master/etc/vta1_explained.png)

The major circularized isoform (exons 2-4; another less prevalent circle appears to include exon 5) is isolated below:


![User example circle](https://github.com/flemingtonlab/SpliceV/blob/master/etc/vta1_circ.png)

## Authors ##
Created by Nathan Ungerleider and Erik Flemington

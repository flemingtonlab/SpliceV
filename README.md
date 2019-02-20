# SpliceV #
Visualize coverage, canonical, and backsplice junctions.

![Example plot](https://github.com/flemingtonlab/SpliceV/blob/master/etc/example.png)

## Documentation ##
See https://splicev.readthedocs.io/en/master/
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

../bin/SpliceV.py -b *bam -gtf *gtf -sj *canon* -bsj *circ* -g VTA1 -f 3
```
This will generate the following plot:
![User example plot](https://github.com/flemingtonlab/SpliceV/blob/master/etc/vta1.png)

The sample data comes from Akata cells (a B Cell Lymphoma line) that were treated with the exonuclease RNase R prior to sequencing. RNase R exclusively digests RNA with a free end, in effect sparing circularized RNA and enriching its relative abundance in the sample. This example plot reveals a prominant circle from exon 2 through exon 4 (evidenced by the back-splice junction reads which are plotted as curves below the axis. This is also demonstrated by the higher level of coverage in exon 2-4, shown by the relative intensity of color).


## Authors ##
Created by Nathan Ungerleider and Erik Flemington

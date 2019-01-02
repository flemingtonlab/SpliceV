# SpliceV #
Visualize coverage, canonical, and backsplice junctions.

![Example plot](https://github.com/flemingtonlab/SpliceV/blob/master/etc/example.png)

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

python ../bin/SpliceV -b example.bam -gtf uba2.gtf -g UBA2 -f 5

```

## Authors ##
Created by Nathan Ungerleider and Erik Flemington

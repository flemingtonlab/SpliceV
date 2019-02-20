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

![User example plot](https://github.com/flemingtonlab/SpliceV/blob/master/etc/vta1.png)

![User example plot] (https://github.com/flemingtonlab/SpliceV/blob/master/etc/vta1.png)

## Authors ##
Created by Nathan Ungerleider and Erik Flemington

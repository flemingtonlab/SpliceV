# spliceV #
Visualize coverage, canonical, and backsplice junctions.

![Example plot](https://github.com/flemingtonlab/spliceV/blob/master/example/example.png)


## Requirements ##
spliceV works with Python 2.7 and 3.0+.
## Dependencies ##
* Matplotlib
* Numpy
## Installation ##
To install spliceV:

```
pip install spliceV
```

Or:

```
git clone https://github.com/flemingtonlab/spliceV.git
```

## Example ##
To run the example dataset:

```
git clone https://github.com/flemingtonlab/spliceV.git 

cd spliceV/example 

wget ftp://ftp.ensembl.org/pub/release-92/gtf/homo_sapiens/Homo_sapiens.GRCh38.92.gtf.gz .

gunzip Homo_sapiens.GRCh38.92.gtf.gz 

python ../bin/spliceV.py --gtf Homo_sapiens.GRCh38.92.gtf 

```

## Authors ##
Created by Nathan Ungerleider and Erik Flemington

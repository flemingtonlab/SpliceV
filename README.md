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

## Example ##
To run the example dataset:

```
git clone https://github.com/flemingtonlab/circVis.git 
cd circVis/example 
wget ftp://ftp.ensembl.org/pub/release-92/gtf/homo_sapiens . 
python build_db.py --gtf Homo_sapiens.GRCh38.92.gtf --wigneg --wigpos --splicejunction --circlejunction --output example
```

## Authors ##
Created by Nathan Ungerleider and Erik Flemington



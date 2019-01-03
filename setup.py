import setuptools

setuptools.setup(
    name='SpliceV',
    version='0.1.3',
    author='Nathan Ungerleider',
    author_email='nungerle@tulane.edu',
    description='Visualize splice junctions, backsplice junctions (circleRNA) and coverage from RNA-Seq datasets',
    url='https://github.com/flemingtonlab/SpliceV',
    packages=setuptools.find_packages(),
    long_description=open('README.md').read(),
    classifiers=('Programming Language :: Python :: 3',
                'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                'Operating System :: OS Independent',
                ),
    install_requires = ['matplotlib', 'numpy', 'pysam', 'pkg_resources'],
    python_requires='>=2.7,>=3',
    scripts=['bin/SpliceV', 'bin/find_circ_convert', 'bin/star_sj_convert'],
)


import setuptools

setuptools.setup(
    name='SpliceV',
    version='0.1.6.8',
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
    install_requires = ['matplotlib', 'numpy', 'pysam'],
    python_requires='>=2.7,>=3',
    scripts=['bin/SpliceV', 'bin/fa.py','bin/find_circ_convert', 'bin/star_sj_convert'],
    include_package_data = True,
)


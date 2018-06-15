import setuptools

setuptools.setup(
    name='circleVis',
    version='0.1.dev4',
    author='Nathan Ungerleider',
    author_email='nungerle@tulane.edu',
    description='Visualize splice junctions, backsplice junctions (circleRNA) and coverage from RNA-Seq datasets',
    url='https://github.com/flemingtonlab/circleVis',
    packages=setuptools.find_packages(),
    long_description=open('README.md').read(),
    classifiers=('Programming Language :: Python :: 3',
                'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                'Operating System :: OS Independent',
                ),
    install_requires = ['matplotlib', 'numpy'],
    python_requires='>=3',
    scripts=['bin/circplot', 'bin/circbuild', 'bin/find_circ_convert', 'bin/star_sj_convert'],
)


from collections import *
import re
import pysam
import numpy as np


gene = 'SPPL2A'
path = '/home/nate/Documents/gtf/gencode.v28.annotation.gtf'

path = '/Users/mac9/circleVis/example/Homo_sapiens.GRCh38.92.gtf'
bampath = '/Volumes/Flemington_lab_experiments/temp/2_Akata_60_induced_cr_RNaseR_FCHL5J5BBXX_L3_RDWHHUMcrmEAAHRAB-209_1.fq_STARAligned.sortedByCoord.out.bam'


def longest_transcript(coordinates):
    ''' Returns longest transcript isoform of a gene '''

    lengths = defaultdict(int)
    for transcript, exons in coordinates.items():
        for exon in exons:
            lengths[transcript] += int(exon[2]) - int(exon[1]) + 1
    
    if len(lengths) == 0:
        sys.exit("%s not found in gtf file.")

    longest = 0
    for transcript, length in lengths.items():
        if length >= longest:
            longest = length
            long_transcript = transcript

    return long_transcript


def exons(path, gene, transcript=False):
    '''Returns exon coordinates from gtf file -> (chromosome, 5prime coord, 3prime coord, strand) for each exon'''

    coordinates = defaultdict(list)

    if transcript:
        prog = re.compile('transcript_id "%s"' % gene)
    else:
        prog = re.compile('gene_name "%s"' % gene)

    with open(path) as gtf:
        Exon = namedtuple('Exon', ['chromosome', 'source', 'feature', 'start', 'stop', 'score', 'strand', 'frame'])
        for line in gtf:

            # Skip header lines
            if line[0] == '#':
                continue

            # Only lines with gene name
            if prog.search(line):
                *line, attributes = line.split('\t')
                exon = Exon(*line)
                if exon.feature != 'exon':
                    continue

                transcript = attributes.split('transcript_id ')[1].split('"')[1]
                coordinates[transcript].append((exon.chromosome, exon.start, exon.stop, exon.strand))

    coords = coordinates[ longest_transcript(coordinates) ]
    return coords


def prep_bam(path):

    sam = pysam.AlignmentFile(path)
    try:
        sam.check_index()
    except ValueError:
        try:
            pysam.index(path)
        except:
            pysam.sort("-o", path.replace('bam', 'sorted.bam'), path)
            pysam.index(path)

    except AttributeError:
        pysam.view('-bho', path.replace('sam', 'bam'), path)
        sam = pysam.AlignmentFile(path.replace('sam', 'bam'))

        try:
            sam.check_index()
        except ValueError:
            try:
                pysam.index(path)
            except:
                pysam.sort("-o", path.replace('bam', 'sorted.bam'), path)

    return sam

def coverage(bam, chromosome, start, stop, strand='+'):
    '''Takes AlignmentFile instance and returns coverage and bases in a dict with position as key'''
    
    try:
        pileup = bam.pileup(chromosome, int(start), int(stop))

    except ValueError:
        if 'chr' not in chromosome:
            pileup = bam.pileup('chr' + chromosome, int(start), int(stop))
        else:
            pileup = bam.pileup(chromosome.replace('chr',''), int(start), int(stop))

    coverage = defaultdict(list)

    for column in pileup:
        for read in column.pileups:         
            if strand == '-' and not read.alignment.is_reverse:
                continue
            if strand == '+' and read.alignment.is_reverse:
                continue
            if not read.is_del and not read.is_refskip:
                base = read.alignment.query_sequence[read.query_position]
                coverage[column.pos].append(base)

    for key in coverage:
        c = Counter(coverage[key])
        sum_c = sum(c.values())
        coverage[key] = (sum_c, c)

    avg = []
    for i in coverage:
        avg.append(coverage[i][0])

    if len(avg)>0:
        return np.mean(avg)
    else:
        return 0


def junctions(bam, chromosome, upstream, downstream, strand):
    try:
        fetched = bam.fetch(chromosome, upstream, downstream)
    except ValueError:   
        if 'chr' not in chromosome:
            bam.fetch('chr' + chromosome, upstream, downstream)
        else:
            try:
                bam.fetch(chromosome.replace('chr', ''), upstream, downstream)
            except ValueError:
                sys.exit("Chromosome %s not found in bam file.." % chromosome)

    if strand == '-':
        stranded = (read for read in fetched if read.is_reverse)
    else:
        stranded = (read for read in fetched if not read.is_reverse)
 
    introns = bam.find_introns(stranded).items()
    
    return [(i[0], i[1], j) for i,j in introns]

cov=[]
bam = prep_bam(bampath)
coords = exons(path, gene)
for chromosome, start, stop, strand in coords:
    cov.append(coverage(bam, chromosome, start, stop, strand))
low_coord = min(int(i[1]) for i in coords)
high_coord = max(int(i[2]) for i in coords)
sj = junctions(bam, chromosome, low_coord, high_coord, strand)
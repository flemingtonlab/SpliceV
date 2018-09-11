from collections import *
import re
import pysam

gene = 'EGFR'
path = '/home/nate/Documents/gtf/gencode.v28.annotation.gtf'



def longest_transcript(coordinates):
    ''' Returns longest transcript isoform of a gene '''

    lengths = defaultdict(int)
    for transcript, exons in coordinates.items():
        for exon in exons:
            lengths[transcript] += int(exon[2]) - int(exon[1]) + 1

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
    return coords.sort(lambda x:x[1])


def prep_bam(path):

    sam = pysam.AlignmentFile(path)
    try:
        sam.check_index()
    except ValueError:
        try:
            pysam.index(path)
        except:
            pysam.sort("-o", path.replace('bam', 'sorted.bam'), path)


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

    pileup = bam.pileup(chromosome, start, stop)
    coverage = defaultdict(list)
    for column in pileup:
       # print("\ncoverage at base %s = %s" % (column.pos,column.n))
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

    
    return np.mean([coverage[i][0] for i in coverage])

def junctions(bam, chromosome, start, stop):

     return [(i[0], i[1], j) for i,j in bam.find_introns(bam.fetch(chromosome, start, stop)).items()]

coords = exons(path, gene)
for chromosome, start, stop, strand in coords:
    cov.append(coverage(bam, chromosome, start, stop, strand))
    
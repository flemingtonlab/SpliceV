from collections import namedtuple, defaultdict
import re

coordinates = defaultdict(list)
gene = 'EGFR'
path = '/Users/mac9/circleVis/example/Homo_sapiens.GRCh38.92.gtf'

with open(path) as gtf:
    Exon = namedtuple('Exon', ['chromosome', 'source', 'feature', 'start', 'stop', 'score', 'strand', 'frame'])
    for line in gtf:
        if line[0] == '#':
            continue
        prog = re.compile('gene_name "%s"'%gene)
        match = prog.search(line)
        if match:
            *line, attributes = line.split('\t')

            exon = Exon(*line)
            if exon.feature != 'exon':
                continue
            transcript = attributes.split('transcript_id ')[1].split('"')[1]

            coordinates[transcript].append((exon.chromosome, exon.start, exon.stop, exon.strand))
            longest = defaultdict(int)
            for transcript, exons in coordinates.items():
                for exon in exons:
                    longest[transcript] += int(exon[2]) - int(exon[1]) + 1
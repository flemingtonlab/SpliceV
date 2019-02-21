import re
import os
import glob
from RNABP import get_rnabp

# Convert ambiguous nucleotides to regex expressions
regex_d = {
    'A':'A',
    'B':'[GTC]',
    'C':'C',
    'D':'[GAT]',
    'G':'G',
    'H':'[ACT]',
    'K':'[GT]',
    'M':'[AC]',
    'T':'T',
    'Y':'[CT]',
    'S':'[GC]',
    'W':'[AT]',
    'V':'[GCA]',
    'N':'[ACTG]',
    'R':'[GA]',
    'U': 'T'
}


rnabp = get_rnabp()

def bp_positions(query, seq, start):

    positions = []
    if query not in rnabp:
        print("RNA binding protein %s not found in list." %query)
        return []
    searchstring = ''.join([regex_d[i] for i in rnabp[query]])
    for i in re.finditer(searchstring, seq):
        mid = sum(i.span()) / 2
        positions.append(mid)
    positions = [i + start for i in positions]
    return positions


def index_fasta(path):
    # NAME	Name of this reference sequence
    # LENGTH	Total length of this reference sequence, in bases
    # OFFSET	Offset in the FASTA/FASTQ file of this sequence's first base
    # LINEBASES	The number of bases on each line
    # LINEWIDTH	The number of bytes in each line, including the newline
        
    fa_ix = []
    f = open(path)
    name = offset = line_len = ini_stripped_line_len = stripped_line_len = 0
    total_seq_len = 0
    line = '_'
    set_name = False
    
    while line:

        if line[0] == '>':
            if name:
                fa_ix.append((name, total_seq_len, offset, ini_stripped_line_len, line_len))
            set_name = True
            line = line.strip()
            name = line.replace('>', '').split()[0]
            offset = f.tell()
            total_seq_len = 0
            print(line, f.tell())
            line = f.readline()
        
        elif set_name:
            line_len = len(line)
            ini_stripped_line_len = len(line.strip())  # Windows adds 2 bytes for newline char, others add 1
            set_name = False

        stripped_line_len = len(line.strip())
        total_seq_len += stripped_line_len

        line = f.readline()

    fa_ix.append((name, total_seq_len, offset, ini_stripped_line_len, line_len))
    seq_map = {}

    with open(path + '.fai', 'w') as index_file:
        for item in fa_ix:
            name = item[0]
            index_file.write('%s\n' % '\t'.join([str(j) for j in item]))
            seq_map[name] = path

    return seq_map

def prep_fasta(paths):
    # TODO: include function to download reference sequences (and ref gtfs?) - maybe have a dictionary of species
    seq_map = {}
    if len(paths) == 1:
        path = paths[0]
        if os.path.isdir(path):
            paths = glob.glob(os.path.join(path, '*.fa'))
            if len(paths) == 0:
                print("No fasta files found in directory '%s'. Make sure fasta files end in '.fa'." % path)
                return

    for path in paths:
        if path.endswith(".fa") and os.path.isfile(path):
            # Single genome fasta file with headers. Make sure indexed or index here
            fasta_index_path = path + '.fai'
            if not os.path.isfile(fasta_index_path):
                print("No fasta index, %s, found. Generating fasta index..." % fasta_index_path)

                temp_seq_map = index_fasta(path)
                seq_map = {**seq_map, **temp_seq_map}
    return seq_map

def rcomp(seq, reverse=True):
    
    dna = dict(zip('ATCG','TAGC'))
    seq = seq.upper()

    complement = ''.join([dna[i] if i in dna else i for i in seq])

    if reverse:
        return complement[::-1]
    
    return complement


def read_fasta(fa_path, chromosome, start, stop, strand):
    '''Requires path to indexed fasta file. Returns the nucleotide sequence'''

    fai_path = fa_path + '.fai'

    if not os.path.exists(fa_path):
        print("Fasta file %s was not found. Please make sure the path is correct and the file exists." % fa_path)
    
    if not os.path.exists(fai_path):
        print("Fasta index for %s was not found (%s). Please make sure the path is correct and the file exists. Indexing now.." % (fa_path, os.path.basename(fai_path)))
        index_fasta(fa_path)
      
    with open(fai_path, 'r') as fai:
        for line in fai:
            name, seq_len, offset, nbases_line, nbytes_line = line.strip().split('\t')

            if name == chromosome:  
                seq_len = int(seq_len)
                offset = int(offset)
                nbases_line = int(nbases_line)
                nbytes_line = int(nbytes_line)
                
                break
            
    start_offset = nbytes_line * start // nbases_line
    stop_offset = nbytes_line * stop // nbases_line
    if strand == '-':
        start_offset, stop_offset = stop_offset, start_offset

    start_offset -= 1
    seek_len = abs(stop_offset - start_offset)

    with open(fa_path, 'r') as fa:
        fa.seek(offset + start_offset)
        seq = fa.read(seek_len).replace('\n', '')
    if strand == '-':
        return rcomp(seq)
    return seq


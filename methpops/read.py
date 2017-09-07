import cigar
from utils import colored

def generate_read_from_bam_line(line):
    return Read(*line.split('\t'))

class Read:
    def __init__(self, readName, flag, chromosome, position, \
                    mappingQuality, CIGAR, mateName, matePosition, templateLength, \
                    seq, quality, *args):
        self.readName = readName
        self.flag = int(flag)
        if self.flag in [0, 99, 147]:
            self.reverse = False
        else:
            self.reverse = True
        self.chromosome = chromosome
        self.position = int(position)
        self.mappingQuality = mappingQuality
        self.CIGAR = CIGAR
        self.mateName = mateName
        self.matePosition = matePosition
        self.templateLength = templateLength
        self.seq = seq
        self.quality = quality
        self.methylation = self._parse_tags(args)
        # print(self.methylation)

    def _parse_tags(self, tags):
        for tag in tags:
            t = tag.split(':')
            if t[0] == 'XM': # XM indicates methylation call string
                return t[2]

    def get_CpGs(self):
        coordinates = cigar.get_genomic_coordinates(self.CIGAR, self.position)
        CpGs = []
        for c, coordinate in zip(self.methylation, coordinates):
            if c == 'z':
                CpGs.append((coordinate - [0, 1][self.reverse], 0))
            elif c == 'Z':
                CpGs.append((coordinate - [0, 1][self.reverse], 1))

        return CpGs

    def get_full_CpG_colored_string(self):
        nucs = []
        for nuc, meth in zip(self.seq, self.methylation):
            if meth == 'z':
                nucs.append(colored(nuc, 'green'))
            elif meth == 'Z':
                nucs.append(colored(nuc, 'red'))
            else:
                nucs.append(nuc)
            
        return ''.join(nucs)

    def get_CpG_colored_string(self):
        nucs = []
        for nuc, meth in zip(self.seq, self.methylation):
            if meth == 'z':
                nucs.append(colored(' ', 'grey'))
            elif meth == 'Z':
                nucs.append(colored(' ', 'red'))
            
        return ''.join(nucs)
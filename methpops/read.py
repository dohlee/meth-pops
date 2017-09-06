import cigar

def generate_read_from_bam_line(line):
	return Read(*line.split('\t'))

class Read:
	def __init__(self, readName, flag, chromosome, position, \
					mappingQuality, CIGAR, mateName, matePosition, templateLength, \
					seq, quality, *args):
		self.readName = readName
		self.flag = flag
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
				CpGs.append((coordinate, 0))
			elif c == 'Z':
				CpGs.append((coordinate, 1))

		return CpGs
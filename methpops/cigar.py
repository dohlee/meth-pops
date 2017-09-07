import re

def get_genomic_coordinates(cigar, start, reverse=False):
    ashes = re.findall(r'(\d+)([A-Z]{1})', cigar)

    curr = start
    coordinates = []
    for length, alignType in ashes:

        # Match
        if alignType == 'M':
            for i in range(int(length)):
                coordinates.append(curr)
                curr += 1

        # Insertion
        elif alignType == 'I':
            for i in range(int(length)):
                # the base in the read does not appear in the reference
                coordinates.append(-1)

        # Deletion
        elif alignType == 'D':
            for i in range(int(length)):
                curr += 1

    return coordinates

if __name__ == '__main__':
    print(get_genomic_coordinates('1M2I12M2D30M', 1))

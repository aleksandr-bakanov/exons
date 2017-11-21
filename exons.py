import csv
import re
import sys

class Exon:
    """Exon represents one line in an output BED file"""

    def __init__(self, exonStart, exonEnd):
        """Constructs new exon.

        Args:
            exonStart (int): Exon start coordinate
            exonEnd (int): Exon end coordinate

        """
        self.exonStart = exonStart
        self.exonEnd = exonEnd

class ChromExons:
    """Contains list of exons plus additional info from a row"""

    """Pattern to extract chromosome ID"""
    chromIdPattern = re.compile('^chr(\d+)')

    def __init__(self, chrom, strand, exons, name, name2):
        """
        Args:
            chrom (str): Chromosome ID
            strand (bool): '+' is True, '-' is False
            exons (list of Exon): Exons produced by processing one row
            name (str): NCBI transcript ID
            name2 (str): HGNC gene ID

        """
        self.chrom = chrom
        self.strand = strand
        self.exons = sorted(exons, key=lambda e: e.exonStart)
        self.name = name
        self.name2 = name2
        # This field will be used for sorting
        self.chromNumber = 0
        if chrom == 'chrY':
            # Y will be last in list
            self.chromNumber = 100
        elif chrom == 'chrX':
            # X will be penultimate in list
            self.chromNumber = 50
        elif chrom.startswith('chr'):
            # Others will be arranged by number follows 'chr'
            match = ChromExons.chromIdPattern.match(chrom)
            self.chromNumber = int(match.group(1)) if match else 9999
        else:
            # Any others will go to the end of a list
            self.chromNumber = 9999

def processNCBILine(row):
    """Function processes one line from input file.

    Args:
        row (dict of str): Dictionary created from line

    Returns:
        ChromExons

    """
    chrom = row['chrom']
    name = row['name']
    name2 = row['name2']
    # if len(e) > 0 check is needed because split(',') may return empty strings
    exonStarts = [int(e) for e in row['exonStarts'].split(',') if len(e) > 0]
    exonEnds = [int(e) for e in row['exonEnds'].split(',') if len(e) > 0]
    strand = row['strand']
    # Assuming that len(exonStarts) == len(exonEnds), need to add check later
    exons = [Exon(exonStarts[i], exonEnds[i]) for i in range(0, len(exonStarts))]
    return ChromExons(chrom, strand == '+', exons, name, name2)

def printChromExons(ce):
    """
    Args:
        ce (ChromExons): ChromExons to be printed
    """
    exonNums = range(1, len(ce.exons) + 1)
    if not ce.strand:
        exonNums = list(reversed(exonNums))
    idx = 0
    for ex in ce.exons:
        print '{}\t{}\t{}\t{}_exon-{}_{}'.format(ce.chrom, str(ex.exonStart),
                str(ex.exonEnd), ce.name2, exonNums[idx], ce.name)
        idx += 1

# Check whether we've been provided by input file name
if (len(sys.argv) < 2):
    print "Usage: python exons.py <NCBI_ref_file>"
    exit()

with open(sys.argv[1]) as csvfile:
    # Read input file
    reader = csv.DictReader(csvfile, delimiter='\t')
    # List of all ChromExons
    chromExonsList = []
    for row in reader:
        chromExonsList.append(processNCBILine(row))
    # Sort by chromosome ID
    chromExonsList.sort(key=lambda ch: ch.chromNumber)
    # Print
    for ch in chromExonsList:
        printChromExons(ch)


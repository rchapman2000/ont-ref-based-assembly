import sys
import os

import argparse as ap
from Bio import SeqIO

def main():

    # Create an argument parser and define arguments for the script.
    parser = ap.ArgumentParser()

    parser.add_argument('-i', '--input', required = True, type=str, \
        help='[Required] - path to a fasta file to calculate coverage of', \
        action = 'store', dest = 'inFile')

    args = parser.parse_args()

    # Check whether the file path provided exists.
    if not os.path.exists(args.inFile):
        # If the input file does not exist, notify the user and exit.
        sys.exit("ERROR: input file {0} does not exist. Please supply an existing file.".format(args.inFile))
    else:
        # If the input file does exist, we can begin parsing it.
        
        # Create empty variables to store the total sequence
        # length of all segments/chromosomes, the total number
        # of Ns in all segments/chromosomes, and the a list
        # of percent coverages of individual segments.
        totalLength = 0
        totalNs = 0
        segmentCovs = []

        # Open the input file as a handle.
        with open(args.inFile) as inputHandle:
            # Parse the handle using biopython and loop over each record.
            for record in SeqIO.parse(inputHandle, "fasta"):

                # Grab the length of the sequence
                seqLen = len(record.seq)

                # Count the number of times N appears in the sequence.
                Ns = record.seq.count("N")

                # Append the length of the sequence
                # and the number of Ns to the total 
                # counters
                totalLength += seqLen
                totalNs += Ns

                # Calculate the coverage of this segment by calculating the
                # percentage of N bases, subtracting this from 1, and multiplying by 100
                segmentCovs.append(round(((1 - (Ns / seqLen)) * 100 ), 4))

        # Calculate the total coverage across all segments by dividing the number
        # of Ns by the total length, subtracting this from one, and multiplying by 100
        totalCov = round(((1 - (totalNs / totalLength)) * 100), 4)

        # Check whether there are multiple segments.
        if len(segmentCovs) == 1:
            # If there is only 1 segment, the total coverage will be the same as the
            # segment coverage. Thus, we only need to print the total coverage.
            print(totalCov)
        else:
            # If there are multiple segments, print the coverage followed by a list
            # of the segment coverages separated by '; ' characters.
            print(str(totalCov) + " [{0}]".format("; ".join(str(cov) for cov in segmentCovs)))

if __name__ == "__main__":
    main()
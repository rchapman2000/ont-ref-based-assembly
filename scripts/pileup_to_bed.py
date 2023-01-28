import sys
import os
import argparse as ap

def main():

    # Create an argument parser and define arguments for the script.
    parser = ap.ArgumentParser()

    parser.add_argument('-i', '--input', required = True, type=str, \
        help='[Required] - path to pileup file', \
        action = 'store', dest = 'infile')
    parser.add_argument('-o', '--output', required = True, type=str, \
        help='[Required] - Output bed file', action='store', dest='outfile')
    parser.add_argument('--minCov', required=False, type=int, \
        help='The minimum coverage required for a position to be output.', \
        action='store', dest='minCov')

    args = parser.parse_args()

    # Create a variant to store the minimum depth of
    # coverage cutoff. The default is value is 0 (meaning
    # that all positions are reported)
    minCov = 0
    if (args.minCov):
        # If the user supplied a minimum coverage value.
        # set that.
        minCov = args.minCov

    f = ""
    # Check whether the input file supplied by the user exists
    if not os.path.exists(args.infile):
        # If not, notify the user and exit.
        sys.exit("ERROR: Input file {0} does not exist. Please supply an existing file.".format(args.infile))
    else:
        # If the file does exist, open it.
        f = open(args.infile, "r")
    
    
    # Open the output file and create it if it does not exist already.
    o = open(args.outfile, "w+")

    # Loop over each line in the file.
    for l in f:
        # Remove any newline characters and split the line 
        # the tab characters.
        split = l.strip().split("\t")
        
        # The third column contains information about the depth at a
        # given position. Check whether the current position is
        # at or above the minimum coverage threshold.
        if int(split[3]) >= minCov:
            # If the position has a depth at or above the minimum coverage
            # threshold, we will output this position to a bed file.
            
            # Grab the chromosome/segment from the first column of
            # the line.
            chrom = split[0]
            
            # Grab the position from the second column of the line.
            pos = int(split[1])

            # Write a line to the bed file for this position. The bed
            # format is:
            # 
            # CHROM	START	END
            #
            # where START is 0-based and END is 1-based.
            #
            o.write("{0}\t{1}\t{2}\n".format(chrom, pos - 1, pos))

    # Close the input and output file streams.
    f.close()
    o.close()

if __name__ == "__main__":
    main()
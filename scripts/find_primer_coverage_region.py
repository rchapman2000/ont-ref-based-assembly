import sys
import os

import argparse as ap

def main():
    
    # Create an argument parser and define arguments for the script.
    parser = ap.ArgumentParser()

    parser.add_argument('-i', '--input', required = True, type=str, \
        help='[Required] - path to primer bed file', \
        action = 'store', dest = 'inFile')
    parser.add_argument('-o', '--output', required = True, type=str, \
        help='[Required] - Output bed file', action='store', dest='outFile')

    args = parser.parse_args()


    # Handle the input file
    i = ""
    # Check whether the file path provided exists.
    if not os.path.exists(args.inFile):
        # If not, notify the user and exit the script.
        sys.exit("ERROR: Input file {0} does not exist. Please supply an existing file.".format(args.inFile))
    else:
        # If the file path does exist. Open it.
        i = open(args.inFile, "r")
    
    # Use the output file prefix provided to create and open
    # a file to store the primer coverage region.
    o = open(args.outFile + "-coverage-region.bed", "w+")

    # We want to determine the region coverage by amplicons for the primer
    # set. Thus, this region would span from the end of the first forward primer
    # to the start of the final reverse primer. By definition, we should never have a 
    # any primer before the first forward primer or a primer after the final reverse
    # primer. 
    # 
    # Thus, we can simply look at the coordinates of the primers
    # and choose the smallest end coordinate value (as the end of
    # the first forward primer) and the largest start position (as
    # the beginning of last reverse primer)
    #
    # Finally, we do not need to correct for 0 or 1-based coordinates
    # because we are switching end/start coordinates. For example, if the
    # first forward primer spans from 1 to 10, it will have a bed entry 
    # of "CHR	0	10", and the amplicon will begin at position 11 in a 1-based system
    # or position 10 in a 0-based system. Thus, to denote the start of hte amplicon region,
    # 10 would be correct start position (as the start coordinate is 0-based in bed format).
    # A similar logic is used for the final reverse primer.

    # Define a dictionary to store the start/end position of the region covered
    # for each chromosome/segment 
    chrPosDict = {}

    # Loop over each line in the primer coordinate bed file.
    for l in i:

        # Remove any newline characters and split the line by
        # tab characters
        split = l.strip("\n").split("\t")

        # Grab the Chromosome/segment (first item in the line),
        # as well as the start and end coordinates (second and third
        # items in the line respectively)
        chr = split[0]
        start = int(split[1])
        end = int(split[2])

        
        # If the chromosome/segment has not already been
        # added to the dictionary, add it and set the start position 
        # equal to the current line's end position and the end position
        # equal to the current line's start position. (As we are going from first 
        # primer's end to last primer's start!)
        if chr not in chrPosDict.keys():
            chrPosDict[chr] = [end,start]
        # If the chromosome/segment has already been added, 
        # then compare the current line's start/end position
        # to the stored positions.
        else:
            # If the current primer ends earlier than 
            # the currently saved start position, replace the
            # current start position with this primer's end position
            if end < chrPosDict[chr][0]:
                chrPosDict[chr][0] = end

            # If the current primer begins after the 
            # currently saved end position, replace the
            # current end position with the primer's
            # start position.
            if start > chrPosDict[chr][1]:
                chrPosDict[chr][1] = start

    # Loop over the dictionary and print
    # each entry in bed format:
    #
    # CHR	START	END
    #
    for chr in chrPosDict.keys():
        start = chrPosDict[chr][0]
        end = chrPosDict[chr][1]

        o.write("{0}\t{1}\t{2}\n".format(chr, start, end))

    # Close the input and output file streams.
    o.close()
    i.close()


if __name__ == "__main__":
    main()
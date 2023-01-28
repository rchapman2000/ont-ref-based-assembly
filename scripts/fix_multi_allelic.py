import sys
import os
import argparse as ap
import vcf

def main():
    
    # Create an argument parser and define arguments for the script.
    parser = ap.ArgumentParser()

    parser.add_argument('-i', '--input', required = True, type=str, \
        help='[Required] - Vcf file', \
        action = 'store', dest = 'infile')
    parser.add_argument('-o', '--output', required = True, type=str, \
        help='[Required] - Output vcf file name', action='store', dest='outfile')

    args = parser.parse_args()
    # Reads the VCF using a vcf reader.
    vcf_reader = vcf.Reader(open(args.infile, 'r'))
    # Opens a vcf output stream and passes the input vcf reader
    # to perserve the heading in the output file.
    vcf_writer = vcf.Writer(open(args.outfile, 'w+'), vcf_reader)

    # Loops over the records in the VCF
    for rec in vcf_reader:
        
        # Creates a copy of the record to write out
        toWrite = rec
        
        # If the record contains more than 1 alternative allele
        # (i.e. it is multiallelic), fix that site.
        if len(rec.ALT) > 1:

            # Find the most prevalent allele by looking at
            # the 'AO' INFO field, which gives
            # abundances of alleles, and choosing the 
            # allele with the largest abundance.
            mostPrevalentAlt = rec.INFO['AO'].index(max(rec.INFO['AO']))

            # Chance the record to be written's alternative 
            # allele to the most prevalent alternative.
            toWrite.ALT = [rec.ALT[mostPrevalentAlt]]

            # Loops over all of the INFO fields 
            # and fixes them to remove any information
            # about the filtered allele.
            newInfo = {}
            for x in rec.INFO.keys():
                
                if isinstance(rec.INFO[x], list) and len(rec.INFO[x]) > 1:
                    # If the INFO field contains more than 1 value,
                    # add only the information for the most prevalent 
                    # allele into the new information.
                    newInfo[x] = rec.INFO[x][mostPrevalentAlt]
                elif x == "NUMALT":
                    newInfo[x] = 1
                else:
                    # If not, addthe information field into the new
                    # information object.
                    newInfo[x] = rec.INFO[x]

            # Change the copy record's information to the 
            # altered info fields.
            toWrite.INFO = newInfo
        
        # Remove any genotyping information from the VCF
        # (As this information is likely changed from multiallelic sites,
        # and not used in my downstream analysis.)
        toWrite.samples = []

        vcf_writer.write_record(toWrite)

    vcf_writer.close()


if __name__ == "__main__":
    main()
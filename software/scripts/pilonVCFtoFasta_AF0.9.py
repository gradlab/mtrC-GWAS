#!/usr/bin/env python

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC,Gapped
import gzip


# check for correct number of inputs
if len(sys.argv) < 2 :
     print("Usage: pilonVCFtoFasta.py <input vcf 1> ... <input vcf n>")
     sys.exit(0)

def read_vcf(inFile):
    """Create strings corresponding to chromosome"""
    with gzip.open(inFile, 'rt') as vcf:
        for line in vcf:
            if ("contig" in line and "length" in line):
                line = line.strip()
                length = line.split("=")[-1].strip(">")
                print(length)
                refID = line.split(",")[0].split("=")[-1]
                chromosome = ["N"] * int(length)
            elif line[0] != "#":
                line = line.strip()
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO = line.split()[0:8]
                FILTERS = FILTER.split(";")
                INFO = INFO.split(";")
                if len(REF) == 1 and len(ALT) == 1:
                    if "Amb" in FILTERS:
                        ALLELE = "N"
                    elif "PASS" in FILTERS:
                        AF = float(INFO[-1].split("=")[-1])
                        if AF != 0 and AF < 0.9:
                            ALLELE = "N"
                        elif ALT == ".":
                            ALLELE = REF
                        else:
                            ALLELE = ALT
                    elif "Del" in FILTERS:
                            ALLELE = "-"
                    else:
                        if FILTER == "LowCov":
                            ALLELE = "N"
                        else:
                            print(FILTER)
                index = int(POS) - 1
                chromosome[index] = ALLELE
    return(chromosome, refID)

def write_fasta(chromosome, RGID, refID):
    """Writes a RGA fasta alignment for each vcf"""
    outFile = RGID + "_pseudogenome.fasta"
    Sample = RGID
    record = SeqRecord(Seq("".join(chromosome), Gapped(IUPAC.ambiguous_dna, '-')), id=Sample, description = "RGA_to_" + refID)
    SeqIO.write(record, outFile, "fasta")

for n in sys.argv[1:(len(sys.argv))] :
    RGID = os.path.basename(n)[0:-7].strip('_pilon')
    chromosome, refID = read_vcf(n)
    write_fasta(chromosome, RGID, refID)

#!/usr/bin/env python

import os,sys
import pandas as pd


def get_options():
    import argparse

    description = 'Script to run Pyseer unitig annotation pipeline followed by custom scripts'
    parser = argparse.ArgumentParser(description=description, prog='annotation_pipeline')

    parser.add_argument('resdir',
                        help='Input / output directory')
    parser.add_argument('resfile',
                        help='Raw unitigs results file from Pyseer GWAS (.results.txt)')
    return parser.parse_args()


def main():

    options = get_options()
    resdir = options.resdir
    resfile = options.resfile
    
    # Format results
    resdat = pd.read_csv('{0}/{1}'.format(resdir, resfile), sep='\t').sort_values('lrt-pvalue')
    resdat.to_csv('{0}/{1}'.format(resdir, resfile.replace('.txt', '.sorted.txt')), sep='\t', index=False)

    # Map unitigs to WHO_F genome (modified to have 1 copy of 23S rRNA) using Pyseer phandango_mapper
    os.system('phandango_mapper {0}/{1} gwas/reference/WHO_F_unique23s.fasta {0}/{2}'
        .format(resdir, resfile.replace('.txt', '.sorted.txt'), resfile.replace('.txt', '.mapped-WHO_F.plot')))

    # Record Bonferroni threshold
    os.system('python software/scripts/pyseer/count_patterns.py {0}/unitig_patterns.txt > {0}/threshold.txt'.format(resdir))
    with open('{0}/threshold.txt'.format(resdir), 'r') as infile:
        for line in infile:
            if line.startswith('Threshold:'):
                threshold = float(line.split()[1])
    resdat[resdat['lrt-pvalue'] <= threshold].to_csv('{0}/{1}'
        .format(resdir, resfile.replace('.txt', '.sorted.significant.txt')), sep='\t', index=False)

    # Annotate significant unitigs using Pyseer pipeline
    os.system('annotate_hits_pyseer {0}/{1} gwas/reference/reference_WHO_F.txt {0}/{2}'
        .format(resdir, resfile.replace('.txt', '.sorted.significant.txt'), resfile.replace('.txt', '.significant.annotated.txt')))
    os.system('rm remaining_kmers*; rm tmp_bed')

    # Further summarization and filtering of Pyseer annotations for final tables
    os.system('python software/scripts/summarize_top_unitigs.py {0} {1}'.format(resdir, resfile.replace('.txt', '.significant.annotated.txt')))


if __name__ == "__main__":
    main()

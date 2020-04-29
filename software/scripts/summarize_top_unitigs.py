#!/usr/bin/env python

import os,sys
import pandas as pd


def get_options():
    import argparse

    description = 'Script to filter and summarise (by grouping in sliding windows) annotated unitigs from Pyseer'
    parser = argparse.ArgumentParser(description=description, prog='summarise_unitigs')

    parser.add_argument('resdir',
                        help='Input / output directory')
    parser.add_argument('resfile',
                        help='Significant annotated unitigs file from Pyseer (.significant.annotated.txt)')
    parser.add_argument("--af-filter",
                        help="Allele frequency filter cutoff (by default keep minor alleles) "
                        "[default=0.5]",
                        default=0.5)
    parser.add_argument("--len-filter",
                        help="Uniting length filter (discard unitigs that are too short) "
                        "[default=20]",
                        default=20)
    parser.add_argument("--window",
                        help="Window size around each unitig to cluster (50 bp = 25 bp from each end) "
                        "[default=50bp]",
                        default=50)
    return parser.parse_args()

 
def main():
    
    options = get_options()
    resdir = options.resdir
    resfile = options.resfile
    n = options.window/2
    len_filter = options.len_filter
    af_filter = options.af_filter
    
    unitigs = []
    top_hits = []

    with open('{0}/{1}'.format(resdir, resfile), 'r') as res:
        for line in res:
            if len(line.rstrip().split('\t')[7].split(';')) != 4: # This removes unitigs that mapped to multiple places in the genome
                continue
            effect = line.rstrip().split('\t')[4]
            af = line.rstrip().split('\t')[1]
            if float(af) > af_filter:
                continue # Skip unitigs that are above the allele frequency filter
            seq = line.rstrip().split('\t')[0]
            if len(seq) < len_filter: 
                continue # Skip short unitigs
            unitigs.append(line.rstrip())

    while unitigs:

        # Starting with most significant unitig, cluster remaining unitigs that fall within specified window
        # Iterate with the next remaining most significant unitig

        hit = unitigs[0]
        count = 0
        l,r = hit.split('\t')[7].split(';')[0].split(':')[1].split('-') # Pull out position that unitig mapped to
        l = int(l)
        r = int(r)
        if l > r:
            sys.exit('Error in interval order')
        new_unitigs = []

        for i,unitig in enumerate(unitigs):
            nl,nr = unitig.split('\t')[7].split(';')[0].split(':')[1].split('-')
            nl = int(nl)
            nr = int(nr)
            if nl > l-n and nl < r+n: # Check to see if either side falls within window
                count += 1
                continue
            elif nr > l-n and nr < r+n:
                count += 1
                continue
            else:
                new_unitigs.append(unitig)
        top_hits.append(hit)
        unitigs = new_unitigs

    with open('{0}/{1}'.format(resdir, resfile.replace('.txt', '.topkmers.txt')), 'w') as out:
        out.write('Gene\tP-value\tFrequency\tBeta\tBeta Std Error\tHeritability\tPosition\tUnitig\n')
        for unitig in top_hits:
            unitig,af,skip,pval,beta,betase,h2,annot = unitig.split('\t')
            pos = annot.split(';')[0].split(':')[1]
            gene = annot.split(';')[2] # Take gene name if unitig maps to a gene
            if gene.startswith('piiC') or gene.startswith('pilE'): # Filter known repetitive genes by annotation
                continue
            if gene == '': # Otherwise, treat this as an intergenic region
                gene = '{0}...{1}'.format(annot.split(';')[1], annot.split(';')[3])
                if 'piiC' in gene or 'pilE' in gene: # Filter known repetitive genes by annotation
                    continue
            gene = gene.replace('mexA', 'mtrC') # Adjust mtrRCDE annotation names
            gene = gene.replace('mexB', 'mtrD')
            gene = gene.replace('acrR', 'mtrR')
            out.write('\t'.join([gene,pval,af,beta,betase,h2,pos,unitig])+'\n')

    # Adjust number of digits to display
    df = pd.read_csv('{0}/{1}'.format(resdir, resfile.replace('.txt', '.topkmers.txt')), sep='\t').sort_values('P-value')
    df['P-value'] = df['P-value'].map('{:0.2e}'.format)
    df['Frequency'] = df['Frequency'].map('{:0.3f}'.format)
    df['Heritability'] = df['Heritability'].map('{:0.3f}'.format)
    df['CI'] = df.apply( lambda row: (round(row['Beta'] - 1.96*row['Beta Std Error'], 3), round(row['Beta'] + 1.96*row['Beta Std Error'], 3)), axis=1)
    df['Beta'] = df['Beta'].map('{:0.2f}'.format)
    df['Beta Std Error'] = df['Beta Std Error'].map('{:0.3f}'.format)
    df = df[[c for c in df if c not in ['Unitig']] + ['Unitig']] # Quick way to resposition
    df.to_csv('{0}/{1}'.format(resdir, resfile.replace('.txt', '.topkmers.txt')), sep='\t', index=False)
    os.system("column -t -s $'\\t' {0}/{1} > {0}/{2}".format(resdir, resfile.replace('.txt', '.topkmers.txt'), resfile.replace('.txt', '.topkmers.aligned.txt')))


if __name__ == "__main__":
    main()

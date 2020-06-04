#!/usr/bin/env python

import sys
import argparse
from datetime import datetime

def get_args():
    parser = argparse.ArgumentParser(description="Create resistance allele matrix and itol files")
    parser.add_argument("vcf", help="VCF of SNPs in pseudogenomes mapped to NCCP11945")
    return parser.parse_args()

color_palettes = {"bugn_2":["#b2e2e2", "#238b45"],
                  "bupu_2":["#b3cde3", "#88419d"],
                  "gnbu_2":["#bae4bc", "#2b8cbe"],
                  "orrd_2":["#fdcc8a", "#d7301f"],
                  "pubu_2":["#bdc9e1", "#0570b0"],
                  "pubugn_2":["#bdc9e1", "#02818a"],
                  "purd_2":["#d7b5d8", "#ce1256"],
                  "rdpu_2":["#fbb4b9", "#ae017e"],
                  "ylgn_2":["#c2e699", "#238443"],
                  "ylgnbu_2":["#a1dab4", "#225ea8"],
                  "ylorbr_2":["#fed98e", "#cc4c02"],
                  "ylorrd_2":["#fecc5c", "#e31a1c"],
                  "bugn_3":["#ccece6", "#66c2a4", "#006d2c"],
                  "bupu_3":["#bfd3e6", "#8c96c6", "#810f7c"],
                  "gnbu_3":["#ccebc5", "#7bccc4", "#0868ac"],
                  "orrd_3":["#fdd49e", "#fc8d59", "#b30000"],
                  "pubu_3":["#d0d1e6", "#74a9cf", "#045a8d"],
                  "purd_3":["#d4b9da", "#df65b0", "#980043"],
                  "rdpu_3":["#fcc5c0", "#f768a1", "#7a0177"],
                  "ylgn_3":["#d9f0a3", "#78c679", "#006837"],
                  "ylgnbu_6":["#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84"],
                  "pubugn_6":["#d0d1e6", "#a6bddb", "#67a9cf", "#3690c0", "#02818a", "#016450"],
                  "ylgnbu_7":["#edf8b1", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#253494", "#081d58"],
                  "pubugn_8":["#fff7fb", "#ece2f0", "#d0d1e6", "#a6bddb", "#67a9cf", "#3690c0", "#02818a", "#016450"]}

def read_vcf(vcf_file):
    res_positions = ["109623", "1524644", "1524645", "1524522", "1524495",
                     "1524494", "2050038", "2050039", "2050040", "2050037",
                     "2050036", "1332981", "1333000", "1051666", "1051667",
                     "1051678", "1051679", "195711", "195708", "195707",
                     "195696", "195695", "2031479", "1307220", "2033259",
                     "2033260"]

    res_alleles = {}
    with open(vcf_file, "r") as infile:
        for line in infile:
            if line[0:2] == "##":
                continue
            elif line[0] == "#":
                line = line.strip().split()
                strains = line[10:]
            else:
                line = line.strip().split()
                POS = line[1]
                if POS in res_positions:
                    nucleotides = [line[3]] + line[4].split(",")
                    alleles = [nucleotides[int(i)] for i in line[10:]]
                    res_alleles[POS] = alleles
    return(strains, res_alleles)

def single_nt_change(locus_name, position, ref_nt, alt_nt, ref_aa, alt_aa, color):
    try:
        # check that position is polymorphic, if not assign reference allele
        # to all strains
        nt = res_alleles[position]
    except KeyError:
        return(dict(zip(strains, [ref_aa]*len(strains))))
    aa = []
    for n in nt:
        if n == "*":
            aa.append("NA")
        elif n == ref_nt:
            aa.append(ref_aa)
        elif n == alt_nt:
            aa.append(alt_aa)
        else:
            print("{0} has multiple alternate alleles in this dataset. Please edit code.")
            sys.exit(0)
    locus_dict = dict(zip(strains,aa))
    with open("itol_{0}.txt".format(locus_name.replace(" ", "_")), "w") as itol_file:
        itol_file.write("DATASET_COLORSTRIP\n\n")
        itol_file.write("SEPARATOR TAB\n\n")
        itol_file.write("DATASET_LABEL\t{0}\nCOLOR\t{1}\n\n".format(locus_name, color[-1]))
        itol_file.write("LEGEND_TITLE\t{0}\nLEGEND_SHAPES\t1\t1\t1\nLEGEND_COLORS\t{1}\t#D3D3D3\nLEGEND_LABELS\t{2}\t{3}\tunknown\n\n".format(locus_name, "\t".join(color), ref_aa, alt_aa))
        itol_file.write("BORDER_WIDTH\t0.25\nBORDER_COLOR\t#CCCCCC\n\n")
        itol_file.write("DATA\n")
        for strain,allele in locus_dict.items():
            if allele == ref_aa:
                c = color[0]
            elif allele == alt_aa:
                c = color[1]
            else:
                c = "#D3D3D3"
            itol_file.write("{0}\t{1}\t{2}\n".format(strain, c, allele))

    return(locus_dict)

def multiple_nt_change(locus_name, position, ref_nt, alt_nt, ref_aa, alt_aa, color):
    if len(set(position) - set(res_alleles.keys())) == 0:
        # check that both positions were actually polymorphic
        nt = ["".join(x) for x in list(zip(*[res_alleles[p] for p in position]))]
    else:
        nt_list = []
        for i,p in enumerate(position):
           if p in res_alleles:
               nt_list.append(res_alleles[p])
           else:
               nt_list.append([ref_nt[0][i]]*len(res_alleles[list(res_alleles.keys())[0]]))
        nt = ["".join(x) for x in list(zip(*nt_list))]
    aa = []
    for n in nt:
        if "*" in n:
            aa.append("NA")
        elif n in ref_nt:
            aa.append(ref_aa)
        elif n in alt_nt:
            aa.append(alt_aa[alt_nt.index(n)])
        else:
            print("{0} is not a listed allele for {1}. Please edit code.".format(n, locus_name))
            sys.exit(0)
    locus_dict = dict(zip(strains,aa))
    with open("itol_{0}.txt".format(locus_name.replace(" ","_")), "w") as itol_file:
        itol_file.write("DATASET_COLORSTRIP\n\n")
        itol_file.write("SEPARATOR TAB\n\n")
        itol_file.write("DATASET_LABEL\t{0}\nCOLOR\t{1}\n\n".format(locus_name, color[-1]))
        itol_file.write("LEGEND_TITLE\t{0}\nLEGEND_SHAPES\t{1}\nLEGEND_COLORS\t{2}\t#D3D3D3\nLEGEND_LABELS\t{3}\t{4}\tunknown\n\n".format(locus_name,
                                                                                                                                          "\t".join(['1']*(len(set(alt_aa)) + 2)),
                                                                                                                                          "\t".join(color),
                                                                                                                                          ref_aa,
                                                                                                                                          "\t".join(list(dict.fromkeys(alt_aa)))))
        itol_file.write("BORDER_WIDTH\t0.25\nBORDER_COLOR\t#CCCCCC\n\n")
        itol_file.write("DATA\n")
        for strain,allele in locus_dict.items():
            if allele == ref_aa:
                c = color[0]
            elif allele in alt_aa:
                c = color[list(dict.fromkeys(alt_aa)).index(allele) + 1]
            else:
                c = "#D3D3D3"
            itol_file.write("{0}\t{1}\t{2}\n".format(strain, c, allele))

    return(locus_dict)

args = get_args()


strains, res_alleles = read_vcf(args.vcf)

res_loci_dict = {}
res_loci_dict["PBP1_421"] = single_nt_change("PBP1 421", "109623", "T", "C", "L", "P",  color_palettes["bugn_2"])
res_loci_dict["PBP2_501"] = multiple_nt_change("PBP2 501", ["1524645", "1524644"], ["CG"], ["TG", "CA"], "A", ["T", "V"], color_palettes["bupu_3"])
res_loci_dict["PBP2_542"] = single_nt_change("PBP2 542", "1524522", "C", "T", "G", "S", color_palettes["gnbu_2"])
res_loci_dict["PBP2_551"] = multiple_nt_change("PBP2 551", ["1524495", "1524494"], ["GG"], ["GA", "AG", "TG"], "P", ["L", "S", "T"], color_palettes["ylgnbu_6"])
res_loci_dict["PorB_120"] = multiple_nt_change("PorB 120",
                                               ["2050040", "2050039", "2050038"],
                                               ["CCG", "CCA"],
                                               ["TTC","TCC","CTG","TTG","TTA","GTC", "TGC", "CTC"],
                                               "G",
                                               ["K","R","D","N","N","Q","C", "E"],
                                               color_palettes["pubugn_8"])
res_loci_dict["PorB_121"] = multiple_nt_change("PorB 121",
                                               ["2050037", "2050036"],
                                               ["CG"],
                                               ["CC","CA","CT","TT"],
                                               "A",
                                               ["G","V","D","N"],
                                               color_palettes["pubugn_6"][:-1])
res_loci_dict["MtrR_39"] = single_nt_change("MtrR 39", "1332981", "G", "A", "A", "T", color_palettes["pubu_2"])
res_loci_dict["MtrR_45"] = single_nt_change("MtrR 45", "1333000", "G", "A", "G", "D", color_palettes["purd_2"])
res_loci_dict["GyrA_91"] = multiple_nt_change("GyrA 91",
                                              ["1051666", "1051667"],
                                              ["TC"],
                                              ["TA", "AC", "TT"],
                                              "S",
                                              ["Y", "T", "F"],
                                              color_palettes["ylgnbu_6"][:-2])
res_loci_dict["GyrA_95"] = multiple_nt_change("GyrA 95",
                                              ["1051678", "1051679"],
                                              ["GA"],
                                              ["AA", "TA", "GG", "GC"],
                                              "D",
                                              ["N", "Y", "G", "A"],
                                              color_palettes["pubugn_6"][:-1])
res_loci_dict["ParC_86"] = single_nt_change("ParC 86", "195711", "C", "T", "D", "N", color_palettes["ylgnbu_2"])
res_loci_dict["ParC_87"] = multiple_nt_change("ParC 87",
                                              ["195708", "195707"],
                                              ["TC"],
                                              ["TT", "GC", "TA", "AC", "AT"],
                                              "S",
                                              ["N", "R", "I", "C", "Y"],
                                              color_palettes["pubugn_6"])
res_loci_dict["ParC_91"] = multiple_nt_change("ParC 91",
                                              ["195696", "195695"],
                                              ["CT"],
                                              ["TT", "CC", "CG", "GT"],
                                              "E",
                                              ["K", "G", "A", "Q"],
                                              color_palettes["ylgnbu_6"][:-1])
res_loci_dict["RpsJ_57"] = single_nt_change("RpsJ 57", "2031479", "G", "A", "V", "M", color_palettes["ylorbr_2"])
res_loci_dict["FolP_229"] = multiple_nt_change("FolP 229",
                                               ["1307220"],
                                               ["G","A"],
                                               ["T"],
                                               "R",
                                               ["S"],
                                               color_palettes["ylorrd_2"])
res_loci_dict["RplD_70"] = multiple_nt_change("RplD_70",
                                              ["2033259", "2033260"],
                                              ["GG"],
                                              ["GC","CG","AG","GA"],
                                              "G",
                                              ["A","R","S","D"],
                                              color_palettes["pubugn_6"])


with open('{0}_gc_resistance_alleles.tsv'.format(datetime.strftime(datetime.now(), '%Y-%m-%d')), "w") as outfile:
    loci = list(res_loci_dict.keys())
    outfile.write("Accession\t" + "\t".join(loci) + "\n")
    for s in strains:
        a = [res_loci_dict[l][s] for l in loci]
        outfile.write("{0}\t{1}\n".format(s, "\t".join(a)))

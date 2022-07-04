#!/usr/bin/python

# Usage: merge_vcf.py vcf_file1.vcf vcf_file2.vcf output.vcf
# Author: Veerendra Gadekar (gpveerendra09@gmail.com)

import csv
import sys

def read_vcf(vcf_file):
    """
    Read a vcf file to list of dict.
    :argument - str vcf_file: Path to a vcf file.
    """
    vcf_dict = []
    with open(vcf_file, 'r') as invcf:
        for line in invcf:
            if line.startswith('#'):
                continue
            line = line.strip().split()
            vcf_dict.append({'CHROMPOS' : line[0]+':'+line[1], 'CHROM':line[0], 'POS' : line[1], 'ID' : line[2], 'REF': line[3], 
                             'ALT':line[4], 'QUAL':line[5], 'FILTER':line[6], 'INFO':line[7], 
                             'FORMAT':line[8], 'Sample1':line[9]})
    return vcf_dict 
    

def map_snp_loc(dict1, dict2):
    ''' check if dict elements map'''
    return dict1['CHROMPOS'] in dict2

def get_vcf_rows(x, y):
    vcf_tup  = tuple(x)
    rows_tup = tuple(y)
    out = [z for i,z in enumerate(x) if y[i]]
    return(out)

def edit_vcf_INFO(x, y):
    x['INFO'] = x['INFO'] + ';' + y
    x.pop('CHROMPOS', None)
    return(x)

if __name__ == "__main__":

    vcf1 = sys.argv[1]
    vcf2 = sys.argv[2]
    output = sys.argv[3]
    
    freebayes = read_vcf(vcf1)
    varscan   = read_vcf(vcf2)    

    freebayes_CHROMPOS = set(x['CHROMPOS'] for x in freebayes)
    varscan_CHROMPOS   = set(x['CHROMPOS'] for x in varscan)    

    # unique and common snp position between the vsf files
    freebayes_only = freebayes_CHROMPOS - varscan_CHROMPOS
    varscan_only   = varscan_CHROMPOS   - freebayes_CHROMPOS
    common         = varscan_CHROMPOS   & freebayes_CHROMPOS

    # get the rows unique and common to two vcf files
    freebayes_snp_loc_common = list(map(lambda p: map_snp_loc(p, common), freebayes))
    varscan_snp_loc_common   = list(map(lambda p: map_snp_loc(p, common), varscan))
    freebayes_snp_loc_unique = list(map(lambda p: map_snp_loc(p, freebayes_only), freebayes))
    varscan_snp_loc_unique   = list(map(lambda p: map_snp_loc(p, varscan_only), varscan))

    fv_common_vcf_dict      = get_vcf_rows(freebayes, freebayes_snp_loc_common)
    vf_common_vcf_dict      = get_vcf_rows(varscan, varscan_snp_loc_common)
    freebayes_only_vcf_dict = get_vcf_rows(freebayes, freebayes_snp_loc_unique)
    varscan_only_vcf_dict   = get_vcf_rows(varscan, varscan_snp_loc_unique)

    # edit the info tags
    fv_common_vcf_dict        = list(map(lambda x: edit_vcf_INFO(x, "calledBy=Freebayes+VarScan"), fv_common_vcf_dict))
    vf_common_vcf_dict        = list(map(lambda x: edit_vcf_INFO(x, "calledBy=VarScan+Freebayes"), vf_common_vcf_dict))
    freebayes_only_vcf_edited = list(map(lambda x: edit_vcf_INFO(x, "calledBy=Freebayes"), freebayes_only_vcf_dict))
    varscan_only_vcf_edited   = list(map(lambda x: edit_vcf_INFO(x, "calledBy=varscan"), varscan_only_vcf_dict))

    # put all the edited dict together
    main = varscan_only_vcf_edited + freebayes_only_vcf_edited + vf_common_vcf_dict


    # header for merged vcf file
    header = """##fileformat=VCF
    ##fileDate=20220704
    ##source=merge_vcf.py
    #CHROM POS ID REF ALT QUAL FILTER INFO
    """

    output_VCF = output
    with open(output_VCF, 'w') as vcf:
        vcf.write(header)
        keys = main[0].keys()
        dict_writer = csv.DictWriter(vcf, keys, delimiter='\t')
        dict_writer.writeheader()
        dict_writer.writerows(main)


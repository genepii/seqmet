#!/usr/bin/env python3
import argparse

##Count the number of minor variants in a target vcf reported as major variant in a reference vcf
#v0.0.8

def list_intersect(containedl, containingl, mode, threshold):
    if mode == 'match':
        return [ x for x in containedl if x in containingl ]
    elif mode == 'missing':
        return [ x for x in containedl if x not in containingl ]

def count_commented(file):
    lines = open(file, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')
    count = 0
    for line in lines:
        if line[0] == "#":
            count += 1
    return count

def list_flatten(toflat_list):
    toflat_list = [item for ilist in toflat_list for item in ilist]
    return toflat_list

def list_bgapos(bga):
    bgapos = []
    for i in range(len(bga)):
        if bga[i][0] not in [ x[0] for x in bgapos ]:
            bgapos.append([bga[i][0], []])
            loop_index = [ x[0] for x in bgapos ].index(bga[i][0])
        else:
            loop_index = [ x[0] for x in bgapos ].index(bga[i][0])
        bgapos[loop_index][1].append([int(bga[i][3]) for x in range(int(bga[i][1]),int(bga[i][2]))])
    return bgapos

def list_bedpos(bed):
    bedpos = []
    for i in range(len(bed)):
        if bed[i][0] not in [ x[0] for x in bedpos ]:
            bedpos.append([bed[i][0], []])
            loop_index = [ x[0] for x in bedpos ].index(bed[i][0])
        else:
            loop_index = [ x[0] for x in bedpos ].index(bed[i][0])
        bedpos[loop_index][1].append([int(x) for x in range(int(bed[i][1]),int(bed[i][2]))])
    return bedpos

def vcf_parse(vcf, fields):
    if len(vcf) == 0:
        return []
    fields_header = [ x.split('=')[0] for x in vcf[0].split('\t')[7].split(';') if '=' in x ]
    fields_index = [ fields_header.index(x) for x in fields.split(',') ]
    var = [ [x.split('\t')[0], int(x.split('\t')[1])-1, x.split('\t')[3], x.split('\t')[4]] for x in vcf]
    for i in range(len(var)):
        for j in range(len(fields_index)):
            var[i].append(vcf[i].split('\t')[7].split(';')[fields_index[j]].split('=')[1])
        if var[i][3][0] == '-':
            var_temp = var[i][2]
            var[i][2] = var[i][2] + var[i][3][1:]
            var[i][3] = var_temp
        if var[i][3][0] == '+':
            var[i][3] = var[i][2] + var[i][3][1:]
    return var

def value_median_index(value_list, lowest_index, highest_index):
    if len(value_list) < 3:
        return 0
    nvalue = highest_index - lowest_index + 1
    nvalue = (nvalue + 1) // 2 - 1
    return nvalue + lowest_index
 
def value_IQR(value_list, lowest_index, highest_index):
    if len(value_list) < 3:
        return 0
    median_index = value_median_index(value_list, lowest_index, highest_index)
    Q1 = value_list[value_median_index(value_list, lowest_index, median_index)]
    Q3 = value_list[value_median_index(value_list, median_index + 1, highest_index)]
    return (Q3 - Q1)

parser = argparse.ArgumentParser(description='Count the number of minor variants in a target vcf reported as major variant in a reference vcf')
debugmode = parser.add_mutually_exclusive_group()
debugmode.add_argument('-v', '--verbose', action='store_true')
debugmode.add_argument('-q', '--quiet', action='store_true')
parser.add_argument('--version', action='version', version='0.0.6')
parser.add_argument('-t', '--target', help='the base')
parser.add_argument('-r', '--reference', help='the base')
parser.add_argument('-e', '--exclusion', help='the base')
parser.add_argument('-m', '--mode', help='the base')
parser.add_argument('-d', '--depth', help='the base')
parser.add_argument('-b', '--bed', help='the base')
parser.add_argument('--min_depth', type=int, default=100, help='the base')
parser.add_argument('--min_freq', type=float, default=0.05, help='the base')
parser.add_argument('-o', '--output', help='output file', default='./')

if __name__ == '__main__':
    args = parser.parse_args()

tvcf = open(args.target, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')[count_commented(args.target):]
rvcf = open(args.reference, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')[count_commented(args.reference):]
bga = [x.split('\t') for x in open(args.depth, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')]
bed = [x.split('\t') for x in open(args.bed, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')]

depth = [ [x[0], list_flatten(x[1])] for x in list_bgapos(bga) ]
region = [ [x[0], list_flatten(x[1])] for x in list_bedpos(bed) ]

#chrom, pos, ref, alt, af, dp
tvar = vcf_parse(tvcf, 'AF')

rvar = vcf_parse(rvcf, 'AF')

depth_chrom = [y[0] for y in depth]
region_chrom = [y[0] for y in region]

tvar_major = [ x for x in tvar if float(x[4])>=0.5 and depth[depth_chrom.index(x[0])][1][x[1]]>args.min_depth and x[1] in region[region_chrom.index(x[0])][1] ]
tvar_minor = [ x for x in tvar if float(x[4])<0.5 and float(x[4])>args.min_freq and depth[depth_chrom.index(x[0])][1][x[1]]>args.min_depth and x[1] in region[region_chrom.index(x[0])][1] ]
rvar_major = [ x for x in rvar if float(x[4])>=0.5 and depth[depth_chrom.index(x[0])][1][x[1]]>args.min_depth ]

tvar_major_profile = [ [ x[0], x[1], x[2], x[3] ] for x in tvar_major ]
tvar_major_profile_reversed = [ [ x[0], x[1], x[3], x[2] ] for x in tvar_major if len(x[2]) == 1 and len(x[3]) == 1]
tvar_minor_profile = [ [ x[0], x[1], x[2], x[3] ] for x in tvar_minor ]
rvar_major_profile = [ [ x[0], x[1], x[2], x[3] ] for x in rvar_major ]
rvar_major_profile_reversed = [ [ x[0], x[1], x[3], x[2] ] for x in rvar_major if len(x[2]) == 1 and len(x[3]) == 1]

expected = []
common = []

if args.mode == 'raw':
    expected += list_intersect( rvar_major_profile, tvar_major_profile, 'missing', 100 )
    expected += list_intersect( tvar_major_profile_reversed, rvar_major_profile_reversed , 'missing', 100 )
if args.mode == 'min':
    expected += list_intersect( rvar_major_profile, tvar_major_profile, 'missing', 100 )
elif args.mode == 'maj':
    expected += rvar_major_profile
    expected += list_intersect( tvar_major_profile_reversed, rvar_major_profile_reversed , 'missing', 100 )

if args.exclusion:
    evcf = open(args.exclusion, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')[count_commented(args.exclusion):]
    evar = vcf_parse(evcf, 'AF')
    evar_major = [ x for x in evar if float(x[4])>=0.5 ]
    evar_major_profile = [ [ x[0], x[1], x[2], x[3] ] for x in evar_major ]
    evar_major_profile_reversed = [ [ x[0], x[1], x[3], x[2] ] for x in evar_major if len(x[2]) == 1 and len(x[3]) == 1]
    expected = list_intersect( expected, evar_major_profile , 'missing', 100 )
    expected = list_intersect( expected, evar_major_profile_reversed , 'missing', 100 )


if args.mode == 'min' or args.mode == 'raw':
    common += list_intersect( tvar_minor_profile , expected, 'match', 100 )
elif args.mode == 'maj':
    common += list_intersect( tvar_minor_profile , expected, 'match', 100 )
    common += list_intersect( tvar_major_profile , expected, 'match', 100 )

common = [list(x) for x in set(tuple(x) for x in common)]
common.sort()

tvar_profile_full = tvar_major_profile + tvar_major_profile_reversed + tvar_minor_profile
tvar_profile_full_af = [ x[4] for x in tvar_major ] + [ x[4] for x in tvar_major if len(x[2]) == 1 and len(x[3]) == 1] + [ x[4] for x in tvar_minor ]
common_af = [ int(float(tvar_profile_full_af[tvar_profile_full.index(x)])*100) for x in common ]

common_length = len(common)
expected_length = len(expected)
common_af.sort()
median = [ common_af[x] if x!=0 else 0 for x in [value_median_index(common_af, 0, len(common_af))] ][0]
IQR = value_IQR(common_af, 0, len(common_af))

refname = '.'.join(args.reference.split('/')[-1].split('.')[:-1])
w = open(args.output, 'a+')

if args.mode == 'raw':
    w.write(refname + "\t" + str(common_length) + "/" + str(expected_length) + "\n")
if args.mode == 'min' or args.mode == 'maj':
    if expected_length > 0:
        w.write(refname + "\t" + str(common) + "\t" + str(expected) + "\t" + str(common_length) + "\t" + str(expected_length) + "\t" + str(round(float(common_length)/float(expected_length), 3)) + "\t" + str(median) + "\t" + str(IQR) + "\n")
    else:
        w.write(refname + "\t" + str(common) + "\t" + str(expected) + "\t" + str(common_length) + "\t" + str(expected_length) + "\t0.000\t" + str(median) + "\t" + str(IQR) + "\n")

w.close()

from __future__ import print_function
import os
import sys
import getopt

##Recalculate scores and metrics of vcf files generated with freebayes
#v0.0.1

def main(argv):
    global inp
    global oup
    global ctype
    global minfreq
    global maxfreq
    inp = ''
    oup = ''
    ctype = ['snp','del','ins', 'mnp', 'complex']
    minfreq = 0.0
    maxfreq = 1.0
    try:
        opts, args = getopt.getopt(argv, 'hi:o:t:n:x:', ['help', 'input', 'output', 'calltype', 'minfreq', 'maxfreq'])
        for opt, arg in opts:
            if opt in ('-h', '--help'):
                usage()
                sys.exit()
            elif opt in ('-i', '--input'):
                inp = arg
            elif opt in ('-o', '--output'):
                oup = arg
            elif opt in ('-t', '--calltype'):
                ctype = []
                for i in range(len(arg.split(','))):
                    ctype.append(arg.split(',')[i])
            elif opt in ('-n', '--minfreq'):
                minfreq = arg
            elif opt in ('-x', '--maxfreq'):
                maxfreq = arg
        if inp == '':
            usage()
            sys.exit()
        if oup == '':
            oup = inp.split("/")[-1].split(".")[0] + '_fbvcf.' + inp.split("/")[-1].split(".")[1]
    except getopt.GetoptError:
        usage()
        sys.exit(2)

def usage():
    print('usage: ' + sys.argv[0] + ' -h --help --vcf [vcf] -o --output [tsv] -t --calltype [snp,del,ins,mnp,complex] -n --minfreq [float] -x --maxfreq [float]')

if __name__ == '__main__':
    main(sys.argv[1:])

def count_commented(file):
    lines = open(file, 'r').read().rstrip('\n').split('\n')
    count = 0
    for line in lines:
        if line[0] == "#":
            count += 1
    return count

def truncate(f, n):
    '''Truncates/pads a float f to n decimal places without rounding'''
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    return '.'.join([i, (d+'0'*n)[:n]])

print(ctype)

flatten = lambda t: [item for sublist in t for item in sublist]

vcf = open(inp, 'r').read().rstrip('\n').split('\n')[count_commented(inp):]
comment = open(inp, 'r').read().rstrip('\n').split('\n')[:count_commented(inp)]

for i in range(len(vcf)):
    vcf[i] = vcf[i].split('\t')
    vcf[i][7] = vcf[i][7].split(';')

for i in range(len(vcf)):
    vcf_info_name = [x.split('=')[0] for x in vcf[i][7]]
    vcf_info_value = [x.split('=')[1] for x in vcf[i][7]]
    
    if vcf_info_name[3] != 'AF' and vcf_info_name[5] != 'AO' and vcf_info_name[7] != 'DP' and vcf_info_name[22] != 'PAO' and vcf_info_name[25] != 'PRO' and vcf_info_name[28] != 'RO' and vcf_info_name[40] != 'TYPE':
        sys.exit('Fields wrongly formatted')
    if vcf_info_value[40] not in ctype:
        vcf_info_value[3] = str(truncate(float(int(vcf_info_value[5])) / float(vcf_info_value[7]), 3))
        
        new_field = []
        for y in range(len(vcf_info_name)):
            new_field.append(vcf_info_name[y] + '=' + vcf_info_value[y])
        vcf[i][7] = ';'.join(new_field)
        continue
    
    if truncate(float(int(vcf_info_value[5])) / float(vcf_info_value[7]), 3) < truncate(minfreq, 3):
        vcf_info_value[3] = truncate(minfreq, 3)
    elif truncate(float(int(vcf_info_value[5])) / float(vcf_info_value[7]), 3) > truncate(maxfreq, 3):
        vcf_info_value[3] = truncate(maxfreq, 3)
    else:
        vcf_info_value[3] = str(truncate(float(int(vcf_info_value[5])) / float(vcf_info_value[7]), 3))
    
    new_field = []
    for y in range(len(vcf_info_name)):
        new_field.append(vcf_info_name[y] + '=' + vcf_info_value[y])
    
    vcf[i][7] = ';'.join(new_field)

w = open(oup, 'w')

for i in range(len(comment)):
    w.write(comment[i] + '\n')

for i in range(len(vcf)):
    w.write('\t'.join(vcf[i]) + '\n')

w.close()

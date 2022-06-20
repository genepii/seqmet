from __future__ import print_function
import os
import sys
import getopt

##Mask lower-case nucleotides in a fasta file
#v0.0.1

def main(argv):
    global inp
    global oup
    inp = ''
    oup = ''
    try:
        opts, args = getopt.getopt(argv, 'hi:o:s:e:', ['help', 'input', 'output'])
        for opt, arg in opts:
            if opt in ('-h', '--help'):
                usage()
                sys.exit()
            elif opt in ('-i', '--input'):
                inp = arg
            elif opt in ('-o', '--output'):
                oup = arg
        if oup == '':
            oup = inp.split("/")[-1].rstrip('.' + inp.split('.')[-1]) + '_masked.' + inp.split(".")[-1]
    except getopt.GetoptError:
        usage()
        sys.exit(2)

def usage():
    print('usage: ' + sys.argv[0] + ' -h --help -o --output [path]')

if __name__ == '__main__':
    main(sys.argv[1:])

seq = [[x.replace('\r\n','\n').split('\n')[0], ''.join(x.replace('\r\n','\n').split('\n')[1:]).replace(' ','')] for x in open(inp, 'r').read().rstrip('\n').split('>')[1:]]

w = open(oup, 'w')

for i in range(len(seq)):
    w.write('>' + seq[i][0] + '\n' + seq[i][1].replace('a','N').replace('c','N').replace('g','N').replace('t','N') + '\n')

w.close()

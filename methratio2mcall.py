#! /usr/bin/env python
import sys, time, os, array, optparse
usage = "usage: %prog [options] BSMAP_MAPPING_FILES"
parser = optparse.OptionParser(usage=usage)
parser.add_option("-o", "--out", dest="outfile", metavar="FILE", help="output methylation ratio file name. [default: STDOUT]", default="")
options, infiles = parser.parse_args()

#################### will be delted after testing
# outfile = './data/moabs/swap/liver.2.shuf.moabs.txt' 
# infiles = ['./data/moabs/swap/liver.2.shuf.txt']
# fout = open(outfile, 'w')
#####################
fin = open(infiles[0], 'r')
fout = open(options.outfile, 'w')
fout.write('#chrom\tstart\tend\tratio\ttotalC\tmethC\tstrand\tnext\tPlus\ttotalC\tmethC\tMinus\ttotalC\tmethC\tlocalSeq\n')
i=0
chrom, start, end, ratio, totalC_both, methC_both, strand, next_N, Plus, totalC_plus, methC_plus, Minus, totalC_minus, methC_minus, localSeq = 'NA',0,0,0,0,0,'NA','G','NA',0,0,'NA',0,0,'CG' 
last_line = ''
for line in fin:
    i=i+1
    if 'strand' in line:
        continue
    if last_line != '':  
        if last_line.split('\t')[2] == '+' and int(line.split('\t')[1]) == (int(last_line.split('\t')[1]) + 1) and last_line != '':
            chrom, start, end, ratio, totalC_both, methC_both, strand, Plus, totalC_plus, methC_plus, Minus, totalC_minus, methC_minus = last_line.split('\t')[0],int(last_line.split('\t')[1]),int(last_line.split('\t')[1])+2,(float(last_line.split('\t')[6]) + float(line.split('\t')[6]))/(float(last_line.split('\t')[5]) + float(line.split('\t')[5])),(float(last_line.split('\t')[5]) + float(line.split('\t')[5])),(float(last_line.split('\t')[6]) + float(line.split('\t')[6])),'B','+',float(last_line.split('\t')[5]),float(last_line.split('\t')[6]),'-',float(line.split('\t')[5]),float(line.split('\t')[6] )  
            fout.write('%s\t%d\t%d\t%.3f\t%.3f\t%d\t%s\t%s\t%s\t%.3f\t%d\t%s\t%.3f\t%d\t%s\n' % (chrom, start, end, ratio, totalC_both, methC_both, strand, next_N, Plus, totalC_plus, methC_plus, Minus, totalC_minus, methC_minus, localSeq))
            last_line = ''
        elif last_line.split('\t')[2] == '+' and int(line.split('\t')[1]) != (int(last_line.split('\t')[1]) + 1) and last_line != '':
            chrom, start, end, ratio, totalC_both, methC_both, strand, Plus, totalC_plus, methC_plus, Minus, totalC_minus, methC_minus = last_line.split('\t')[0],int(last_line.split('\t')[1]),int(last_line.split('\t')[1])+2,float(last_line.split('\t')[6])/float(last_line.split('\t')[5]),float(last_line.split('\t')[5]),float(last_line.split('\t')[6]),'+','+',float(last_line.split('\t')[5]),float(last_line.split('\t')[6]),'-',0,0 
            fout.write('%s\t%d\t%d\t%.3f\t%.3f\t%d\t%s\t%s\t%s\t%.3f\t%d\t%s\t%.3f\t%d\t%s\n' % (chrom, start, end, ratio, totalC_both, methC_both, strand, next_N, Plus, totalC_plus, methC_plus, Minus, totalC_minus, methC_minus, localSeq))
            last_line = line
        else:
            chrom, start, end, ratio, totalC_both, methC_both, strand, Plus, totalC_plus, methC_plus, Minus, totalC_minus, methC_minus = last_line.split('\t')[0],int(last_line.split('\t')[1])-2,int(last_line.split('\t')[1]),float(last_line.split('\t')[6])/float(last_line.split('\t')[5]),float(last_line.split('\t')[5]),float(last_line.split('\t')[6]),'-','+',0,0,'-',float(last_line.split('\t')[5]),float(last_line.split('\t')[6])
            fout.write('%s\t%d\t%d\t%.3f\t%.3f\t%d\t%s\t%s\t%s\t%.3f\t%d\t%s\t%.3f\t%d\t%s\n' % (chrom, start, end, ratio, totalC_both, methC_both, strand, next_N, Plus, totalC_plus, methC_plus, Minus, totalC_minus, methC_minus, localSeq))
            last_line = line
    else:
        last_line = line

# if last_line.split('\t')[2] == '+' and int(line.split('\t')[1]) == (int(last_line.split('\t')[1]) + 1) and last_line != '':
#     chrom, start, end, ratio, totalC_both, methC_both, strand, Plus, totalC_plus, methC_plus, Minus, totalC_minus, methC_minus = last_line.split('\t')[0],int(last_line.split('\t')[1]),int(last_line.split('\t')[1])+2,(float(last_line.split('\t')[6]) + float(line.split('\t')[6]))/(float(last_line.split('\t')[5]) + float(line.split('\t')[5])),(float(last_line.split('\t')[5]) + float(line.split('\t')[5])),(float(last_line.split('\t')[6]) + float(line.split('\t')[6])),'B','+',float(last_line.split('\t')[5]),float(last_line.split('\t')[6]),'-',float(line.split('\t')[5]),float(line.split('\t')[6] )  
#     fout.write('%s\t%d\t%d\t%.3f\t%.3f\t%d\t%s\t%s\t%s\t%.3f\t%d\t%s\t%.3f\t%d\t%s\n' % (chrom, start, end, ratio, totalC_both, methC_both, strand, next_N, Plus, totalC_plus, methC_plus, Minus, totalC_minus, methC_minus, localSeq))
#     last_line = ''
# elif last_line.split('\t')[2] == '+' and int(line.split('\t')[1]) != (int(last_line.split('\t')[1]) + 1) and last_line != '':
#     chrom, start, end, ratio, totalC_both, methC_both, strand, Plus, totalC_plus, methC_plus, Minus, totalC_minus, methC_minus = last_line.split('\t')[0],int(last_line.split('\t')[1]),int(last_line.split('\t')[1])+2,float(last_line.split('\t')[6])/float(last_line.split('\t')[5]),float(last_line.split('\t')[5]),float(last_line.split('\t')[6]),'+','+',float(last_line.split('\t')[5]),float(last_line.split('\t')[6]),'-',0,0 
#     fout.write('%s\t%d\t%d\t%.3f\t%.3f\t%d\t%s\t%s\t%s\t%.3f\t%d\t%s\t%.3f\t%d\t%s\n' % (chrom, start, end, ratio, totalC_both, methC_both, strand, next_N, Plus, totalC_plus, methC_plus, Minus, totalC_minus, methC_minus, localSeq))
#     last_line = line
# else:
#     chrom, start, end, ratio, totalC_both, methC_both, strand, Plus, totalC_plus, methC_plus, Minus, totalC_minus, methC_minus = last_line.split('\t')[0],int(last_line.split('\t')[1])-2,int(last_line.split('\t')[1]),float(last_line.split('\t')[6])/float(last_line.split('\t')[5]),float(last_line.split('\t')[5]),float(last_line.split('\t')[6]),'-','+',0,0,'-',float(last_line.split('\t')[5]),float(last_line.split('\t')[6])
#     fout.write('%s\t%d\t%d\t%.3f\t%.3f\t%d\t%s\t%s\t%s\t%.3f\t%d\t%s\t%.3f\t%d\t%s\n' % (chrom, start, end, ratio, totalC_both, methC_both, strand, next_N, Plus, totalC_plus, methC_plus, Minus, totalC_minus, methC_minus, localSeq))
#     last_line = line

fin.close()
fout.close()

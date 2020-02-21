#!/usr/bin/python3
'''=================================================================================================
This script is a "hack" on the original infer_experiment.py script provided by RSeQC.  
The output format is not very easy to parse, so we create an alternate version that directly reports
the findings.

Note that it uses a threshold value for determining strandedness.

Everything up to the end is the same as the original script.
================================================================================================='''

#import built-in modules
import os,sys
if sys.version_info[0] != 3:
	print("\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " This verion of RSeQC needs python3!\n", file=sys.stderr)
	sys.exit()	

import re
import string
from optparse import OptionParser
import warnings
import string
import collections
import math
from time import strftime

#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *
from scipy.stats import binom_test

#import my own modules
from qcmodule import SAM
#changes to the paths

#changing history to this module


__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="3.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


def printlog (mesg):
    '''print progress into stderr and log file'''
    mesg="@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
    LOG=open('class.log','a')
    print(mesg, file=sys.stderr)
    print(mesg, file=LOG)


def main():
    usage="%prog [options]" + "\n"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Input alignment file in SAM or BAM format")
    parser.add_option("-r","--refgene",action="store",type="string",dest="refgene_bed",help="Reference gene model in bed fomat.")
    parser.add_option("-s","--sample-size",action="store",type="int",dest="sample_size",default=200000, help="Number of reads sampled from SAM/BAM file. default=%default")	
    parser.add_option("-q","--mapq",action="store",type="int",dest="map_qual",default=30,help="Minimum mapping quality (phred scaled) for an alignment to be considered as \"uniquely mapped\". default=%default")
    # extra options added by me:
    parser.add_option("-p", "--pval", action="store", type="float", dest="pval_threshold", default=1e-5, help="Binomial p-value for rejecting null hypothesis that experiment was unstranded. default=%default")
    parser.add_option("-o", "--outfile", action="store", type="string", dest="output_file", help="Name of the output file to write the result to.")

    (options,args)=parser.parse_args()

    if not (options.input_file and options.refgene_bed):
        parser.print_help()
        print('\n\n' + __doc__, file=sys.stderr)
        sys.exit(0)
    for f in (options.input_file,options.refgene_bed):
        if not os.path.exists(f):
            print('\n\n' + f + " does NOT exists." + '\n', file=sys.stderr)
            sys.exit(0)
    if options.sample_size <1000:
        print("Warn: Sample Size too small to give a accurate estimation", file=sys.stderr)
    obj = SAM.ParseBAM(options.input_file)
    (protocol,sp1,sp2,other)=obj.configure_experiment(refbed=options.refgene_bed, sample_size = options.sample_size, q_cut = options.map_qual)
    if other <0: 
        other=0.0
    '''
    #Below is original "status" message that gets printed to the console
	if protocol == "PairEnd":
		print("\n\nThis is PairEnd Data")
		print("Fraction of reads failed to determine: %.4f" % other)
		print("Fraction of reads explained by \"1++,1--,2+-,2-+\": %.4f" % sp1)
		print("Fraction of reads explained by \"1+-,1-+,2++,2--\": %.4f" % sp2)
		
	elif protocol == "SingleEnd":
		print("\n\nThis is SingleEnd Data")
		print("Fraction of reads failed to determine: %.4f" % other)
		print("Fraction of reads explained by \"++,--\": %.4f" % sp1)
		print("Fraction of reads explained by \"+-,-+\": %.4f" % sp2)
		
	else:
		print("Unknown Data type")
	#print mesg
    '''

    '''
    sp1 and sp2 are floats giving the fraction of reads.
    The "safest" option is to assume that the experiment is NOT stranded, which then makes tools like 
    featureCounts skip the ambiguous areas.  However, if we have strong enough evidence to support a stranded
    protocol, that would likely be more accurate for quantification.

    To that end, we perform a binomial test.  We know the number of trials (m) and the three fractions: other, sp1, sp2
    From those, we can get rough estimates on the number of reads assigned to each group ( floor(m*other), floor(m*sp1), floor(m*sp2))
    Then, we will ignore the unassigned reads and use N=floor(m*sp1) + floor(m*sp2) as an estimate of the total trials.  From there, 
    a binomial test with N trials and m*sp1 successes will give us a p-value for the null hypothesis that it is unstranded, which assumes
    reads have an equal probability from each strand. (p=0.5).  It is usually very obvious if the protocol is stranded, e.g. sp1=0.05, sp2=0.95
    so the p-value will be VERY stringent
    '''
    fout = open(options.output_file, 'w')
    m = options.sample_size
    n1 = int(m*sp1)
    n2 = int(m*sp2)
    N = n1 + n2 # total reads that were assigned to one "style" or the other
    pval = binom_test(n1, N, p=0.5)

    header = 'sp1_fraction,sp2_fraction,total_sampled,total_assigned,n1,n2,pval_threshold,pval,strand_option\n'
    fout.write(header)
    # if we reject the null hypothesis by the chosen threshold:
    if pval < options.pval_threshold:
        if sp2 > sp1:
            # this corresponds to "reverse stranded" in featureCounts parlance.  Created by dUTP, for instance.
            # The option for featurecounts is -s2
            strand_option = 2
        else:
            # this corresponds to the other stranded protocols
            strand_option = 1
    else:
        # not rejecting null hypothesis of unstranded
            strand_option = 0
    data='%.4f,%4f,%d,%d,%d,%d,%.4E,%.4E,%d' % (sp1,sp2,m,N,n1,n2,options.pval_threshold, pval, strand_option)
    fout.write(data)
    fout.close()


if __name__ == '__main__':
	main()

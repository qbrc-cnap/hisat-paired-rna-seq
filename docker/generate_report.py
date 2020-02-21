#!/usr/bin/python3

import subprocess as sp
import json
import argparse
import pandas as pd
import sys
import os
from jinja2 import Environment, FileSystemLoader

# some variables for common reference.
# many refer to keys set in the json file of
# config variables
TEMPLATE = 'template'
OUTPUT = 'output'
ANNOTATIONS = 'annotations_file'
DESEQ_OUTPUT = 'deseq2_output_filename'
CFG = 'config_vars'
R1 = 'r1_files'
R2 = 'r2_files'
GIT_REPO = 'git_repo'
GIT_COMMIT = 'git_commit'
GENOME = 'genome'
NC_FILE_SUFFIX = 'normalized_counts_file_suffix'
VERSUS_SEP = 'versus_sep'
ADJ_PVAL = 'adj_pval'


class InputDisplay(object):
    '''
    A simple object to carry info to the markdown report.
    '''
    def __init__(self, sample_name, r1, r2):
        self.sample_name = sample_name
        self.r1 = r1
        self.r2 = r2


class AnnotationDisplay(object):
    '''
    A simple object to carry info to the markdown report.
    '''
    def __init__(self, sample_name, condition):
        self.name = sample_name
        self.condition = condition


class ContrastDisplay(object):
    '''
    A simple object for handling relevant info about an experimental 
    contrast in differential expression
    '''
    def __init__(self, 
        base_condition, 
        exp_condition, 
        up_counts, 
        down_counts,
        contrast_name):
        self.base_condition = base_condition
        self.exp_condition = exp_condition
        self.up_counts = up_counts
        self.down_counts = down_counts
        self.contrast_name = contrast_name


def get_jinja_template(template_path):
    '''
    Returns a jinja template to be filled-in
    '''
    template_dir = os.path.dirname(template_path)
    env = Environment(loader=FileSystemLoader(template_dir), lstrip_blocks=True, trim_blocks=True)
    return env.get_template(
        os.path.basename(template_path)
    )

def create_annotation_objects(annotation_file):
    df = pd.read_csv(
        annotation_file, 
        sep='\t', 
        header=None, 
        names=['sample','condition']
    )
    return [AnnotationDisplay(r['sample'], r['condition']) for i,r in df.iterrows()]


def run_cmd(cmd, return_stderr=False):
    '''
    Runs a command through the shell
    '''
    p = sp.Popen(cmd, shell=True, stderr=sp.PIPE, stdout=sp.PIPE)
    stdout, stderr = p.communicate()
    if return_stderr:
        return stderr.decode('utf-8')
    return stdout.decode('utf-8')


def get_versions():
    '''
    Runs the command to get the versions
    '''
    hisat_version_str = run_cmd('hisat2 --version')
    samtools_version_str = run_cmd('samtools --version')
    fc_version_str = run_cmd('featureCounts -v', return_stderr=True)
    multiqc_version_str = run_cmd('multiqc --version')
    fastqc_version_str = run_cmd('fastqc --version')
    rseqc_version_str = run_cmd('pip3 freeze | grep RSeQC')
    picard_md_str = run_cmd('java -jar /opt/software/picard/picard.jar MarkDuplicates --version', return_stderr=True)
    
    # post-process these to extract just the version number:
    hisat_version = hisat_version_str.split('\n')[0].split(' ')[-1]
    samtools_version = samtools_version_str.split('\n')[0].split(' ')[1]
    fc_version = fc_version_str.strip().split(' ')[-1]
    multiqc_version = multiqc_version_str.strip().split(',')[-1].split(' ')[-1]
    fastqc_version = fastqc_version_str.strip().split(" ")[-1]
    rseqc_version = rseqc_version_str.strip().split('==')[1]
    picard_md_version = picard_md_str.strip().split(' ')[-1]

    d = {}
    d['hisat_version'] = hisat_version
    d['samtools_version'] = samtools_version
    d['featurecounts_version'] = fc_version
    d['multiqc_version'] = multiqc_version
    d['fastqc_version'] = fastqc_version
    d['rseqc_version'] = rseqc_version
    d['picard_mark_duplicates_version'] = picard_md_version
    return d

def summarize_contrasts(deseq_outputs, versus_sep, adj_pval_threshold):
    '''
    Parses through the contrasts and extracts summary info
    '''
    summary_list = []
    for f in deseq_outputs:
        # f is the filename, e.g. A_versus_B.deseq.tsv (or something like that)
        # the base and experimental condition are separated by `versus_sep`, which
        # we show as '_versus_' above.  The first (A) is considered the experimental
        # condition.
        contrast_name = os.path.basename(f).split('.')[0]
        experimental, base = contrast_name.split(versus_sep)
        df = pd.read_csv(f, sep='\t')
        sig = df.loc[df['padj'] <= adj_pval_threshold]
        up_count = pd.np.sum(sig['log2FoldChange'] > 0)
        down_count = pd.np.sum(sig['log2FoldChange'] < 0)
        c = ContrastDisplay(base, experimental, up_count, down_count, contrast_name)
        summary_list.append(c)
    return summary_list
    
def get_r_environment():
    '''
    Gets the sessionInfo for all bioc packages
    '''
    cmd = 'R -e "library(\'DESeq2\'); sessionInfo();"'
    return run_cmd(cmd)


def parse_input():
    '''
    Parses the commandline input, returns a dict
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', required=True, dest=ANNOTATIONS)
    parser.add_argument('-d', required=True, dest=DESEQ_OUTPUT, nargs='+')
    parser.add_argument('-t', required=True, dest=TEMPLATE)
    parser.add_argument('-o', required=True, dest=OUTPUT)
    parser.add_argument('-j', required=True, dest=CFG)
    parser.add_argument('-r1', required=True, dest=R1, nargs='+')
    parser.add_argument('-r2', required=True, dest=R2, nargs='+')

    args = parser.parse_args()
    return vars(args)


def fill_template(context, template_path, output):
    if os.path.isfile(template_path):
        template = get_jinja_template(template_path)
        with open(output, 'w') as fout:
            fout.write(template.render(context))
    else:
        print('The report template was not valid: %s' % template_path)
        sys.exit(1)


if __name__ == '__main__':

    # parse commandline args and separate the in/output from the
    # context variables
    arg_dict = parse_input()
    output_file = arg_dict.pop(OUTPUT)
    input_template_path = arg_dict.pop(TEMPLATE)

    # get all the software versions:
    versions_dict = get_versions()

    # get the info about R/Bioc we used:
    session_info = get_r_environment()

    # get a list of objects for showing the sample annotations:
    annotations_object_list = create_annotation_objects(arg_dict[ANNOTATIONS])

    # parse the json file which has additional variables
    j = json.load(open(arg_dict[CFG]))

    # alter how the files are displayed:
    r1_files = arg_dict[R1]
    r2_files = arg_dict[R2]
    samples = [os.path.basename(x)[:-len('_R1.fastq.gz')] for x in r1_files]
    file_display = []
    for r1, r2, s in zip(r1_files, r2_files, samples):
        ipd = InputDisplay(s, r1, r2)
        file_display.append(ipd)

    # get info about the contrasts performed:
    contrast_display_list = summarize_contrasts(
        arg_dict[DESEQ_OUTPUT], 
        j[VERSUS_SEP], 
        float(j[ADJ_PVAL])
    )

    # get the suffix for the deseq files:
    deseq_suffix = '.'.join(arg_dict[DESEQ_OUTPUT][0].split('.')[1:])

    # make the context dictionary
    context = {}
    context.update(versions_dict)
    context.update(arg_dict)
    context.update(j)
    
    context[ANNOTATIONS] = os.path.basename(context[ANNOTATIONS])
    context.update({'session_info': session_info})
    context.update({'annotation_objs': annotations_object_list})
    context.update({'file_display': file_display})
    context.update({'contrast_display': contrast_display_list})
    context.update({'deseq2_output_file_suffix': deseq_suffix})

    # fill and write the completed report:
    fill_template(context, input_template_path, output_file)

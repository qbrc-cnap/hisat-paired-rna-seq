import pandas as pd
import argparse

SAMPLE = 'sample'
CONDITION = 'condition'
COL_NAMES = [SAMPLE, CONDITION]
ANNOTATIONS = 'annotations'
OUTPUT_ANN = 'ann_output'

def get_commandline_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', required=True, dest=ANNOTATIONS)
    parser.add_argument('-o', required=True, dest=OUTPUT_ANN)
    args = parser.parse_args()
    return vars(args)


def read_annotations(annotation_filepath):
    '''
    Tries to parse the annotation file.  If it cannot, issue return some sensible
    message
    '''
    generic_problem_message = '''A problem occurred when trying to parse your annotation 
        file, which was inferred to be in %s format.  
        Please ensure it follows our expected formatting.'''
    file_extension = annotation_filepath.split('.')[-1].lower()
    if file_extension == 'tsv':
        df = pd.read_csv(annotation_filepath, header=None, sep='\t')
    elif file_extension == 'csv':
        df = pd.read_csv(annotation_filepath, header=None, sep=',')
    elif ((file_extension == 'xlsx') or (file_extension == 'xls')):    
        df = pd.read_excel(annotation_filepath, header=None)
    else:
        return None

    # now that we have successfully parsed something.  
    if df.shape[1] < 2:
        return (None, 
                ['The file extension of the annotation file was %s, but the' 
                ' file reader parsed %d column(s).  Please check your annotation file.  Could it have the wrong file extension?' % (file_extension, df.shape[1])])

    # in case the client put extra columns that are blank, just keep the first two
    if df.shape[1] > 2:
        df = df.ix[:,[0,1]]

    # drop any completely empty rows:
    df = df.dropna(how='all')

    # check for NAs in partially filled rows:
    if df.dropna().shape != df.shape:
        return (None, ['There were missing inputs in the annotation table.  Look for blank cells in particular.'])

    # we have two cols and had no NAs.  Now name the cols:
    df.columns = COL_NAMES
    return df


if __name__ == '__main__':
    arg_dict = get_commandline_args()
    df = read_annotations(arg_dict[ANNOTATIONS])
    df.to_csv(arg_dict[OUTPUT_ANN], sep='\t', index=False, header=False)
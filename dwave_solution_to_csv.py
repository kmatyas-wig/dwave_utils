#!/usr/bin/env python3

import sys
import pickle
from io import StringIO
from argparse import ArgumentParser
import pandas as pd
import dimod

version = "0.1"
myname = sys.argv[0].replace('./', '')
parser = ArgumentParser(
    description='Extract sample in csv from picked DWave solution.',
    prog=myname)

parser.add_argument('dwavepickle', nargs=1,
                    type=str,
                    help = \
                    'Pickled dwave result, e.g. the output of solve_dwave_aws_QUBO_from_csv.py')

parser.add_argument('-o', '--outfile', type=str,
                    help='File to save the csv into. STDOUT if not specified')

parser.add_argument('-n', '--no-chain-break', action='store_true',
                    help='Do not include chain_break_fraction')


args = parser.parse_args()

if args.outfile is not None:
    output = args.outfile
else:
    output = StringIO()

with open(args.dwavepickle[0], 'rb') as infile:
    sample_dict = pickle.load(infile)
    infile.close()


df = dimod.SampleSet.from_serializable(sample_dict).to_pandas_dataframe()

try:
    df=df[['energy','num_occurrences','chain_break_fraction']+list(
        df.keys()[0:-3])].sort_values('energy')
except:
    df=df[['energy','num_occurrences']+list(
        df.keys()[0:-2])].sort_values('energy')

if args.no_chain_break:
    df.drop(columns=['chain_break_fraction'], inplace=True)

df.to_csv(output, index=False)

if args.outfile is None:
    print(output.getvalue())

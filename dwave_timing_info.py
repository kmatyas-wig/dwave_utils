#!/usr/bin/env python3

import pickle
import sys
import csv
from argparse import ArgumentParser

version = "0.1"
myname = sys.argv[0].replace('./', '')
parser = ArgumentParser(
    description='''Solve QUBO from a csv with header on DWave.
See also https://docs.dwavesys.com/docs/latest/c_qpu_timing.html
Times are in microseconds.
''',
    prog=myname)

parser.add_argument('dwavepickle', nargs=1,
                    type=str,
                    help = \
                    'Dwave result pickle')

parser.add_argument('--noheader', action="store_true",
                    help='Do not write header')

parser.add_argument('--fields', type=str,
                    help='''Comma-separated list of fields to include, subset of 
qpu_delay_time_per_sample,post_processing_overhead_time,
qpu_anneal_time_per_sample,minimum,
qpu_access_overhead_time,qpu_sampling_time,total_post_processing_time,
qpu_programming_time,qpu_access_time,qpu_readout_time_per_sample''')


ARGS = parser.parse_args()

with open(ARGS.dwavepickle[0], 'rb') as ifi:
    data = pickle.load(ifi)
    ifi.close()

output = data['info']['timing']
output['minimum'] = sorted(data['vectors']['energy']['data'])[0]
if ARGS.fields is None:
    fields = fieldnames=output.keys()
    finalOutput = output
else:
    fields = ARGS.fields.split(',')
    finalOutput = { thekey: output[thekey] for thekey in fields}
writer = csv.DictWriter(sys.stdout, fieldnames=fields)

if not ARGS.noheader:
    writer.writeheader()
writer.writerow(finalOutput)


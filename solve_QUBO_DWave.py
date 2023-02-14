#!/usr/bin/env python3

import sys
import os
import pickle
from argparse import ArgumentParser
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

version = "0.1"
myname = sys.argv[0].replace('./', '')
parser = ArgumentParser(
    description='Solve QUBO from a csv with header on DWave.',
    prog=myname)

parser.add_argument('csvfile', nargs=1,
                    type=str,
                    help = \
                    'Csv file with the sparse QUBO matrix Q. Symmetric, upper, or lower triangular.')

parser.add_argument('--outfile', type=str,
                    help='File to save DWave result in',
                    default="last_dwave_result.pickle")

parser.add_argument('--dry', action="store_true",
        help='Do not solve the problem, just test reading.')

parser.add_argument('--noautoscale', action="store_true",
                    help='Disable autoscale.')

parser.add_argument('--num-reads', type=int,
                    help='Number of reads, default=1000',
                    default=1000)

parser.add_argument('--annealing-time', type=float,
                    help='Annealing time for each sample in microseconds, deafult=20.0',
                    default=20.0)

args = parser.parse_args()

Q={}
with open(args.csvfile[0], 'rt') as infile:
    for rowraw in infile.readlines()[1:]:
        row = rowraw.strip()
        (i_str, j_str, Qij_str) = row.split(',')
        i = int(i_str)
        j = int(j_str)
        #Sparse format.
        #Here we make it upper triangular
        #Can also be symmetric, this is how
        #sample_qubo would deal with it anyway
        if i == j:
            Q[(i,j)] = float(Qij_str)
        elif i<j:
            try:
                Q[(i,j)] += float(Qij_str)
            except:
                Q[(i,j)] = float(Qij_str)
        else:
            try:
                Q[(j,i)] += float(Qij_str)
            except:
                Q[(j,i)] = float(Qij_str)

autoscale = not args.noautoscale
if not args.dry:
    dwave_key=os.getenv('DWAVETOKEN')
    if dwave_key is None:
        raise ValueError('Error: invalid DWave token, set DWAVETOKEN properly')
        exit(-1)
    sampler = EmbeddingComposite(DWaveSampler(endpoint='https://cloud.dwavesys.com/sapi',
                                            token=dwave_key))
    sampleset = sampler.sample_qubo(Q,
                                    num_reads=args.num_reads,
                                    annealing_time=args.annealing_time,
                                    auto_scale=autoscale)
    sdf = sampleset.to_serializable()
    with open(args.outfile, "wb") as ofi:
        pickle.dump(sdf, ofi)
        ofi.close()
else:
    print("Dry run, echoing Q")
    print(Q)
    print("num_reads=", args.num_reads)
    print("annealing_time=",args.annealing_time)

print("Problem has concluded.")

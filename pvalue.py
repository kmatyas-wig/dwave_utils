#!/usr/bin/env python3
#
# python implementation of
#   https://github.com/iitis/StatisticalVerificationOfIsingEnergies
# please cite
#  Statistical quality assessment of Ising-based annealer outputs
#  by K Domino, M Koniorczyk, Z Puchała
# Quantum Information Processing 21 (8), 1-19
# https://link.springer.com/article/10.1007/s11128-022-03623-5

description = """
python implementation of
   https://github.com/iitis/StatisticalVerificationOfIsingEnergies
 Please cite
  Statistical quality assessment of Ising-based annealer outputs
  by K Domino, M Koniorczyk, Z Puchała
 Quantum Information Processing 21 (8), 1-19
 https://link.springer.com/article/10.1007/s11128-022-03623-5
"""
import sys
import pickle
from argparse import ArgumentParser
import pandas as pd
import numpy as np
import scipy.stats as st
import random
from statsmodels.distributions.empirical_distribution import ECDF
import dimod

def estimate_ground_state_energy(energies, alpha=0.19):
    """
    estimate gropund state enetguy from energies of Ising output and alpha perameter
    """
    nu = np.mean(energies)
    sigma = st.tstd(energies)
    eta = st.skew(energies, bias = True)

    return nu - (alpha + 2)/(alpha + 1)*sigma/eta

def bootstrap_hists_of_mins(energies, alpha, s):
    """ returns s estimated ground state energies from bootstraping """
    l = len(energies)
    ret = np.zeros(s)
    for i in range(s):
        e_random = random.choices(energies, k=l)
        ret[i] = estimate_ground_state_energy(e_random, alpha)
    return ret


def bootstrap_get_pvalue(energies, alpha=0.19, s=1000):
    """
    Estimates p-value that the ground state is in Ising based energies.
    """
    y = bootstrap_hists_of_mins(energies, alpha, s)
    ecdf = ECDF(y)
    Hmin = np.min(energies)
    return 1.0 - ecdf(Hmin)

if __name__ == "__main__":

    version = "0.1"
    myname = sys.argv[0].replace('./', '')
    parser = ArgumentParser(
        description=description,
        prog=myname)

    parser.add_argument('infilename', nargs=1,
                        type=str,
                        help=\
                        'Csv sample with "energy" and "num_occurrences fields", or pickled output of sampleset.to_serializable()')

    parser.add_argument('--alpha',
                        type=float,
                        help='alpha value, default: 0.19, see paper',
                        default=0.19,
                        required=False)

    parser.add_argument('-s', '--bootstrap-samples',
                        type=int,
                        help='number of bootstrap samples, default 1000',
                        default=1000,
                        required=False)

    parser.add_argument('--estimate-energy',
                        help='Provide an energy estimate from the raw sample',
                        action='store_true')



    args = parser.parse_args()

    have_data = False
    with open(args.infilename[0], 'rb') as infile:
        try:
            sample_dict = pickle.load(infile)
            have_data = True
            df = dimod.SampleSet.from_serializable(sample_dict).to_pandas_dataframe()
        except:
            have_data = False
        infile.close()
    if not have_data:
        df = pd.read_csv(args.infilename[0])

    energies = []
    for index, row in df.iterrows():
        for _ in range(int(row['num_occurrences'])):
            energies.append(row['energy'])

    p = bootstrap_get_pvalue(energies, alpha=args.alpha, s=args.bootstrap_samples)
    print(p)

    if args.estimate_energy:
        gs_energy = estimate_ground_state_energy(energies, args.alpha)
        print("H_min=%f"%gs_energy)

#!/usr/bin/env python3

import pickle
import sys
import morejson as json

with open(sys.argv[1], 'rb') as ifi:
    data = pickle.load(ifi)
    ifi.close()

print(json.dumps(data))

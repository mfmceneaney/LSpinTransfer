import pandas as pd
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.dates as md
from datetime import datetime

if len(sys.argv)==1:
    print("Usage: "+os.path.abspath(sys.argv[0])+" /path/to/efficiencies.txt")
    sys.exit(0)

# Set input filename
filename  = sys.argv[1]
delimiter = " " #TODO: Make this an argument...
verbose   = True

# Read file
df = pd.read_csv(filename,delimiter=delimiter)
ks = df.keys()


for key in ks:
    if key != 'cpu_used': continue

    # Set data
    x = df[key].to_list()
    print("DEBUGGING: x = ",x)
    print("DEBBUGGING: type(x) = ",type(x))
    print("DEBBUGGING: type(x[0]) = ",type(x[0]))
    for idx, el in enumerate(x):
        if '-' in el:
            t = el.split('-')[1]
            d = str(int(el.split('-')[0])+1) #NOTE: Add an extra day
            new_el = '-'.join([d,t])
            print("DEBUGGING: new_el = ",new_el)
            x[idx] = datetime.strptime(new_el, "%d-%H:%M:%S")
            print(idx,x[idx])
        else:
            x[idx] = datetime.strptime(el, "%H:%M:%S")
            if idx<10: print(idx,x[idx])

    # Plot data
    bins = 100
    figsize = (16,10)
    f = plt.figure(figsize=figsize)
    h = plt.hist(x, bins=bins, label=key)
    plt.xlabel(key+" + 1day")
    plt.ylabel('Counts')
    # minor_steps = 1
    # plt.gca().xaxis.set_major_locator(plt.MultipleLocator(5*minor_steps))
    # plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(minor_steps))
    plt.gca().xaxis.set_major_formatter(md.DateFormatter('%d %H:%M:%S'))
    f.savefig(key+'.pdf')

if verbose: plt.show()
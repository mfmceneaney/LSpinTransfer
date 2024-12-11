
#--------------------------------------------------#
# Script to plot CLAS12 data from TGraphErrors
# saved to ROOT files.

import uproot as ur
import numpy as np
import matplotlib.pyplot as plt

def offset_graph_x(g, offset):
    npoints = len(g[0])
    for idx in range(npoints):
        g[0][idx] += offset

def get_plots(
    method = 'HB',
    inpath1 = 'HB_CT1/HB_costheta1_mass_ppim_0.000_100.000_sgasym_0.00_bgasym_0.00.root',
    inpath2 = 'HB_CT2/HB_costheta2_mass_ppim_0.000_100.000_sgasym_0.00_bgasym_0.00.root',
    injected_asym = 0.00,
    ):

    file1 = ur.open(inpath1)
    file2 = ur.open(inpath2)

    g1 = [np.array(file1["Graph"].member('fX')), file1["Graph"].member('fY'), file1["Graph"].member('fEX'), file1["Graph"].member('fEY')]
    g2 = [np.array(file2["Graph"].member('fX')), file2["Graph"].member('fY'), file2["Graph"].member('fEX'), file2["Graph"].member('fEY')]

    # Get costheta 1 results
    asym_ct1 = g1[1][0]
    asym_ct1_err = g1[3][0]

    # Get costheta 2 results
    asym_ct2 = g2[1][0]
    asym_ct2_err = g2[3][0]

    ret_str_1 = f"{injected_asym:.2f} & "+("" if asym_ct1<0 else " ")+f"{asym_ct1:.5f} $\pm$ {asym_ct1_err:.5f} \\\\"
    ret_str_2 = f"{injected_asym:.2f} & "+("" if asym_ct1<0 else " ")+f"{asym_ct1:.5f} $\pm$ {asym_ct1_err:.5f} \\\\"

    return [ret_str_1,ret_str_2]

if __name__=="__main__":

    verbose = True
    injected_asyms = [-0.10,-0.01,0.00,0.01,0.10]
    # asym_str = f'sgasym_{injected_asym:.2f}_bgasym_0.00'
    method = 'HB'

    packages = []
    for injected_asym in injected_asyms:
        asym_str = f'sgasym_{injected_asym:.2f}_bgasym_0.00'
        folder = 'mc_asym_injection_'+(f"{injected_asym:.2f}" if np.abs(injected_asym)<=0.01 else f"{injected_asym:.1f}")+'/'
        packages.append({
            'method' : method,
            'inpath1' : folder+method+'_CT1/'+method+'_costheta1_mass_ppim_0.000_100.000_'+asym_str+'.root',
            'inpath2' : folder+method+'_CT2/'+method+'_costheta2_mass_ppim_0.000_100.000_'+asym_str+'.root',
            'injected_asym' : injected_asym,
        })

    strs1, strs2 = [], []
    for pack in packages:
        str1, str2 = get_plots(**pack)

        strs1.append(str1)
        strs2.append(str2)

    # Show table costheta 1
    print("Cos theta 1 table")
    for el in strs1:
        print(el)
    print()

    # Show table costheta 2
    print("Cos theta 2 table")
    for el in strs2:
        print(el)
    print()

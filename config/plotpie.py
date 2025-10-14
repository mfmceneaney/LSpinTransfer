import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from cycler import cycler
plt.rcParams['text.usetex'] = True #NOTE: Force use of LaTeX in text rendering

labels = ['String', '$\Sigma^{0}$', '$\Sigma^{* 0}$', '$\Sigma^{* +}$', '$\Sigma^{-}$']#, 'Other'] #NOTE: LAST ENTRY IS TOO SMALL AND TOO CLOSE TO PREVIOUS AND MAKES THE PLOT TOO CLUTTERED SO OMIT FOR NOW.
sizes = [66.7, 19.0, 7.2, 5.6, 1.3]#, 0.2]

# Set font sizes
plt.rc('font', size=25) #controls default text size
plt.rc('axes', titlesize=50) #fontsize of the title
plt.rc('axes', labelsize=50) #fontsize of the x and y labels
plt.rc('xtick', labelsize=25) #fontsize of the x tick labels
plt.rc('ytick', labelsize=25) #fontsize of the y tick labels
plt.rc('legend', fontsize=25) #fontsize of the legend

# Get some nicer plot settings
plt.rcParams['font.family'] = 'serif'
plt.rcParams['figure.autolayout'] = True


# Create custom cycler
name = 'Pastel1'
cmap=mpl.colormaps[name]
custom_cycler = (cycler(color=[cmap.colors[i] if i>1 else cmap.colors[1] if i==0 else cmap.colors[0] for i in range(len(cmap.colors))]) +
                 cycler(lw=[i+1 for i in range(len(cmap.colors))]))

# Get some nicer plot settings
plt.rcParams['font.family'] = 'serif'
plt.rcParams['figure.autolayout'] = True

fig, ax = plt.subplots(figsize=(10,10))
ax.set_prop_cycle(custom_cycler)

plt.title('$\Lambda$ parent $x_{F}>0$',usetex=True,pad=20)
ax.pie(sizes, labels=labels, autopct='%1.1f%%',radius=1.0,pctdistance=0.8,labeldistance=1.1)
fig.savefig('pieplot.pdf')
plt.show()

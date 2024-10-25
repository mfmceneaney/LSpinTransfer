# $\Lambda$ Spin Transfer Analysis

Analysis for $\Lambda$ hyperon longitudinal spin transfer in SIDIS at CLAS12.

Use [CLAS12-Analysis](https://github.com/mfmceneaney/CLAS12-Analysis.git) or your own software to produce the input ROOT trees with event by event kinematics selecting all $e^{-}p\pi^{-}$ combinations.

# Build the project

This is a CMake project so you can build wherever but this is probably the simplest way to go:
```bash
mkdir build
cd build
cmake ..
make
```
You should now have several executables in your `build` directory.

# Compute binned results

Put your binning arguments in a yaml file.  An example is included at [args.yaml](args.yaml).
You may add arbitrary 1D binning in any kinematic variable.  Details are all in [analysis.cpp](analysis.cpp) and [include/analysis.h](include/analysis.h).

Also since this is a $\Lambda$ analysis you have to isolate the mass peak.
Hence, there is a mass signal fit done for each kinematic bin.  Details are all in [include/massfit.h](include/massfit.h).


# Make nice plots

Finally, ROOT is a bit finicky about plots so I use matplotlib instead for the final $D^{\Lambda}_{LL'}$ results.
This plotting is run with the [aggregate_results.py](aggregate_results.py) script, although input file paths may need to be changed.

# Plot Kinematics
The ROOT macros [PlotComparisons.C](PlotComparisons.C) and [PlotBeforeAndAfterCuts.C](PlotBeforeAndAfterCuts.C) are used to plot the 1D and 2D kinematics distributions comparing data and MC and the effects of different cuts.

#

Contact: matthew.mceneaney@duke.edu

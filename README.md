# $$\Lambda$$ Spin Transfer Analysis

Analysis for $$\Lambda$$ hyperon longitudinal spin transfer in SIDIS at CLAS12.

Use [CLAS12-Analysis](https://github.com/mfmceneaney/CLAS12-Analysis.git) or your own software to produce the input ROOT trees with event by event kinematics selecting all $$e^{-}p\pi^{-}$$ combinations.

# Build the project

This is a CMake project so you can build wherever but this is probably the simplest way to go:
```bash
mkdir build
cd build
cmake ..
make
```
You should now have an executable called `analysis` in your `build` directory.

# Compute binned results

Put your binning arguments in a yaml file.  An example is included at `args.yaml`.
You may add arbitrary 1D binning in any kinematic variable.  Details are all in `analysis.cpp` and `include/analysis.h`.

Also since this is a $$\Lambda$$ analysis you have to isolate the mass peak.
Hence, there is a mass signal fit done for each kinematic bin.  Details are all in `include/massfit.h`.


# Make nice plots

Finally, ROOT is a bit finicky about plots so I use matplotlib instead for the final $$D^{\Lambda}_{LL'}$$ results.
This plotting is run with the `convert_to_python_plot.py` script, although input file paths may need to be changed.

#

Contact: matthew.mceneaney@duke.edu

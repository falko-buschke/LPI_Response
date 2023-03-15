# R-code to simulate the Living Planet Index for fluctuating populations

This code accompanies the following manuscript:

* Buschke, F.T., Hagan, J.G. Santini, L. & Coetzee, B.W.T. *The relative prevalence of additive and multiplicative stochasticity in nature is distinct from their effects on the Living Planet Index*.

The code was accurate as of 15 March 2023. For enquiries, contact `falko.buschk@gmail.com`

## General comments

The code presented here  relies on the dedicated `rlpi (v0.1.0)` package for calculating the *Living Planet Index*. This code is not on the official CRAN repository, so it must be accessed and installed directly from the [rlpi GitHub repository](https://github.com/Zoological-Society-of-London/rlpi), which also requires the `devtools (v2.4.5.)` package. The code needed to install these packages is included in the R-scripts. The code also uses the `stats (v4.1.3)` package.

## Guidance for use

The code is self-contained and produces Figure 1 as presented in the aforementioned manuscript. To function, take note of the following:

1. The working directory must contain an empty folder titled `LPI_files`. This folder is used by the 'rlpi' package to store working files while calculating the Living Planet Index. 

2. The code simulates 5,000 populations for 16 combinations of additve and multiplicative population strochasticity, and calculates the Living Planet Index for each of these. This whole process can take more than an hour (possibly longer, depending on your computer specs) becasue it has to fit 80,000 generalised additive models. If you merely want to test the code, reduce the number of species' populations being simulated (parameter `S`, line 20) to 100.

## File structure

The repository contains the following:

* An empty folder titled `LPI_files` to store the preliminary files from the 'LPIMain` function.
* The script titled `LPI_toymodel`, containing all the self-contained R-code.
* The figure titled `FigureR1` which is the final output of the script.


# OSFD_fork
A third party fork of the open source feature detector published in https://doi.org/10.1002/rcm.9206

I have added a command line interface using the optparse package for ease of use and a conda environment definition file.

To set up the conda environment run:
`conda env create -n osfd_env -f env.yml`

To print a help message with the command line parameters run:
`conda acivate osfd_env`
`Rscript peakpicking.R --help`

This simulation was used in the following article:
  Neymotin SA, Dura-Bernal S, Lakatos P, Sanger TD, Lytton WW.
  Multitarget Multiscale Simulation for Pharmacological Treatment of Dystonia in Motor Cortex
  Frontiers in Pharmacology 7:157 2016
  http://journal.frontiersin.org/article/10.3389/fphar.2016.00157/full

The code in this folder generates a raster plot from a dystonia simulation.

The simulation was tested/developed on LINUX systems, but may run on Microsoft Windows or Mac OS.

To run the demo, you will need the NEURON simulator (version 7.4 and up; available at http://www.neuron.yale.edu)
compiled with python and MPI enabled. You will need Matplotlib to draw the output ( http://matplotlib.org/ ).

Instructions to setup the model:
  unzip the file
  cd dystdemo
  nrnivmodl 

The nrnivmodl command will produce an architecture-dependent folder with a script called special.
On 64 bit systems the folder is x86_64.

-----------------------------------------------------------------------------------------------------
Run a network simulation from the terminal with the following command:
mpiexec -np 8 nrniv -python -mpi mpisim.py netcfg.cfg

This will run an individual simulation and save its output data.

Next, to draw the raster plots from the simulation run:
python
execfile("simdat.py")

That will load the data and draw an example raster plot from
the set of simulations in the paper (this one has strong beta
in layer 5).

This simulation uses MPI for parallelization. The simulation can take a long time to run,
depending on your hardware setup. The simulation saves output data to the data subdirectory. 
simdat.py loads the data and draws the output.

Note: The -np 8 above specifies that mpiexec should use 8 cores. You can change this number depending
on your hardware setup. If you change the number of cores (-np 8), make sure to update the defnCPU
parameter on line 47 of simdat.py; this is because mpisim.py saves 1 output file per core.
The simulation takes ~8-9 minutes on Linux with 8 cores (Intel(R) Core(TM) i7-4940MX CPU @ 3.10GHz)
and runs in ~4 minutes on Linux with 24  cores (Intel(R) Xeon(R) CPU E5-4610 0 @ 2.40GHz).
-----------------------------------------------------------------------------------------------------

For questions/comments email:  
  samn at neurosim dot downstate dot edu

20160915 This updated version from the Lytton lab allows their models
which contain misc.mod and misc.h to compile on the mac.

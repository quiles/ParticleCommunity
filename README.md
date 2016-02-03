# ParticleCommunity

Instructions to compile:

   1) Download the SNAP/C++ package: https://snap.stanford.edu/releases/Snap-2.4.zip
          [https://snap.stanford.edu] 

   2) Download the source code: https://github.com/quiles/ParticleCommunity/archive/master.zip

   3) Copy the folder into the subfolder ./examples of the Snap package.
          [i.e. ./Snap-2.4/examples/Particle]
   4) Compile the source code: $ make

   Obs) This source code was tested using the following Operating Systems/GCC
      a) OSX El Capitan 10.11.1 / Snap 2.4 / GCC 4.2.1 and
      b) CentOS Linux 6.5 / Snap 2.4 / GCC 4.4.7 

Instructions to run:

1) Using an specific set of parameters

    a) running with default parameters
         $ ./Particle filename.dat 

    b) choosing the parameters
         $ ./Particle filename.dat -b 0.1 -tr 0.1 -v
         This command will run the Particle Model with beta=0.1, \theta_R=0.1 (stop condition), and in verbose mode       

    Help: type $ ./Particle


2) Run for several values of beta

    a) Evaluating the model with disting values of beta (possible hierarchical community detection)
         $ ./Particle filename.dat -sb 0.1 0.8 0.01
         This command will evaluate the network by varying beta from 0.1 to 0.8 with step 0.01. 


3) Run on time-varying networks

    a) To run with time-varying data, you must compile a file with the names of each network state (one file per line)
         $ ./Particle -dynamic filename.dat [plus regular options....]



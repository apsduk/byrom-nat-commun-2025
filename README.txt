This folder contains the the following:
    2 Scripts (.m)
    1 Data file (.mat)
    13 Functions (.m)

Here, we outline briefly the purpose of each file. Further details of their
workings are provided within the files themselves.

The 2 scripts can be run, and do the following:
    CLpATXopt.m:
        This file runs a multi-objective optimisation for the CLpATX circuit,
        aiming to simultaneously maximise initial output P0, half-life T50, and
        time to taken to fall outside a window of +-10% of P0, T+-10, by varying
        three circuit parameters: wA (maximal transcription rate), bA (ribosome
        binding rate) and kA (transcription factor threshold parameter, or 'control
        strength').
    CLpATXvsOL.m:
        Given an initial output P0, this file plots time-series population-wide outputs P
        for an open-loop system and a close-to-optimal closed-loop system which produce
        that initial output.
We include these 2 scripts of example pieces of code which demonstrate the workings
of our model and how we have used it throughout the paper.

The 1 data file contains the following:
    CLpATXOptimisations.mat
        This file shows results of the multi-optimisation defined in CLpATXopt.m, which
        are used in the CLpATXvsOL.m script.

The 13 functions do the following:
    findbijectiveregion.m: Splits a multi-valued function into upper and lower portions
    generatemutrates.m: Generates a transition matrix between competing subpopulations
    gfpCLpATXJbatch.m: Defines a Jacobian for the CLpATX circuit.
    gfpCLpATXODE.m: Defines an ODE model for the CLpATX circuit.
    gfpcost.m: Defines a cost function for the multi-objective optimisation.
    gfpOLJbatch.m: Defines a Jacobian for the open-loop circuit.
    gfpOLODE.m: Defines an ODE model for the open-loop circuit.
    hostonlyODE: Defines an ODE model for a host with no gene circuit.
    loadhostparameters: Defines default parameters and initial conditions for the host model.
    mknewdir.m: Makes a new directory with a given name.
    repeatbatch.m: Performs repeated batch simulations where nutrients and population size are reset each day.
    xlifecalculator.m: Calculates the time taken for a vector to fall below a specified percentage of its original value.
    xlifecalculatoreitherside.m: Calculates the time taken for a vector to fall below or above a specified percentage of its original value.

        
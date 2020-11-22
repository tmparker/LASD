# LASD

Replication files for "Loss aversion and the welfare ranking of policy
interventions"

## Contents
The main directory contains all-R files that illustrate how the loss
averse-sensitive dominance test statistics are computed (although simulations
depended on a C++ version for speed), and a few files that generate figures used
in the appendix.

The subdirectories contain code that was used to conduct the empirical
illustration and to run simulations:

* ``empirical`` contains the data extract that was used along with the files
 that are needed to reproduce the tests illustrated in the section of the paper
 about Jobs First and AFDC and in the appendix.
* ``simulation`` contains the files that produce the results of three simulation
  experiments detailed in the appendix.  See the other readme file
  ``/simulation/READ_about_computation.md`` for notes about how to replicate
  results.

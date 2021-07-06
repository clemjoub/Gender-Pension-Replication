
This repository contains the codes and data sets to replicate the results of the paper:  Gender Pension Gaps
in a Private Retirement Accounts System. (Joubert & Todd)

## Folder architecture:

- Current: contains the solution code and input files
- Dofiles: contain the dofiles used to create the estimation moments and initial conditions
- DATA: Raw data and working datasets should be saved here
- routines: minor routines 

## Computing the estimation moments

The raw data (EPS surveys 2002, 2004, 2006, 2009 and linked pension account balances) should be requested from the Subsecretaria de prevision social de Chile at:
https://www.previsionsocial.gob.cl/sps/biblioteca/encuesta-de-proteccion-social/bases-de-datos-eps/

.\Dofiles\Master.do can be run to create ./Current/inputs/initialconditions.txt and  ./dataformoments.txt. 

Compile and run createmdata.f90 and msm.f90 to create the data moments from dataformoments.txt

IN PROGRESS 4/13/21: subdofiles that assemble the raw data into intermediate datasets are provided in  .\Dofiles\Dataconstruction\upstream but they have not been ported to run in this folder or on the latest release of the raw data yet


## Solving the Model: 

The fortran code solves and simulates a dynamic model of saving and labor supply decisions. 
It returns the simulated method of moments criterion to be minimized to estimate the model's parameters.

The list of source files is contained in the file routine/makeprog. The executable has the following arguments:

 - infile.asc                      file that contains the subset of model parameters that are estimated
  
 - outfile.asc                     file in which the value of the SMM criterion is returned
 
 - 1                               (HOPSPACK argument, not necessary otherwise)
 

The programs also needs to read the following files:

 - ./inputs/input.txt              contains the rest of the model's parameters
  
 - ./inputs/Initialconditions.txt  contains the intial conditions from which simulations start
  
 - ./inputs/mdata.txt              contains the data moments that the model simulations approximate
  
 - ./inputs/momdef.txt             contains the moment definitions  
  
as well as other minor input files:   ./inputs/CNU.txt, ./inputs/dp.txt, ./inputs/Lifetable_men.txt, ./inputs/Lifetable_women.txt
 
The code can be compiled and run using the file: ./RUN, option 1. 
Options can be set within that file to solve the model under the baseline or counterfactual rules and 
also to compute the parameter gradient to get the standard errors.

  
## Estimating the model's parameters:

The code is set up for optimization using HOPSPACK, but it should be easy to adapt to another algorithm.
The optimization routine (e.g. Matlab global optimization toolbox) would have to pass each new parameter
vector to be evaluated to a file named infile.asc which is then read by the fortran code. 
The first line of the file is the total number of parameters being estimated and each of the 
following lines contains one element of the parameter vector to be evaluated. 

The fortran code computes the value of the criterion function to be minimized and writes it to the file
 outfile.asc. The criterion value must then be read and passed to the optimization routine.

# Parameter estimates are stored in:

Current/inputs/input.txt

## Computing standard errors
REPLICABILITY IN PROGRESS 4/13/2021

## Simulating policy experiments
[Ivan please add description of procedure here]
REPLICABILITY IN PROGRESS 4/13/2021

## Producing tables and figures
[Ivan please improve on the description of the procedures here]
### Analysing Reform impacts
After producing simulations under the old and the new systems, run Dofiles/Analysis/ReformImpacts/reform-imp.do
### Analysing Counterfactuals
After producing simulations under different counterfactuals, run Dofiles/Analysis/Counterfactuals/counterfactuals2.do

REPLICABILITY IN PROGRESS 4/13/2021

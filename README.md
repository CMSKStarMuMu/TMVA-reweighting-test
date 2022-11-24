# TMVA-reweighting-test
test program for using  TMVA BDT in the MC reweighting procedure
Tested with rel CMSSW_10_4_0

compile with make

Data and MC files point to the the soft links in eos (you can change to your needs). In this example the files for different "parity" are the same and arity selections are defined in the program.

Be aware that TMVA Multithread is active and the requested #treads is 40!  

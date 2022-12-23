# TMVA-reweighting-test
test program for using  TMVA BDT in the MC reweighting procedure
Tested with rel CMSSW_10_4_0

compile with make

Data and MC files point to the the soft links in eos (you can change to your needs). In this example the files for different "parity" are
the same and parity selections are defined in the program.

Be aware that TMVA Multithread is active and the requested #treads is 40! 

 
example: ./tmva-test-allYears 2016 bdt  ; if you want to train a new bdt & plots

         ./tmva-test-allYears 2016	; if you want just produce the  plots

         ./tmva-test-allYears 2016 save ; if you want save the MC weights in a root file

seed = 123
seqfile=../MySeq.txt
Imapfile=../MyImap.txt
outfile=out.txt
mcmcfile=smcmc.txt
speciesdelimitation=1 1 2 1
speciestree=1
species&tree=3 A B C
               10 10 10 
               ((A,B),C);
usedata=1
nloci=13
locusrate=0
cleandata=0
thetaprior=3 0.02 e
tauprior=3 0.06
finetune=1: .01 .0001 .005 .0005 .2 .01 .01 .01
print=1 0 0 0
burnin=20000
sampfreq=1
nsample=20000
threads=12

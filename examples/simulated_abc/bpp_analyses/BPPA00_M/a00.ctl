seed = 123
seqfile=../MySeq.txt
Imapfile=../MyImap.txt
outfile=out.txt
mcmcfile=smcmc.txt
speciesdelimitation=0
speciestree=0
species&tree=3 A B C
               10 10 10 
               ((A,B),C);
usedata=1
locusrate=0
cleandata=0
thetaprior=3 0.04 e
tauprior=3 0.04
finetune=1: .01 .0001 .005 .0005 .2 .01 .01 .01
print=1 0 0 0
burnin=2000
sampfreq=1
nsample=2000
threads=12

migprior=40 20
migration = 2
	A B
	B A

# BackgroundStudies

This repositery is dedicated to HSCP analysis, for the study of the background estimation method

## Run the first part of the code: produce histograms, with given bins & given cuts (regions definitions)

Configure configFile.txt respecting dedicated columns
```bash
root file   Ias cut>    pt cut>     Ih cut>     momentum cut<       eta num. of bins    Ih num. of bins       p num. of bins      mass num. of bins     invIso  invMET 
```
where you choose the cut values on Ias, pt and Ih. You can also cut on momentum (avoid overflow problem). 
invIso and invMET correspond to a boolean which set inverse cuts on isolation or MET (in order to be sure to not have any signal). 
then run,
```bash
root -l -q -b HscpCandidates.C
```

## Run the second part: read histograms & produce mass plots 

Configure configFile_readHist.txt respecting dedicated columns
```bash
root file   rebin   rebin eta   rebin Ih    rebin p     rebin mass  variable bins mass      variable bins p     threshold mass      threshold p
```
where rebin is a boolean in order to know if we want rebin the different distributions. Then we give the values of rebinning.
For a variable rebinning, we give a threshold and we need to set variable bins mass (or p) to 1 (it is a boolean). 

then run,
```bash
root -l -q -b readHisto.C
```

## Make plots

Use the script makePlots.py in order to produce automatically some plots.

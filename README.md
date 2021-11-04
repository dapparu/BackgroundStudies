# BackgroundStudies

This repositery is dedicated to HSCP analysis, for the study of the background estimation method

## Run the first part of the code: produce histograms, with given bins & given cuts (regions definitions)

Configure configFile.txt respecting dedicated columns
```bash
root file       Ias cut>    pt cut>     Ih cut>     momentum cut<       eta num. of bins    Ih num. of bins       p num. of bins      mass num. of bins     invIso  invMET 
```
then run,
```bash
`root -q -b HscpCandidates.C
```

## Run the second part: read histograms & produce mass plots 

Configure configFile_readHist.txt respecting dedicated columns
```bash
root file       rebin   rebin eta   rebin Ih    rebin p     rebin mass  variable bins mass      variable bins p     threshold mass      threshold p
```
rebin is a boolean in order to know if we want rebin the different distributions. Then we give the values of rebinning.
For a variable rebinning, we give a threshold.

then run,
```bash
`root -q -b readHisto.C
```

## Make plots

Use the script makePlots.py in order to produce automatically some plots

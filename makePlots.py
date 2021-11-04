import subprocess
import os,subprocess,sys

listFiles = ["outfile_ias25_pt60_ih0_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed"
        ,"outfile_ias25_pt60_ih29_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed"
        ,"outfile_ias25_pt60_ih34_p-1_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed"
        ,"outfile_ias25_pt60_ih34_p4000_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed"
        ,"outfile_ias25_pt60_ih34_p4000_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_pThresholdVarBins80_analysed"
        ,"outfile_ias25_pt60_ih34_p4000_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_massThresholdVarBins500_analysed"
        ,"outfile_ias25_pt60_ih34_p4000_etabins120_ihbins1000_pbins2000_massbins2000_invIso1_invMET0_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed"
        ,"outfile_ias25_pt60_ih34_p4000_etabins120_ihbins1000_pbins2000_massbins2000_invIso0_invMET1_rebinEta2_rebinIh5_rebinP10_rebinMass25_analysed"]

for f in listFiles:
    command = ["rootprint -d dir_"+f+" -D colz "+f+".root:ih_p_all"+" &"]
    process = subprocess.Popen(command,shell=True)

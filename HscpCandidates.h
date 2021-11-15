//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep 15 15:27:19 2021 by ROOT version 6.14/09
// from TTree HscpCandidates/HscpCandidates
// found on file: /opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/resultSingleMu_TkOnly_UL2017C_v2-4_tree/histoSingleMu_TkOnly_UL2017C_v2-4_tree.root
//////////////////////////////////////////////////////////

#ifndef HscpCandidates_h
#define HscpCandidates_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TRandom3.h"

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

#include <iostream>
#include <fstream>


// Values determined by Dylan 
float K(2.37), C(2.93);

// Function returning the MassErr as function momentum, dEdx, and errors on momentum and dEdx
// Not take into account any erros coming from K&C factors because this function is used to see the impact of binning in p and dedx on mass error
double GetMassErr (double P, double PErr, double dEdx, double dEdxErr, double M, double dEdxK, double dEdxC)
{
   if (M < 0) return -1;

   double Criteria = dEdx - dEdxC;
   double MassErr = sqrt(pow(P*dEdxErr,2)/(4*dEdxK*Criteria)+(Criteria/dEdxK)*pow(PErr,2));

   if (std::isnan(MassErr) || std::isinf(MassErr)) MassErr = -1;

   return MassErr/M;
}

// Function returning an 1D-histogram which contains the ratio between two histograms   
TH1F* ratioHist(TH1F* h1, TH1F* h2)
{
    TH1F* res = (TH1F*) h1->Clone();
    for(int i=0;i<h1->GetNbinsX();i++)
    {
        double ratio = h2->GetBinContent(i)>0 ? h1->GetBinContent(i)/h2->GetBinContent(i) : 0;
        res->SetBinContent(i,ratio);
    }
    res->Divide(h2);
    res->Sumw2();
    return res;
}

// Return the mass as a function of momentum, dEdx, K and C. 
// It corresponds to the Bethe-Bloch parametrisation used in the Hscp analysis
float GetMass(float& p, float& ih, float& k, float& c)
{
    return (ih-C)<0?-1:sqrt((ih-c)/k)*p;
}

// Scale the 1D-histogram given to the unit 
void scale(TH1F* h)
{
    h->Scale(1./h->Integral(0,h->GetNbinsX()+1));
}

// Inverse the scaling ; scale the 1D-histogram to its number of entries
void invScale(TH1F* h)
{
    h->Sumw2();
    h->Scale(h->GetEntries());
}

// Function returning the ratio of right integer (from x to infty) for two 1D-histograms
// This function is used in the Hscp data-driven background estimate to test the mass shape prediction
// The argument to use this type of ratio is that we're in case of cut & count experiment 
TH1F* ratioIntegral(TH1F* h1, TH1F* h2)
{
    TH1F* res = (TH1F*) h1->Clone(); res->Reset();
    for(int i=0;i<h1->GetNbinsX()+1;i++)
    {   
        double Perr=0, Derr=0;
        double P=h1->IntegralAndError(i,h1->GetNbinsX()+1,Perr); if(P<=0) continue;
        double D=h2->IntegralAndError(i,h2->GetNbinsX()+1,Derr);
        res->SetBinContent(i,D/P);
        res->SetBinError(i,sqrt(pow(Derr*P,2)+pow(Perr*D,2))/pow(P,2));
    }
    return res;
}

// Function returning chi2/ndof compatibility test for two 1D-histograms
// As a reference, we've access to the number of degrees of freedom 
float chi2test(TH1F* h1, TH1F* h2,int& dof)
{
    float res=0;
    int ndof=0;
    for(int i=1;i<h1->GetNbinsX();i++)
    {
        res += h2->GetBinContent(i)>0 ? pow((h1->GetBinContent(i)-h2->GetBinContent(i)),2)/h2->GetBinContent(i) : 0 ;
        if(h2->GetBinContent(i)>0) ndof++;
    }
    dof=ndof;
    return res/ndof;
}

// Function returning a float number corresponding to a linear interpolation for a certain value x
// FIXME WIP 
float invExpoInterpolate(TH1F* h1, float x, std::vector<float> vect)
{
    Int_t xbin = h1->FindBin(x);
    Double_t x0,x1,y0,y1;
           
    if(x<=h1->GetBinCenter(1)) {
        return h1->GetBinContent(1);
    } else if(x>=h1->GetBinCenter(h1->GetNbinsX())) {
        return h1->GetBinContent(h1->GetNbinsX());
    } else {
        /*if(x<=h1->GetBinCenter(xbin)) {
            y0 = h1->RetrieveBinContent(xbin-1);
            x0 = h1->GetBinCenter(xbin-1);
            y1 = h1->RetrieveBinContent(xbin);
            x1 = h1->GetBinCenter(xbin);
        } else {
            y0 = h1->RetrieveBinContent(xbin);
            x0 = h1->GetBinCenter(xbin);
            y1 = h1->RetrieveBinContent(xbin+1);
            x1 = h1->GetBinCenter(xbin+1);
        }*/
        y0 = h1->GetBinContent(xbin-2);
        x0 = vect[xbin-2];
        y1 = h1->GetBinContent(xbin-1);
        x1 = vect[xbin-1];
        
        return y0 + (x-x0)*((y1-y0)/(x1-x0)) > 0 ? y0 + (x-x0)*((y1-y0)/(x1-x0)) : 0;
    }
}

// Function returning a 1D-histogram corresponding to interpolation of 1D-histo
// FIXME WIP 
TH1F* histoInterpolated(TH1F* h1)
{
    TH1F* h2 = (TH1F*) h1->Clone(); h2->Reset();
    for(int i=0;i<h1->GetNbinsX()+1;i++)
    {
        h2->SetBinContent(i,h1->Interpolate(h1->GetBinCenter(i)));
    }
    return h2;
}

// Function returning a canvas divide in two windows
// The first one contains the two 1D-histograms, given as arguments, superposed.
// There also is a legend associated to this window where the names are defined as arguments.
// The second window contains the ratio of these 1D-histograms or the ratio of right integers of them. 
// We define which kind of ratio we want with tha 'ratioSimple' boolean.
// The 'name' given corresponds to the name of the canvas 
TCanvas* plotting(TH1F* h1, TH1F* h2, bool ratioSimple=true, std::string name="", std::string leg1="", std::string leg2="")
{
    h1->Sumw2(); h2->Sumw2();
    TCanvas* c1 = new TCanvas(("plotting_"+name).c_str(),"");
    c1->Divide(1,2);
    gStyle->SetOptStat(0);
    c1->cd(1);
    TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(h1,leg1.c_str(),"lep");
    leg->AddEntry(h2,leg2.c_str(),"lep");
    h1->Draw();
    h2->SetLineColor(2);
    h2->Draw("esame");
    leg->Draw("same");
    c1->SetLogy();
    c1->cd(2);
    TH1F* tmp = (TH1F*) h1->Clone(); tmp->Reset();
    if(ratioSimple)
    {
        h1->Scale(1./h1->Integral());
        h2->Scale(1./h2->Integral());
        tmp = (TH1F*)h1->Clone();
        tmp->Divide(h2);
        tmp->GetYaxis()->SetTitle("#frac{N_{obs}}{N_{pred}}");
        tmp->GetYaxis()->SetTitleSize(0.06);
    }
    else
    {
        h1->Scale(1./h1->Integral());
        h2->Scale(1./h2->Integral());
        tmp=ratioIntegral(h2,h1);
        tmp->GetYaxis()->SetTitle("#int_{M}^{#infty} dm_{obs} / #int_{M}^{#infty} dm_{pred}");
        tmp->GetYaxis()->SetTitleSize(0.06);
    }

    int dof;
    float chi2=chi2test(h1,h2,dof);
    tmp->GetYaxis()->SetRangeUser(0,2);
    tmp->Draw();
    return c1;
}

// Function returning a 1D-histogram with variable binning.
// To determine it, we set a minimum value 'min' corresponding to the number of entries wanted in each bins of the distribution.
// If a given bin has a number of entries below than 'min' then this bin is merged with the next one and then there is another test, etc. 
// Finally, as a reference argument, the function gives a vector composed by the binning previously determined 
TH1F* rebinning(TH1F* h1, float min, std::vector<double>& vect)
{
    std::vector<double> v_val;
    int i=0;
    while(i<h1->GetNbinsX()+1)
    {
        float cont=0;
        while(cont<min && i<h1->GetNbinsX()+1)
        {
            cont += h1->GetBinContent(i);
            i++;
        }
        v_val.push_back(h1->GetBinLowEdge(i));
    }

    TH1F* h2 = new TH1F("Rebinned","",v_val.size()-1,v_val.data());
    for(int j=0;j<h1->GetNbinsX()+1;j++)
    {
        for(int c=0;c<h1->GetBinContent(j);c++)
        {
            h2->Fill(h1->GetBinLowEdge(j));
        }
    }
    h2->Sumw2();
    vect = v_val;
    return h2;
}

// Function doing the eta reweighing between two 2D-histograms as done in the Hscp background estimate method,
// because of the correlation between variables (momentum & transverse momentum). 
// The first given 2D-histogram is weighted in respect to the second one. 
void etaReweighingP(TH2F* eta_p_1, TH2F* eta_p_2)
{
    TH1F* eta1 = (TH1F*) eta_p_1->ProjectionY(); eta1->Scale(1./eta1->Integral());
    TH1F* eta2 = (TH1F*) eta_p_2->ProjectionY(); eta2->Scale(1./eta2->Integral());
    eta2->Divide(eta1);
    for(int i=0;i<eta_p_1->GetNbinsX()+1;i++)
    {
        for(int j=0;j<eta_p_1->GetNbinsY()+1;j++)
        {
            float val_ij = eta_p_1->GetBinContent(i,j);
            float err_ij = eta_p_1->GetBinError(i,j);
            eta_p_1->SetBinContent(i,j,val_ij*eta2->GetBinContent(j));
            eta_p_1->SetBinError(i,j,err_ij*eta2->GetBinContent(j));
        }
    }
    eta_p_1->Sumw2();
}

//FIXME Still necessary ? 
void predMassEtaBinning(TH1F* res, TH2F* eta_p, TH2F* ih_eta, float norm=0)
{
    TH1F* eta = (TH1F*) ih_eta->ProjectionX();
    for(int i=1;i<eta->GetNbinsX();i++)
    {
        TH1F* p = (TH1F*) eta_p->ProjectionX("proj_p",i,i);
        TH1F* ih = (TH1F*) ih_eta->ProjectionY("proj_ih",i,i);
        scale(p);
        for(int j=1;j<p->GetNbinsX()+1;j++)
        {
            for(int k=1;k<ih->GetNbinsX()+1;k++)
            {
                if(p->GetBinContent(j)<0) continue;
                if(ih->GetBinContent(k)<0) continue;
                float mom = p->GetBinCenter(j);
                float dedx = ih->GetBinCenter(k);
                float prob = p->GetBinContent(j) * ih->GetBinContent(k);
                if(prob>=0)
                {
                    res->Fill(GetMass(mom,dedx,K,C),prob*norm);
                }
            }
        }
        delete p;
        delete ih;
    }
    res->Sumw2();
    //scale(massFrom1DTemplatesEtaBinning);
}

void crossHistos(TH2F* res, TH1F* h1, TH1F* h2)
{
    scale(h1); 
    for(int i=0;i<h1->GetNbinsX()+1;i++)
    {
        for(int j=0;j<h2->GetNbinsX()+1;j++)
        {   

            float mom = h1->GetBinCenter(i);
            float dedx = h2->GetBinCenter(j);
            double prob = h1->GetBinContent(i)*h2->GetBinContent(j);
            if(prob>=0)
            {
                res->Fill(mom,dedx,prob);
                /*for(int l=0;l<prob;l++)
                {
                    res->Fill(mom,dedx);
                }*/
            }
        }
    }
    res->Sumw2();
}

void crossHistosEtaBinning(TH2F* res, TH2F* eta_p, TH2F* ih_eta)
{
    TH1F* eta = (TH1F*) ih_eta->ProjectionX();
    for(int i=0;i<eta->GetNbinsX()+1;i++)
    {
        TH1F* p = (TH1F*) eta_p->ProjectionX("proj_p",i,i+1);
        TH1F* ih = (TH1F*) ih_eta->ProjectionY("proj_ih",i,i+1);
        scale(p);
        for(int j=0;j<p->GetNbinsX()+1;j++)
        {
            for(int k=0;k<ih->GetNbinsX()+1;k++)
            {
                float mom = p->GetBinCenter(j);
                float dedx = ih->GetBinCenter(k);
                float prob = p->GetBinContent(j) * ih->GetBinContent(k);
                float err_prob = prob*sqrt((1./(ih->GetBinContent(k)))+(1./(p->GetBinContent(j)*ih->Integral())));
                if(prob>=0)
                {
                    //res->Fill(mom,dedx,prob);
                    /*for(int l=0;l<prob;l++)
                    {   
                        res->Fill(mom,dedx);
                    }*/
                    res->SetBinContent(j,k,res->GetBinContent(j,k)+prob);
                    res->SetBinError(j,k,sqrt(pow(res->GetBinError(j,k),2)+pow(err_prob,2)));

                }
                //std::cout << p->GetBinContent(j) << std::endl;
                //res->SetBinContent(i,j,prob);
            }
        }
        delete p;
        delete ih;
    }
    res->Sumw2();
}

void mapOfDifferences(TH2F* res, TH2F* h1, TH2F* h2)
{
    for(int i=0;i<h1->GetNbinsX();i++)
    {
        for(int j=0;j<h1->GetNbinsY();j++)
        {
            //double diff = h1->GetBinContent(i,j)>0 ? (h2->GetBinContent(i,j))/h1->GetBinContent(i,j) : 0;
            //
            double err=h1->GetBinError(i,j);
            //double diff = h2->GetBinContent(i,j)-h1->GetBinContent(i,j);
            double diff=h1->GetBinContent(i,j)>0?abs((h2->GetBinContent(i,j)-h1->GetBinContent(i,j))/h1->GetBinContent(i,j)):0;
            if(diff>=0) res->SetBinContent(i,j,diff);
        }
    }
}

class region{

    public:
        int nbins;
        float* xbins;
        int np;
        double* xp;
        float plow;
        float pup;
        int npt;
        float ptlow;
        float ptup;
        int nih;
        float ihlow;
        float ihup;
        int nias;
        float iaslow;
        float iasup;
        int neta;
        float etalow;
        float etaup;
        int nmass;
        float masslow;
        float massup;
        std::vector<float> vect;
        region();
        region(std::string suffix,int& etabins,int& ihbins,int& pbins,int& massbins);
        region(std::string suffix,int nbins, float* xbins, std::vector<double> v_pbins, std::vector<float> vect);
        ~region();
        void initHisto();
        void initHisto(int& etabins,int& ihbins,int& pbins,int& massbins);
        void rebinEtaP(std::vector<double> vect);
        void fill(float& eta, float& nhits, float&p, float& pt, float& ih, float& ias, float& m);
        void fillStdDev();
        void fillQuantile();
        void fillMassFrom1DTemplates();
        void fillMassFrom1DTemplatesEtaBinning();
        void plotMass();
        void cross1D();
        void write();
        void setBins(int nb,float* b);
        void rebinQuantiles(int rebin);
        std::string suffix_;
        TCanvas* c;
        TH2F* ih_pt;
        TH2F* ias_pt;
        TH2F* ih_ias;
        TH2F* ih_nhits;
        TH2F* ias_nhits;
        TH2F* eta_pt;
        TH2F* eta_p;
        TH2F* nhits_pt;
        TH2F* eta_nhits;
        TH2F* ih_eta;
        TH2F* ih_p;
        TH2F* ias_p;
        TH1F* stdDevIh_p;
        TH1F* stdDevIas_pt;
        TH1F* stdDevIas_p;
        TH1F* stdDevIas_p_y;
        TH1F* quantile99Ih_p;
        TH1F* quantile90Ih_p;
        TH1F* quantile70Ih_p;
        TH1F* quantile50Ih_p;
        TH1F* quantile30Ih_p;
        TH1F* quantile10Ih_p;
        TH1F* quantile01Ih_p;
        TH1F* quantile99Ias_pt;
        TH1F* quantile90Ias_pt;
        TH1F* quantile70Ias_pt;
        TH1F* quantile50Ias_pt;
        TH1F* quantile30Ias_pt;
        TH1F* quantile10Ias_pt;
        TH1F* quantile01Ias_pt;
        TH1F* quantile99Ias_p;
        TH1F* quantile90Ias_p;
        TH1F* quantile70Ias_p;
        TH1F* quantile50Ias_p;
        TH1F* quantile30Ias_p;
        TH1F* quantile10Ias_p;
        TH1F* quantile01Ias_p;
        TH1F* mass;
        TH1F* massFrom1DTemplates;
        TH1F* massFrom1DTemplatesEtaBinning;
        TH2F* eta_p_rebinned;
        TH3F* ih_p_eta;
        std::vector<double> VectOfBins_P_;
        TH1F* errMass;
        TH2F* Mass_errMass;
        TH2F* cross1Dtemplates;
        TH1F* ih_used;
        TH2F* mapM800;

};

region::region()
{
    ih_pt       = 0;
    ias_pt      = 0;
    ih_ias      = 0;
    ih_nhits    = 0;
    ias_nhits   = 0;
    eta_pt      = 0;
    eta_p       = 0;
    nhits_pt    = 0;
    eta_nhits   = 0;
    ih_eta      = 0;
    ih_p        = 0;
    ias_p       = 0;
    stdDevIh_p  = 0;
    stdDevIas_pt= 0;
    stdDevIas_p = 0;
    stdDevIas_p_y = 0;
    mass        = 0;
    massFrom1DTemplates = 0;
    massFrom1DTemplatesEtaBinning = 0;
    eta_p_rebinned = 0;
    quantile99Ih_p      = 0;
    quantile90Ih_p      = 0;
    quantile70Ih_p      = 0;
    quantile50Ih_p      = 0;
    quantile30Ih_p      = 0;
    quantile10Ih_p      = 0;
    quantile01Ih_p      = 0;
    quantile99Ias_pt    = 0;
    quantile90Ias_pt    = 0;
    quantile70Ias_pt    = 0;
    quantile50Ias_pt    = 0;
    quantile30Ias_pt    = 0;
    quantile10Ias_pt    = 0;
    quantile01Ias_pt    = 0;
    quantile99Ias_p    = 0;
    quantile90Ias_p    = 0;
    quantile70Ias_p    = 0;
    quantile50Ias_p    = 0;
    quantile30Ias_p    = 0;
    quantile10Ias_p    = 0;
    quantile01Ias_p    = 0;

    ih_p_eta            = 0;

    Mass_errMass    = 0;
    cross1Dtemplates= 0;
    ih_used         = 0;
    mapM800         = 0;
    initHisto();
} 

region::region(std::string suffix,int& etabins,int& ihbins,int& pbins,int& massbins)
{

    suffix_ = suffix;

    ih_pt       = 0;
    ias_pt      = 0;
    ih_ias      = 0;
    ih_nhits    = 0;
    ias_nhits   = 0;
    eta_pt      = 0;
    eta_p       = 0;
    nhits_pt    = 0;
    eta_nhits   = 0;
    ih_eta      = 0;
    ih_p        = 0;
    ias_p       = 0;
    stdDevIh_p  = 0;
    stdDevIas_pt= 0;
    stdDevIas_p = 0;
    stdDevIas_p_y = 0;
    mass        = 0;
    massFrom1DTemplates = 0;
    massFrom1DTemplatesEtaBinning = 0;
    eta_p_rebinned = 0;
    quantile99Ih_p      = 0;
    quantile90Ih_p      = 0;
    quantile70Ih_p      = 0;
    quantile50Ih_p      = 0;
    quantile30Ih_p      = 0;
    quantile10Ih_p      = 0;
    quantile01Ih_p      = 0;
    quantile99Ias_pt    = 0;
    quantile90Ias_pt    = 0;
    quantile70Ias_pt    = 0;
    quantile50Ias_pt    = 0;
    quantile30Ias_pt    = 0;
    quantile10Ias_pt    = 0;
    quantile01Ias_pt    = 0;
    quantile99Ias_p    = 0;
    quantile90Ias_p    = 0;
    quantile70Ias_p    = 0;
    quantile50Ias_p    = 0;
    quantile30Ias_p    = 0;
    quantile10Ias_p    = 0;
    quantile01Ias_p    = 0;

    ih_p_eta            = 0;

    Mass_errMass    = 0;
    cross1Dtemplates= 0;
    ih_used         = 0;
    mapM800         = 0;

    initHisto(etabins,ihbins,pbins,massbins);
} 

region::region(std::string suffix,int nbins_,float* xbins_,std::vector<double> v_pbins,std::vector<float> vect_)
{
    nbins = nbins_;
    xbins = xbins_;


    //np = 2000;
    plow = 0;
    pup = 4000;

    //npt = 2000;
    ptlow = 0;
    ptup = 4000;

    //nih = 500;
    ihlow = 0;
    ihup = 10;

    //nias = 500;
    iaslow = 0;
    iasup = 1;

    //neta = 120;
    etalow = -3;
    etaup = 3;

    //nmass = 2000;
    masslow = 0;
    massup = 4000;


    suffix_ = suffix;

    vect = vect_;

    c = new TCanvas(suffix.c_str(),"");

    ih_p_eta = new TH3F(("ih_p_eta"+suffix).c_str(),";#eta;p [GeV];I_{h} [MeV/cm]",neta,etalow,etaup,np,plow,pup,nih,ihlow,ihup);

    ih_pt = new TH2F(("ih_pt"+suffix).c_str(),";pt [GeV];I_{h} [MeV/cm]",npt,ptlow,ptup,nih,ihlow,ihup);
    ias_pt = new TH2F(("ias_pt"+suffix).c_str(),";pt [GeV];I_{as}",npt,ptlow,ptup,nias,iaslow,iasup);
    ih_ias = new TH2F(("ias_ih"+suffix).c_str(),";I_{as};I_{h} [MeV/cm]",nias,iaslow,iasup,nih,ihlow,ihup);
    ih_nhits = new TH2F(("ih_nhits"+suffix).c_str(),";nhits;I_{h} [MeV/cm]",20,0,20,nih,ihlow,ihup);
    ias_nhits = new TH2F(("ias_nhits"+suffix).c_str(),";nhits;I_{as}",20,0,20,nias,iaslow,iasup);
    eta_pt = new TH2F(("eta_pt"+suffix).c_str(),";pt [GeV];#eta",npt,ptlow,ptup,neta,etalow,etaup);
    eta_p = new TH2F(("eta_p"+suffix).c_str(),";p [GeV];#eta",np,plow,pup,neta,etalow,etaup);
    eta_p_rebinned = 0;
    //eta_p = new TH2F(("eta_p"+suffix).c_str(),";p [GeV];#eta",v_pbins.size()-1,v_pbins.data(),60,-3,3);
    nhits_pt = new TH2F(("nhits_pt"+suffix).c_str(),";pt [GeV];nhits",npt,ptlow,ptup,20,0,20);
    eta_nhits = new TH2F(("eta_nhits"+suffix).c_str(),";nhits;#eta",20,0,20,neta,etalow,etaup);
    ih_eta = new TH2F(("ih_eta"+suffix).c_str(),";#eta;I_{h} [MeV/cm]",neta,etalow,etaup,nih,ihlow,ihup);
    ih_p = new TH2F(("ih_p"+suffix).c_str(),";p [GeV];I_{h} [MeV/cm]",np,plow,pup,nih,ihlow,ihup);
    //ih_p = new TH2F(("ih_p"+suffix).c_str(),";p [GeV];I_{h} [MeV/cm]",v_pbins.size()-1,v_pbins.data(),100,0,10);
    ias_p = new TH2F(("ias_p"+suffix).c_str(),";p [GeV];I_{as}",np,plow,pup,nias,iaslow,iasup);
    //ias_p = new TH2F(("ias_p"+suffix).c_str(),";p [GeV];I_{as}",v_pbins.size()-1,v_pbins.data(),100,0,1);
    stdDevIh_p = new TH1F(("stdDevIh_p"+suffix).c_str(),";p [GeV];StdDev Ih [MeV/cm]",np,plow,pup);
    //stdDevIh_p = new TH1F(("stdDevIh_p"+suffix).c_str(),";p [GeV];StdDev Ih [MeV/cm]",v_pbins.size()-1,v_pbins.data());
    stdDevIas_pt = new TH1F(("stdDevIas_pt"+suffix).c_str(),";pt [GeV];StdDev I_{as}",npt,ptlow,ptup);
    stdDevIas_p  = new TH1F(("stdDevIas_p"+suffix).c_str(),";p [GeV];StdDev I_{as}",np,plow,pup);
    stdDevIas_p_y  = new TH1F(("stdDevIas_p_y"+suffix).c_str(),";I_{as};StdDev p [GeV]",nias,iaslow,iasup);
    
    /*quantile99Ih_p = new TH1F(("quantile99Ih_p"+suffix).c_str(),";p [GeV];0.99-quantile Ih [MeV/cm]",np,xp);
    quantile90Ih_p = new TH1F(("quantile90Ih_p"+suffix).c_str(),";p [GeV];0.90-quantile Ih [MeV/cm]",np,xp);
    quantile70Ih_p = new TH1F(("quantile70Ih_p"+suffix).c_str(),";p [GeV];0.70-quantile Ih [MeV/cm]",np,xp);
    quantile50Ih_p = new TH1F(("quantile50Ih_p"+suffix).c_str(),";p [GeV];0.50-quantile Ih [MeV/cm]",np,xp);
    quantile30Ih_p = new TH1F(("quantile30Ih_p"+suffix).c_str(),";p [GeV];0.30-quantile Ih [MeV/cm]",np,xp);
    quantile01Ih_p = new TH1F(("quantile01Ih_p"+suffix).c_str(),";p [GeV];0.01-quantile Ih [MeV/cm]",np,xp);*/

    quantile99Ih_p = new TH1F(("quantile99Ih_p"+suffix).c_str(),";p [GeV];0.99-quantile Ih [MeV/cm]",np,plow,pup);
    quantile90Ih_p = new TH1F(("quantile90Ih_p"+suffix).c_str(),";p [GeV];0.90-quantile Ih [MeV/cm]",np,plow,pup);
    quantile70Ih_p = new TH1F(("quantile70Ih_p"+suffix).c_str(),";p [GeV];0.70-quantile Ih [MeV/cm]",np,plow,pup);
    quantile50Ih_p = new TH1F(("quantile50Ih_p"+suffix).c_str(),";p [GeV];0.50-quantile Ih [MeV/cm]",np,plow,pup);
    quantile30Ih_p = new TH1F(("quantile30Ih_p"+suffix).c_str(),";p [GeV];0.30-quantile Ih [MeV/cm]",np,plow,pup);
    quantile10Ih_p = new TH1F(("quantile10Ih_p"+suffix).c_str(),";p [GeV];0.10-quantile Ih [MeV/cm]",np,plow,pup);
    quantile01Ih_p = new TH1F(("quantile01Ih_p"+suffix).c_str(),";p [GeV];0.01-quantile Ih [MeV/cm]",np,plow,pup);

    quantile99Ias_pt = new TH1F(("quantile99Ias_pt"+suffix).c_str(),";pt [GeV];0.99-quantile I_{as}",npt,ptlow,ptup);
    quantile90Ias_pt = new TH1F(("quantile90Ias_pt"+suffix).c_str(),";pt [GeV];0.90-quantile I_{as}",npt,ptlow,ptup);
    quantile70Ias_pt = new TH1F(("quantile70Ias_pt"+suffix).c_str(),";pt [GeV];0.70-quantile I_{as}",npt,ptlow,ptup);
    quantile50Ias_pt = new TH1F(("quantile50Ias_pt"+suffix).c_str(),";pt [GeV];0.50-quantile I_{as}",npt,ptlow,ptup);
    quantile30Ias_pt = new TH1F(("quantile30Ias_pt"+suffix).c_str(),";pt [GeV];0.30-quantile I_{as}",npt,ptlow,ptup);
    quantile10Ias_pt = new TH1F(("quantile10Ias_pt"+suffix).c_str(),";pt [GeV];0.10-quantile I_{as}",npt,ptlow,ptup);
    quantile01Ias_pt = new TH1F(("quantile01Ias_pt"+suffix).c_str(),";pt [GeV];0.01-quantile I_{as}",npt,ptlow,ptup);

    quantile99Ias_p = new TH1F(("quantile99Ias_p"+suffix).c_str(),";p [GeV];0.99-quantile I_{as}",np,plow,pup);
    quantile90Ias_p = new TH1F(("quantile90Ias_p"+suffix).c_str(),";p [GeV];0.90-quantile I_{as}",np,plow,pup);
    quantile70Ias_p = new TH1F(("quantile70Ias_p"+suffix).c_str(),";p [GeV];0.70-quantile I_{as}",np,plow,pup);
    quantile50Ias_p = new TH1F(("quantile50Ias_p"+suffix).c_str(),";p [GeV];0.50-quantile I_{as}",np,plow,pup);
    quantile30Ias_p = new TH1F(("quantile30Ias_p"+suffix).c_str(),";p [GeV];0.30-quantile I_{as}",np,plow,pup);
    quantile10Ias_p = new TH1F(("quantile10Ias_p"+suffix).c_str(),";p [GeV];0.10-quantile I_{as}",np,plow,pup);
    quantile01Ias_p = new TH1F(("quantile01Ias_p"+suffix).c_str(),";p [GeV];0.01-quantile I_{as}",np,plow,pup);

    
    mass = new TH1F(("massFromTree"+suffix).c_str(),";Mass [GeV]",nmass,masslow,massup);
    massFrom1DTemplates = new TH1F(("massFrom1DTemplates"+suffix).c_str(),";Mass [GeV]",nmass,masslow,massup);
    massFrom1DTemplatesEtaBinning = new TH1F(("massFrom1DTemplatesEtaBinning"+suffix).c_str(),";Mass [GeV]",nmass,masslow,massup);


    errMass = new TH1F(("errMass"+suffix).c_str(),";Mass error [GeV]",nmass,masslow,massup);

    Mass_errMass = new TH2F(("Mass_errMass"+suffix).c_str(),";Mass [GeV];Mass error [GeV]",nmass,masslow,massup,nmass,masslow,massup);
    
    cross1Dtemplates = new TH2F(("cross1Dtemplates_ih_p_"+suffix).c_str(),";p [GeV];I_{h} [MeV/cm]",np,plow,pup,nih,ihlow,ihup);

    ih_used         = new TH1F(("ih_used"+suffix).c_str(),";I_{h} [MeV/cm];",nih,ihlow,ihup);
    mapM800         = new TH2F(("mapM800_ih_p_"+suffix).c_str(),";p [GeV];I_{h} [MeV/cm]",np,plow,pup,nih,ihlow,ihup);

    //mass = new TH1F(("massFromTree"+suffix).c_str(),";Mass [GeV]",nbins,xbins);
    //massFrom1DTemplates = new TH1F(("massFrom1DTemplates"+suffix).c_str(),";Mass [GeV]",nbins,xbins);
    //massFrom1DTemplatesEtaBinning = new TH1F(("massFrom1DTemplatesEtaBinning"+suffix).c_str(),";Mass [GeV]",nbins,xbins);

}   

region::~region()
{
    /*delete ih_p_eta;
    delete ih_pt;
    delete ias_pt;
    delete ih_ias;
    delete ih_nhits;
    delete ias_nhits;
    delete eta_pt;
    delete eta_p;
    delete nhits_pt;
    delete eta_nhits;
    delete ih_eta;
    delete ih_p;
    delete ias_p;
    delete stdDevIh_p;
    delete stdDevIas_pt;
    delete stdDevIas_p;
    delete stdDevIas_p_y;
    delete quantile01Ih_p;
    delete quantile10Ih_p;
    delete quantile30Ih_p;
    delete quantile50Ih_p;
    delete quantile70Ih_p;
    delete quantile90Ih_p;
    delete quantile99Ih_p;
    delete quantile01Ias_pt;
    delete quantile10Ias_pt;
    delete quantile30Ias_pt;
    delete quantile50Ias_pt;
    delete quantile70Ias_pt;
    delete quantile90Ias_pt;
    delete quantile99Ias_pt;
    delete quantile01Ias_p;
    delete quantile10Ias_p;
    delete quantile30Ias_p;
    delete quantile50Ias_p;
    delete quantile70Ias_p;
    delete quantile90Ias_p;
    delete quantile99Ias_p;

    delete mass;
    delete massFrom1DTemplates;
    delete massFrom1DTemplatesEtaBinning;
    delete Mass_errMass;
    delete cross1Dtemplates;
    delete c;*/
}

void region::initHisto()
{
    np = 2000;
    plow = 0;
    pup = 4000;

    npt = 2000;
    ptlow = 0;
    ptup = 4000;

    nih = 1000;
    ihlow = 0;
    ihup = 20;

    nias = 500;
    iaslow = 0;
    iasup = 1;

    neta = 120;
    etalow = -3;
    etaup = 3;

    nmass = 2000;
    masslow = 0;
    massup = 4000;


    std::string suffix = suffix_;

    c = new TCanvas(suffix.c_str(),"");

    ih_p_eta = new TH3F(("ih_p_eta"+suffix).c_str(),";#eta;p [GeV];I_{h} [MeV/cm]",neta,etalow,etaup,np,plow,pup,nih,ihlow,ihup);

    ih_pt = new TH2F(("ih_pt"+suffix).c_str(),";pt [GeV];I_{h} [MeV/cm]",npt,ptlow,ptup,nih,ihlow,ihup);
    ias_pt = new TH2F(("ias_pt"+suffix).c_str(),";pt [GeV];I_{as}",npt,ptlow,ptup,nias,iaslow,iasup);
    ih_ias = new TH2F(("ias_ih"+suffix).c_str(),";I_{as};I_{h} [MeV/cm]",nias,iaslow,iasup,nih,ihlow,ihup);
    ih_nhits = new TH2F(("ih_nhits"+suffix).c_str(),";nhits;I_{h} [MeV/cm]",20,0,20,nih,ihlow,ihup);
    ias_nhits = new TH2F(("ias_nhits"+suffix).c_str(),";nhits;I_{as}",20,0,20,nias,iaslow,iasup);
    eta_pt = new TH2F(("eta_pt"+suffix).c_str(),";pt [GeV];#eta",npt,ptlow,ptup,neta,etalow,etaup);
    eta_p = new TH2F(("eta_p"+suffix).c_str(),";p [GeV];#eta",np,plow,pup,neta,etalow,etaup);
    eta_p_rebinned = 0;
    //eta_p = new TH2F(("eta_p"+suffix).c_str(),";p [GeV];#eta",v_pbins.size()-1,v_pbins.data(),60,-3,3);
    nhits_pt = new TH2F(("nhits_pt"+suffix).c_str(),";pt [GeV];nhits",npt,ptlow,ptup,20,0,20);
    eta_nhits = new TH2F(("eta_nhits"+suffix).c_str(),";nhits;#eta",20,0,20,neta,etalow,etaup);
    ih_eta = new TH2F(("ih_eta"+suffix).c_str(),";#eta;I_{h} [MeV/cm]",neta,etalow,etaup,nih,ihlow,ihup);
    ih_p = new TH2F(("ih_p"+suffix).c_str(),";p [GeV];I_{h} [MeV/cm]",np,plow,pup,nih,ihlow,ihup);
    //ih_p = new TH2F(("ih_p"+suffix).c_str(),";p [GeV];I_{h} [MeV/cm]",v_pbins.size()-1,v_pbins.data(),100,0,10);
    ias_p = new TH2F(("ias_p"+suffix).c_str(),";p [GeV];I_{as}",np,plow,pup,nias,iaslow,iasup);
    //ias_p = new TH2F(("ias_p"+suffix).c_str(),";p [GeV];I_{as}",v_pbins.size()-1,v_pbins.data(),100,0,1);
    stdDevIh_p = new TH1F(("stdDevIh_p"+suffix).c_str(),";p [GeV];StdDev Ih [MeV/cm]",np,plow,pup);
    //stdDevIh_p = new TH1F(("stdDevIh_p"+suffix).c_str(),";p [GeV];StdDev Ih [MeV/cm]",v_pbins.size()-1,v_pbins.data());
    stdDevIas_pt = new TH1F(("stdDevIas_pt"+suffix).c_str(),";pt [GeV];StdDev I_{as}",npt,ptlow,ptup);
    stdDevIas_p  = new TH1F(("stdDevIas_p"+suffix).c_str(),";p [GeV];StdDev I_{as}",np,plow,pup);
    stdDevIas_p_y  = new TH1F(("stdDevIas_p_y"+suffix).c_str(),";I_{as};StdDev p [GeV]",nias,iaslow,iasup);
    
    /*quantile99Ih_p = new TH1F(("quantile99Ih_p"+suffix).c_str(),";p [GeV];0.99-quantile Ih [MeV/cm]",np,xp);
    quantile90Ih_p = new TH1F(("quantile90Ih_p"+suffix).c_str(),";p [GeV];0.90-quantile Ih [MeV/cm]",np,xp);
    quantile70Ih_p = new TH1F(("quantile70Ih_p"+suffix).c_str(),";p [GeV];0.70-quantile Ih [MeV/cm]",np,xp);
    quantile50Ih_p = new TH1F(("quantile50Ih_p"+suffix).c_str(),";p [GeV];0.50-quantile Ih [MeV/cm]",np,xp);
    quantile30Ih_p = new TH1F(("quantile30Ih_p"+suffix).c_str(),";p [GeV];0.30-quantile Ih [MeV/cm]",np,xp);
    quantile01Ih_p = new TH1F(("quantile01Ih_p"+suffix).c_str(),";p [GeV];0.01-quantile Ih [MeV/cm]",np,xp);*/

    quantile99Ih_p = new TH1F(("quantile99Ih_p"+suffix).c_str(),";p [GeV];0.99-quantile Ih [MeV/cm]",np,plow,pup);
    quantile90Ih_p = new TH1F(("quantile90Ih_p"+suffix).c_str(),";p [GeV];0.90-quantile Ih [MeV/cm]",np,plow,pup);
    quantile70Ih_p = new TH1F(("quantile70Ih_p"+suffix).c_str(),";p [GeV];0.70-quantile Ih [MeV/cm]",np,plow,pup);
    quantile50Ih_p = new TH1F(("quantile50Ih_p"+suffix).c_str(),";p [GeV];0.50-quantile Ih [MeV/cm]",np,plow,pup);
    quantile30Ih_p = new TH1F(("quantile30Ih_p"+suffix).c_str(),";p [GeV];0.30-quantile Ih [MeV/cm]",np,plow,pup);
    quantile10Ih_p = new TH1F(("quantile10Ih_p"+suffix).c_str(),";p [GeV];0.10-quantile Ih [MeV/cm]",np,plow,pup);
    quantile01Ih_p = new TH1F(("quantile01Ih_p"+suffix).c_str(),";p [GeV];0.01-quantile Ih [MeV/cm]",np,plow,pup);

    quantile99Ias_pt = new TH1F(("quantile99Ias_pt"+suffix).c_str(),";pt [GeV];0.99-quantile I_{as}",npt,ptlow,ptup);
    quantile90Ias_pt = new TH1F(("quantile90Ias_pt"+suffix).c_str(),";pt [GeV];0.90-quantile I_{as}",npt,ptlow,ptup);
    quantile70Ias_pt = new TH1F(("quantile70Ias_pt"+suffix).c_str(),";pt [GeV];0.70-quantile I_{as}",npt,ptlow,ptup);
    quantile50Ias_pt = new TH1F(("quantile50Ias_pt"+suffix).c_str(),";pt [GeV];0.50-quantile I_{as}",npt,ptlow,ptup);
    quantile30Ias_pt = new TH1F(("quantile30Ias_pt"+suffix).c_str(),";pt [GeV];0.30-quantile I_{as}",npt,ptlow,ptup);
    quantile10Ias_pt = new TH1F(("quantile10Ias_pt"+suffix).c_str(),";pt [GeV];0.10-quantile I_{as}",npt,ptlow,ptup);
    quantile01Ias_pt = new TH1F(("quantile01Ias_pt"+suffix).c_str(),";pt [GeV];0.01-quantile I_{as}",npt,ptlow,ptup);

    quantile99Ias_p = new TH1F(("quantile99Ias_p"+suffix).c_str(),";p [GeV];0.99-quantile I_{as}",np,plow,pup);
    quantile90Ias_p = new TH1F(("quantile90Ias_p"+suffix).c_str(),";p [GeV];0.90-quantile I_{as}",np,plow,pup);
    quantile70Ias_p = new TH1F(("quantile70Ias_p"+suffix).c_str(),";p [GeV];0.70-quantile I_{as}",np,plow,pup);
    quantile50Ias_p = new TH1F(("quantile50Ias_p"+suffix).c_str(),";p [GeV];0.50-quantile I_{as}",np,plow,pup);
    quantile30Ias_p = new TH1F(("quantile30Ias_p"+suffix).c_str(),";p [GeV];0.30-quantile I_{as}",np,plow,pup);
    quantile10Ias_p = new TH1F(("quantile10Ias_p"+suffix).c_str(),";p [GeV];0.10-quantile I_{as}",np,plow,pup);
    quantile01Ias_p = new TH1F(("quantile01Ias_p"+suffix).c_str(),";p [GeV];0.01-quantile I_{as}",np,plow,pup);

    
    mass = new TH1F(("massFromTree"+suffix).c_str(),";Mass [GeV]",nmass,masslow,massup);
    massFrom1DTemplates = new TH1F(("massFrom1DTemplates"+suffix).c_str(),";Mass [GeV]",nmass,masslow,massup);
    massFrom1DTemplatesEtaBinning = new TH1F(("massFrom1DTemplatesEtaBinning"+suffix).c_str(),";Mass [GeV]",nmass,masslow,massup);
    

    errMass = new TH1F(("errMass"+suffix).c_str(),";Mass error [GeV]",nmass,masslow,massup);
    Mass_errMass = new TH2F(("Mass_errMass"+suffix).c_str(),";Mass [GeV];Mass error [GeV]",nmass,masslow,massup,nmass,masslow,massup);
    cross1Dtemplates = new TH2F(("cross1Dtemplates_ih_p_"+suffix).c_str(),";p [GeV];I_{h} [MeV/cm]",np,plow,pup,nih,ihlow,ihup);
    
    ih_used         = new TH1F(("ih_used"+suffix).c_str(),";I_{h} [MeV/cm];",nih,ihlow,ihup);
    mapM800         = new TH2F(("mapM800_ih_p_"+suffix).c_str(),";p [GeV];I_{h} [MeV/cm]",np,plow,pup,nih,ihlow,ihup);

    //mass = new TH1F(("massFromTree"+suffix).c_str(),";Mass [GeV]",nbins,xbins);
    //massFrom1DTemplates = new TH1F(("massFrom1DTemplates"+suffix).c_str(),";Mass [GeV]",nbins,xbins);
    //massFrom1DTemplatesEtaBinning = new TH1F(("massFrom1DTemplatesEtaBinning"+suffix).c_str(),";Mass [GeV]",nbins,xbins);

}

void region::initHisto(int& etabins,int& ihbins,int& pbins,int& massbins)
{
    np = pbins;
    plow = 0;
    pup = 4000;

    npt = pbins;
    ptlow = 0;
    ptup = 4000;

    nih = ihbins;
    ihlow = 0;
    ihup = 20;

    nias = ihbins;
    iaslow = 0;
    iasup = 1;

    neta = etabins;
    etalow = -3;
    etaup = 3;

    nmass = massbins;
    masslow = 0;
    massup = 4000;


    std::string suffix = suffix_;

    c = new TCanvas(suffix.c_str(),"");

    ih_p_eta = new TH3F(("ih_p_eta"+suffix).c_str(),";#eta;p [GeV];I_{h} [MeV/cm]",neta,etalow,etaup,np,plow,pup,nih,ihlow,ihup);

    ih_pt = new TH2F(("ih_pt"+suffix).c_str(),";pt [GeV];I_{h} [MeV/cm]",npt,ptlow,ptup,nih,ihlow,ihup);
    ias_pt = new TH2F(("ias_pt"+suffix).c_str(),";pt [GeV];I_{as}",npt,ptlow,ptup,nias,iaslow,iasup);
    ih_ias = new TH2F(("ias_ih"+suffix).c_str(),";I_{as};I_{h} [MeV/cm]",nias,iaslow,iasup,nih,ihlow,ihup);
    ih_nhits = new TH2F(("ih_nhits"+suffix).c_str(),";nhits;I_{h} [MeV/cm]",20,0,20,nih,ihlow,ihup);
    ias_nhits = new TH2F(("ias_nhits"+suffix).c_str(),";nhits;I_{as}",20,0,20,nias,iaslow,iasup);
    eta_pt = new TH2F(("eta_pt"+suffix).c_str(),";pt [GeV];#eta",npt,ptlow,ptup,neta,etalow,etaup);
    eta_p = new TH2F(("eta_p"+suffix).c_str(),";p [GeV];#eta",np,plow,pup,neta,etalow,etaup);
    eta_p_rebinned = 0;
    //eta_p = new TH2F(("eta_p"+suffix).c_str(),";p [GeV];#eta",v_pbins.size()-1,v_pbins.data(),60,-3,3);
    nhits_pt = new TH2F(("nhits_pt"+suffix).c_str(),";pt [GeV];nhits",npt,ptlow,ptup,20,0,20);
    eta_nhits = new TH2F(("eta_nhits"+suffix).c_str(),";nhits;#eta",20,0,20,neta,etalow,etaup);
    ih_eta = new TH2F(("ih_eta"+suffix).c_str(),";#eta;I_{h} [MeV/cm]",neta,etalow,etaup,nih,ihlow,ihup);
    ih_p = new TH2F(("ih_p"+suffix).c_str(),";p [GeV];I_{h} [MeV/cm]",np,plow,pup,nih,ihlow,ihup);
    //ih_p = new TH2F(("ih_p"+suffix).c_str(),";p [GeV];I_{h} [MeV/cm]",v_pbins.size()-1,v_pbins.data(),100,0,10);
    ias_p = new TH2F(("ias_p"+suffix).c_str(),";p [GeV];I_{as}",np,plow,pup,nias,iaslow,iasup);
    //ias_p = new TH2F(("ias_p"+suffix).c_str(),";p [GeV];I_{as}",v_pbins.size()-1,v_pbins.data(),100,0,1);
    stdDevIh_p = new TH1F(("stdDevIh_p"+suffix).c_str(),";p [GeV];StdDev Ih [MeV/cm]",np,plow,pup);
    //stdDevIh_p = new TH1F(("stdDevIh_p"+suffix).c_str(),";p [GeV];StdDev Ih [MeV/cm]",v_pbins.size()-1,v_pbins.data());
    stdDevIas_pt = new TH1F(("stdDevIas_pt"+suffix).c_str(),";pt [GeV];StdDev I_{as}",npt,ptlow,ptup);
    stdDevIas_p  = new TH1F(("stdDevIas_p"+suffix).c_str(),";p [GeV];StdDev I_{as}",np,plow,pup);
    stdDevIas_p_y  = new TH1F(("stdDevIas_p_y"+suffix).c_str(),";I_{as};StdDev p [GeV]",nias,iaslow,iasup);
    
    /*quantile99Ih_p = new TH1F(("quantile99Ih_p"+suffix).c_str(),";p [GeV];0.99-quantile Ih [MeV/cm]",np,xp);
    quantile90Ih_p = new TH1F(("quantile90Ih_p"+suffix).c_str(),";p [GeV];0.90-quantile Ih [MeV/cm]",np,xp);
    quantile70Ih_p = new TH1F(("quantile70Ih_p"+suffix).c_str(),";p [GeV];0.70-quantile Ih [MeV/cm]",np,xp);
    quantile50Ih_p = new TH1F(("quantile50Ih_p"+suffix).c_str(),";p [GeV];0.50-quantile Ih [MeV/cm]",np,xp);
    quantile30Ih_p = new TH1F(("quantile30Ih_p"+suffix).c_str(),";p [GeV];0.30-quantile Ih [MeV/cm]",np,xp);
    quantile01Ih_p = new TH1F(("quantile01Ih_p"+suffix).c_str(),";p [GeV];0.01-quantile Ih [MeV/cm]",np,xp);*/

    quantile99Ih_p = new TH1F(("quantile99Ih_p"+suffix).c_str(),";p [GeV];0.99-quantile Ih [MeV/cm]",np,plow,pup);
    quantile90Ih_p = new TH1F(("quantile90Ih_p"+suffix).c_str(),";p [GeV];0.90-quantile Ih [MeV/cm]",np,plow,pup);
    quantile70Ih_p = new TH1F(("quantile70Ih_p"+suffix).c_str(),";p [GeV];0.70-quantile Ih [MeV/cm]",np,plow,pup);
    quantile50Ih_p = new TH1F(("quantile50Ih_p"+suffix).c_str(),";p [GeV];0.50-quantile Ih [MeV/cm]",np,plow,pup);
    quantile30Ih_p = new TH1F(("quantile30Ih_p"+suffix).c_str(),";p [GeV];0.30-quantile Ih [MeV/cm]",np,plow,pup);
    quantile10Ih_p = new TH1F(("quantile10Ih_p"+suffix).c_str(),";p [GeV];0.10-quantile Ih [MeV/cm]",np,plow,pup);
    quantile01Ih_p = new TH1F(("quantile01Ih_p"+suffix).c_str(),";p [GeV];0.01-quantile Ih [MeV/cm]",np,plow,pup);

    quantile99Ias_pt = new TH1F(("quantile99Ias_pt"+suffix).c_str(),";pt [GeV];0.99-quantile I_{as}",npt,ptlow,ptup);
    quantile90Ias_pt = new TH1F(("quantile90Ias_pt"+suffix).c_str(),";pt [GeV];0.90-quantile I_{as}",npt,ptlow,ptup);
    quantile70Ias_pt = new TH1F(("quantile70Ias_pt"+suffix).c_str(),";pt [GeV];0.70-quantile I_{as}",npt,ptlow,ptup);
    quantile50Ias_pt = new TH1F(("quantile50Ias_pt"+suffix).c_str(),";pt [GeV];0.50-quantile I_{as}",npt,ptlow,ptup);
    quantile30Ias_pt = new TH1F(("quantile30Ias_pt"+suffix).c_str(),";pt [GeV];0.30-quantile I_{as}",npt,ptlow,ptup);
    quantile10Ias_pt = new TH1F(("quantile10Ias_pt"+suffix).c_str(),";pt [GeV];0.10-quantile I_{as}",npt,ptlow,ptup);
    quantile01Ias_pt = new TH1F(("quantile01Ias_pt"+suffix).c_str(),";pt [GeV];0.01-quantile I_{as}",npt,ptlow,ptup);

    quantile99Ias_p = new TH1F(("quantile99Ias_p"+suffix).c_str(),";p [GeV];0.99-quantile I_{as}",np,plow,pup);
    quantile90Ias_p = new TH1F(("quantile90Ias_p"+suffix).c_str(),";p [GeV];0.90-quantile I_{as}",np,plow,pup);
    quantile70Ias_p = new TH1F(("quantile70Ias_p"+suffix).c_str(),";p [GeV];0.70-quantile I_{as}",np,plow,pup);
    quantile50Ias_p = new TH1F(("quantile50Ias_p"+suffix).c_str(),";p [GeV];0.50-quantile I_{as}",np,plow,pup);
    quantile30Ias_p = new TH1F(("quantile30Ias_p"+suffix).c_str(),";p [GeV];0.30-quantile I_{as}",np,plow,pup);
    quantile10Ias_p = new TH1F(("quantile10Ias_p"+suffix).c_str(),";p [GeV];0.10-quantile I_{as}",np,plow,pup);
    quantile01Ias_p = new TH1F(("quantile01Ias_p"+suffix).c_str(),";p [GeV];0.01-quantile I_{as}",np,plow,pup);

    
    mass = new TH1F(("massFromTree"+suffix).c_str(),";Mass [GeV]",nmass,masslow,massup);
    massFrom1DTemplates = new TH1F(("massFrom1DTemplates"+suffix).c_str(),";Mass [GeV]",nmass,masslow,massup);
    massFrom1DTemplatesEtaBinning = new TH1F(("massFrom1DTemplatesEtaBinning"+suffix).c_str(),";Mass [GeV]",nmass,masslow,massup);
    

    errMass = new TH1F(("errMass"+suffix).c_str(),";Mass error [GeV]",nmass,masslow,massup);
    Mass_errMass = new TH2F(("Mass_errMass"+suffix).c_str(),";Mass [GeV];Mass error [GeV]",nmass,masslow,massup,nmass,masslow,massup);
    cross1Dtemplates = new TH2F(("cross1Dtemplates_ih_p_"+suffix).c_str(),";p [GeV];I_{h} [MeV/cm]",np,plow,pup,nih,ihlow,ihup);

    ih_used         = new TH1F(("ih_used"+suffix).c_str(),";I_{h} [MeV/cm];",nih,ihlow,ihup);
    mapM800         = new TH2F(("mapM800_ih_p_"+suffix).c_str(),";p [GeV];I_{h} [MeV/cm]",np,plow,pup,nih,ihlow,ihup);

    //mass = new TH1F(("massFromTree"+suffix).c_str(),";Mass [GeV]",nbins,xbins);
    //massFrom1DTemplates = new TH1F(("massFrom1DTemplates"+suffix).c_str(),";Mass [GeV]",nbins,xbins);
    //massFrom1DTemplatesEtaBinning = new TH1F(("massFrom1DTemplatesEtaBinning"+suffix).c_str(),";Mass [GeV]",nbins,xbins);

}


void region::setBins(int nb, float* b)
{
    nbins = nb;
    xbins = b;
}

void region::rebinEtaP(std::vector<double> vect)
{
    //eta_p_rebinned = new TH2F(("eta_p_rebinned"+suffix_).c_str(),";p [GeV];#eta",vect.size()-1,vect.data(),60,-3,3);
    eta_p_rebinned = (TH2F*) eta_p->Clone();
    eta_p_rebinned->Rebin(vect.size()-1,"",vect.data());
    eta_p_rebinned->Sumw2();
}

void region::fill(float& eta, float& nhits, float& p, float& pt, float& ih, float& ias, float& m)
{
   ih_p_eta->Fill(eta,p,ih);
   ih_pt->Fill(pt,ih);
   ias_pt->Fill(pt,ias);
   ih_ias->Fill(ias,ih);
   ih_nhits->Fill(nhits,ih);
   ias_nhits->Fill(nhits,ias);
   eta_pt->Fill(pt,eta);
   eta_p->Fill(p,eta);
   nhits_pt->Fill(pt,nhits);
   eta_nhits->Fill(nhits,eta);
   ih_eta->Fill(eta,ih);
   ih_p->Fill(p,ih);
   ias_p->Fill(p,ias);
   mass->Fill(m);
}

void region::fillStdDev()
{
    for(int i=0;i<ih_p->GetNbinsX();i++)
    {
        float stddev=ih_p->ProjectionY("",i,i+1)->GetStdDev();
        float stddeverr=ih_p->ProjectionY("",i,i+1)->GetStdDevError();
        stdDevIh_p->SetBinContent(i,stddev);
        stdDevIh_p->SetBinError(i,stddeverr);
    }
    //stdDevIh_p->Rebin(20);
    for(int i=0;i<ias_pt->GetNbinsX();i++)
    {
        float stddev=ias_pt->ProjectionY("",i,i+1)->GetStdDev();
        float stddeverr=ias_pt->ProjectionY("",i,i+1)->GetStdDevError();
        stdDevIas_pt->SetBinContent(i,stddev);
        stdDevIas_pt->SetBinError(i,stddeverr);
    }
    //stdDevIas_pt->Rebin(20);
    for(int i=0;i<ias_p->GetNbinsX();i++)
    {
        float stddev=ias_p->ProjectionY("",i,i+1)->GetStdDev();
        float stddeverr=ias_p->ProjectionY("",i,i+1)->GetStdDevError();
        stdDevIas_p->SetBinContent(i,stddev);
        stdDevIas_p->SetBinError(i,stddeverr);
    }
    //stdDevIas_p->Rebin(20);
    for(int j=0;j<ias_p->GetNbinsY();j++)
    {
        float stddev=ias_p->ProjectionX("",j,j+1)->GetStdDev();
        float stddeverr=ias_p->ProjectionX("",j,j+1)->GetStdDevError();
        stdDevIas_p_y->SetBinContent(j,stddev);
        stdDevIas_p_y->SetBinError(j,stddeverr);
    }
    //stdDevIas_p_y->Rebin(5);
}

void region::fillQuantile()
{
    int n_quan=7;
    double p[7]={0.01,0.10,0.30,0.50,0.70,0.90,0.99};
    double q[7];
    for(int i=0;i<ih_p->GetNbinsX();i++)
    {
        ih_p->ProjectionY("",i,i+1)->GetQuantiles(7,q,p);
        quantile01Ih_p->SetBinContent(i,q[0]);
        quantile10Ih_p->SetBinContent(i,q[1]);
        quantile30Ih_p->SetBinContent(i,q[2]);
        quantile50Ih_p->SetBinContent(i,q[3]);
        quantile70Ih_p->SetBinContent(i,q[4]);
        quantile90Ih_p->SetBinContent(i,q[5]);
        quantile99Ih_p->SetBinContent(i,q[6]);
    }
    for(int i=0;i<ias_pt->GetNbinsX();i++)
    {
        ias_pt->ProjectionY("",i,i+1)->GetQuantiles(7,q,p);
        quantile01Ias_pt->SetBinContent(i,q[0]);
        quantile10Ias_pt->SetBinContent(i,q[1]);
        quantile30Ias_pt->SetBinContent(i,q[2]);
        quantile50Ias_pt->SetBinContent(i,q[3]);
        quantile70Ias_pt->SetBinContent(i,q[4]);
        quantile90Ias_pt->SetBinContent(i,q[5]);
        quantile99Ias_pt->SetBinContent(i,q[6]);
    }
    for(int i=0;i<ias_p->GetNbinsX();i++)
    {
        ias_pt->ProjectionY("",i,i+1)->GetQuantiles(7,q,p);
        quantile01Ias_p->SetBinContent(i,q[0]);
        quantile10Ias_p->SetBinContent(i,q[1]);
        quantile30Ias_p->SetBinContent(i,q[2]);
        quantile50Ias_p->SetBinContent(i,q[3]);
        quantile70Ias_p->SetBinContent(i,q[4]);
        quantile90Ias_p->SetBinContent(i,q[5]);
        quantile99Ias_p->SetBinContent(i,q[6]);
    }
    
}

void region::fillMassFrom1DTemplates()
{
    TH1F* p = (TH1F*) ih_p->ProjectionX();
    TH1F* ih = (TH1F*) ih_p->ProjectionY();
    scale(p); scale(ih);
    for(int i=1;i<p->GetNbinsX()+1;i++)
    {
        for(int j=1;j<ih->GetNbinsX()+1;j++)
        {
            float mom = p->GetBinCenter(i);
            float dedx = ih->GetBinCenter(j);
            float prob = p->GetBinContent(i) * ih->GetBinContent(j);
            massFrom1DTemplates->Fill(GetMass(mom,dedx,K,C),prob*p->Integral());
        }
    }

    massFrom1DTemplates->Sumw2();
    delete p;
    delete ih;
}

// in order to compute properly the uncertainties we use the methods SetBinContent SetBinError instead of Fill
// as several couples of bins in (p,ih) can provide the same mass estimate we need to properly sum the entries and errors
// for a couple of bins in (p,ih) where the bin content were (N_p,N_ih) the associated quantities should be 
// content: (N_p * N_ih) / N_total, where N_total represents the total number of events in the region (integral of p, ih & mass distributions)
// error: content * sqrt( 1 / N_p + 1 / N_ih ) where we assume Poisson uncertainties in both distributions (independent distributions) and we neglect the uncertainty on N_total
// While combining the input for several couples leading to the same mass: 
// contents are added 
// errors: the sqrt of the squared uncertainties are added
// 

void region::fillMassFrom1DTemplatesEtaBinning() 
{
    errMass = new TH1F(("errMass"+suffix_).c_str(),";Mass error",200,0,2000);
    TH1F* eta = (TH1F*) ih_eta->ProjectionX();
    //ih_p_eta->GetYaxis()->SetRange(ih_p_eta->GetYaxis()->FindBin(0.),ih_p_eta->GetYaxis()->FindBin(200.));
    for(int i=1;i<eta->GetNbinsX();i++)
    {
        TH1F* p = (TH1F*) eta_p->ProjectionX("proj_p",i,i+1);
        if(VectOfBins_P_.size()>1) p = (TH1F*)p->Rebin(VectOfBins_P_.size()-1,"",VectOfBins_P_.data());
        /*TH1F* ih = (TH1F*)((TH2F*)ih_p_eta->Project3D("zx"))->ProjectionY("proj_ih_pcut",i,i);
        for(int x=0;x<ih->GetNbinsX()+1;x++)
        {
            ih_used->Fill(ih->GetBinCenter(x),ih->GetBinContent(x));
        }*/
        TH1F* ih = (TH1F*) ih_eta->ProjectionY("proj_ih",i,i+1);
        scale(p); //
        
        //scale(ih);
        for(int j=1;j<p->GetNbinsX();j++)
        {
            for(int k=1;k<ih->GetNbinsX();k++)
            {
                if(p->GetBinContent(j)<=0) continue;
                if(ih->GetBinContent(k)<=0) continue;
                //p=histoInterpolated(p);
                float mom = p->GetBinCenter(j);
                //if(mom>600) continue;
                //float mom = vect[j];
                //std::cout << "vetc: " << vect[j] << " center: " << p->GetBinCenter(j) << std::endl;
                float dedx = ih->GetBinCenter(k);

                //float prob = p->GetBinContent(j) * ih->GetBinContent(k);
                //std::cout << "interpolate: " << p->Interpolate(mom) << " invExpo: " << invExpoInterpolate(p,mom,vect) << std::endl;
                //float probP = mom<=1000?p->Interpolate(mom):invExpoInterpolate(p,mom,vect);
                
                //float probP = p->Interpolate(mom);
                //float probP = p->GetBinContent(j);
                float prob = p->GetBinContent(j) * ih->GetBinContent(k);
                float weight = prob*p->Integral();
                float err_weight = weight*sqrt((1./(ih->GetBinContent(k)))+(1./(p->GetBinContent(j)*ih->Integral())));
                float mass = GetMass(mom,dedx,K,C);
                int bin_mass = massFrom1DTemplatesEtaBinning->FindBin(mass);
                float mass_err = mass*GetMassErr(mom,p->GetBinWidth(j)/sqrt(12),dedx,ih->GetBinWidth(k)/sqrt(12),mass,K,C); 
                if(prob>=0)
                {
                    // first version : wrong
                    //massFrom1DTemplatesEtaBinning->Fill(GetMass(mom,dedx,K,C),prob*p->Integral());
                    massFrom1DTemplatesEtaBinning->SetBinContent(bin_mass,massFrom1DTemplatesEtaBinning->GetBinContent(bin_mass)+weight);
                    massFrom1DTemplatesEtaBinning->SetBinError(bin_mass,sqrt(pow(massFrom1DTemplatesEtaBinning->GetBinError(bin_mass),2)+pow(err_weight,2)));
                    errMass->Fill(mass_err);
                    Mass_errMass->Fill(mass,mass_err);
                    if(mass>800) mapM800->Fill(mom,dedx,prob);
                }
            }
        }
        delete p;
        delete ih;
    }
    massFrom1DTemplatesEtaBinning->Sumw2();
    //scale(massFrom1DTemplatesEtaBinning);
}

void region::cross1D()
{
    crossHistosEtaBinning(cross1Dtemplates,eta_p,ih_eta);
}

void region::plotMass()
{
    c->Divide(1,2);
    c->cd(1);
    scale(mass); scale(massFrom1DTemplates); scale(massFrom1DTemplatesEtaBinning);
    mass->Draw();
    //massFrom1DTemplates->Draw("same");
    //massFrom1DTemplates->SetLineColor(2);
    massFrom1DTemplatesEtaBinning->Draw("same");
    massFrom1DTemplatesEtaBinning->SetLineColor(2);
    c->SetLogy();
    c->cd(2);
    //TH1F* rat = ratioIntegral(massFrom1DTemplates,mass);
    TH1F* rat2 = ratioIntegral(massFrom1DTemplatesEtaBinning,mass);
    //TH1F* rat3 = ratioHist(massFrom1DTemplatesEtaBinning,mass);
    //rat->Draw();
    rat2->Draw();
    //rat2->Draw("same");
    //rat3->Draw("same");

    //printf(("Kolmogorov test: %f Chi2: %f dof: %f "+suffix_+"\n").c_str(),mass->KolmogorovTest(massFrom1DTemplates),chi2test(mass,massFrom1DTemplates),mass->GetNbinsX());
}

void region::write()
{
    //plotMass();
    //cross1D();
    ih_p_eta->Write();
    ih_pt->Write();
    ias_pt->Write();
    ih_ias->Write();
    ih_nhits->Write();
    ias_nhits->Write();
    eta_pt->Write();
    eta_p->Write();
    nhits_pt->Write();
    eta_nhits->Write();
    ih_eta->Write();
    ih_p->Write();
    cross1Dtemplates->Write();
    ias_p->Write();
    stdDevIh_p->Write();
    stdDevIas_pt->Write();
    stdDevIas_p->Write();
    stdDevIas_p_y->Write();
    quantile01Ih_p->Write();
    quantile10Ih_p->Write();
    quantile30Ih_p->Write();
    quantile50Ih_p->Write();
    quantile70Ih_p->Write();
    quantile90Ih_p->Write();
    quantile99Ih_p->Write();
    quantile01Ias_pt->Write();
    quantile10Ias_pt->Write();
    quantile30Ias_pt->Write();
    quantile50Ias_pt->Write();
    quantile70Ias_pt->Write();
    quantile90Ias_pt->Write();
    quantile99Ias_pt->Write();
    quantile01Ias_p->Write();
    quantile10Ias_p->Write();
    quantile30Ias_p->Write();
    quantile50Ias_p->Write();
    quantile70Ias_p->Write();
    quantile90Ias_p->Write();
    quantile99Ias_p->Write();

    mass->Write();
    massFrom1DTemplates->Write();
    massFrom1DTemplatesEtaBinning->Write();
    errMass->Write();
    Mass_errMass->Write();

    ih_used->Write();
    mapM800->Write();

    c->Write();

    ih_pt->ProjectionY()->Write();

    ih_p->ProfileX()->Write();
    ias_p->ProfileX()->Write();
    ih_pt->ProfileX()->Write();
    ias_pt->ProfileX()->Write();

    
}

void region::rebinQuantiles(int rebin)
{
    quantile01Ias_p->RebinX(rebin);
    quantile10Ias_p->RebinX(rebin);
    quantile30Ias_p->RebinX(rebin);
    quantile50Ias_p->RebinX(rebin);
    quantile70Ias_p->RebinX(rebin);
    quantile90Ias_p->RebinX(rebin);
    quantile99Ias_p->RebinX(rebin);
}
class HscpCandidates {
public :        
    
    TH1F* hA;
    TH1F* hB;
    TH1F* hC;
    TH1F* hD;

    TH2F* regD_ih;
    TH2F* regD_p; 
    TH2F* regD_mass;

    float iascut_;
    float ptcut_;
    float ihcut_;
    float pcut_;
    int etabins_;
    int ihbins_;
    int pbins_;
    int massbins_;
    bool invIso_;
    bool invMET_;
    float etacut_;

    std::string outfilename_;
   
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          Trig;
   UInt_t          Run;
   ULong64_t       Event;
   UInt_t          Lumi;
   UInt_t          PileUp;
   UInt_t          nofVtx;
   UInt_t          Hscp;
   Float_t         Weight;
   Bool_t          HLT_Mu50;
   Bool_t          HLT_PFMET120_PFMHT120_IDTight;
   Bool_t          HLT_PFHT500_PFMET100_PFMHT100_IDTight;
   Bool_t          HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
   Bool_t          HLT_MET105_IsoTrk50;
   Float_t         CaloMET;
   Float_t         RecoPFMET;
   Float_t         RecoPFMHT;
   Float_t         HLTPFMET;
   Float_t         HLTPFMHT;
   vector<bool>    *passCutPt55;
   vector<bool>    *passPreselection_noIsolation_noIh;
   vector<bool>    *passPreselection;
   vector<bool>    *passSelection;
   vector<float>   *Charge;
   vector<float>   *Pt;
   vector<float>   *PtErr;
   vector<float>   *Ias;
   vector<float>   *Ih;
   vector<float>   *Ick;
   vector<float>   *Fmip;
   vector<float>   *ProbQ;
   vector<float>   *TOF;
   vector<float>   *TOFErr;
   vector<unsigned int> *TOF_ndof;
   vector<float>   *DTTOF;
   vector<float>   *DTTOFErr;
   vector<unsigned int> *DTTOF_ndof;
   vector<float>   *CSCTOF;
   vector<float>   *CSCTOFErr;
   vector<unsigned int> *CSCTOF_ndof;
   vector<float>   *Mass;
   vector<float>   *MassErr;
   vector<float>   *dZ;
   vector<float>   *dXY;
   vector<float>   *dR;
   vector<float>   *eta;
   vector<float>   *phi;
   vector<unsigned int> *NOH;
   vector<unsigned int> *NOPH;
   vector<float>   *FOVH;
   vector<unsigned int> *NOMH;
   vector<float>   *FOVHD;
   vector<unsigned int> *NOM;
   vector<float>   *iso_TK;
   vector<float>   *iso_ECAL;
   vector<float>   *iso_HCAL;
   vector<float>   *MuonPFIsolationR03_sumChargedHadronPt;
   vector<float>   *MuonPFIsolationR03_sumNeutralHadronPt;
   vector<float>   *MuonPFIsolationR03_sumPhotonPt;
   vector<float>   *MuonPFIsolationR03_sumPUPt;
   vector<float>   *Ih_noL1;
   vector<float>   *Ih_15drop;
   vector<float>   *Ih_StripOnly;
   vector<float>   *Ih_StripOnly_15drop;
   vector<float>   *Ih_SaturationCorrectionFromFits;
   vector<vector<float> > *clust_charge;
   vector<vector<float> > *clust_pathlength;
   vector<vector<bool> > *clust_ClusterCleaning;
   vector<vector<unsigned int> > *clust_nstrip;
   vector<vector<bool> > *clust_sat254;
   vector<vector<bool> > *clust_sat255;
   vector<vector<unsigned int> > *clust_detid;
   vector<vector<bool> > *clust_isStrip;
   vector<vector<bool> > *clust_isPixel;
   vector<float>   *GenId;
   vector<float>   *GenCharge;
   vector<float>   *GenMass;
   vector<float>   *GenPt;
   vector<float>   *GenEta;
   vector<float>   *GenPhi;

   // List of branches
   TBranch        *b_Trig;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_Lumi;   //!
   TBranch        *b_PileUp;   //!
   TBranch        *b_nofVtx;   //!
   TBranch        *b_Hscp;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_HLT_Mu50;   //!
   TBranch        *b_HLT_PFMET120_PFMHT120_IDTight;   //!
   TBranch        *b_HLT_PFHT500_PFMET100_PFMHT100_IDTight;   //!
   TBranch        *b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;   //!
   TBranch        *b_HLT_MET105_IsoTrk50;   //!
   TBranch        *b_CaloMET;   //!
   TBranch        *b_RecoPFMET;   //!
   TBranch        *b_RecoPFMHT;   //!
   TBranch        *b_HLTPFMET;   //!
   TBranch        *b_HLTPFMHT;   //!
   TBranch        *b_passCutPt55;   //!
   TBranch        *b_passPreselection_noIsolation_noIh;   //!
   TBranch        *b_passPreselection;   //!
   TBranch        *b_passSelection;   //!
   TBranch        *b_Charge;   //!
   TBranch        *b_Pt;   //!
   TBranch        *b_PtErr;   //!
   TBranch        *b_Ias;   //!
   TBranch        *b_Ih;   //!
   TBranch        *b_Ick;   //!
   TBranch        *b_Fmip;   //!
   TBranch        *b_ProbQ;   //!
   TBranch        *b_TOF;   //!
   TBranch        *b_TOFErr;   //!
   TBranch        *b_TOF_ndof;   //!
   TBranch        *b_DTTOF;   //!
   TBranch        *b_DTTOFErr;   //!
   TBranch        *b_DTTOF_ndof;   //!
   TBranch        *b_CSCTOF;   //!
   TBranch        *b_CSCTOFErr;   //!
   TBranch        *b_CSCTOF_ndof;   //!
   TBranch        *b_Mass;   //!
   TBranch        *b_MassErr;   //!
   TBranch        *b_dZ;   //!
   TBranch        *b_dXY;   //!
   TBranch        *b_dR;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_NOH;   //!
   TBranch        *b_NOPH;   //!
   TBranch        *b_FOVH;   //!
   TBranch        *b_NOMH;   //!
   TBranch        *b_FOVHD;   //!
   TBranch        *b_NOM;   //!
   TBranch        *b_iso_TK;   //!
   TBranch        *b_iso_ECAL;   //!
   TBranch        *b_iso_HCAL;   //!
   TBranch        *b_MuonPFIsolationR03_sumChargedHadronPt;   //!
   TBranch        *b_MuonPFIsolationR03_sumNeutralHadronPt;   //!
   TBranch        *b_MuonPFIsolationR03_sumPhotonPt;   //!
   TBranch        *b_MuonPFIsolationR03_sumPUPt;   //!
   TBranch        *b_Ih_noL1;   //!
   TBranch        *b_Ih_15drop;   //!
   TBranch        *b_Ih_StripOnly;   //!
   TBranch        *b_Ih_StripOnly_15drop;   //!
   TBranch        *b_Ih_SaturationCorrectionFromFits;   //!
   TBranch        *b_clust_charge;   //!
   TBranch        *b_clust_pathlength;   //!
   TBranch        *b_clust_ClusterCleaning;   //!
   TBranch        *b_clust_nstrip;   //!
   TBranch        *b_clust_sat254;   //!
   TBranch        *b_clust_sat255;   //!
   TBranch        *b_clust_detid;   //!
   TBranch        *b_clust_isStrip;   //!
   TBranch        *b_clust_isPixel;   //!
   TBranch        *b_GenId;   //!
   TBranch        *b_GenCharge;   //!
   TBranch        *b_GenMass;   //!
   TBranch        *b_GenPt;   //!
   TBranch        *b_GenEta;   //!
   TBranch        *b_GenPhi;   //!

   HscpCandidates(TTree *tree=0);
   //HscpCandidates(TTree *tree=0,float iascut,float ptcut,float ihcut,float pcut,float etacut,int etabins,int ihbins,int pbins,int massbins,bool invIso,bool invMET);
   virtual ~HscpCandidates();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef HscpCandidates_cxx



HscpCandidates::HscpCandidates(TTree *tree) : fChain(0) 
{
    ifstream ifile;
    ifile.open("./configFile.txt");
    std::string line;
    std::string filename;
    float iascut,ptcut,ihcut,pcut,etacut;
    int etabins,ihbins,pbins,massbins;
    bool invIso, invMET;

    while(std::getline(ifile,line))
    {
        if(std::strncmp(line.c_str(),"#",1)==0) continue;
        std::cout << line << std::endl;
        std::stringstream ss(line);
        ss >> filename >> iascut >> ptcut >> ihcut >> pcut >> etabins >> ihbins >> pbins >> massbins >> invIso >> invMET >> etacut;
    }
    iascut_ = iascut;
    ptcut_ = ptcut;
    ihcut_ = ihcut;
    pcut_ = pcut;
    etabins_ = etabins;
    ihbins_ = ihbins;
    pbins_ = pbins;
    massbins_ = massbins;
    invIso_ = invIso;
    invMET_ = invMET;
    etacut_ = etacut;


    outfilename_ = "outfile_ias"+to_string((int)(1000*iascut_))+"_pt"+to_string((int)ptcut_)+"_ih"+to_string((int)(10*ihcut_))+"_p"+to_string((int)pcut_)+"_etabins"+to_string(etabins_)+"_ihbins"+to_string(ihbins_)+"_pbins"+to_string(pbins_)+"_massbins"+to_string(massbins_)+"_invIso"+to_string(invIso_)+"_invMET"+to_string(invMET_)+"_eta"+to_string((int)(10*etacut_));

    std::cout << outfilename_ << std::endl;

// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("histoMET_TkOnly_UL2017C_v2-4_pt50.root");
      //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("histoSingleMu_TkOnly_UL2017C_v2-4_tree.root");
      //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("histoSingleMu_TkOnly_UL2017C_v2-5_wTree.root");
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename.c_str());
      if (!f || !f->IsOpen()) {
         //f = new TFile("histoMET_TkOnly_UL2017C_v2-4_pt50.root");
         //f = new TFile("histoSingleMu_TkOnly_UL2017C_v2-4_tree.root");
         //f = new TFile("histoSingleMu_TkOnly_UL2017C_v2-5_wTree.root");
         f = new TFile(filename.c_str());
      }
      //TDirectory * dir = (TDirectory*)f->Get("histoMET_TkOnly_UL2017C_v2-4_pt50.root:/analyzer/Data_2017_UL");
      //TDirectory * dir = (TDirectory*)f->Get("histoSingleMu_TkOnly_UL2017C_v2-4_tree.root:/analyzer/Data_2017_UL");
      //TDirectory * dir = (TDirectory*)f->Get("histoSingleMu_TkOnly_UL2017C_v2-5_wTree.root:/analyzer/Data_2017_UL");
      TDirectory * dir = (TDirectory*)f->Get((filename+":/analyzer/Data_2017_UL").c_str());
      dir->GetObject("HscpCandidates",tree);


      regD_ih = (TH2F*) f->Get((filename+":/analyzer/Data_2017_UL/RegionD_I").c_str());
      regD_p = (TH2F*) f->Get((filename+":/analyzer/Data_2017_UL/RegionD_P").c_str());
      regD_mass = (TH2F*) f->Get((filename+":/analyzer/Data_2017_UL/Mass").c_str());
      
      //regD_ih = (TH2F*) f->Get("histoSingleMu_TkOnly_UL2017C_v2-5_wTree.root:/analyzer/Data_2017_UL/RegionD_I");
      //regD_ih = (TH2F*) f->Get("histoMET_TkOnly_UL2017C_v2-4_pt50.root:/analyzer/Data_2017_UL/RegionD_I");
      //regD_p = (TH2F*) f->Get("histoSingleMu_TkOnly_UL2017C_v2-5_wTree.root:/analyzer/Data_2017_UL/RegionD_P");
      //regD_p = (TH2F*) f->Get("histoMET_TkOnly_UL2017C_v2-4_pt50.root:/analyzer/Data_2017_UL/RegionD_P");
      //regD_mass = (TH2F*) f->Get("histoSingleMu_TkOnly_UL2017C_v2-5_wTree.root:/analyzer/Data_2017_UL/Mass");
      //regD_mass = (TH2F*) f->Get("histoMET_TkOnly_UL2017C_v2-4_pt50.root:/analyzer/Data_2017_UL/Mass");

      //dir->GetObject("RegionD_I",regD_ih);
      //dir->GetObject("RegionD_P",regD_p);
      //dir->GetObject("Mass",regD_mass);

   }
   Init(tree);
   Loop();
}

//HscpCandidates(TTree *tree=0,float iascut,float ptcut,float ihcut,float pcut,float etacut,int etabins,intihbins,int pbins,int massbins,bool invIso,bool invMET);

/*HscpCandidates::HscpCandidates(TTree *tree,float iascut,float ptcut,float ihcut,float pcut,float etacut,int etabins,int ihbins,int pbins,int massbins,bool invIso,bool invMET) : fChain(0) 
{
    iascut_ = iascut;
    ptcut_ = ptcut;
    ihcut_ = ihcut;
    pcut_ = pcut;
    etabins_ = etabins;
    ihbins_ = ihbins;
    pbins_ = pbins;
    massbins_ = massbins;
    invIso_ = invIso;
    invMET_ = invMET;
    etacut_ = etacut;

    outfilename_ = "outfile_ias"+to_string((int)(1000*iascut_))+"_pt"+to_string((int)ptcut_)+"_ih"+to_string((int)(10*ihcut_))+"_p"+to_string((int)pcut_)+"_eta"+to_string((int)(10*etacut_))+"_etabins"+to_string(etabins_)+"_ihbins"+to_string(ihbins_)+"_pbins"+to_string(pbins_)+"_massbins"+to_string(massbins_)+"_invIso"+to_string(invIso_)+"_invMET"+to_string(invMET_);

    std::cout << outfilename_ << std::endl;

// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename.c_str());
      if (!f || !f->IsOpen()) {
         f = new TFile(filename.c_str());
      }
      TDirectory * dir = (TDirectory*)f->Get((filename+":/analyzer/Data_2017_UL").c_str());
      dir->GetObject("HscpCandidates",tree);


      regD_ih = (TH2F*) f->Get((filename+":/analyzer/Data_2017_UL/RegionD_I").c_str());
      regD_p = (TH2F*) f->Get((filename+":/analyzer/Data_2017_UL/RegionD_P").c_str());
      regD_mass = (TH2F*) f->Get((filename+":/analyzer/Data_2017_UL/Mass").c_str());
      
   }
   Init(tree);
   Loop();


}*/

HscpCandidates::~HscpCandidates()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t HscpCandidates::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HscpCandidates::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void HscpCandidates::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   passCutPt55 = 0;
   passPreselection_noIsolation_noIh = 0;
   passPreselection = 0;
   passSelection = 0;
   Charge = 0;
   Pt = 0;
   PtErr = 0;
   Ias = 0;
   Ih = 0;
   Ick = 0;
   Fmip = 0;
   ProbQ = 0;
   TOF = 0;
   TOFErr = 0;
   TOF_ndof = 0;
   DTTOF = 0;
   DTTOFErr = 0;
   DTTOF_ndof = 0;
   CSCTOF = 0;
   CSCTOFErr = 0;
   CSCTOF_ndof = 0;
   Mass = 0;
   MassErr = 0;
   dZ = 0;
   dXY = 0;
   dR = 0;
   eta = 0;
   phi = 0;
   NOH = 0;
   NOPH = 0;
   FOVH = 0;
   NOMH = 0;
   FOVHD = 0;
   NOM = 0;
   iso_TK = 0;
   iso_ECAL = 0;
   iso_HCAL = 0;
   MuonPFIsolationR03_sumChargedHadronPt = 0;
   MuonPFIsolationR03_sumNeutralHadronPt = 0;
   MuonPFIsolationR03_sumPhotonPt = 0;
   MuonPFIsolationR03_sumPUPt = 0;
   Ih_noL1 = 0;
   Ih_15drop = 0;
   Ih_StripOnly = 0;
   Ih_StripOnly_15drop = 0;
   Ih_SaturationCorrectionFromFits = 0;
   clust_charge = 0;
   clust_pathlength = 0;
   clust_ClusterCleaning = 0;
   clust_nstrip = 0;
   clust_sat254 = 0;
   clust_sat255 = 0;
   clust_detid = 0;
   clust_isStrip = 0;
   clust_isPixel = 0;
   GenId = 0;
   GenCharge = 0;
   GenMass = 0;
   GenPt = 0;
   GenEta = 0;
   GenPhi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Trig", &Trig, &b_Trig);
   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("Lumi", &Lumi, &b_Lumi);
   fChain->SetBranchAddress("PileUp", &PileUp, &b_PileUp);
   fChain->SetBranchAddress("nofVtx", &nofVtx, &b_nofVtx);
   fChain->SetBranchAddress("Hscp", &Hscp, &b_Hscp);
/*   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("HLT_Mu50", &HLT_Mu50, &b_HLT_Mu50);
   fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight", &HLT_PFMET120_PFMHT120_IDTight, &b_HLT_PFMET120_PFMHT120_IDTight);
   fChain->SetBranchAddress("HLT_PFHT500_PFMET100_PFMHT100_IDTight", &HLT_PFHT500_PFMET100_PFMHT100_IDTight, &b_HLT_PFHT500_PFMET100_PFMHT100_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_MET105_IsoTrk50", &HLT_MET105_IsoTrk50, &b_HLT_MET105_IsoTrk50);
   fChain->SetBranchAddress("CaloMET", &CaloMET, &b_CaloMET);
   fChain->SetBranchAddress("RecoPFMET", &RecoPFMET, &b_RecoPFMET);
   fChain->SetBranchAddress("RecoPFMHT", &RecoPFMHT, &b_RecoPFMHT);
   fChain->SetBranchAddress("HLTPFMET", &HLTPFMET, &b_HLTPFMET);
   fChain->SetBranchAddress("HLTPFMHT", &HLTPFMHT, &b_HLTPFMHT);*/
   fChain->SetBranchAddress("passCutPt55", &passCutPt55, &b_passCutPt55);
   fChain->SetBranchAddress("passPreselection_noIsolation_noIh", &passPreselection_noIsolation_noIh, &b_passPreselection_noIsolation_noIh);
   fChain->SetBranchAddress("passPreselection", &passPreselection, &b_passPreselection);
   fChain->SetBranchAddress("passSelection", &passSelection, &b_passSelection);
   fChain->SetBranchAddress("Charge", &Charge, &b_Charge);
   fChain->SetBranchAddress("Pt", &Pt, &b_Pt);
   fChain->SetBranchAddress("PtErr", &PtErr, &b_PtErr);
   fChain->SetBranchAddress("Ias", &Ias, &b_Ias);
   fChain->SetBranchAddress("Ih", &Ih, &b_Ih);
   fChain->SetBranchAddress("Ick", &Ick, &b_Ick);
   fChain->SetBranchAddress("Fmip", &Fmip, &b_Fmip);
   fChain->SetBranchAddress("ProbQ", &ProbQ, &b_ProbQ);
   fChain->SetBranchAddress("TOF", &TOF, &b_TOF);
/*   fChain->SetBranchAddress("TOFErr", &TOFErr, &b_TOFErr);
   fChain->SetBranchAddress("TOF_ndof", &TOF_ndof, &b_TOF_ndof);
   fChain->SetBranchAddress("DTTOF", &DTTOF, &b_DTTOF);
   fChain->SetBranchAddress("DTTOFErr", &DTTOFErr, &b_DTTOFErr);
   fChain->SetBranchAddress("DTTOF_ndof", &DTTOF_ndof, &b_DTTOF_ndof);
   fChain->SetBranchAddress("CSCTOF", &CSCTOF, &b_CSCTOF);
   fChain->SetBranchAddress("CSCTOFErr", &CSCTOFErr, &b_CSCTOFErr);
   fChain->SetBranchAddress("CSCTOF_ndof", &CSCTOF_ndof, &b_CSCTOF_ndof);*/
   fChain->SetBranchAddress("Mass", &Mass, &b_Mass);
   fChain->SetBranchAddress("MassErr", &MassErr, &b_MassErr);
   fChain->SetBranchAddress("dZ", &dZ, &b_dZ);
   fChain->SetBranchAddress("dXY", &dXY, &b_dXY);
   fChain->SetBranchAddress("dR", &dR, &b_dR);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("NOH", &NOH, &b_NOH);
   fChain->SetBranchAddress("NOPH", &NOPH, &b_NOPH);
   fChain->SetBranchAddress("FOVH", &FOVH, &b_FOVH);
   fChain->SetBranchAddress("NOMH", &NOMH, &b_NOMH);
   fChain->SetBranchAddress("FOVHD", &FOVHD, &b_FOVHD);
   fChain->SetBranchAddress("NOM", &NOM, &b_NOM);
   fChain->SetBranchAddress("iso_TK", &iso_TK, &b_iso_TK);
   fChain->SetBranchAddress("iso_ECAL", &iso_ECAL, &b_iso_ECAL);
   fChain->SetBranchAddress("iso_HCAL", &iso_HCAL, &b_iso_HCAL);
   fChain->SetBranchAddress("MuonPFIsolationR03_sumChargedHadronPt", &MuonPFIsolationR03_sumChargedHadronPt, &b_MuonPFIsolationR03_sumChargedHadronPt);
   fChain->SetBranchAddress("MuonPFIsolationR03_sumNeutralHadronPt", &MuonPFIsolationR03_sumNeutralHadronPt, &b_MuonPFIsolationR03_sumNeutralHadronPt);
   fChain->SetBranchAddress("MuonPFIsolationR03_sumPhotonPt", &MuonPFIsolationR03_sumPhotonPt, &b_MuonPFIsolationR03_sumPhotonPt);
   fChain->SetBranchAddress("MuonPFIsolationR03_sumPUPt", &MuonPFIsolationR03_sumPUPt, &b_MuonPFIsolationR03_sumPUPt);
/*   fChain->SetBranchAddress("Ih_noL1", &Ih_noL1, &b_Ih_noL1);
   fChain->SetBranchAddress("Ih_15drop", &Ih_15drop, &b_Ih_15drop);
   fChain->SetBranchAddress("Ih_StripOnly", &Ih_StripOnly, &b_Ih_StripOnly);
   fChain->SetBranchAddress("Ih_StripOnly_15drop", &Ih_StripOnly_15drop, &b_Ih_StripOnly_15drop);
   fChain->SetBranchAddress("Ih_SaturationCorrectionFromFits", &Ih_SaturationCorrectionFromFits, &b_Ih_SaturationCorrectionFromFits);
   fChain->SetBranchAddress("clust_charge", &clust_charge, &b_clust_charge);
   fChain->SetBranchAddress("clust_pathlength", &clust_pathlength, &b_clust_pathlength);
   fChain->SetBranchAddress("clust_ClusterCleaning", &clust_ClusterCleaning, &b_clust_ClusterCleaning);
   fChain->SetBranchAddress("clust_nstrip", &clust_nstrip, &b_clust_nstrip);
   fChain->SetBranchAddress("clust_sat254", &clust_sat254, &b_clust_sat254);
   fChain->SetBranchAddress("clust_sat255", &clust_sat255, &b_clust_sat255);
   fChain->SetBranchAddress("clust_detid", &clust_detid, &b_clust_detid);
   fChain->SetBranchAddress("clust_isStrip", &clust_isStrip, &b_clust_isStrip);
   fChain->SetBranchAddress("clust_isPixel", &clust_isPixel, &b_clust_isPixel);
   fChain->SetBranchAddress("GenId", &GenId, &b_GenId);
   fChain->SetBranchAddress("GenCharge", &GenCharge, &b_GenCharge);
   fChain->SetBranchAddress("GenMass", &GenMass, &b_GenMass);
   fChain->SetBranchAddress("GenPt", &GenPt, &b_GenPt);
   fChain->SetBranchAddress("GenEta", &GenEta, &b_GenEta);
   fChain->SetBranchAddress("GenPhi", &GenPhi, &b_GenPhi);*/
   Notify();
}

Bool_t HscpCandidates::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HscpCandidates::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HscpCandidates::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef HscpCandidates_cxx

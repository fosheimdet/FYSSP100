/*
 * Project:        Exercises 11.1-11.3
 * File:           Walkthrough_skeleton.C
 * Author:         Ivo van Vulpen, Aart Heijboer
 * Version (date): 1.0 (23.06.2013)
 *
 * Copyright (C) 2013, Ivo van Vulpen, Aart Heijboer
 * All rights reserved.
 *
 * Description:
 * A code skeleton for the searches part.
 *
 * This code is distributed with the solution manual to the book
 *
 * Data Analysis in High Energy Physics: A Practical Guide to Statistical Methods,
 * Wiley-VCH (2013),
 * O. Behnke, K. Kroeninger, G. Schott, Th. Schoerner-Sadenius (editors)
 */

#include"TLatex.h"

#include "TBox.h"
#include "TArrow.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "Math/QuantFuncMathCore.h"
#include "Math/ProbFuncMathCore.h" // for ROOT::Math::gaussian_cdf

#include <iostream>
using namespace std;

//-------------------------------------------------------------------------------------------------------------------------
//-- full functions
TH1D * GetMassDistribution(int Itype = 1, double scalefactor = 1.00);
void MassPlot(int Irebin = 20);

//-- skeleton functions
void SideBandFit(int Irebin = 10);
double Get_TestStatistic(TH1D *h_mass_dataset, TH1D *h_template_bgr, TH1D *h_template_sig);
TH1D * GenerateToyDataSet(TH1D *h_mass_template);

//-- some fusefull functions
double IntegratePoissonFromRight(double mu, int N_obs);
double IntegrateFromRight(TH1D * h_X_bgr, double X_value);
vector<double> Get_Quantiles( TH1D* hist );
void AddText( Double_t txt_x = 0.50, Double_t txt_y = 0.50, const char * txt = "dummy", Double_t txt_size = 0.045,
	      Double_t txt_angle = 0., const char * Alignment = "left", Int_t UseNormalizedSize = 1, Int_t txt_color =1 );
//-------------------------------------------------------------------------------------------------------------------------



//========================================
// S O M E   F I N A L   F U N C T I O N S
//========================================


//========================================================
TH1D * GetMassDistribution(int Itype, double scalefactor){
//========================================================
 //----------------------------------------------------------
 // Goal: return the histogram of the 4-lepton invariant mass
 //       for given type with an optional scale factor
 //
 //       Itype 1 = ZZ SM background
 //             2 = data
 //           125 = Higgs 125
 //           200 = Higgs 200
 //
 //      scalefactor: histograms will be scaled with this number
 //
 //  Note: Histograms have ~200 MeV bins, so need to rebin
 //---------------------------------------------------------

  //-- [1] Get histogram from the file
  TH1D *h_mass = 0;
  TDirectory* dir = gDirectory;
  TFile *file = new TFile("Histograms_fake.root", "READ");
  dir->cd();

  //-- Higgs 125
  if(Itype == 125){
    h_mass  = (TH1D*) file->Get("h_m4l_Higgs125_fake")->Clone("h_mass");
  }
  //-- Higgs 200
  if(Itype == 200){
    h_mass  = (TH1D*) file->Get("h_m4l_Higgs200_fake")->Clone("h_mass");
  }
  //-- ZZ SM background
  if(Itype == 1){
    h_mass  = (TH1D*) file->Get("h_m4l_ZZ_fake")->Clone("h_mass");
  }
  //-- data
  if(Itype == 2){
    h_mass  = (TH1D*) file->Get("h_m4l_data_fake")->Clone("h_mass");
  }

  //-- [2] scale histograms
  int Nbins = h_mass->GetNbinsX();
  for (int i_bin = 1; i_bin < Nbins; i_bin++){
    double mu_bin = h_mass->GetBinContent(i_bin);
    h_mass -> SetBinContent( i_bin, scalefactor * mu_bin);
  }


  file->Close();
  //-- [3] return histogram
  return h_mass;

  //===========================
} // end GetMassDistribution()
  //===========================




//========================
void MassPlot(int Irebin){
//========================
  // ------------------------------------------
  // Goal: produce SM+Higgs+data plot
  //       Note: rebinning is only for plotting
  // Rebinning joins bins together. E.g. h_mass->Rebin(5) will merge 5 bins to 1. I.e., the "granularity" of the histogram is reduced.
  // ------------------------------------------

  //------------------------------------
  //-- Standard stuff and prepare canvas
  //------------------------------------
  gROOT->Clear();
  gROOT->Delete();

  //-- Prepare canvas and plot histograms
  TCanvas * canvas1 = new TCanvas("canvas1","Standard Canvas",600,400);
//  canvas1->SetLeftMargin(0.125);
  canvas1->SetLeftMargin(0.125);
  canvas1->SetBottomMargin(0.125);
	canvas1->SetGrid();
  canvas1->cd();

  //------------------------------------------------------------------
  //-- [1] Prepare histograms
  //--     o Get histograms from the files (signal, background and data)
  //--     o Make cumulative histograms (for signal and background)
  //------------------------------------------------------------------

  //-- Get histograms from the files (higgs, zz and data)
  TH1D *h_sig, *h_bgr, *h_data;
  h_sig  = GetMassDistribution(125);
//  h_sig  = GetMassDistribution(200);
  h_bgr  = GetMassDistribution(1);
  h_data = GetMassDistribution(2);

  //-----------------------------------
  //-- [2] Plot histograms and make gif
  //--     o rebin histograms
  //--     o prepare cumulative histogram
  //--     o make plot + opsmuk + gif
  //-----------------------------------

  //-- Rebin histograms (only for plotting)
  h_sig->Rebin(Irebin);
  h_bgr->Rebin(Irebin);
  h_data->Rebin(Irebin);

  //-- Prepare cumulative histogram for signal + background
  TH1D *h_sig_plus_bgr = (TH1D* ) h_bgr->Clone("h_sig_plus_bgr");
  h_sig_plus_bgr->Reset();
  for (int i_bin = 1; i_bin < h_bgr->GetNbinsX(); i_bin++){
       h_sig_plus_bgr->SetBinContent( i_bin, h_sig->GetBinContent(i_bin) + h_bgr->GetBinContent(i_bin));
       printf("  REBINNED HISTOGRAM:  bin %d, Ndata = %d\n",i_bin,(int)h_data->GetBinContent(i_bin));
  }

  //-- prepare histograms and plot them on canvas
  double Data_max = h_data->GetBinContent(h_data->GetMaximumBin());
  double Ymax_plot = 1.10* (Data_max + TMath::Sqrt(Data_max));
  h_sig_plus_bgr->SetFillColor(7);
  h_sig_plus_bgr->SetAxisRange(0.,Ymax_plot,"Y");
  h_sig_plus_bgr->SetAxisRange(0.,400.,"X");
  h_bgr->SetFillColor(2);
  h_sig_plus_bgr->Draw("hist");
  h_bgr->Draw("same");
  h_bgr->Draw("axis same");
  h_data->Draw("e same");

  //-- some nice axes and add legend
  AddText( 0.900, 0.035, "4-lepton invariant mass [GeV]",0.060, 0.,"right");                             // X-axis
  AddText( 0.040, 0.900, Form("Number of events / %3.1f GeV",h_bgr->GetBinWidth(1)) ,0.060,90.,"right"); // Y-axis
  TLegend *leg1 = new TLegend(0.61,0.61,0.85,0.8);
  leg1->SetBorderSize(0); leg1->SetFillColor(0);
  TLegendEntry *leg1a = leg1->AddEntry(h_bgr,          " SM(ZZ)", "f");  leg1a->SetTextSize(0.038);
  TLegendEntry *leg1b = leg1->AddEntry(h_sig_plus_bgr, " Higgs" , "f");  leg1b->SetTextSize(0.038);
  leg1->Draw();

  //-- prepare gif
  canvas1->Print(Form("./MassPlot_rebin%d.gif",Irebin));

  return;

   //===============
 } // end MassPlot()
   //===============

//=========================================================================================================================
//=========================================================================================================================
//=========================================================================================================================
//=========================================================================================================================
//=========================================================================================================================
//=========================================================================================================================
//=========================================================================================================================


//===========================================================================================================
double Get_TestStatistic(TH1D *h_mass_dataset, TH1D *h_template_bgr, TH1D *h_template_sig){
//===========================================================================================================
  //printf(" dummy = %d\n",h_mass_dataset->GetNbinsX() + h_template_bgr->GetNbinsX() + h_template_sig->GetNbinsX());

  double test_statistic = 0.;

  //-- do likelihood fit
	int Nbins = h_mass_dataset->GetXaxis()->GetNbins();
	int nObs;
	double s,b;

	double loglik_bgr = 0;
	double loglik_sig_plus_bgr = 0;
  for (int i_bin = 1; i_bin <= Nbins; i_bin++){
		//Get number of observed events, background and signal for the given bin:
		  nObs = h_mass_dataset->GetBinContent(i_bin);
			b = h_template_bgr->GetBinContent(i_bin);
			s = h_template_sig->GetBinContent(i_bin);

      loglik_bgr += TMath::Log(TMath::Poisson(nObs,b));  // likelihood for the mu=0 (b-only) scenario
      loglik_sig_plus_bgr += TMath::Log(TMath::Poisson(nObs, s+b));  // likelihood for the mu=1 (s+b) scenario
  } // end loop over bins

  //-- compute test statistic
  //double test_statistic  =  // find this out yourself
	test_statistic = 2*loglik_bgr - 2*loglik_sig_plus_bgr;

  //-- return test_statistic
  return test_statistic;

} // end Get_TestStatistic()

double observedTestStatistic(){

 TH1D* h_bgr = GetMassDistribution(1);
 TH1D* h_data = GetMassDistribution(2);
 TH1D* h_sig = GetMassDistribution(125);

 double Xobs = Get_TestStatistic(h_data, h_bgr, h_sig);

 return Xobs;

}



//===============================================
TH1D * GenerateToyDataSet(TH1D *h_mass_template){
//===============================================
  //-------------------------------------------------
  // Goal: Generate Toy data set from input histogram
  // How:  dumb way -> draw random Poisson in each bin
  //-------------------------------------------------
  TRandom3 *R = new TRandom3(0);

  //-- Create new histogram for the data-set
  TH1D *h_mass_toydataset = (TH1D*) h_mass_template->Clone("h_mass_toydataset"); h_mass_toydataset->Reset();


  //-- Loop over bins and draw Poisson number of event in each bin
	for (int i = 1; i<= h_mass_toydataset->GetNbinsX(); i++){
		h_mass_toydataset->SetBinContent(i,R->Poisson(h_mass_template->GetBinContent(i))); //For each bin, pick a random variable distributed
		//according to a Poisson with mean equal to the bin content.
	}

  //-- return histogram of toy data-set
  return h_mass_toydataset;

} // end GenerateToyDataSet()


tuple<double,double,double,double> X_dist_hypotheses(int Ntoys, double mu, bool plotAndPrint){

	TH1D* h_bgr = GetMassDistribution(1);
	TH1D* h_data = GetMassDistribution(2); //Our "mass template"
	TH1D* h_sig = GetMassDistribution(125);

	double X;

	TH1D* Xdist_b = new TH1D("Xdist_b", "", 100, -33,22);  //Histogram to hold the values of X for all the background MC simulations
	TH1D* Xdist_bs = new TH1D("Xdist_bs", "", 100, -33, 22); //Histogram to hold the values of X for all the background+signal MC simulations

	TH1D* h_bs = (TH1D*) h_bgr->Clone("h_bs");

	h_bs->Reset();
	for(int i = 0; i<h_bgr->GetNbinsX(); i++){
		h_bs->SetBinContent(i,h_bgr->GetBinContent(i)+mu*h_sig->GetBinContent(i));
	}

	for(int i = 0; i<Ntoys; i++){
	 Xdist_b->Fill(Get_TestStatistic(GenerateToyDataSet(h_bgr), h_bgr, h_sig)); //Under the bgr hypothesis, the bin count should vary according to P(n;b)
	 Xdist_bs->Fill(Get_TestStatistic(GenerateToyDataSet(h_bs), h_bgr, h_sig)); //under the bgr+sig hypothesis, the bin count should vary according to P(n; s+b)
	}
	double X_obs= Get_TestStatistic(h_data, h_bgr, h_sig);
//=====================================================================================
//Make histograms for the 68 and 95% quantiles of Xdist_b
//=====================================================================================
	vector<double> a = Get_Quantiles(Xdist_b); //Quantiles for the background distribution
	vector<double> k = Get_Quantiles(Xdist_bs); //Quantiles for the signal+background distribution
	TH1D* green = (TH1D*) Xdist_b->Clone("green");
	green->Reset();
	TH1D* yellow = (TH1D*) Xdist_b->Clone("yellow");
	yellow->Reset();

  double Xdist_b_count = Xdist_b->Integral();
	double Xdist_bs_count = Xdist_bs->Integral();

	for(int i = 0; i<Xdist_b->GetNbinsX(); i++){
     yellow->SetBinContent(i,0);
		 green->SetBinContent(i,0);
		 double x = Xdist_b->GetBinCenter(i);
		 double Xb = Xdist_b->GetBinContent(i);
		if(a[0]<x && x<a[4]){yellow->SetBinContent(i,Xb);}
		if(a[2]<x && x<a[3]){green->SetBinContent(i,Xb);}
	}
	//cout<< "a[0]:  " << a[0] << "  " << "a[1]:  " << a[1] << "a[2]:  "<<a[2] << "a[3]:  " << a[3] << "a[4]:  " << a[4] << endl;
green->SetFillColorAlpha(kGreen, 0.4);
yellow->SetFillColorAlpha(kYellow,0.4);
//=====================================================================================
//																		Plotting
//=====================================================================================
if(plotAndPrint){


TCanvas* c_x = new TCanvas("Xdist", "", 600, 400);

c_x->cd();
c_x->SetGrid();

Xdist_b->Draw("hist");
Xdist_bs->SetLineColor(kRed);


yellow->Draw("hist same");
green->Draw("hist same");
Xdist_bs->Draw("hist same");

AddText( 0.45, 0.57, "t_{obs}",0.060, 0.,"right");

TArrow* X_obs_arrow = new TArrow(X_obs, Xdist_b->GetMaximum()/2, X_obs, 0.3,0.02, "|>");

X_obs_arrow->SetLineWidth(2);
AddText( 0.900, 0.035, "Test statistic, t(n)",0.060, 0.,"right"); // X-axis
AddText( 0.040, 0.900, "#MC simulations" ,0.060,90.,"right");   // Y-axis

TLegend *leg = new TLegend(0.7, 0.6, 0.88, 0.88);
leg->AddEntry(Xdist_b, "H_0/b distribution", "l");
leg->AddEntry(Xdist_bs, "H_1/b+s distribution", "l");
// leg->AddEntry(Xobs,"Expected p-value", "f");
//   leg->AddEntry(fit, "Fit function", "l");
leg->Draw("same");

X_obs_arrow->Draw("same |>");

} // End bool for plot
//=====================================================================================
//=====================================================================================

double X_b_med = a[2];
double X_bs_med = k[2];

//======================Calculate p-values (i.e. 1-CL_b) ==============================
if(plotAndPrint)cout << "===================p-values and significances=======================" << endl;
double count_obs = Xdist_b->Integral(0,Xdist_b->FindBin(X_obs)); //Count number of entries in all bins from bin zero to the bin corresponding to the specified t-value
double count_b = Xdist_b->Integral(0,Xdist_b->FindBin(X_b_med));
double count_bs = Xdist_b->Integral(0, Xdist_b->FindBin(X_bs_med));

double tot_b = Xdist_b->Integral(); //Total bin-count for background distribution

double p_obs =count_obs/double(tot_b);
double p_b = count_b/tot_b;
double p_bs = count_bs/tot_b;
if(plotAndPrint == true){
	cout << "p_obs:  "<<p_obs << endl;
	cout << "p_b_med: "<<p_b << endl;
	cout << "p_bs_med:  "<<p_bs << endl;
}


double Z_obs = ROOT::Math::gaussian_quantile_c(p_obs,1); //Observed significance
double Z_b = ROOT::Math::gaussian_quantile_c(p_b, 1);
double Z_bs = ROOT::Math::gaussian_quantile_c(p_bs, 1); //Expected significance

if(plotAndPrint == true){
	cout << "Z_obs:  " << Z_obs << endl;
	cout << "Z_b:  " << Z_b << endl;
	cout << "Z_bs:  " << Z_bs << endl;
}



//======================Calculate CL_s+b========================================
if(plotAndPrint)cout << "===================CL_{s+b}'s=======================" << endl;


int Nbins_bs = Xdist_bs->GetXaxis()->GetNbins();

double count_CL_obs = Xdist_bs->Integral(Xdist_bs->FindBin(X_obs), Nbins_bs);
double count_CL_b = Xdist_bs->Integral(Xdist_bs->FindBin(a[2]), Nbins_bs);
double count_CL_bs = Xdist_bs->Integral(Xdist_bs->FindBin(k[2]), Nbins_bs);


double tot_bs = Xdist_bs->Integral();


double CL_obs = count_CL_obs/tot_bs; //CL_s+b t_obs
double CL_b = count_CL_b/tot_bs;     //CL_s+b for median of b dist
double CL_bs = count_CL_bs/tot_bs;    //CL_s+b for median of b+s dist

if(plotAndPrint == true){
	cout << "bin number, median bs: " << Xdist_bs->FindBin(k[2]) << endl;
	cout << "total number of bins, bs: " << Nbins_bs << endl;
	cout << "tot_bs:  " << tot_bs << endl;
	cout << "CL_obs:   " << CL_obs << endl;
	cout << "CL_b:  " << CL_b << endl;
	cout << "CL_bs:  " << CL_bs << endl;
}


tuple<double,double,double, double> medians = make_tuple(X_obs, X_b_med, X_bs_med, CL_obs);
return medians;

}

//====================================================================================================
//===================== Confidence level vs signal cross section scalefactor =========================
//====================================================================================================
void CL_vs_mu(int Ntoys, int Nmus){
// int Nmus = 20;
TH1D* mus = new TH1D("mus", "",Nmus, 0,4);
//TH1D* mus = new TH1D("mus", "",Nmus, 2.40,3.50);
TH1D* p = (TH1D*) mus->Clone("p"); p->Reset();
tuple<double,double,double,double> MCrunRes;

double current_mu; //To be used in for loop
double current_CL;

//Loop over various signal cross section scale factors and calculate CL_sb_obs each time (via X_dist_hypotheses())
for(int i = 0; i<=Nmus; i++){
	current_mu = mus->GetBinCenter(i);
	MCrunRes = X_dist_hypotheses(Ntoys,current_mu, false);
	current_CL = get<3>(MCrunRes);
	mus->SetBinContent(i,current_CL);
	p->SetBinContent(i,0.05);
}

TCanvas* c = new TCanvas("c", "", 600, 400);
c->SetBottomMargin(0.12);
c->cd();
c->SetGrid();

mus->Draw();
p->SetLineColor(kRed);
p->Draw("same");
AddText( 0.900, 0.035, "Scalefactor for signal cross-section",0.060, 0.,"right"); // X-axis
AddText( 0.2, 0.95, "CL_{s+b}[t_{obs}]" ,0.060,0.,"right");   // Y-axis

}







//=========================================================================================================================
//=========================================================================================================================
//=========================================================================================================================
//=========================================================================================================================
//=========================================================================================================================
//=========================================================================================================================
//=========================================================================================================================

   //================================================
   // S O M E   U S E F U L L   F U N C T I O N S
   //================================================





   //====================================================
   double IntegratePoissonFromRight(double mu, int N_obs){
   //====================================================
   // --------------------------------------------------------------------
   // Compute p-value for case zero background uncertainty, i.e.
   //         just integrate Poisson from the right from N_obs to infinity
   // --------------------------------------------------------------------

     double integral = 1.;
     for(int i_obs = 0; i_obs < N_obs; i_obs++){
       integral -= TMath::Poisson(i_obs,mu);
     }

   //  if (integral == 1.){integral = 1. - 1e-9;}
     if (integral == 1.){integral = 0.5;}
     return integral;

   } // end IntegratePoissonFromRight()


   //========================================================
   double IntegrateFromRight(TH1D * h_X_bgr, double X_value){
   //========================================================
   // --------------------------------------------------------------------
   // Compute p-value: integrate number of events from X_value to infinity
   // --------------------------------------------------------------------

     //-- Integrate distributions
     int Nbins = h_X_bgr->GetNbinsX();
     int X_bin = h_X_bgr->FindBin(X_value);

     //-- Compute integral from X-value to infinity
     double pvalue = h_X_bgr->Integral(X_bin,Nbins) / h_X_bgr->Integral();

     return pvalue;

   } // end IntegrateFrom Right()




   //=========================================
   vector<double> Get_Quantiles( TH1D* hist ){
   //=========================================
   // Quantiles returns a vector<double> with 5 entries.
   // Entries 0 and 4 are the values on the histogram x-axis
   // so that 95% of the content lies between these values.
   // Entries 1 and 3 bracket 68% in the same way.
   // Entry 2 is the median of the histogram.

     //-- define quantiles
     double fraction_1sigma = ROOT::Math::gaussian_cdf(-1.,1.,0.); // 15.8655 %
     double fraction_2sigma = ROOT::Math::gaussian_cdf(-2.,1.,0.); //  2.2750 %
     double probs[5] = {fraction_2sigma, fraction_1sigma, 0.50, 1.00-fraction_1sigma, 1.00-fraction_2sigma };

     //-- output of the quantiles
     double Xvalues[5];

     //-- extract quantiles
     hist->GetQuantiles( 5, Xvalues, probs );

     vector<double> Xvalues_output(5);
     for (int i=0; i<5; i++)
       {
         Xvalues_output[i] = Xvalues[i];
       }

     return Xvalues_output;
   } // end Get_Quantiles()




   //=======================================================================================================================
   void AddText( Double_t txt_x, Double_t txt_y, const char * txt, Double_t txt_size,
                 Double_t txt_angle, const char * Alignment, Int_t UseNormalizedSize, Int_t txt_color)
   //=======================================================================================================================

   {
     Int_t txt_align = 12;
     if ( !strcmp(Alignment, "left"))   { txt_align = 12; } // left
     if ( !strcmp(Alignment, "right"))  { txt_align = 32; } // right
     if ( !strcmp(Alignment, "center")) { txt_align = 22; } // center

     TLatex* t1 = new TLatex( txt_x, txt_y, txt);
     if(UseNormalizedSize) {t1->SetNDC(kTRUE);} // <- use NDC coordinate
     t1->SetTextSize(txt_size);
     t1->SetTextAlign(txt_align);
     t1->SetTextAngle(txt_angle);
     t1->SetTextColor(txt_color);
     t1->Draw();

   } // end AddText()

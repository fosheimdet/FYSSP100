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
	// canvas1->SetTicksX();
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
  TLegend *leg1 = new TLegend(0.7,0.7,0.80,0.80);
  leg1->SetBorderSize(0); leg1->SetFillColor(0);
  TLegendEntry *leg1a = leg1->AddEntry(h_bgr,          " SM(ZZ)", "f");  leg1a->SetTextSize(0.04);
  TLegendEntry *leg1b = leg1->AddEntry(h_sig_plus_bgr, " Higgs" , "f");  leg1b->SetTextSize(0.04);
  leg1->Draw();

  //-- prepare gif
  canvas1->Print(Form("./MassPlot_rebin%d.gif",Irebin));

  return;

   //===============
 } // end MassPlot()
   //===============




//===============================================
// S O M E   S K E L E T O N    F U N C T I O N S
//===============================================


//=============================================================
tuple<double,double,double,double> Significance_Optimization(double Lumi_scalefactor, double bgr_scalefactor, bool plot){
//=============================================================

  printf("\n Significance_Optimization()\n\n");

  //------------------------------------------------------------------
  //-- [1] Prepare histograms
  //--     o Get histograms from the files (signal, background and data)
  //--     o scale to correct luminosity
  //------------------------------------------------------------------

  //-- Get histograms from the files (higgs, zz and data)
  TH1D *h_sig, *h_bgr, *h_data;
  printf ("\n  INFO: Mass distribution in the 4 lepton channel\n");
  h_sig  = GetMassDistribution(125, Lumi_scalefactor);
  h_bgr  = GetMassDistribution(1, Lumi_scalefactor);
  h_data = GetMassDistribution(2, Lumi_scalefactor);

  //-------------------------------------------------------
  //-- [2] Compute significance for various mass windows
  //--     o try various options for window (use histogram)
  //--     o compute expected and observed significance
  //-------------------------------------------------------

 // -- Define histogram that defines mass windows to try


  TH1D *h_masswindow          = new TH1D("h_masswindow","",250,0.,25.);           // make a mass window - full width between 0 and 25 GeV
  TH1D *h_masswindow_expected = new TH1D("h_masswindow_expected","",250,0.,25.);  // histogram to hold results for expected
  TH1D *h_masswindow_observed = new TH1D("h_masswindow_observed","",250,0.,25.);  // histogram to hold results for observed

  //---------------------------------
  //-- Loop over various mass windows
  //---------------------------------

  double largest_exp_sig = -10;
  double largest_obs_sig = -10;
  double optimal_masswindow_sig = 0;
  double optimal_masswindow_obs = 0;


  double binCount_OptWindow_sigPlusbgr = 0; // Get the number of predicted events in the optimal mass window for s+b
	double binCount_OptWindow_bgr = 0;
	int binCount_optWindow_obs = 0;    // Number of observed data events in the optimal mass window for

  for (int i_bin = 0; i_bin<=h_masswindow->GetNbinsX(); i_bin++ ){

    //-- get full width of mass window (bin center of our histogram) and the number of events in mass window for each event type
    double masswindow_fullwidth = h_masswindow->GetBinCenter(i_bin);


    //-- [a] determine the number of events in the mass window for each event type
    //       Ndata_win, Nbgr_win and Nsig_win
      Int_t binNumber_125 = h_bgr->FindBin(125);

      Double_t x_min = 125-masswindow_fullwidth/2;
      Double_t x_max = 125+masswindow_fullwidth/2;

      Int_t binMin = h_bgr->FindBin(x_min);
      Int_t binMax = h_bgr->FindBin(x_max);
//    Int_t binNumber_100 = h_bgr->FindBin(100);
//    Int_t binNumber_150 = h_bgr->FindBin(125+25);


    Double_t error_sig, error_bgr, error_data;
    double nEvents_bgr = h_bgr->IntegralAndError(binMin,binMax, error_bgr);
    double nEvents_sig = h_sig->IntegralAndError(binMin,binMax, error_sig);
    Int_t nEvents_data = h_data->IntegralAndError(binMin,binMax, error_data);


		//Scale by found background scalefactor
		nEvents_bgr = bgr_scalefactor*nEvents_bgr;


//    cout <<"Bin number at which m_ll = 125GeV:  " << binNumber_125 << endl;



    //-- [b] compute EXPECTED significance and save in histogram
    double pvalue_expected    =  IntegratePoissonFromRight(nEvents_bgr,nEvents_bgr+nEvents_sig); //IntegratePoissonFromRight(\mu,n_obs)

     // you need to do this yourself

    double significance_expected = ROOT::Math::gaussian_quantile_c(pvalue_expected,1);
    h_masswindow_expected->SetBinContent(i_bin, significance_expected);

//	double significance_expected = TMath::Sqrt(4);
//	cout << significance_expected<<endl;

    //-- [c] compute OBSERVED significance and save in histogram
    double pvalue_observed       = IntegratePoissonFromRight(nEvents_bgr, nEvents_data); // you need to do this yourself


    double significance_observed = ROOT::Math::gaussian_quantile_c(pvalue_observed,1);
    h_masswindow_observed->SetBinContent(i_bin, significance_observed);

//   printf("Mass window, full width: %5.2f \n", masswindow_fullwidth);
//   printf("   Trying as mass window: %5.2f GeV ->  ", masswindow_fullwidth);

//printf("Expected significance: %5.2f	||	Observed significance: %5.2f \n",significance_expected, significance_observed);



  if(largest_exp_sig  < significance_expected){
  largest_exp_sig = significance_expected;
  optimal_masswindow_sig = masswindow_fullwidth;
  }

    if(largest_obs_sig  < significance_observed){
  largest_obs_sig = significance_observed;
  optimal_masswindow_obs = masswindow_fullwidth;
  }


  } // end loop over width mass window

  //-- print optimum to the screen
  //



///////////////////////////////////////////////////////////////////////////
  cout << "largest expected significance:  " << largest_exp_sig << endl;
  cout << "Optimal masswindow fullwidth, sig: " << optimal_masswindow_sig << endl;
    cout << "largest observed significance:  " << largest_obs_sig << endl;
  cout << "Optimal masswindow fullwidth, obs: " << optimal_masswindow_obs<< endl;
///////////////////////////////////////////////////////////////////////////



//  //----------------------------------
//  //-- [3] Plot histogram and make gif
//  //----------------------------------
if(plot==1){

  TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);
  canvas1->SetLeftMargin(0.125);
  canvas1->SetBottomMargin(0.125);
	canvas1->SetTopMargin(0.125);
	canvas1->SetGrid();
  canvas1->cd();

  h_masswindow_expected->SetLineColor(1);
  h_masswindow_expected->SetLineWidth(2);
  h_masswindow_observed->SetLineColor(4);
  h_masswindow_observed->SetLineWidth(2);

  h_masswindow_expected->SetAxisRange(-1.,6.,"Y");
  h_masswindow_expected->Draw("l");
  if(fabs(Lumi_scalefactor-1.00)<0.01){
    h_masswindow_observed->Draw("l same");
  }

	EColor sigCol = kBlack;
	if(Lumi_scalefactor>1){
		sigCol = kBlue;
	}

	TLine* optWinSig = new TLine(optimal_masswindow_sig,6, optimal_masswindow_sig,-1.0);
	optWinSig->SetLineWidth(2);
	optWinSig->SetLineColor(sigCol);
	optWinSig->SetLineStyle(2);
	optWinSig->Draw("same");

	TLine* optSigSig = new TLine(0,largest_exp_sig, 25,largest_exp_sig);
	optSigSig->SetLineWidth(2);
	optSigSig->SetLineColor(sigCol);
	optSigSig->SetLineStyle(2);
	optSigSig->Draw("same");

  //-- axes
	//====================================================================================================
	// Axes and legend
	//====================================================================================================
  AddText( 0.900, 0.035, "Mass window[GeV]",0.060, 0.,"right"); // X-axis
  AddText( 0.25, 0.950, "Significance" ,0.060,0.,"right");   // Y-axis
  AddText( 0.45, 0.825, Form("Luminosity scalefactor = %5.1f",Lumi_scalefactor),0.0450, 0.,"left");
	//====================================================================================================

	//====================================================================================================
	//Add largest expected significance and the corresponding mass window as numbers
	//====================================================================================================
	  AddText( 0.3, 0.930, Form("%5.2f",optimal_masswindow_sig),0.0450, 0.,"left"); //Optimal exp. window
		AddText( 0.9, 0.5, Form("%5.2f",largest_exp_sig),0.0450, 0.,"left"); //largest exp. sig
	//====================================================================================================



	//====================================================================================================
	//If luminosity scalefactor larger than 1, don't add text and numbers for observed values:
	//===================================================================================================
  if(fabs(Lumi_scalefactor-1.00)<0.01){
    AddText( 0.700, 0.300, "Observed significance",0.050, 0.,"right",1,4);
		AddText( 0.700, 0.200, "Expected significance",0.050, 0.,"right",1,1);
		AddText( 0.9, 0.7, Form("%5.2f",largest_obs_sig),0.0450, 0.,"left"); //Optimal exp. window
		AddText( 0.2, 0.930, Form("%5.2f",optimal_masswindow_obs),0.0450, 0.,"left"); //largest exp. sig

		TLine* optWinObs = new TLine(optimal_masswindow_obs,6, optimal_masswindow_obs,-1.);
		optWinObs->SetLineWidth(2);
		optWinObs->SetLineColor(kBlue);
		optWinObs->SetLineStyle(2);
		optWinObs->Draw("same");

		TLine* optSigObs = new TLine(0,largest_obs_sig, 25,largest_obs_sig);
		optSigObs->SetLineWidth(2);
		optSigObs->SetLineColor(kBlue);
		optSigObs->SetLineStyle(2);
		optSigObs->Draw("same");
  }
	//====================================================================================================
	//====================================================================================================

  //-- prepare gif
  canvas1->Print(Form("./Significance_Optimization_lumiscalefactor%d.gif",int(Lumi_scalefactor)));
}// end plot condition


	tuple<double,double,double,double> expSig_res = make_tuple(optimal_masswindow_sig, largest_exp_sig, optimal_masswindow_obs,largest_obs_sig);
  return expSig_res;

  //================================
} // end Significance_Optimization()
  //================================




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

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
  TLegend *leg1 = new TLegend(0.65,0.65,0.90,0.85);
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





//====================================================================================================================================
//====================================================================================================================================
//========================= Side Band Fit ============================================================================================
//====================================================================================================================================
//====================================================================================================================================
tuple<double,double,double> SideBandFit(int Irebin, double lowerMass, double upperMass, double massWindow,bool plot){

  printf("\n SideBandFit()\n\n");

  //-------------------------
  //-- [1] Prepare histograms
  //-------------------------
  TH1D *h_bgr, *h_data;
  h_bgr  = GetMassDistribution(1);
  h_data = GetMassDistribution(2);

  //-- rebin histograms if necessary
  h_bgr->Rebin(Irebin);
  h_data->Rebin(Irebin);
  printf(" INFO: Rebinning the histograms with a factor %d. Binwidth is now %5.2f GeV\n", Irebin, h_data->GetBinWidth(1));

  //-----------------------------------------
  //-- [2] Loop over scale factor (alpha_bgr)
  //-----------------------------------------
  int nBins_alpha = 1e4;
  TH1D *h_scalefactor_bgr = new TH1D("h_scalefactor_bgr","",nBins_alpha,0.1,3.1); // you'll want to put more here I guess
  double scalefactor_bgr;

  //----------------------------------------------
  // [2a] Loop 1: loop over scalefactors in alpha
  //----------------------------------------------

  double max_log= 0;
  double alpha_opt;
	double alpha_opt2 = 0;
	double max_lik2 = -1e5;

  int lower_sideband_bin = h_bgr->GetXaxis()->FindBin(lowerMass);
	int upper_sideband_bin = h_bgr->GetXaxis()->FindBin(upperMass);

	cout << "lower mass, sideband: " <<h_bgr->GetXaxis()->GetBinCenter(lower_sideband_bin)<<endl;
	cout << "upper mass, sideband: " <<h_bgr->GetXaxis()->GetBinCenter(upper_sideband_bin)<<endl;

  for (int i_bin_sf = 1; i_bin_sf <=h_scalefactor_bgr->GetNbinsX(); i_bin_sf++){

    //-- determine the scale factor for the background
    scalefactor_bgr = h_scalefactor_bgr->GetBinCenter(i_bin_sf);
//    printf(" Loop 1: I am now trying alpha = %5.2f\n",scalefactor_bgr);

    //-----------------------------------------------------------------------------------
    // [2b] Loop 2: loop over bins in the histogram, compute loglik and save in histogram
    //-----------------------------------------------------------------------------------
    double loglik = 1e-10;
    for (int i_bin = lower_sideband_bin; i_bin <=  upper_sideband_bin; i_bin++){
//      printf("        bin %d, mass = %5.2f, Ndata = %5.2f, Nbgr = %5.2f\n",i_bin, h_bgr->GetBinCenter(i_bin),h_bgr->GetBinContent(i_bin),h_data->GetBinContent(i_bin));
      double data_i = h_data->GetBinContent(i_bin);
      double bgr_i = h_bgr->GetBinContent(i_bin);

      loglik += TMath::Log(TMath::Poisson(data_i, scalefactor_bgr*bgr_i)); //Use alpha*b as expectation value and find the log probability to observe data_i events in this bin

//      loglik += TMath::Poisson(h_bgr->GetBinContent(i_bin) ,scalefactor_bgr*h_bgr->GetBinContent(i_bin)); //First argument of Poisson: its mean/expectation value. Second argument is equivalent to n

    } // end loop over bins

    h_scalefactor_bgr->SetBinContent(i_bin_sf,loglik);
		// max_lik2 = h_scalefactor_bgr->GetXaxis()->GetBinCenter(0);
		// cout << max_lik2 << "    " << loglik << endl;
		if(max_lik2<loglik){
			max_lik2 = loglik;
			alpha_opt2 = scalefactor_bgr;
		}


  } // end loop over scale factors for the background

	cout << "alpha_opt2: "<< alpha_opt2 << endl;

    //----------------------------------------------------
    //-- [3] Interpret the -2Log (Likelihood distribution)
    //----------------------------------------------------


			// Set histogram for max_log-log

    //-- Find minimum
    // todo
//================== Find alpha corresponding to the maximum of the likelihood function ======================
		int maxbin_i = h_scalefactor_bgr->GetMaximumBin();

		alpha_opt = h_scalefactor_bgr->GetBinCenter(maxbin_i);
		max_log = h_scalefactor_bgr->GetBinContent(maxbin_i);
//============================================================================================================

//====================Create histograms for plots and uncertainty estimation==================================
		TH1D* h_logMax_minusLog  = (TH1D*) h_scalefactor_bgr->Clone("h_logMax_minusLog"); //Will hold the "rescaled" likelihood function. I chose to rescaled it by max_log-loglik.
		//This means one can find the +- 1sigma uncertainties where this scaled likelihood is equal to 1/2.
    TH1D* h_scalefactor_bgr_rescaled = (TH1D*) h_bgr->Clone("h_scalefactr_bgr_rescaled");

		TH1D* h_ones = (TH1D*) h_scalefactor_bgr->Clone("h_ones");  //Used to plot a straight line at 1/2 to read off uncertainties. Might change to a more suitable name later.

	    for(int i = 0; i<h_scalefactor_bgr->GetNbinsX()+1; i++){
			h_scalefactor_bgr_rescaled->SetBinContent(i,alpha_opt*h_scalefactor_bgr_rescaled->GetBinContent(i));
	    h_logMax_minusLog->SetBinContent(i,max_log-h_scalefactor_bgr->GetBinContent(i));
			h_ones->SetBinContent(i,0.5);
		 }
//============================================================================================================


//========================Find uncertainties on alpha=========================================================

//////////////////////////////////////////////////////
			int lower_uncertainty_bin = maxbin_i;
			double val_uncertainty_lower = h_logMax_minusLog->GetBinContent(maxbin_i);

			while(val_uncertainty_lower<=0.5){
			 lower_uncertainty_bin--;
			 val_uncertainty_lower = h_logMax_minusLog->GetBinContent(lower_uncertainty_bin);
			}
			double sigma_lower = h_logMax_minusLog->GetBinCenter(lower_uncertainty_bin);
//////////////////////////////////////////////////////
			int upper_uncertainty_bin = maxbin_i;
			double val_uncertainty_upper = h_logMax_minusLog->GetBinContent(maxbin_i);

			while(val_uncertainty_upper<=0.5){
			 upper_uncertainty_bin++;
			 val_uncertainty_upper = h_logMax_minusLog->GetBinContent(upper_uncertainty_bin);
			}

			double sigma_upper = h_logMax_minusLog->GetBinCenter(upper_uncertainty_bin);
//////////////////////////////////////////////////////
//===========================================================================================================

//===========================================================================================================
//====== Count number of background events in signal region for the rescaled histogram ======================
//===========================================================================================================

TH1D *h_bgr2, *h_data2, *h_sig2;
h_bgr2  = GetMassDistribution(1);
h_sig2  = GetMassDistribution(125); //Excess signal events?
h_data2 = GetMassDistribution(2);
// h_data2 = GetMassDistribution(2);

double windowHalfWidth = massWindow/2;
cout << "half width: 	" << windowHalfWidth << endl;


int binLower_count = h_bgr2->FindBin(125-windowHalfWidth);
int binUpper_count = h_bgr2->FindBin(125+windowHalfWidth);

// int binLower_count = h_scalefactor_bgr_rescaled->FindBin(125-windowHalfWidth);
// int binUpper_count = h_scalefactor_bgr_rescaled->FindBin(125+windowHalfWidth);

cout <<" Lower bound for mass window:   "<< h_scalefactor_bgr_rescaled->GetXaxis()->GetBinCenter(binLower_count) << endl;
cout << "Upper bound for mass window:   "<< h_scalefactor_bgr_rescaled->GetXaxis()->GetBinCenter(binUpper_count) << endl;

double background_estimate_unscaled = h_bgr2->Integral(binLower_count,binUpper_count);
// double background_estimate_rescaled = h_scalefactor_bgr_rescaled->Integral(binLower_count,binUpper_count);
double background_estimate_rescaled = background_estimate_unscaled*alpha_opt;

 // Find s for this new background scaling (will have decreased). This doesn't make much sense, as the signal to background
 //ratio should remain constant. This will be accounted for when fitting the signal.
double new_signal_estimate = h_sig2->Integral(binLower_count,binUpper_count);
double data_in_window = h_data2->Integral(binLower_count,binUpper_count);


double bCount_uncertMinus_rescaled = background_estimate_unscaled*(alpha_opt-sigma_lower);
double bCount_uncertPlus_rescaled = background_estimate_unscaled*(sigma_upper-alpha_opt);

// double bCount_uncertMinus_unscaled = background_estimate_unscaled*((alpha_opt-sigma_lower)/alpha_opt);
// double bCount_uncertPlus_unscaled = background_estimate_unscaled*(sigma_upper-alpha_opt)/alpha_opt;

cout << "Lower bin in optimal mass window (significance), rescaled bgr: " <<binLower_count <<endl;
cout << "Upper bin in optimal mass window (significance), rescaled bgr: " <<binUpper_count <<endl;
cout << "b_r:   " << background_estimate_rescaled << endl;
cout << "s_r:   " << new_signal_estimate << endl;
cout << "n_obs  " << data_in_window << endl;
cout << "lower uncertainty on b count: " << bCount_uncertMinus_rescaled << endl;
cout << "upper uncertainty on b count: " << bCount_uncertPlus_rescaled << endl;
cout << "b_u:  " << background_estimate_unscaled << endl;




//===========================================================================================================
//======================================Plotting=============================================================

if(plot ==true){
	//Create a TCanvas with two panels. Upper panel for likelihood plot, lower for
	//scaled likelihood for "finding" uncertainties.
	TCanvas* c1 = new TCanvas("c_scale", "Likelihood vs scalefactor", 600, 400);
	c1->cd();
	c1->SetGrid();
	// double leftMarg = 0.1;
	// double bottomMarg = 0.12;
	// c1->SetGrid(1,1);
	// c1->SetLeftMargin(leftMarg);
	// c1->SetBottomMargin(bottomMarg);
  // c1->SetRightMargin(0.1);

  double leftMarg = 0.1;
  double bottomMarg = 0.12;
  c1->SetLeftMargin(leftMarg+0.05);
  c1->SetBottomMargin(bottomMarg);
	// c_scale->Divide(1,2);

  ////////////////////////////////////////////////////////
	//Preparing upper panel plot:
  ////////////////////////////////////////////////////////
	// c_scale->cd(1);


     TAxis* scaleFac_x = h_scalefactor_bgr->GetXaxis();
     TAxis* scaleFac_y = h_scalefactor_bgr->GetYaxis();
     // scaleFac_x->SetTitleOffset(0.7);
     // scaleFac_y->SetTitleOffset(0.6);


		 int Nbins_scaleFac = scaleFac_x->GetNbins();
		 scaleFac_x->SetRange(1,Nbins_scaleFac-5);


     h_scalefactor_bgr->SetLineColor(kRed);
     h_scalefactor_bgr->SetTitle("Log likelihood vs scalefactor");


     int a = scaleFac_y->FindBin(max_log);
    	h_scalefactor_bgr->Draw();
			//-- some nice axes and add legend
      AddText( leftMarg+0.07,0.950, "ln[L(#alpha)]",0.060,0.,"right"); // Y-axis
			AddText( 0.900, 0.035, "Scalefactor, #alpha",0.060, 0.,"right"); // X-axis
//h_scalefactor_bgr->GetBinContent(0)
      TLine* max = new TLine(alpha_opt, -695, alpha_opt, -190);
      max ->SetLineColor(kBlue);
      max->SetLineWidth(2);
      max->SetLineStyle(2);
      max->Draw("same");

      // TLine* vertical_upper = new TLine(sigma_upper, 0, sigma_upper, 1.01);
      // vertical_upper->SetLineColor(kBlue);
      // vertical_upper->SetLineWidth(2);
      // vertical_upper->SetLineStyle(2);
      // vertical_upper->Draw("same");

      // AddText(0.3, 0.930, "1.05",0.0450, 0.,"right"); // X-axis
      AddText(0.5, 0.930, "1.11",0.0450, 0.,"right"); // X-axis
      //AddText(0.7, 0.930, "1.18",0.0450, 0.,"right"); // X-axis

     ////////////////////////////////////////////////////////
		 //Preparting lower panel plot:
		 ////////////////////////////////////////////////////////

     // c_scale->cd(2);
		 TCanvas* c2 = new TCanvas("c2", "", 600, 400);
		 c2->cd();
     c2->SetGrid();
		 double leftMarg2 = 0.1;
		 double bottomMarg2 = 0.12;
		 c2->SetLeftMargin(leftMarg2+0.05);
		 c2->SetBottomMargin(bottomMarg2);
     // c2->SetRightMargin(0.12);

     TAxis* rescaledX = h_logMax_minusLog->GetXaxis();
		 // rescaledX->SetTitle("Scalefactor, #alpha");
		 // rescaledX->SetTitleOffset(0.8);

     TAxis* rescaledY = h_logMax_minusLog->GetYaxis();
     // rescaledY->SetTitle("ln(#alpha_{max})-ln(#alpha)");
     // rescaledY->SetTitleOffset(0.6);



     // double C_L = 95/300;
     // double C_U = 115/300;
		 //
     // int Lb = int(C_L*nBins_alpha);
     // int Ub = int(C_U*nBins_alpha);
     rescaledX->SetRange(lower_uncertainty_bin-80,upper_uncertainty_bin+80);
//      rescaledX->SetRange(3000,4000);
//     rescaledX->SetNdivisions(20,0,0,0);
     TAxis* onesX = h_ones->GetXaxis();
     TAxis* onesY = h_ones->GetYaxis();
     onesX->SetNdivisions(-414);


     onesY->SetTitle("lnmax");
		 //Set colors:
		 h_logMax_minusLog->SetLineColor(kRed);
     h_ones->SetLineColor(kBlack);

		 //Draw lower panel histograms:
		 h_logMax_minusLog->Draw();
	   h_ones->Draw("same");

		 AddText( 0.900, 0.035, "Scalefactor, #alpha",0.060, 0.,"right"); // X-axis
		 AddText( leftMarg2+0.11,0.950, "ln(#alpha_{max})-ln(#alpha)",0.060,0.,"right"); // Y-axis

//======================= Adding lines and numbers =============================
     // TLine* horizontal = new TLine(0.98,mu_max, 1.25, mu_max);
     // horizontal->SetLineColor(kBlack);
     // horizontal->SetLineWidth(2);
     // horizontal->SetLineStyle(2);
     // horizontal->Draw("same");

     TLine* vertical_lower = new TLine(sigma_lower, 0, sigma_lower, 1.01);
     vertical_lower->SetLineColor(kBlue);
     vertical_lower->SetLineWidth(2);
     vertical_lower->SetLineStyle(2);
     vertical_lower->Draw("same");

     TLine* vertical_middle = new TLine(alpha_opt, 0, alpha_opt, 1.01);
     vertical_middle ->SetLineColor(kBlue);
     vertical_middle->SetLineWidth(2);
     vertical_middle->SetLineStyle(2);
     vertical_middle->Draw("same");

     TLine* vertical_upper = new TLine(sigma_upper, 0, sigma_upper, 1.01);
     vertical_upper->SetLineColor(kBlue);
     vertical_upper->SetLineWidth(2);
     vertical_upper->SetLineStyle(2);
     vertical_upper->Draw("same");

     AddText(0.3, 0.930, "1.05",0.0450, 0.,"right"); // X-axis
     AddText(0.5, 0.930, "1.11",0.0450, 0.,"right"); // X-axis
     AddText(0.7, 0.930, "1.18",0.0450, 0.,"right"); // X-axis

     //
     // AddText(0.9, 0.8, "1.28",0.0450, 0.,"right"); // X-axis

     cout << sigma_lower << endl;
     cout << sigma_upper << endl;

	 } //End If-statement for plotting
//===========================================================================================================
//========================Print and return optimal alpha and uncertainties===================================
				  cout << "Optimal alpha: " << alpha_opt << endl;
					cout << "Lower uncertainty: " << alpha_opt-sigma_lower << endl;
					cout << "Upper uncertainty: " << sigma_upper-alpha_opt<< endl;

	tuple<double,double,double> results = make_tuple(alpha_opt,alpha_opt-sigma_lower,sigma_upper-alpha_opt);
  return results;
	//===========================================================================================================

  //==================
} // end SideBandFit()
  //==================



//============================================================================================================
//===================================alpha and mu optimization================================================
 tuple<double, double, double> alpha_mu(int rescale, int Nbins){
using namespace TMath;
//
 //TH2D* h2 = new TH2D("h2", "", Nbins, 0.5,2., Nbins, 0.,5.); //Alphas along x-axis, mu's along y
TH2D* h2 = new TH2D("h2", "", Nbins, 0.98,1.25, Nbins, 0, 3); //Alphas along x-axis, mu's along y
//TH2D* h2 = new TH2D("h2", "", Nbins, 1.0,1.3, Nbins, 0, 2.5); //Alphas along x-axis, mu's along y

TH1D* h_bgr = GetMassDistribution(1); h_bgr->Rebin(10);
TH1D* h_data = GetMassDistribution(2); h_data->Rebin(10);
TH1D* h_sig = GetMassDistribution(125); h_sig->Rebin(10);

double lnlik = 1e-9;
double NbinsAlpha = h2->GetXaxis()->GetNbins();
double NbinsMu = h2->GetYaxis()->GetNbins();
double Nbins_bgr = h_bgr->GetXaxis()->GetNbins();

double n, a, u, bgr, sig;
double maxlik = -1e6;

for(int i = 0; i<=NbinsAlpha; i++){
	for(int j = 0; j<=NbinsMu; j++){
		for(int k = 0; k<=Nbins_bgr; k++){
			a = h2->GetXaxis()->GetBinCenter(i);
			u = h2->GetYaxis()->GetBinCenter(j);
			n = h_data->GetBinContent(k);
			bgr = h_bgr->GetBinContent(k);
			sig = h_sig->GetBinContent(k);
			lnlik += Log(Poisson(n, a*bgr + u*sig));
		}
		h2->SetBinContent(i,j,lnlik);
		if(maxlik<lnlik){maxlik = lnlik;}
		lnlik = 1e-9;
	}
}

cout << "maxlik: " << maxlik << endl;

int alpha_max_bin, mu_max_bin, zmax, we;



we = h2->GetMaximumBin(alpha_max_bin, mu_max_bin, zmax);

double alpha_max = h2->GetXaxis()->GetBinCenter(alpha_max_bin);
double mu_max = h2->GetYaxis()->GetBinCenter(mu_max_bin);
double lnlikMax = h2->GetBinContent(alpha_max, mu_max);

cout << "alpha_max:  " << alpha_max << endl;
cout << "mu_max:  " << mu_max << endl;
cout << "lnlikMax:  " << lnlikMax << endl;


   // Int_t bla = h2->GetMinimumBin();
   // Int_t x,y,z;
   // h2->GetBinXYZ(bla, x, y, z);


   // printf("The bin having the maximum value is (%d,%d)\n",x,y);
	 // double lnlikMax = h2->GetBinContent(x,y);
	 // cout << lnlikMax << endl;

TH2D* h2_rescaled = (TH2D*) h2->Clone("h2_rescaled");
h2_rescaled->Reset();
//
//
// int alpha_max_bin_p, mu_max_bin_p, zmax_p;
// wep = h2_rescaled->GetMaximumBin(alpha_max_bin_p, mu_max_bin_p, zmax_p);
//
// double alpha_max_p = h2_rescaled->GetXaxis()->GetBinCenter(alpha_max_bin_p);
// double mu_max_p = h2_rescaled->GetYaxis()->GetBinCenter(mu_max_bin_p);
// double lnlikMax_p = h2_rescaled->GetBinContent(alpha_max_p, mu_max_p);

//Rescale histogram to get lnL_max-lnlik contours
for(int i = 0; i<=NbinsAlpha; i++){
	for(int j = 0; j<=NbinsAlpha; j++){
		double binContent = maxlik - h2->GetBinContent(i,j);
		h2_rescaled->SetBinContent(i,j,binContent);
	}
}

TH2D* h2_sigma = (TH2D*) h2->Clone("h2_sigma");
h2_sigma->Reset();

for(int i = 0; i<=NbinsAlpha; i++){
	for(int j = 0; j<=NbinsAlpha; j++){
		if((h2_rescaled->GetBinContent(i,j))<1.){
			// cout << h2_rescaled->GetBinContent(i,j) << endl;
			// h2_sigma->SetBinContent(i,j,lnlikMax-h2->GetBinContent(i,j));
			h2_sigma->SetBinContent(i,j,h2_rescaled->GetBinContent(i,j));
		}
	}
}
//========================Find uncertainties on alpha=========================================================

double alpha_lower;
double alpha_upper;

//////////////////////////////////////////////////////
			int lower_uncertainty_bin_a = alpha_max_bin;
			double val_uncertainty_lower_a = h2_rescaled->GetBinContent(alpha_max_bin, mu_max_bin);

			while(val_uncertainty_lower_a<=0.5){
			 lower_uncertainty_bin_a--;
			 val_uncertainty_lower_a = h2_rescaled->GetBinContent(lower_uncertainty_bin_a,mu_max_bin);
			}
			 alpha_lower = h2_rescaled->GetXaxis()->GetBinCenter(lower_uncertainty_bin_a);
//////////////////////////////////////////////////////
			int upper_uncertainty_bin_a = alpha_max_bin;
			double val_uncertainty_upper_a = h2_rescaled->GetBinContent(alpha_max_bin, mu_max_bin);

			while(val_uncertainty_upper_a<=0.5){
			 upper_uncertainty_bin_a++;
			 val_uncertainty_upper_a = h2_rescaled->GetBinContent(upper_uncertainty_bin_a, mu_max_bin);
			}

			alpha_upper = h2_rescaled->GetXaxis()->GetBinCenter(upper_uncertainty_bin_a);
//////////////////////////////////////////////////////

cout << "DeltaPlus_alpha:  " << alpha_upper-alpha_max << endl;
cout << "DeltaMinus_alpha: " << alpha_max - alpha_lower << endl;
//===========================================================================================================

//========================Find uncertainties on alpha=========================================================
double mu_lower;
double mu_upper;
//////////////////////////////////////////////////////
			int lower_uncertainty_bin_u = mu_max_bin;
			double val_uncertainty_lower_u = h2_rescaled->GetBinContent(alpha_max_bin, mu_max_bin);

			while(val_uncertainty_lower_u<=0.5){
			 lower_uncertainty_bin_u--;
			 val_uncertainty_lower_u = h2_rescaled->GetBinContent(alpha_max_bin,lower_uncertainty_bin_u);
			}
			 mu_lower = h2_rescaled->GetYaxis()->GetBinCenter(lower_uncertainty_bin_u);
//////////////////////////////////////////////////////
			int upper_uncertainty_bin_u = mu_max_bin;
			double val_uncertainty_upper_u = h2_rescaled->GetBinContent(alpha_max_bin, mu_max_bin);

			while(val_uncertainty_upper_u<=0.5){
			 upper_uncertainty_bin_u++;
			 val_uncertainty_upper_u = h2_rescaled->GetBinContent(alpha_max_bin, upper_uncertainty_bin_u);
			}

			mu_upper = h2_rescaled->GetYaxis()->GetBinCenter(upper_uncertainty_bin_u);
//////////////////////////////////////////////////////

cout << "DeltaPlus_mu:  " << mu_upper-mu_max << endl;
cout << "DeltaMinus_mu: " << mu_max - mu_lower << endl;
//===========================================================================================================


//==========================================================================================
//============================= 3D surface plot=============================================
//==========================================================================================
TCanvas* C3d = new TCanvas("C3d", "", 600, 400);
C3d->cd();
C3d->SetGrid();
// gStyle->Reset();
// TH2D* h2_rescaled2 = (TH2D*) h2_rescaled->Clone("h_rescaled2");

h2->Draw("surf3");
//
//
// //====================================================================================
// //==============================Contour Plotting==============================================
// //====================================================================================
//
// // int palette[8];
// // palette[0] = kRed;
// // palette[1] = 20;
// // palette[2] = 23;
// // palette[3] = 30;
// // palette[4] = 32;
// // palette[5] = 35;
// // // palette[6] = 40;
// // palette[6] = 40;
// // palette[7] = 40;
// // // palette[8] = 45;
//
//
// //
// // int palette[17];
// // palette[0] = kRed+0;
// // palette[1] = kRed+1;
// // palette[2] = kRed+2;
// // palette[3] = kRed+3;
// // palette[4] = kRed+4;
// // palette[5] = kMagenta+3;
// // // palette[6] = 40;
// // palette[6] = kMagenta+2;
// // palette[7] = kMagenta+1;
// // palette[8] = kMagenta+0;
// // palette[9] = kMagenta-1;
// // palette[10] = kBlue-1;
// // palette[11] = kBlue-2;
// // palette[12] = kBlue-3;
// // palette[13] = kBlue-4;
// // palette[14] = kBlue+1;
// // palette[15] = kBlue+2;
// // palette[16] = kBlue+3;
//
// // int palette[17];
// // EColor first = kRed; EColor second = kMagenta; EColor third = kBlue;
// // // first.SetAlpha(1);
// // palette[0] = first+0;
// // palette[1] = first+1;
// // palette[2] = first+2;
// // palette[3] = first+3;
// // palette[4] = first+4;
// // palette[5] = second+3;
// // palette[6] = second+2;
// // palette[7] = second+1;
// // palette[8] = second+0;
// // palette[9] = second-1;
// // palette[10] = third-1;
// // palette[11] = third-2;
// // palette[12] = third-3;
// // palette[13] = third-4;
// // palette[14] = third+1;
// // palette[15] = third+2;
// // palette[16] = third+3;
//
// int palette[17];
// EColor onlyCol = kBlue;
// palette[0] = onlyCol;
// palette[1] = onlyCol;
// palette[2] = onlyCol;
// palette[3] = onlyCol;
// palette[4] = onlyCol;
// palette[5] = onlyCol;
// palette[6] = onlyCol;
// palette[7] = onlyCol;
// palette[8] = onlyCol;
// palette[9] = onlyCol;
// palette[10] = onlyCol;
// palette[11] = onlyCol;
// palette[12] = onlyCol;
// palette[13] = onlyCol;
// palette[14] = onlyCol;
// palette[15] = onlyCol;
// palette[16] = onlyCol;
//
// double levels[17] = {0.0,0.5,1.0,1.5,2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0};
// TH1D* rescaled2 = (TH1D*) h2_rescaled->Clone("rescaled2");
// h2_rescaled->SetContour(17,levels);
// // h2->SetContourLevel(0,-64);
// // h2->SetContourLevel(1,-65);
//
// gStyle->SetPalette(17,palette);
//
//
// // h2_rescaled->Draw("CONT1");
// // h2_rescaled->GetZaxis()->SetTickSize(0.01);
// // h2_rescaled->GetZaxis()->SetLabelOffset(0.01);
// h2_rescaled->GetZaxis()->SetRange(0,5);
// // rescaled2->Draw("cont4z");
// //h2_rescaled->Draw("same cont1, CJUST");
//
// TCanvas* C = new TCanvas("c", "", 600, 400);
// C->SetBottomMargin(0.125);
// C->SetGrid();
//
// h2_rescaled->Draw("CONT1");
//
//
// // TMarker* m = new TMarker(alpha_max,mu_max,47);
// // m->SetMarkerSize(1.3);
// // m->Draw("same L");
//
// TLine* horizontal = new TLine(0.98,mu_max, 1.25, mu_max);
// horizontal->SetLineColor(kBlack);
// horizontal->SetLineWidth(2);
// horizontal->SetLineStyle(2);
// horizontal->Draw("same");
//
// TLine* vertical = new TLine(alpha_max, 0, alpha_max, 3);
// vertical->SetLineColor(kBlack);
// vertical->SetLineWidth(2);
// vertical->SetLineStyle(2);
// vertical->Draw("same");
//
// TLine* alpha_low = new TLine(alpha_lower,0, alpha_lower, 3);
// alpha_low->SetLineColor(kRed);
// alpha_low->SetLineWidth(3);
// alpha_low->SetLineStyle(2);
// alpha_low->Draw("same");
//
// TLine* alpha_up = new TLine(alpha_upper,0, alpha_upper, 3);
// alpha_up->SetLineColor(kRed);
// alpha_up->SetLineWidth(3);
// alpha_up->SetLineStyle(2);
// alpha_up->Draw("same");
//
// TLine* mu_low = new TLine(0.98,mu_lower, 1.25, mu_lower);
// mu_low->SetLineColor(kGreen);
// mu_low->SetLineWidth(3);
// mu_low->SetLineStyle(2);
// mu_low->Draw("same");
//
// TLine* mu_up = new TLine(0.98,mu_upper,1.25, mu_upper);
// mu_up->SetLineColor(kGreen);
// mu_up->SetLineWidth(3);
// mu_up->SetLineStyle(2);
// mu_up->Draw("same");
//
// cout << "alpha_lower: "<<alpha_lower << endl;
// cout << "alpha_upper: "<<alpha_upper << endl;
//
// cout << "mu_lower: " << mu_lower << endl;
// cout << "mu_upper: " << mu_upper << endl;


AddText( 0.900, 0.035, "Background scalefactor, #alpha",0.060, 0.,"right"); // X-axis
AddText( 0.33, 0.95, "Signal scalefactor, #mu" ,0.060,0.,"right");   // Y-axis
AddText( 0.33, 0.95, "ln L(#alpha,#mu;n)" ,0.060,0.,"right");   // Y-axis

// AddText(0.3, 0.930, "1.05",0.0450, 0.,"right"); // X-axis
// AddText(0.5, 0.930, "1.11",0.0450, 0.,"right"); // X-axis
// AddText(0.7, 0.930, "1.17",0.0450, 0.,"right"); // X-axis
//
//
// AddText(0.9, 0.5, "0.74",0.0450, 0.,"right"); // X-axis
// AddText(0.9, 0.7, "1.28",0.0450, 0.,"right"); // X-axis
// AddText(0.9, 0.9, "1.92",0.0450, 0.,"right"); // X-axis
//AddText(0.9, 0.8, "3.92",0.0450, 0.,"right"); // X-axis
//h2->Draw("surf4z");
// h2_rescaled->Draw("surf4z");


// h2->Draw("surf2");
// h2_rescaled->Draw("cont1z");
// h2_sigma->SetFillColor(kGreen);
// h2_sigma->Draw("same L");

// h2->Draw("surf4");


tuple<double, double, double> results = make_tuple(1,1,1);
return results;

}






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


// void AddText( Double_t txt_x = 0.50, Double_t txt_y = 0.50, const char * txt = "dummy", Double_t txt_size = 0.045,
// 	      Double_t txt_angle = 0., const char * Alignment = "left", Int_t UseNormalizedSize = 1, Int_t txt_color =1 );
#include "TArrow.h"
void func(){

//  gROOT->Clear();
//  gROOT->Delete();
//
//   auto c2 = new TCanvas("c2","c2",200,10,600,400);
//
//   c2->SetFillColor(42);
//   c2->SetGrid();



double delta_L0 = 0.023;
int n0 = 500;

double C = delta_L0*n0;
// we want DeltaL * n = DeltaL_0 * n_0
int n = 100;
double delta_L = C/double(n);  //Find the increment in L which results in a good plot for the chosen n

double largest_sigs[100];
double luminosities[100];

     // gStyle->SetTitleFontSize(0.03);
//     gStyle->SetLabelFont(50);



TH1D* fiveSig = new TH1D("fiveSig", "", n, 0, n*delta_L);
TH1D* sigVSlum = new TH1D("sigVSlum","", n, 0, n*delta_L);
TH1D* sqrtLum = new TH1D("sqrtLum", "", n , 0, n*delta_L);

tuple<double,double,double,double> results;

for(int i = 0; i<=n; i++){

results = Significance_Optimization(delta_L*i,1, false);
// largest_sigs[i] = Significance_Optimization(delta_L*i);
largest_sigs[i] = get<1>(results);
luminosities[i] = delta_L*double(i);



sigVSlum->SetBinContent(i,largest_sigs[i]);
fiveSig->SetBinContent(i,5);
sqrtLum->SetBinContent(i,sqrt(delta_L*i));


}

   // TGraph* gr = new TGraph(n, luminosities, largest_sigs);
   //
   // gr->SetTitle("Lum vs sig; Luminosity; Significance");

 //  TAxis* m_xaxis = gr->GetXaxis();
//   m_xaxis->SetTitle("Luminosity");
//   m_xaxis->SetTitleFontSize(0.5);



//   gr->GetYaxis()->SetTitle("Significance");
// gROOT->Clear();
// gROOT->Delete();
   TCanvas* c = new TCanvas("c","",600,400);
   c->SetLeftMargin(0.09);
   c->SetBottomMargin(0.1);
   c->SetGrid(1,1);
   c->SetGridx(1);
   c->SetGridy(1);
   // gr->GetYaxis()->SetTitleOffset(1.0);
   // gr->GetXaxis()->SetTitleOffset(1.0);

 //  AddText( 0.900, 0.535, "Mass window GeV",0.060, 0.,"right"); // X-axis

   // c->Update();
   c->cd();
   sigVSlum->Draw("l");
   fiveSig->SetLineColor(kGreen);
   fiveSig->Draw("l same");

   TLine* reqlum = new TLine(5.36,7.65,5.39,0.00);
   reqlum->SetLineStyle(2);
   reqlum->SetLineWidth(2);
   reqlum->SetLineColor(kBlue);
   reqlum->Draw("same");

   sqrtLum->SetLineStyle(9);
   sqrtLum->Draw("same"); 

   // TArrow* alpha_max_arrow = new TArrow(5.36, 7,5.36,0.1,0.02, "|>");
   //
   // alpha_max_arrow->SetLineWidth(2);
   //
   // alpha_max_arrow->Draw("same |>");

   TLegend *leg = new TLegend(0.7, 0.6-0.3, 0.88, 0.88-(0.3+0.1));
   leg->AddEntry(sigVSlum, "Z(L)", "l");
   leg->AddEntry(fiveSig, "5#sigma discovery requirement", "l");

// leg->AddEntry(alpha_max_arrow,"Expected p-value", "f");
//   leg->AddEntry(fit, "Fit function", "l");
leg->Draw("same");

   AddText( 0.950, 0.025, "Integrated Luminosity scalefactor",0.060, 0.,"right"); // X-axis
   AddText( 0.35, 0.950, "Expected Significance" ,0.060,0.,"right");   // Y-axis
//   AddText( 0.60, 0.8, "Discovery int. lum." ,0.060,0.,"right");   // Y-axis

   AddText(0.3, 0.930, "5.38",0.0450, 0.,"right"); // X-axis
// AddText(0.9, 0.8, "4.79",0.0450, 0.,"right"); // X-axis


   // gr->Draw("al");

//   gr = new TGraph(n,x,y);
//   gr->SetLineColor(2);
//   gr->SetLineWidth(4);
//   gr->SetMarkerColor(4);
//   gr->SetMarkerSize(1.5);
//   gr->SetMarkerStyle(21);
//   gr->SetTitle("Option ACP example");
//   gr->GetXaxis()->SetTitle("X title");
//   gr->GetYaxis()->SetTitle("Y title");
//   gr->Draw("ACP");
//
//      // TCanvas::Update() draws the frame, after which one can change it
//   c2->Update();
//   c2->GetFrame()->SetFillColor(21);
//   c2->GetFrame()->SetBorderSize(12);
//   c2->Modified();

}

// //=======================================================================================================================
// void AddText( Double_t txt_x, Double_t txt_y, const char * txt, Double_t txt_size,
//               Double_t txt_angle, const char * Alignment, Int_t UseNormalizedSize, Int_t txt_color)
// //=======================================================================================================================
//
// {
//   Int_t txt_align = 12;
//   if ( !strcmp(Alignment, "left"))   { txt_align = 12; } // left
//   if ( !strcmp(Alignment, "right"))  { txt_align = 32; } // right
//   if ( !strcmp(Alignment, "center")) { txt_align = 22; } // center
//
//   TLatex* t1 = new TLatex( txt_x, txt_y, txt);
//   if(UseNormalizedSize) {t1->SetNDC(kTRUE);} // <- use NDC coordinate
//   t1->SetTextSize(txt_size);
//   t1->SetTextAlign(txt_align);
//   t1->SetTextAngle(txt_angle);
//   t1->SetTextColor(txt_color);
//   t1->Draw();

//} // end AddText()

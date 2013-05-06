#include "TFile.h"
#include "TChain.h"
#include "THStack.h"
#include "TF1.h"
#include "TH1.h"
#include "TH3.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TList.h"
#include "TLatex.h"
#include "TLine.h"
#include "TObject.h"
#include "TDirectory.h"
#include "TGraphAsymmErrors.h"
#include "TKey.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <exception>
#include <cmath> 
#include <iomanip>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sstream>

#include "Math/QuantFuncMathCore.h"
#include "Math/QuantFuncMathMore.h"
#include "TMath.h"

using namespace std;

const int Nsample_all = 31;
const int Nsample_comb = 6;
double lumi_err = 0.044;
double rate_err_ttbb = 0.5;
double frac_rate_err_ttbb_2 = rate_err_ttbb * rate_err_ttbb;
int lepSel = 2;
const double alpha = 1 - 0.6827;

// Here set the location of your input histogram files:

TString baseDir = "/data2/ttH/";
TString basePrefix = "Trees_withoutCVSreshape4/histos/yggdrasil_treeReader_53x_newCSV_v14_3rd40_topPtWgt_using_v7_trees";

void PrintMessage();
TFile* openTheFiles(int lepSelection);
void closeTheFiles();

class MySample{
public:
  TString sampleName;
  Color_t sampleColor;
  double QCDscale_unc;
  double PDF_unc;
  TString process;
  MySample();
  MySample(TString,Color_t);
  void set_up(int lepSelection, int index);
  void combine();
  
};

class MyPlot{
public:
  MyPlot();
  MyPlot(string);
  string plotName;
  double plotYMax;
  std::vector<std::string> sys_cat_labels;
  TH1D *TotalBKG_rateUnc;                   // Total BKG with bin-by-bin systematic uncertainties (eg. PDF)
  TH1D *h_bkg_mc_stat;                      // Total BKG with statistical uncertainty only
  TH1D *h_bkg_err_1sig;                     // Total BKG with Stat+Sys error (hatched band)
  vector<TH1D*> h_sum_bkg_sys;              // Individual systematic uncertainties (samples combined)
  TH1D *h_dataMC_ratio;                     // Ratio histogram
  TH1D *h_dataMC_ratio_1sig;                // Ratio error band (green band)
  TGraphAsymmErrors *asymmRatio;            // Ratio with asymmetric errors
  TGraphAsymmErrors * g;                    // Asymmetric errors for data

  TH1D *h_ewk;
  TH1D *h_ttbar_bbbar;
  TH1D *h_ttbar_ccbar;
  TH1D *h_ttbar_lf;
  TH1D *h_singlet;
  TH1D *h_ttV;
  TH1D *h_data;
  TH1D *h_ttH;

  void drawShapeComparison(vector<string> samplesToCompare);
  void drawStack(bool drawData, bool setLog, bool drawLegend);
  vector<double> getDefaultXRange(string name);
  int getRebin(string name);
  string getCategoryLabel();
  string getCategoryDir();
  void getTotalBKG_rateUnc(TH1D *hists[]);
  void getSysShapes(int numSys);
  void getDataMCRatio();
  void setPlotMax(double newMax);
  TCanvas* setUpCanvas(bool withRatio);
};

MySample::MySample(){
};
MyPlot::MyPlot(){
};

TFile* myFiles[Nsample_all];
MySample *allSamples = new MySample[Nsample_all];

MyPlot::MyPlot(string plotname){
  
  h_bkg_err_1sig = 0;
  TotalBKG_rateUnc = 0;
  h_bkg_mc_stat = 0;
  h_ttbar_bbbar = 0;
  h_ttbar_ccbar = 0;
  h_ttbar_lf = 0;
  h_singlet = 0;
  h_ttV = 0;
  h_ewk = 0;
  h_ttH = 0;
  h_data = 0;
  
  plotYMax = -999.;

  plotName=plotname;

  sys_cat_labels.clear();
  sys_cat_labels.push_back("");
  sys_cat_labels.push_back("_lepIdSFUp");
  sys_cat_labels.push_back("_lepIdSFDown");
  sys_cat_labels.push_back("_PUUp");
  sys_cat_labels.push_back("_PUDown");
  sys_cat_labels.push_back("_hfSFUp");
  sys_cat_labels.push_back("_hfSFDown");
  sys_cat_labels.push_back("_lfSFUp");
  sys_cat_labels.push_back("_lfSFDown");
  sys_cat_labels.push_back("_JERUp");
  sys_cat_labels.push_back("_JERDown");
  sys_cat_labels.push_back("_JESUp");
  sys_cat_labels.push_back("_JESDown");
  sys_cat_labels.push_back("_Q2scaleUp");
  sys_cat_labels.push_back("_Q2scaleDown");
  sys_cat_labels.push_back("_HTwgtUp");
  sys_cat_labels.push_back("_HTwgtDown");

  int NumSysCat = int(sys_cat_labels.size());

  // FILL SOME TH1D OBJECTS

  TH1D* hist_sys[Nsample_all][NumSysCat];

  for( int j=0; j<Nsample_all; j++ ){
    for( int b=0; b<NumSysCat; b++ ){
      string tempname = this->plotName;
      if( b!=0 ) tempname += sys_cat_labels[b];
      TH1D* hist_sys_temp = (TH1D*)myFiles[j]->Get(tempname.c_str());
      hist_sys[j][b] = (TH1D*)hist_sys_temp->Clone("temp");
      hist_sys[j][b]->Rebin(this->getRebin(tempname));
    }
  }
  
  TH1D* hist[Nsample_all];
  for( int j=0; j<Nsample_all; j++ ) {
        hist[j] = (TH1D*)hist_sys[j][0]->Clone(Form("hist_data_cat_%d",j));
  } 
  if( hist[0]->Integral()<1 ) cout << "WARNING: Data integral less than zero for " << this->plotName << endl;
  
  int nbins = hist[0]->GetNbinsX();
  int first_bin = 0;
  int last_bin  = nbins+1;

  
  h_bkg_err_1sig = (TH1D*)hist[0]->Clone();
  h_bkg_mc_stat = (TH1D*)hist[0]->Clone();                      
  h_bkg_err_1sig->Reset();
  h_bkg_mc_stat->Reset();
  

  // HISTO: BIN BY BIN ADDING THE NOMINAL BKG CONTENTS AND HARD-CODED ERRORS

  this->getTotalBKG_rateUnc(hist);
  TH1::SetDefaultSumw2();

  // GROUP TOGETHER THE PROCESSES

  h_ttbar_bbbar = (TH1D*)hist[0]->Clone("h_ttbar_bbbar");
  h_ttbar_ccbar = (TH1D*)hist[0]->Clone("h_ttbar_ccbar");
  h_ttbar_lf = (TH1D*)hist[0]->Clone("h_ttbar_lf");
  h_singlet = (TH1D*)hist[0]->Clone("h_singlet");
  h_ttV = (TH1D*)hist[0]->Clone("h_ttV");
  h_ewk = (TH1D*)hist[0]->Clone("h_ewk");
  h_ttH = (TH1D*)hist[0]->Clone("h_ttH");
  h_data = (TH1D*)hist[0]->Clone("h_data");
  
  h_ttbar_bbbar->Reset();
  h_ttbar_ccbar->Reset();
  h_ttbar_lf->Reset();
  h_singlet->Reset();
  h_ttV->Reset();
  h_ewk->Reset();
  h_ttH->Reset();
  h_data->Reset();
  
  for( int iSample=0; iSample<Nsample_all; iSample++ ){
    if(allSamples[iSample].process=="EWK"){
      h_ewk->Add(hist[iSample]);
      h_ewk->SetFillColor(allSamples[iSample].sampleColor);
      h_ewk->SetLineColor(allSamples[iSample].sampleColor);
    }
    if(allSamples[iSample].process=="singleT"){
      h_singlet->Add(hist[iSample]);
      h_singlet->SetFillColor(allSamples[iSample].sampleColor);
      h_singlet->SetLineColor(allSamples[iSample].sampleColor);
    }
    if(allSamples[iSample].process=="TTV"){
      h_ttV->Add(hist[iSample]);
      h_ttV->SetFillColor(allSamples[iSample].sampleColor);
      h_ttV->SetLineColor(allSamples[iSample].sampleColor);
    }    
    if(allSamples[iSample].process=="TTBB"){
      h_ttbar_bbbar->Add(hist[iSample]);
      h_ttbar_bbbar->SetFillColor(allSamples[iSample].sampleColor);
      h_ttbar_bbbar->SetLineColor(allSamples[iSample].sampleColor);
    }
    if(allSamples[iSample].process=="TTCC"){
      h_ttbar_ccbar->Add(hist[iSample]);
      h_ttbar_ccbar->SetFillColor(allSamples[iSample].sampleColor);
      h_ttbar_ccbar->SetLineColor(allSamples[iSample].sampleColor);
    }
    if(allSamples[iSample].process=="TTLF"){
      h_ttbar_lf->Add(hist[iSample]);
      h_ttbar_lf->SetFillColor(allSamples[iSample].sampleColor);
      h_ttbar_lf->SetLineColor(allSamples[iSample].sampleColor);
    }
    if(allSamples[iSample].process=="signal"){
      h_ttH->Add(hist[iSample]);
      h_ttH->SetLineColor(allSamples[iSample].sampleColor);
    }
    if(allSamples[iSample].process=="data"){
      h_data->Add(hist[iSample]);
      h_data->SetFillColor(allSamples[iSample].sampleColor);
      h_data->SetLineColor(allSamples[iSample].sampleColor);
    }
  }
  
  
  // MAKE STATS-ONLY HISTOGRAM FOR SUM-MC

  h_bkg_mc_stat = (TH1D*)h_ewk->Clone("h_bkg_mc_stat");
  h_bkg_mc_stat->Add(h_ttV);
  h_bkg_mc_stat->Add(h_singlet);
  h_bkg_mc_stat->Add(h_ttbar_bbbar);
  h_bkg_mc_stat->Add(h_ttbar_ccbar);
  h_bkg_mc_stat->Add(h_ttbar_lf);
  
  
  // MAKE HISTOGRAM ADDING ALL THE SYSTEMATICS OF THE BACKGROUNDS TOGETHER
  
  h_sum_bkg_sys.clear();

  for( int b=0; b<NumSysCat; b++ ){
    bool firstSample = true;
    TH1D* temphisto = (TH1D*)hist_sys[1][b]->Clone(Form("hist_sys_%d",b));;
    for( int iSample=1; iSample<Nsample_all; iSample++ ){
      bool SkipSample = false;
      if(allSamples[iSample].process!="signal"){
	if( firstSample )firstSample = false;
	else temphisto->Add(hist_sys[iSample][b]);
      }
    }
    h_sum_bkg_sys.push_back(temphisto);
  }
  
  // MAKING THE HATCHED HISTOGRAM FOR BKG ERROR

  this->getSysShapes(NumSysCat);
  

};

void PrintMessage(){

  std::cout << "hello world" << std::endl;


}

TFile* openTheFiles(int lepSelection){

  std::cout << "hello world" << std::endl;

  TString lepType = "";
  if( lepSelection==0 )      lepType = "mu_sel";
  else if( lepSelection==1 ) lepType = "ele_sel";
  else                       lepType = "lep_sel";

  TString baseSuffix = "_" + lepType + "_histo_1.root";

  for( int iSample=0; iSample<Nsample_all; iSample++ ){
    TString typePrefix = ( allSamples[iSample].sampleName.Contains("Single") ) ? "_data_" : "_mc_";
    myFiles[iSample] = new TFile(baseDir+basePrefix+typePrefix+allSamples[iSample].sampleName+baseSuffix);
  }

  return *myFiles;

}

void closeTheFiles(){

  for( int j=0; j<Nsample_all; j++ ){
    myFiles[j]->Close();
  }
}


void MyPlot::drawShapeComparison(vector<string> samplesToCompare){

  TCanvas *c1 = this->setUpCanvas(false);
  
  
  vector<TH1D*> samplehistos;
  double max_tot = 0.;

  for(int k = 0; k<int(samplesToCompare.size()); k++){

    cout << "sample name = " << samplesToCompare[k] << endl;
    if(samplesToCompare[k]=="data") samplehistos[k] = (TH1D*)h_data->Clone("samplehistos");
    else if(samplesToCompare[k]=="ttbb") samplehistos.push_back(h_ttbar_bbbar);
    else if(samplesToCompare[k]=="ttcc")  samplehistos.push_back(h_ttbar_ccbar);
    else if(samplesToCompare[k]=="ttlf")  samplehistos.push_back(h_ttbar_lf);
    else if(samplesToCompare[k]=="ttV")  samplehistos.push_back(h_ttV);
    else if(samplesToCompare[k]=="ewk")  samplehistos.push_back(h_ewk);
    else if(samplesToCompare[k]=="ttH")  samplehistos.push_back(h_ttH);
    else cout << "Sample "<< samplesToCompare[k] << " not found." << endl;
  
    samplehistos[k]->Scale(1./samplehistos[k]->Integral());

    samplehistos[k]->SetFillColor(0);
    samplehistos[k]->SetLineWidth(2);
    
    if (k==0){
      samplehistos[k]->SetStats(0);
      samplehistos[k]->GetYaxis()->SetTitle("");
      samplehistos[k]->GetYaxis()->SetTitleSize(0.05);
      samplehistos[k]->GetYaxis()->SetTitleOffset(1.);
      max_tot = samplehistos[k]->GetBinContent(samplehistos[k]->GetMaximumBin());
      vector<double> range = this->getDefaultXRange(this->plotName);    
      if((range[1]-range[0])>0.0001) samplehistos[k]->GetXaxis()->SetRangeUser(range[0],range[1]);

    }    
    else max_tot = std::max(max_tot,samplehistos[k]->GetBinContent(samplehistos[k]->GetMaximumBin()));
  }
  if(plotYMax>0.)samplehistos[0]->SetMaximum(this->plotYMax);
  else samplehistos[0]->SetMaximum(max_tot*1.4);
  
  for(int k = 0; k<int(samplesToCompare.size()); k++){ 
    if(k==0) samplehistos[k]->Draw("hist");
    else samplehistos[k]->Draw("same hist");
  }

  string categorydir = this->getCategoryDir(); 
  TString outputHistoname = this->plotName.c_str();
  c1->Print("Images/"+categorydir+outputHistoname+"_compare.pdf");
  
  delete c1;



}

void MyPlot::drawStack(bool drawData, bool setLog, bool drawLegend){

  cout << "drawing..." << endl;
  
  char *countsname = new char[200];
  char *sumname = new char[200];

  int nbins = this->h_data->GetNbinsX();
  int first_bin = 0;
  int last_bin  = nbins+1;

  cout << this->h_data->Integral() << endl;

  // GET TOTAL ERROR FOR YIELDS PRINTOUT (TAKE LARGEST, NOT AVERAGE)

  double sum_up2 = 0, sum_down2 = 0;
  double sum_previous_diff = 0.;
  for( int b=1; b<sys_cat_labels.size(); b++ ){
    double diff = this->h_sum_bkg_sys[b]->Integral(first_bin,last_bin) - this->h_sum_bkg_sys[0]->Integral(first_bin,last_bin);
    double diff_2 = diff*diff;
          if( b!=1 && b%2==0 && (diff*sum_previous_diff)>0 ){
            if( fabs(diff)>fabs(sum_previous_diff) ) diff_2 = diff*diff - sum_previous_diff*sum_previous_diff;
          } 
          if( diff>0 ) sum_up2 += diff_2;
          else if( diff<0. ) sum_down2 += diff_2;
          sum_previous_diff = diff;
  }
  
  
  double sum_sys_err_up   = sqrt( sum_up2 );
  double sum_sys_err_down = sqrt( sum_down2 );
  double sum_sys_err = ( sum_sys_err_up>sum_sys_err_down ) ? sum_sys_err_up : sum_sys_err_down;
  


  // FIND APPROPRIATE HISTO MAXIMUMS
  
  int max_bin_data = this->h_data->GetMaximumBin();
  double max_data = this->h_data->GetBinContent(max_bin_data) + this->h_data->GetBinError(max_bin_data);
  int max_bin_bkg = TotalBKG_rateUnc->GetMaximumBin();
  double max_bkg  = TotalBKG_rateUnc->GetBinContent(max_bin_bkg) + TotalBKG_rateUnc->GetBinError(max_bin_bkg);
  
  int min_bin_data = this->h_data->GetMinimumBin();
  double min_data = this->h_data->GetBinContent(min_bin_data) + this->h_data->GetBinError(min_bin_data);
  int min_bin_bkg = TotalBKG_rateUnc->GetMinimumBin();
  double min_bkg  = TotalBKG_rateUnc->GetBinContent(min_bin_bkg) + TotalBKG_rateUnc->GetBinError(min_bin_bkg);

  double max_tot = max_bkg;
  if(drawData) max_tot = std::max(max_data,max_bkg);
  double min_tot = std::min(min_data,min_bkg);
  max_tot = std::max( max_tot, this->h_ttH->GetMaximum() );
  min_tot = std::max(min_tot,1E-1);

  // SET UP THE PLOT HEADING TEXT 

  std::string cmsinfo = "CMS                         #sqrt{s} = 8 TeV, L = 19.4 fb^{-1}";
  TLatex CMSInfoLatex(0.18, 0.91, cmsinfo.c_str());
  CMSInfoLatex.SetNDC();
  CMSInfoLatex.SetTextFont(42);
  CMSInfoLatex.SetTextSize(0.055);

  std::string str_lepton = "";
  if( lepSel==0 )      str_lepton = "#mu";
  else if( lepSel==1 ) str_lepton = "e";
  else                 str_lepton = "Lepton";
  string str_suffix = this->getCategoryLabel();
  std::string selectioninfo = str_lepton + " + " + str_suffix;
    
  TLatex SELECTIONInfoLatex(0.48, 0.84, selectioninfo.c_str());
  SELECTIONInfoLatex.SetNDC();
  SELECTIONInfoLatex.SetTextFont(42);
  SELECTIONInfoLatex.SetTextSize(0.05);
  
  TLatex SELECTIONInfoLatex_noRat(0.14, 0.93, selectioninfo.c_str());
  SELECTIONInfoLatex_noRat.SetNDC();
  SELECTIONInfoLatex_noRat.SetTextFont(42);
  SELECTIONInfoLatex_noRat.SetTextSize(0.025);


  // YIELD PRINTOUTS

  double nbkg_total = TotalBKG_rateUnc->Integral(first_bin,last_bin);
  sprintf(countsname,"N_{data} = %d",int(this->h_data->Integral(first_bin,last_bin)+0.001));

  TLatex DATA_INTEGRALInfoLatex(0.39, 0.87, countsname);
  DATA_INTEGRALInfoLatex.SetNDC();
  DATA_INTEGRALInfoLatex.SetTextFont(42);
  DATA_INTEGRALInfoLatex.SetTextSize(0.03);
  
  TLatex DATA_INTEGRALInfoLatex_noRat(0.39, 0.89, countsname);
  DATA_INTEGRALInfoLatex_noRat.SetNDC();
  DATA_INTEGRALInfoLatex_noRat.SetTextFont(42);
  DATA_INTEGRALInfoLatex_noRat.SetTextSize(0.025);

  double sum_mc_err_tot_2 = 0.;
  for( int bin=0; bin<nbins; bin++ ){
    double mcstat = this->h_bkg_mc_stat->GetBinError(bin+1);
    sum_mc_err_tot_2 += mcstat*mcstat;
  }
  sprintf(countsname,"N_{bkg}  = %4.1f #pm%4.1f",nbkg_total,sqrt(sum_sys_err*sum_sys_err + sum_mc_err_tot_2));
  if( nbkg_total>999.9 ) sprintf(countsname,"N_{bkg}  = %4.0f #pm%4.0f",nbkg_total,sqrt(sum_sys_err*sum_sys_err + sum_mc_err_tot_2));
  
  TLatex BKG_INTEGRALInfoLatex(0.39, 0.83, countsname);
  BKG_INTEGRALInfoLatex.SetNDC();
  BKG_INTEGRALInfoLatex.SetTextFont(42);
  BKG_INTEGRALInfoLatex.SetTextSize(0.03);
  
  TLatex BKG_INTEGRALInfoLatex_noRat(0.39, 0.85, countsname);
  BKG_INTEGRALInfoLatex_noRat.SetNDC();
  BKG_INTEGRALInfoLatex_noRat.SetTextFont(42);
  BKG_INTEGRALInfoLatex_noRat.SetTextSize(0.025);
  
  double int_data = this->h_data->Integral(first_bin,last_bin);
  double int_bkg  = TotalBKG_rateUnc->Integral(first_bin,last_bin);
  double int_ttlf = this->h_ttbar_lf->Integral(first_bin,last_bin);
  double int_ttcc = this->h_ttbar_ccbar->Integral(first_bin,last_bin);
  double int_ttbb = this->h_ttbar_bbbar->Integral(first_bin,last_bin);
  double int_singlet = this->h_singlet->Integral(first_bin,last_bin);
  double int_ttbarV  = this->h_ttV->Integral(first_bin,last_bin);
  double int_ewk     = this->h_ewk->Integral(first_bin,last_bin);
  double int_ttH  =  this->h_ttH->Integral(first_bin,last_bin);
  
  double scale_ttH = ( int_ttH!=0 ) ? 30. : 1.;
  this->h_ttH->Scale(scale_ttH);

  sprintf(sumname,"%d",int(int_data+0.1));
  std::string str_sum_data = std::string(sumname);
  sprintf(sumname,"%4.1f",int_bkg);
  std::string str_sum_bkg = std::string(sumname);
  sprintf(sumname,"%4.1f",int_ttlf);
  std::string str_sum_ttlf = std::string(sumname);
  sprintf(sumname,"%4.1f",int_ttcc);
  std::string str_sum_ttcc = std::string(sumname);
  sprintf(sumname,"%4.1f",int_ttbb);
  std::string str_sum_ttbb = std::string(sumname);
  sprintf(sumname,"%4.1f",int_singlet);
  std::string str_sum_singlet = std::string(sumname);
  sprintf(sumname,"%4.1f",int_ttbarV);
  std::string str_sum_ttbarV = std::string(sumname);
  sprintf(sumname,"%4.1f",int_ewk);
  std::string str_sum_ewk = std::string(sumname);
  sprintf(sumname,"%4.1f",int_ttH);
  std::string str_sum_ttH = std::string(sumname);

  std::ostringstream strs_scale_H;
  strs_scale_H << int(scale_ttH + 0.500000001);
  std::string str_scale_H = strs_scale_H.str();
  
  std::string ttH_label1 = "t#bar{t}H125";
  std::string ttH_label2 = "(#times " + str_scale_H + ")";
  

  // SET UP LEGEND
  
  TLegend *legend = new TLegend(0.155,0.75,0.92,0.89);
  
  legend->SetFillColor(kWhite);
  legend->SetLineColor(kWhite);
  legend->SetShadowColor(kWhite);
  legend->SetTextFont(42);
  legend->SetTextSize(0.035);
  legend->SetNColumns(3);

  legend->AddEntry(this->h_ttbar_lf,std::string("t#bar{t} + lf (" + str_sum_ttlf + ")").c_str(),"f");
  legend->AddEntry(this->h_ttbar_ccbar,std::string("t#bar{t} + c#bar{c} (" + str_sum_ttcc + ")").c_str(),"f");
  legend->AddEntry(this->h_ttbar_bbbar,std::string("t#bar{t} + b#bar{b} (" + str_sum_ttbb + ")").c_str(),"f");
  legend->AddEntry(this->h_singlet,std::string("single t (" + str_sum_singlet + ")").c_str(),"f");
  legend->AddEntry(this->h_ttV,std::string("t#bar{t} + W,Z (" + str_sum_ttbarV + ")").c_str(),"f");
  legend->AddEntry(this->h_ewk,std::string("EWK (" + str_sum_ewk + ")").c_str(),"f");
  


  // MAKE THE ASYMMETRIC ERRORS FOR DATA

  TGraphAsymmErrors *g = new TGraphAsymmErrors(this->h_data);
  for (int iBin = 0; iBin < g->GetN(); ++iBin) {
    int N = g->GetY()[iBin];
    double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
    double U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) ;
    g->SetPointEYlow(iBin, N-L);
    g->SetPointEYhigh(iBin, U-N);
  }
  

  // MAKE STACK; FORMAT AND DRAW EVERYTHING

   
  THStack *hs = new THStack("hs","");
  hs->Add(this->h_ewk);
  hs->Add(this->h_ttV);
  hs->Add(this->h_singlet);
  hs->Add(this->h_ttbar_bbbar);
  hs->Add(this->h_ttbar_ccbar);
  hs->Add(this->h_ttbar_lf);
  
  this->h_bkg_err_1sig->SetFillStyle(3654);
  this->h_bkg_err_1sig->SetFillColor(kBlack);


  TCanvas *c1 = this->setUpCanvas(drawData);

  g->SetLineColor(kBlack);
  g->SetMarkerStyle(20);
  g->SetLineWidth(2);

  this->h_ttH->SetLineWidth(4);
  
  this->h_data->SetMarkerStyle(20);
  this->h_data->SetMarkerSize(0.75);
  this->h_data->SetLineWidth(2);


  string withoutDataLabel = "";
  string logPlotLabel = "";
  string withoutLegendLabel = "";


  if(drawData){

    c1->cd(1);

    this->h_data->SetStats(0);
    this->h_data->GetYaxis()->SetTitle("Events");
    this->h_data->GetYaxis()->SetTitleSize(0.05);
    this->h_data->GetYaxis()->SetTitleOffset(1.);
    if(setLog){
      if(plotYMax>0.)this->h_data->SetMaximum(this->plotYMax);
      else this->h_data->SetMaximum(max_tot*100.);
      this->h_data->SetMinimum(1.5E-1);
    }
    else{
      if(plotYMax>0.)this->h_data->SetMaximum(this->plotYMax);
      else this->h_data->SetMaximum(max_tot*1.4);
    }
    vector<double> range = this->getDefaultXRange(this->plotName);
    if((range[1]-range[0])>0.0001) this->h_data->GetXaxis()->SetRangeUser(range[0],range[1]);

    this->h_data->Draw("pe1");
    hs->Draw("hist same");
    this->h_bkg_err_1sig->Draw("e2same");
    this->h_ttH->Draw("histsame");
    g->Draw("psame");

    legend->AddEntry(this->h_data,std::string("Data (" + str_sum_data + ")").c_str(),"lpe");
    
    c1->GetPad(1)->RedrawAxis();
    if(setLog){
      c1->GetPad(1)->SetLogy(1);
      logPlotLabel = "_log";
    }  
    else{
      c1->GetPad(1)->SetLogy(0);
    }
  }


  else{

    this->h_ttH->SetStats(0);
    this->h_ttH->GetYaxis()->SetTitle("Events");
    this->h_ttH->GetYaxis()->SetTitleSize(0.05);
    this->h_ttH->GetYaxis()->SetTitleOffset(1.);
    if(setLog){
      if(plotYMax>0.)this->h_ttH->SetMaximum(this->plotYMax);
      else this->h_ttH->SetMaximum(max_tot*100.);
      this->h_ttH->SetMinimum(1.5E-1);
    }
    else{
      if(plotYMax>0.)this->h_ttH->SetMaximum(this->plotYMax);
      else this->h_ttH->SetMaximum(max_tot*1.4);
    }
    vector<double> range = this->getDefaultXRange(this->plotName);    
    if((range[1]-range[0])>0.0001) this->h_ttH->GetXaxis()->SetRangeUser(range[0],range[1]);

    this->h_ttH->Draw("hist");
    hs->Draw("hist same");
    this->h_bkg_err_1sig->Draw("e2same");
    this->h_ttH->Draw("histsame");
    withoutDataLabel = "_noData";

    c1->GetPad(0)->RedrawAxis();
    if(setLog){
      c1->GetPad(0)->SetLogy(1);
      logPlotLabel = "_log";
    }  
    else{
      c1->GetPad(0)->SetLogy(0);
    }

  }
  legend->AddEntry(this->h_bkg_err_1sig,std::string("Sum MC (" + str_sum_bkg + ")").c_str(),"f");
  legend->AddEntry(this->h_ttH,std::string("t#bar{t}H125 (" + str_sum_ttH + "#times" + str_scale_H + ")").c_str(),"l");
  
  CMSInfoLatex.Draw();
  
  if(drawLegend)legend->Draw();
  else{
    SELECTIONInfoLatex.Draw();
    withoutLegendLabel = "_noLegend";
  }

  // DRAW RATIO PLOT FOR PLOTS WITH DATA

  if(drawData){
    
    this->getDataMCRatio();

    c1->cd(2);
    gPad->SetTopMargin(1.e-5);
    gPad->SetTickx();
    gPad->Modified();
    
    h_dataMC_ratio->Reset();
    h_dataMC_ratio->Draw("pe1");
    h_dataMC_ratio_1sig->Draw("e2same");
    asymmRatio->Draw("p0same");
    
    TLine* myLine;
    myLine = new TLine(h_data->GetBinLowEdge(1), 1, h_data->GetBinLowEdge(nbins) + h_data->GetBinWidth(nbins), 1);
    myLine->Draw();
    
    
  }
  
  string categorydir = this->getCategoryDir(); 
  TString outputHistoname = this->plotName.c_str();
  c1->Print("Images/"+categorydir+outputHistoname+withoutDataLabel+logPlotLabel+withoutLegendLabel+".pdf");
  
  delete c1;
  
 
}

TCanvas* MyPlot::setUpCanvas(bool withRatio){

  TCanvas *tempCanvas;

  if(withRatio){
    tempCanvas = new TCanvas("tempCanvas", "tempCanvas", 600,700);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    Float_t small = 1.e-5;
    tempCanvas->Divide(1,2,small,small);
    const float padding=1e-5; const float ydivide=0.3;
    tempCanvas->GetPad(1)->SetPad( padding, ydivide + padding , 1-padding, 1-padding);
    tempCanvas->GetPad(2)->SetPad( padding, padding, 1-padding, ydivide-padding);
    tempCanvas->GetPad(1)->SetLeftMargin(.11);
    tempCanvas->GetPad(2)->SetLeftMargin(.11);
    tempCanvas->GetPad(1)->SetRightMargin(.05);
    tempCanvas->GetPad(2)->SetRightMargin(.05);
    tempCanvas->GetPad(1)->SetBottomMargin(.3);
    tempCanvas->GetPad(2)->SetBottomMargin(.3);
    tempCanvas->GetPad(1)->Modified();
    tempCanvas->GetPad(2)->Modified();
    tempCanvas->cd(1);
    gPad->SetBottomMargin(small);
    gPad->Modified();
    tempCanvas->cd(2);
    gPad->SetTopMargin(small);
    gPad->SetTickx();
    gPad->Modified();
    
  }
  else{
    tempCanvas = new TCanvas("tempCanvas","tempCanvas",600,500);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    Float_t small = 1.e-5;
    const float padding=1e-5; const float ydivide=0.0;
    tempCanvas->GetPad(0)->SetPad( padding, ydivide + padding , 1-padding, 1-padding);
    tempCanvas->GetPad(0)->SetLeftMargin(.11);
    tempCanvas->GetPad(0)->SetRightMargin(.05);
    tempCanvas->GetPad(0)->SetBottomMargin(.11);
    tempCanvas->GetPad(0)->SetTickx();
    tempCanvas->GetPad(0)->Modified();
    TGaxis::SetMaxDigits(3);  
  }

  return tempCanvas;

}

void MyPlot::setPlotMax(double newMax){

  plotYMax = newMax;

}

void MyPlot::getDataMCRatio(){
  
  TH1D* data_forbinning = (TH1D*)h_sum_bkg_sys[0]->Clone();
  int nbins = data_forbinning->GetNbinsX();
  
  TH1D* myRatio = new TH1D("ratio", "", nbins,data_forbinning->GetBinLowEdge(1), data_forbinning->GetBinLowEdge(nbins) + data_forbinning->GetBinWidth(nbins));
  TH1D* myRatio_1sig = new TH1D("ratio_1sig", "", nbins, data_forbinning->GetBinLowEdge(1), data_forbinning->GetBinLowEdge(nbins) + data_forbinning->GetBinWidth(nbins));
  
  
  myRatio->Sumw2();
  
  double ratioMax = 2.3;
  double ratioMin = 0.0;

  
  if( (this->plotName.find("h_numTag")!=std::string::npos) ){
    myRatio->GetXaxis()->SetBinLabel(1,"2");
    myRatio->GetXaxis()->SetBinLabel(2,"3");
    myRatio->GetXaxis()->SetBinLabel(3,"4");
    
    if( (this->plotName.find("full")!=std::string::npos) ){
            myRatio->GetXaxis()->SetBinLabel(4,"5");
    }
    
    myRatio->GetXaxis()->SetLabelSize(0.15);
  }
  if( (this->plotName.find("h_numJet")!=std::string::npos) ){
    myRatio->GetXaxis()->SetBinLabel(1,"4");
    myRatio->GetXaxis()->SetBinLabel(2,"5");
    myRatio->GetXaxis()->SetBinLabel(3,"6");
    myRatio->GetXaxis()->SetBinLabel(4,"7");
    myRatio->GetXaxis()->SetBinLabel(5,"8");
    myRatio->GetXaxis()->SetBinLabel(6,"9");
    myRatio->GetXaxis()->SetBinLabel(7,"10"); 
    myRatio->GetXaxis()->SetLabelSize(0.15);
  }
   

  myRatio->SetStats(0);
  myRatio->Sumw2();
  myRatio->SetLineColor(kBlack);
  myRatio->SetMarkerColor(kBlack);

  for( int bin=0; bin<nbins; bin++ ){
    double bkg  = h_bkg_err_1sig->GetBinContent(bin+1);
    double bkg_1sig  = h_bkg_err_1sig->GetBinError(bin+1);
    double data = h_data->GetBinContent(bin+1);
    double ratio = ( bkg>0. ) ? data/bkg : 0.;
    double ratio_err = ( bkg>0. ) ? sqrt(data)/bkg : 0.;
    
    double bkg_noshift  = TotalBKG_rateUnc->GetBinContent(bin+1);

    double up_err = bkg + bkg_1sig;
    double down_err = bkg - bkg_1sig;
    
    if( bkg>0. && bkg_noshift>0. ){
      myRatio->SetBinContent(bin+1,ratio);
      myRatio->SetBinError(bin+1,ratio_err);
            
      up_err *= 1./bkg_noshift;
      down_err *= 1./bkg_noshift;
      
      double new_ave = 0.5 * ( up_err + down_err );
      
      myRatio_1sig->SetBinContent(bin+1,new_ave);
      myRatio_1sig->SetBinError(bin+1,up_err - new_ave);
      
    }
    
    if( (ratio>ratioMax) && ((ratio - ratio_err)<ratioMax) ){
      double minner = ratio - ratio_err;
      myRatio->SetBinContent(bin+1,ratioMax-0.0001);
      myRatio->SetBinError(bin+1,ratioMax-0.0001-minner);
    }
    
  }
  
  TGraphAsymmErrors *gRatio = new TGraphAsymmErrors(this->h_data);
  
  for (int iBin = 0; iBin < gRatio->GetN(); ++iBin) {
    double xPoint = h_data->GetBinCenter(iBin+1);
    double xWidth = 0.5*h_data->GetBinWidth(iBin+1);
    
    double yG = gRatio->GetY()[iBin];
    double yG_low  = gRatio->GetEYlow()[iBin];
    double yG_high = gRatio->GetEYhigh()[iBin];
    double yData = h_data->GetBinContent(iBin+1);
    
    double yBkg = TotalBKG_rateUnc->GetBinContent(iBin+1);
    
    double yG_ratio = ( yBkg>0. ) ? yG/yBkg : 0.;
    double yG_ratio_low = ( yBkg>0. ) ? yG_low/yBkg : 0.;
    double yG_ratio_high = ( yBkg>0. ) ? yG_high/yBkg : 0.;
    
    if( yData>0 ){
      gRatio->SetPoint(iBin, xPoint, yG_ratio );
      gRatio->SetPointEYlow(iBin, yG_ratio_low);
      gRatio->SetPointEYhigh(iBin, yG_ratio_high);
      
      gRatio->SetPointEXlow(iBin, xWidth);
      gRatio->SetPointEXhigh(iBin, xWidth);
      
      //if( (yG_ratio>ratioMax) && ((yG_ratio - yG_ratio_low)<ratioMax) ){
      //  double minner = yG_ratio_low - (yG_ratio - ratioMax-0.0001);
      //  asymmRatio->SetPoint(iBin, xPoint, ratioMax-0.0001 );
      //  asymmRatio->SetPointEYlow(iBin, minner);
      //	}
    }
    
  }
  gRatio->SetLineColor(kBlack);
  gRatio->SetLineWidth(2);
  
  myRatio->SetMinimum(ratioMin);
  myRatio->SetMaximum(ratioMax);
  myRatio->GetYaxis()->SetNdivisions(50000+404);
  myRatio->GetYaxis()->SetLabelSize(0.1); //make y label bigger
  myRatio->GetYaxis()->SetTitle("Data/MC");
  myRatio->GetYaxis()->SetTitleSize(0.1);
  myRatio->GetYaxis()->SetTitleOffset(.45);
  myRatio->GetYaxis()->CenterTitle();

  //myRatio->GetXaxis()->SetLabelSize(0.1); //make x label bigger
  myRatio->GetXaxis()->SetTitle(h_data->GetXaxis()->GetTitle()); //make y label bigger
  //if( (temp.find("numJet")!=std::string::npos) ) myRatio->GetXaxis()->SetTitle("Number of jets");
  //if( (temp.find("numTag")!=std::string::npos) ) myRatio->GetXaxis()->SetTitle("Number of tags");
  myRatio->GetXaxis()->SetLabelSize(0.12);
  myRatio->GetXaxis()->SetLabelOffset(0.02);
  myRatio->GetXaxis()->SetTitleOffset(0.90);
  myRatio->GetXaxis()->SetTitleSize(0.14);

  myRatio_1sig->SetMarkerColor(kGreen);
  myRatio_1sig->SetFillColor(kGreen);
  myRatio->SetLineWidth(2);
  
  h_dataMC_ratio = (TH1D*)myRatio->Clone("h_dataMC_ratio");
  h_dataMC_ratio_1sig = (TH1D*)myRatio_1sig->Clone("h_dataMC_ratio_1sig");
  asymmRatio = (TGraphAsymmErrors*)gRatio->Clone("asymmRatio");

  vector<double> range = this->getDefaultXRange(this->plotName);
  if((range[1]-range[0])>0.0001) this->h_dataMC_ratio->GetXaxis()->SetRangeUser(range[0],range[1]);

  delete gRatio;
  delete myRatio;
  delete myRatio_1sig;

}

void MyPlot::getSysShapes(int numSys){
  
  // MAKING THE HATCHED HISTOGRAM FOR BKG ERROR
  TH1D* data_forbinning = (TH1D*)h_sum_bkg_sys[0]->Clone();
  int nbins = data_forbinning->GetNbinsX();

  double sum_mc_err_tot_2 = 0.;

  // Check the bin-by-bin and added sys histos of the background sum has the same bin contents
 
  for( int bin=0; bin<nbins; bin++ ){
    if( fabs(TotalBKG_rateUnc->GetBinContent(bin+1) - h_sum_bkg_sys[0]->GetBinContent(bin+1))>0.0001 ){
            std::cout << " FREAK OUT !!! TotalBKG_rateUnc bin = " << TotalBKG_rateUnc->GetBinContent(bin+1) << ", h_sum_bkg_sys[0] bin = " << h_sum_bkg_sys[0]->GetBinContent(bin+1) << std::endl;
    }
      
    double nom = h_sum_bkg_sys[0]->GetBinContent(bin+1);
    double stat = TotalBKG_rateUnc->GetBinError(bin+1);
    double mcstat = h_bkg_mc_stat->GetBinError(bin+1);
    sum_mc_err_tot_2 += mcstat*mcstat;

    double up2 = 0, down2 = 0;
    double previous_diff = 0;

    //Get the systematics one-by-one

    for( int b=1; b<numSys; b++ ){
      double diff = h_sum_bkg_sys[b]->GetBinContent(bin+1) - nom;
      double diff_2 = diff*diff;
      if( b!=1 && b%2==0 && (diff*previous_diff)>0 ){ // DOWNWARD SYS
	if( fabs(diff)>fabs(previous_diff) ) diff_2 = diff*diff - previous_diff*previous_diff;
      } 
      if( diff>0 ) up2 += diff_2;
      else if( diff<0. ) down2 += diff_2;
      previous_diff = diff;
    }

    up2 += stat*stat + mcstat*mcstat;
    down2 += stat*stat + mcstat*mcstat;
    
    double up_err   = nom + sqrt(up2);
    double down_err = nom - sqrt(down2);
    
    double ave = 0.5 * ( up_err + down_err );
    h_bkg_err_1sig->SetBinContent(bin+1,ave);
    h_bkg_err_1sig->SetBinError(bin+1,up_err-ave);
  }


}


void MyPlot::getTotalBKG_rateUnc(TH1D *hists[]){
  
  TH1D* data_forbinning = (TH1D*)hists[0]->Clone();
  int nbins = data_forbinning->GetNbinsX();

  TH1D* h_bkg_err_temp = new TH1D("h_bkg_err_temp","", nbins, data_forbinning->GetBinLowEdge(1), data_forbinning->GetBinLowEdge(nbins) + data_forbinning->GetBinWidth(nbins) );

  int minbin_nonzero = 10000, maxbin_nonzero = 0;  

  for( int iBin=0; iBin<nbins; iBin++ ){    
    double num_bkg = 0;
    double err_pdf_2 = 0;
    double err_scale_2 = 0;
    
    for( int iSample=1; iSample<Nsample_all; iSample++ ){
      if(allSamples[iSample].process!="signal"){
	double bkg_int = hists[iSample]->GetBinContent(iBin+1);
	num_bkg += bkg_int;
	err_pdf_2   += bkg_int*bkg_int*allSamples[iSample].PDF_unc*allSamples[iSample].PDF_unc;
	err_scale_2 += bkg_int*bkg_int*allSamples[iSample].QCDscale_unc*allSamples[iSample].QCDscale_unc;
      }   
    }

    double lumi_err_2 = num_bkg * num_bkg * lumi_err * lumi_err;
    double total_err_2 = lumi_err_2 + err_pdf_2 + err_scale_2;    
    double total_err = sqrt( total_err_2 );
    
    h_bkg_err_temp->SetBinContent(iBin+1,num_bkg);
    h_bkg_err_temp->SetBinError(iBin+1,total_err);
  }

  this->TotalBKG_rateUnc = (TH1D*)h_bkg_err_temp->Clone("this->TotalBKG_rateUnc");

  delete h_bkg_err_temp;

  return;
}



string MyPlot::getCategoryLabel(){

  string suffix;
  string tempname = this->plotName;

  if( (tempname.find("4j1t")!=std::string::npos) )      suffix = "  4 jets +   1 b-tags";
  else if( (tempname.find("4j2t")!=std::string::npos) ) suffix = "  4 jets +   2 b-tags";
  else if( (tempname.find("4j3t")!=std::string::npos) ) suffix = "  4 jets +   3 b-tags";
  else if( (tempname.find("4j4t")!=std::string::npos) ) suffix = "  4 jets + #geq4 b-tags";
  else if( (tempname.find("5j1t")!=std::string::npos) ) suffix = "  5 jets +   1 b-tags";
  else if( (tempname.find("5j2t")!=std::string::npos) ) suffix = "  5 jets +   2 b-tags";
  else if( (tempname.find("5j3t")!=std::string::npos) ) suffix = "  5 jets +   3 b-tags";
  else if( (tempname.find("5j4t")!=std::string::npos) ) suffix = "  5 jets + #geq4 b-tags";
  else if( (tempname.find("6j1t")!=std::string::npos) ) suffix = "#geq6 jets +   1 b-tags";
  else if( (tempname.find("6j2t")!=std::string::npos) ) suffix = "#geq6 jets +   2 b-tags";
  else if( (tempname.find("6j3t")!=std::string::npos) ) suffix = "#geq6 jets +   3 b-tags";
  else if( (tempname.find("6j4t")!=std::string::npos) ) suffix = "#geq6 jets + #geq4 b-tags";
  else suffix = "#geq4 jets + #geq2 b-tags";

  return suffix;
}  

string MyPlot::getCategoryDir(){

  string dir;
  string tempname = this->plotName;

  if( (tempname.find("4j1t")!=std::string::npos) )      dir = "4j1t/";
  else if( (tempname.find("4j2t")!=std::string::npos) ) dir = "4j2t/";
  else if( (tempname.find("4j3t")!=std::string::npos) ) dir = "4j3t/";
  else if( (tempname.find("4j4t")!=std::string::npos) ) dir = "4j4t/";
  else if( (tempname.find("5j1t")!=std::string::npos) ) dir = "5j1t/";
  else if( (tempname.find("5j2t")!=std::string::npos) ) dir = "5j2t/";
  else if( (tempname.find("5j3t")!=std::string::npos) ) dir = "5j3t/";
  else if( (tempname.find("5j4t")!=std::string::npos) ) dir = "5j4t/";
  else if( (tempname.find("6j1t")!=std::string::npos) ) dir = "6j1t/";
  else if( (tempname.find("6j2t")!=std::string::npos) ) dir = "6j2t/";
  else if( (tempname.find("6j3t")!=std::string::npos) ) dir = "6j3t/";
  else if( (tempname.find("6j4t")!=std::string::npos) ) dir = "6j4t/";
  else dir = "";

  return dir;
}  

int MyPlot::getRebin(string name){
  
  string tempname = name;
  int rebin = 1;
  if( (tempname.find("disc")!=std::string::npos) && (tempname.find("CFMlpANN")!=std::string::npos) ){
    if( (tempname.find("4j1t")!=std::string::npos) )      rebin = 4;
    else if( (tempname.find("4j2t")!=std::string::npos) ) rebin = 5;
    else if( (tempname.find("4j3t")!=std::string::npos) ) rebin = 5;
    else if( (tempname.find("4j4t")!=std::string::npos) ) rebin = 10;
    else if( (tempname.find("5j1t")!=std::string::npos) ) rebin = 4;
    else if( (tempname.find("5j2t")!=std::string::npos) ) rebin = 5;
    else if( (tempname.find("5j3t")!=std::string::npos) ) rebin = 5;
    else if( (tempname.find("5j4t")!=std::string::npos) ) rebin = 10;
    else if( (tempname.find("6j1t")!=std::string::npos) ) rebin = 5;
    else if( (tempname.find("6j2t")!=std::string::npos) ) rebin = 5;
    else if( (tempname.find("6j3t")!=std::string::npos) ) rebin = 5;
    else if( (tempname.find("6j4t")!=std::string::npos) ) rebin = 10;
  }
              
  if( (tempname.find("avg_btag_disc_btags")!=std::string::npos) ) rebin = 10;
              
  if( (tempname.find("avg_btag_disc_non_btags")!=std::string::npos) ) rebin = 5;
  if( (tempname.find("lowest_btag")!=std::string::npos) ) rebin = 25;
  if( (tempname.find("sphericity")!=std::string::npos) ) rebin = 5;
  if( (tempname.find("aplanarity")!=std::string::npos) ) rebin = 2;
              
  if( (tempname.find("h_dev_from")!=std::string::npos) ) rebin = 5;
              
  if( (tempname.find("h_met_pt")!=std::string::npos) ) rebin = 2;
              
  if( (tempname.find("h_ave_mass")!=std::string::npos) ||  (tempname.find("h_closest_")!=std::string::npos) ){
    if( (tempname.find("4j")!=std::string::npos) ) rebin = 2;
    else                                       rebin = 5;
  }
              
  if( (tempname.find("_csv")!=std::string::npos) ) rebin = 5;
              
  if( (tempname.find("h_h0")!=std::string::npos) ||
      (tempname.find("h_h1")!=std::string::npos) ||
      (tempname.find("h_h2")!=std::string::npos) ||
      (tempname.find("h_h3")!=std::string::npos) ||
      (tempname.find("h_h4")!=std::string::npos) )  rebin = 5;
              
  if( ( (tempname.find("leptonic")!=std::string::npos) || (tempname.find("hadronic")!=std::string::npos) ) &&
      (tempname.find("mass")!=std::string::npos) ) rebin = 2;
              
  if( (tempname.find("minChi2_getTopSystem")!=std::string::npos) )  rebin = 10;
  if( (tempname.find("lepT_hadT_deltaR")!=std::string::npos) )      rebin = 2;
              
  if( (tempname.find("M3")!=std::string::npos) )  rebin = 5;
  if( (tempname.find("Mlb")!=std::string::npos) ) rebin = 2;
              
  if( (tempname.find("j3t")!=std::string::npos) ){
    if( (tempname.find("avg_btag_disc_btags")!=std::string::npos) ) rebin = 20;
    if( (tempname.find("invariant_mass")!=std::string::npos) )      rebin = 2;
    if( (tempname.find("avg_btag_disc_btags")!=std::string::npos) ) rebin = 20;
    if( (tempname.find("all_sum_pt")!=std::string::npos) )          rebin = 2;
    if( (tempname.find("min_dR_tag_tag")!=std::string::npos) )      rebin = 2;
    if( (tempname.find("jet_1_pt")!=std::string::npos) )            rebin = 2;
    if( (tempname.find("jet_2_pt")!=std::string::npos) )            rebin = 2;
    if( (tempname.find("jet_3_pt")!=std::string::npos) )            rebin = 2;
    if( (tempname.find("jet_4_pt")!=std::string::npos) )            rebin = 2;
    if( (tempname.find("jet_1_tag_pt")!=std::string::npos) )        rebin = 2;
    if( (tempname.find("jet_2_tag_pt")!=std::string::npos) )        rebin = 2;
    if( (tempname.find("jet_3_tag_pt")!=std::string::npos) )        rebin = 2;
    if( (tempname.find("jet_4_tag_pt")!=std::string::npos) )        rebin = 2;
    if( (tempname.find("jet_1_untag_pt")!=std::string::npos) )      rebin = 2;
    if( (tempname.find("jet_2_untag_pt")!=std::string::npos) )      rebin = 2;
    if( (tempname.find("jet_3_untag_pt")!=std::string::npos) )      rebin = 2;
    if( (tempname.find("jet_4_untag_pt")!=std::string::npos) )      rebin = 2;
    if( (tempname.find("jet_tag_1_pt")!=std::string::npos) )        rebin = 2;
    if( (tempname.find("jet_tag_2_pt")!=std::string::npos) )        rebin = 2;
    if( (tempname.find("jet_tag_3_pt")!=std::string::npos) )        rebin = 2;
    if( (tempname.find("jet_tag_4_pt")!=std::string::npos) )        rebin = 2;
    if( (tempname.find("jet_untag_1_pt")!=std::string::npos) )      rebin = 2;
    if( (tempname.find("jet_untag_2_pt")!=std::string::npos) )      rebin = 2;
    if( (tempname.find("jet_untag_3_pt")!=std::string::npos) )      rebin = 2;
    if( (tempname.find("jet_untag_4_pt")!=std::string::npos) )      rebin = 2;
  }
              
  if( (tempname.find("j4t")!=std::string::npos) ){
    if( (tempname.find("mass")!=std::string::npos) )                rebin = 5;
    if( (tempname.find("invariant")!=std::string::npos) )           rebin = 2;
    if( (tempname.find("all_sum_pt")!=std::string::npos)  )         rebin = 2;
    if( (tempname.find("avg_btag_disc_btags")!=std::string::npos) ) rebin = 20;
    if( (tempname.find("min_dR_tag_tag")!=std::string::npos) )      rebin = 2;
    if( (tempname.find("jet_1_pt")!=std::string::npos) )            rebin = 2;
    if( (tempname.find("jet_2_pt")!=std::string::npos) )            rebin = 2;
    if( (tempname.find("jet_3_pt")!=std::string::npos) )            rebin = 2;
    if( (tempname.find("jet_4_pt")!=std::string::npos) )            rebin = 2;
    if( (tempname.find("jet_1_tag_pt")!=std::string::npos) )        rebin = 2;
    if( (tempname.find("jet_2_tag_pt")!=std::string::npos) )        rebin = 2;
    if( (tempname.find("jet_3_tag_pt")!=std::string::npos) )        rebin = 2;
    if( (tempname.find("jet_4_tag_pt")!=std::string::npos) )        rebin = 2;
    if( (tempname.find("jet_1_untag_pt")!=std::string::npos) )      rebin = 2;
    if( (tempname.find("jet_2_untag_pt")!=std::string::npos) )      rebin = 2;
    if( (tempname.find("jet_3_untag_pt")!=std::string::npos) )      rebin = 2;
    if( (tempname.find("jet_4_untag_pt")!=std::string::npos) )      rebin = 2;
    if( (tempname.find("jet_tag_1_pt")!=std::string::npos) )        rebin = 2;
    if( (tempname.find("jet_tag_2_pt")!=std::string::npos) )        rebin = 2;
    if( (tempname.find("jet_tag_3_pt")!=std::string::npos) )        rebin = 2;
    if( (tempname.find("jet_tag_4_pt")!=std::string::npos) )        rebin = 2;
    if( (tempname.find("jet_untag_1_pt")!=std::string::npos) )      rebin = 2;
    if( (tempname.find("jet_untag_2_pt")!=std::string::npos) )      rebin = 2;
    if( (tempname.find("jet_untag_3_pt")!=std::string::npos) )      rebin = 2;
    if( (tempname.find("jet_untag_4_pt")!=std::string::npos) )      rebin = 2;
  }
              
  if( (tempname.find("chargedHadronEnergyFraction")!=std::string::npos) ) rebin = 5;
  if( (tempname.find("neutralHadronEnergyFraction")!=std::string::npos) ) rebin = 5;
  if( (tempname.find("chargedEmEnergyFraction")!=std::string::npos) ) rebin = 5;
  if( (tempname.find("neutralEmEnergyFraction")!=std::string::npos) ) rebin = 5;
  if( (tempname.find("chargedMultiplicity")!=std::string::npos) ) rebin = 2;
  if( (tempname.find("neutralMultiplicity")!=std::string::npos) ) rebin = 2;
  if( (tempname.find("nconstituents")!=std::string::npos) ) rebin = 2;
              
  if( (tempname.find("best_higgs_mass")!=std::string::npos) ) rebin = 10;
              
  if( (tempname.find("dR2Mean")!=std::string::npos) ) rebin = 5;
  if( (tempname.find("dRMean")!=std::string::npos) )  rebin = 5;
  if( (tempname.find("frac01")!=std::string::npos) )  rebin = 5;
  if( (tempname.find("frac02")!=std::string::npos) )  rebin = 5;
  if( (tempname.find("frac03")!=std::string::npos) )  rebin = 5;
  if( (tempname.find("beta")!=std::string::npos) )  rebin = 5;
  if( (tempname.find("betaStar")!=std::string::npos) )  rebin = 5;
  if( (tempname.find("leadCandDistFromPV")!=std::string::npos) )  rebin = 5;
              
  if( (tempname.find("_maxDEta_")!=std::string::npos) ){
    if( (tempname.find("2t")!=std::string::npos) )  rebin = 2;
    else                                        rebin = 4;
  }
              
  if( (tempname.find("_aveJetEta")!=std::string::npos) || (tempname.find("_aveTagEta")!=std::string::npos) ){
    if( (tempname.find("2t")!=std::string::npos) )  rebin = 2;
    else                                        rebin = 5;
  }
  if( (tempname.find("_pt_eta")!=std::string::npos) ){
    if( (tempname.find("eta2p1to2p4")!=std::string::npos) )  rebin = 2;
    //else                                        rebin = 4;
  }
              
  if( (tempname.find("_NPV")!=std::string::npos) ) rebin = 2;
              
  return rebin;
}


vector<double> MyPlot::getDefaultXRange(string name){
  
  vector<double> range(2);
  double xmin = -1.;
  double xmax = -1.;
  string tempname = name;

  if( (tempname.find("_flavour")!=std::string::npos) ){ xmin = -6.0; xmax = 23.0; cout << "flavour plot" << endl;}
  if( (tempname.find("min_dR_lep_jet")!=std::string::npos) ){ xmin = 0.; xmax = 5.0; }
  if( (tempname.find("_csv")!=std::string::npos) ){ xmin = -0.04; xmax = 1.048; }
  if( (tempname.find("_csv")!=std::string::npos)&&(tempname.find("_tag")!=std::string::npos) ){ xmin = 0.64; xmax = 1.048; }
  if( (tempname.find("avg_btag_disc_btags")!=std::string::npos) ){ xmin = 0.65; xmax = 1.009; }
  if( (tempname.find("avg_btag_disc_non_btags")!=std::string::npos) ){ xmin = -1.01; xmax = 1.01; }
  if( (tempname.find("h_numJet")!=std::string::npos) ){ xmin = 3.5; xmax = 10.5; }
  if( (tempname.find("h_numTag")!=std::string::npos) ){ xmin = 1.5; xmax = 4.5; }
  if( (tempname.find("h_numTag_full")!=std::string::npos) ){ xmin = 1.5; xmax = 5.49; }
  if( (tempname.find("_1_csv")!=std::string::npos) ){ xmin = 0.6; xmax = 1.048; }
  if( (tempname.find("_2_csv")!=std::string::npos) &&  !(tempname.find("j1t")!=std::string::npos) ){ xmin = 0.6; xmax = 1.048; }
  if( (tempname.find("_jet_3_pt")!=std::string::npos) ){ xmin = 0.0; xmax = 239; }
  if( (tempname.find("_jet_4_pt")!=std::string::npos) ){ xmin = 0.0; xmax = 239; }
  if( (tempname.find("_jet_3_tag_pt")!=std::string::npos) ){ xmin = 0.0; xmax = 239; }
  if( (tempname.find("_jet_4_tag_pt")!=std::string::npos) ){ xmin = 0.0; xmax = 239; }
  if( (tempname.find("_jet_3_untag_pt")!=std::string::npos) ){ xmin = 0.0; xmax = 239; }
  if( (tempname.find("_jet_4_untag_pt")!=std::string::npos) ){ xmin = 0.0; xmax = 239; }
  if( (tempname.find("_jet_tag_3_pt")!=std::string::npos) ){ xmin = 0.0; xmax = 239; }
  if( (tempname.find("_jet_tag_4_pt")!=std::string::npos) ){ xmin = 0.0; xmax = 239; }
  if( (tempname.find("_jet_untag_3_pt")!=std::string::npos) ){ xmin = 0.0; xmax = 239; }
  if( (tempname.find("_jet_untag_4_pt")!=std::string::npos) ){ xmin = 0.0; xmax = 239; }
  if( (tempname.find("_sumTop_pT")!=std::string::npos) ){ xmin = 0.0; xmax = 299.9; }
  if( (tempname.find("_maxDEta")!=std::string::npos) ){ xmin = 0.0; xmax = 4.; }
  
  // 6jet2tag
  if( (tempname.find("_h0_6j2t")!=std::string::npos) ){ xmin = 0.1; xmax = 0.5; }
  
  // 4jet3tag
  if( (tempname.find("_dev_from_avg_disc_btags_4j3t")!=std::string::npos) ){ xmin = 0.0; xmax = 0.02385; }
  if( (tempname.find("_avg_btag_disc_btags_4j3t")!=std::string::npos) ){ xmin = 0.71; xmax = 1.01; }
  if( (tempname.find("_jet_3_pt_4j3t")!=std::string::npos) || (tempname.find("_jet_4_pt_4j3t")!=std::string::npos) ){ xmin = 0.0; xmax = 199.9; }
  if( (tempname.find("_jet_2_pt_4j3t")!=std::string::npos) ){ xmin = 0.0; xmax = 299.9; }
  if( (tempname.find("_all_sum_pt_with_met_4j3t")!=std::string::npos) ){ xmin = 0.0; xmax = 1599.9; }
  if( (tempname.find("_jet_2_csv_4j3t")!=std::string::npos) ){ xmin = 0.7; xmax = 0.999; }
  
  // 5jet3tag
  if( (tempname.find("_dev_from_avg_disc_btags_5j3t")!=std::string::npos) ){ xmin = 0.0; xmax = 0.02385; }
  if( (tempname.find("_avg_btag_disc_btags_5j3t")!=std::string::npos) ){ xmin = 0.71; xmax = 1.01; }
  if( (tempname.find("_jet_3_pt_5j3t")!=std::string::npos) || (tempname.find("_jet_4_pt_5j3t")!=std::string::npos) ){ xmin = 0.0; xmax = 199.9; }
  if( (tempname.find("_jet_2_pt_5j3t")!=std::string::npos) ){ xmin = 0.0; xmax = 299.9; }
  if( (tempname.find("_min_dR_tag_tag_5j3t")!=std::string::npos) ){ xmin = 0.0; xmax = 3.999; }
  if( (tempname.find("_all_sum_pt_with_met_5j3t")!=std::string::npos) ){ xmin = 0.0; xmax = 1499.9; }
  if( (tempname.find("_jet_2_csv_5j3t")!=std::string::npos) ){ xmin = 0.7; xmax = 0.999; }
  
  // 6jet3tag
  if( (tempname.find("_dev_from_avg_disc_btags_6j3t")!=std::string::npos) ){ xmin = 0.0; xmax = 0.02385; }
  if( (tempname.find("_ave_dR_tag_tag_6j3t")!=std::string::npos) ){ xmin = 0.21; xmax = 4.399; }
  if( (tempname.find("_avg_btag_disc_btags_6j3t")!=std::string::npos) ){ xmin = 0.75; xmax = 1.01; }
  if( (tempname.find("_jet_3_pt_6j3t")!=std::string::npos) || (tempname.find("_jet_4_pt_6j3t")!=std::string::npos) ){ xmin = 0.0; xmax = 199.9; }
  if( (tempname.find("_jet_2_pt_6j3t")!=std::string::npos) ){ xmin = 0.0; xmax = 299.9; }
  if( (tempname.find("_jet_2_csv_6j3t")!=std::string::npos) ){ xmin = 0.7; xmax = 0.999; }
  
  // 4jet4tag
  if( (tempname.find("_dev_from_avg_disc_btags_4j4t")!=std::string::npos) ){ xmin = 0.0; xmax = 0.02385; }
  if( (tempname.find("_ave_dR_tag_tag_4j4t")!=std::string::npos) ){ xmin = 0.61; xmax = 3.599; }
  if( (tempname.find("_avg_btag_disc_btags_4j4t")!=std::string::npos) ){ xmin = 0.75; xmax = 1.01; }
  if( (tempname.find("_jet_3_pt_4j4t")!=std::string::npos) || (tempname.find("_jet_4_pt_4j4t")!=std::string::npos) ){ xmin = 0.0; xmax = 199.9; }
  if( (tempname.find("_jet_2_pt_4j4t")!=std::string::npos) ){ xmin = 0.0; xmax = 299.9; }
  if( (tempname.find("_jet_1_pt_4j4t")!=std::string::npos) ){ xmin = 0.0; xmax = 299.9; }
  if( (tempname.find("_invariant_mass_of_everything_4j4t")!=std::string::npos) ){ xmin = 0.0; xmax = 1499.9; }
  if( (tempname.find("_all_sum_pt_with_met_4j4t")!=std::string::npos) ){ xmin = 0.0; xmax = 999.9; }
  if( (tempname.find("_jet_2_csv_4j4t")!=std::string::npos) ){ xmin = 0.7; xmax = 0.999; }
  
  // 5jet4tag
  if( (tempname.find("_dev_from_avg_disc_btags_5j4t")!=std::string::npos) ){ xmin = 0.0; xmax = 0.02385; }
  if( (tempname.find("_ave_dR_tag_tag_5j4t")!=std::string::npos) ){ xmin = 0.61; xmax = 3.999; }
  if( (tempname.find("_min_dR_tag_tag_5j4t")!=std::string::npos) ){ xmin = 0.11; xmax = 2.599; }
  if( (tempname.find("_min_dR_lep_jet_5j4t")!=std::string::npos) ){ xmin = 0.02; xmax = 2.599; }
  if( (tempname.find("_avg_btag_disc_btags_5j4t")!=std::string::npos) ){ xmin = 0.75; xmax = 1.01; }
  if( (tempname.find("_jet_3_pt_5j4t")!=std::string::npos) || (tempname.find("_jet_4_pt_5j4t")!=std::string::npos) ){ xmin = 0.0; xmax = 199.9; }
  if( (tempname.find("_jet_2_pt_5j4t")!=std::string::npos) ){ xmin = 0.0; xmax = 299.9; }
  if( (tempname.find("_jet_1_pt_5j4t")!=std::string::npos) ){ xmin = 0.0; xmax = 299.9; }
  if( (tempname.find("_invariant_mass_of_everything_5j4t")!=std::string::npos) ){ xmin = 0.0; xmax = 1499.9; }
  if( (tempname.find("_jet_1_csv_5j4t")!=std::string::npos) ){ xmin = 0.7; xmax = 0.999; }
  if( (tempname.find("_jet_2_csv_5j4t")!=std::string::npos) ){ xmin = 0.7; xmax = 0.999; }
  if( (tempname.find("_all_sum_pt_with_met_5j4t")!=std::string::npos) ){ xmin = 0.0; xmax = 1399.9; }
  
  // 6jet4tag
  if( (tempname.find("_dev_from_avg_disc_btags_6j4t")!=std::string::npos) ){ xmin = 0.0; xmax = 0.02385; }
  if( (tempname.find("_ave_dR_tag_tag_6j4t")!=std::string::npos) ){ xmin = 0.81; xmax = 3.999; }
  if( (tempname.find("_min_dR_tag_tag_6j4t")!=std::string::npos) ){ xmin = 0.11; xmax = 2.999; }
  if( (tempname.find("_min_dR_lep_jet_6j4t")!=std::string::npos) ){ xmin = 0.02; xmax = 2.599; }
  if( (tempname.find("_avg_btag_disc_btags_6j4t")!=std::string::npos) ){ xmin = 0.75; xmax = 1.01; }
  if( (tempname.find("_jet_3_pt_6j4t")!=std::string::npos) || (tempname.find("_jet_4_pt_6j4t")!=std::string::npos) ){ xmin = 0.0; xmax = 199.9; }
  if( (tempname.find("_jet_2_pt_6j4t")!=std::string::npos) ){ xmin = 0.0; xmax = 299.9; }
  if( (tempname.find("_jet_1_pt_6j4t")!=std::string::npos) ){ xmin = 0.0; xmax = 299.9; }
  if( (tempname.find("_invariant_mass_of_everything_6j4t")!=std::string::npos) ){ xmin = 0.0; xmax = 1499.9; }
  if( (tempname.find("_jet_2_csv_6j4t")!=std::string::npos) ){ xmin = 0.7; xmax = 0.999; }
  if( (tempname.find("_closest_tagged_dijet_mass_6j4t")!=std::string::npos) ){ xmin = 0.0; xmax = 349.9; }
   if( (tempname.find("chargedMultiplicity")!=std::string::npos) ){ xmin = 0.0; xmax = 50.; }
  if( (tempname.find("neutralMultiplicity")!=std::string::npos) ){ xmin = 0.0; xmax = 50.; }
  if( (tempname.find("nconstituents")!=std::string::npos) ){ xmin = 0.0; xmax = 100.; }
  
  range[0]=xmin;
  range[1]=xmax;

  return range;
}

void MySample::set_up(int lepSelection, int index){

  this->sampleName = "";
  
  if(index==0) this->sampleName  = "SingleMu_2012ABCD_BEAN_53xOn53x";
  if(index==1) this->sampleName  = "DYJetsToLL_M10To50_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==2) this->sampleName  = "DY1JetsToLL_M50_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==3) this->sampleName  = "DY2JetsToLL_M50_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==4) this->sampleName  = "DY3JetsToLL_M50_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==5) this->sampleName  = "DY4JetsToLL_M50_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==6) this->sampleName  = "W1JetsToLNu_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==7) this->sampleName  = "W2JetsToLNu_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==8) this->sampleName  = "W3JetsToLNu_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==9) this->sampleName  = "W4JetsToLNu_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==10) this->sampleName = "WW_TuneZ2Star_8TeV_pythia6Tauola_Summer12_53xOn53x";
  if(index==11) this->sampleName = "WZ_TuneZ2Star_8TeV_pythia6Tauola_Summer12_53xOn53x";
  if(index==12) this->sampleName = "ZZ_TuneZ2Star_8TeV_pythia6Tauola_Summer12_53xOn53x";
  if(index==13) this->sampleName = "T_s_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==14) this->sampleName = "T_tW_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==15) this->sampleName = "T_t_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==16) this->sampleName = "Tbar_s_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==17) this->sampleName = "Tbar_tW_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==18) this->sampleName = "Tbar_t_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==19) this->sampleName = "TTWJets_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==20) this->sampleName = "TTZJets_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==21) this->sampleName = "TTJetsLF_HadronicMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==22) this->sampleName = "TTJetsLF_SemiLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==23) this->sampleName = "TTJetsLF_FullLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==24) this->sampleName = "TTJetsCC_HadronicMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==25) this->sampleName = "TTJetsCC_SemiLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==26) this->sampleName = "TTJetsCC_FullLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==27) this->sampleName = "TTJetsBB_HadronicMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==28) this->sampleName = "TTJetsBB_SemiLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==29) this->sampleName = "TTJetsBB_FullLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==30) this->sampleName = "TTH_Inclusive_M_125_8TeV_53xOn53x";

  if(index==0){
    if( lepSelection==0 )       this->sampleName  = "SingleMu_2012ABCD_BEAN_53xOn53x";
    else if( lepSelection==1 )  this->sampleName  = "SingleElectron_2012ABCD_BEAN_53xOn53x";
    else                        this->sampleName  = "SingleLepton_2012ABCD_BEAN_53xOn53x";
  }

  if(index==0) this->QCDscale_unc  = 0; //"SingleMu_2012ABCD_BEAN_53xOn53x";
  if(index==1) this->QCDscale_unc  = 34./3048; //"DYJetsToLL_M10To50_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==2) this->QCDscale_unc  = 34./3048; //"DY1JetsToLL_M50_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==3) this->QCDscale_unc  = 34./3048; //"DY2JetsToLL_M50_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==4) this->QCDscale_unc  = 34./3048; //"DY3JetsToLL_M50_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==5) this->QCDscale_unc  = 34./3048; //"DY4JetsToLL_M50_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==6) this->QCDscale_unc  = 407./31314; //"W1JetsToLNu_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==7) this->QCDscale_unc  = 407./31314; //"W2JetsToLNu_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==8) this->QCDscale_unc = 407./31314; //"W3JetsToLNu_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==9) this->QCDscale_unc = 407./31314; //"W4JetsToLNu_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==10) this->QCDscale_unc = 0; //"WW_TuneZ2Star_8TeV_pythia6Tauola_Summer12_53xOn53x";
  if(index==11) this->QCDscale_unc = 0; //"WZ_TuneZ2Star_8TeV_pythia6Tauola_Summer12_53xOn53x";
  if(index==12) this->QCDscale_unc = 0; //"ZZ_TuneZ2Star_8TeV_pythia6Tauola_Summer12_53xOn53x";
  if(index==13) this->QCDscale_unc = 0; //"T_s_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==14) this->QCDscale_unc = 0; //"T_tW_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==15) this->QCDscale_unc = 0; //"T_t_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==16) this->QCDscale_unc = 0; //"Tbar_s_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==17) this->QCDscale_unc = 0; //"Tbar_tW_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==18) this->QCDscale_unc = 0; //"Tbar_t_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==19) this->QCDscale_unc = 0; //"TTWJets_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==20) this->QCDscale_unc = 0; //"TTZJets_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==21) this->QCDscale_unc = 19.5/157.5; //"TTJetsLF_HadronicMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==22) this->QCDscale_unc = 19.5/157.5; //"TTJetsLF_SemiLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==23) this->QCDscale_unc = 19.5/157.5; //"TTJetsLF_FullLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==24) this->QCDscale_unc = 19.5/157.5; //"TTJetsCC_HadronicMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==25) this->QCDscale_unc = 19.5/157.5; //"TTJetsCC_SemiLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==26) this->QCDscale_unc = 19.5/157.5; //"TTJetsCC_FullLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==27) this->QCDscale_unc = 19.5/157.5; //"TTJetsBB_HadronicMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==28) this->QCDscale_unc = 19.5/157.5; //"TTJetsBB_SemiLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==29) this->QCDscale_unc = 19.5/157.5; //"TTJetsBB_FullLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==30) this->QCDscale_unc = 0.; //TTH_Inclusive_M_125_8TeV_53xOn53x";

  if(index==0) this->PDF_unc  = 0; //"SingleMu_2012ABCD_BEAN_53xOn53x";
  if(index==1) this->PDF_unc  = 128./3048; //"DYJetsToLL_M10To50_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==2) this->PDF_unc  = 128./3048; //"DY1JetsToLL_M50_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==3) this->PDF_unc  = 128./3048; //"DY2JetsToLL_M50_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==4) this->PDF_unc  = 128./3048; //"DY3JetsToLL_M50_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==5) this->PDF_unc  = 128./3048; //"DY4JetsToLL_M50_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==6) this->PDF_unc  = 1504./31314; //"W1JetsToLNu_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==7) this->PDF_unc  = 1504./31314; //"W2JetsToLNu_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==8) this->PDF_unc  = 1504./31314; //"W3JetsToLNu_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==9) this->PDF_unc  = 1504./31314; //"W4JetsToLNu_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==10) this->PDF_unc = 1.5/43.;  //"WW_TuneZ2Star_8TeV_pythia6Tauola_Summer12_53xOn53x";
  if(index==11) this->PDF_unc = 0.7/18.2; //"WZ_TuneZ2Star_8TeV_pythia6Tauola_Summer12_53xOn53x";
  if(index==12) this->PDF_unc = 0.15/5.9; //"ZZ_TuneZ2Star_8TeV_pythia6Tauola_Summer12_53xOn53x";
  if(index==13) this->PDF_unc = 0.16/4.21; //"T_s_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==14) this->PDF_unc = 0.8/10.6;  //"T_tW_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==15) this->PDF_unc = 3.0/64.6;  //"T_t_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==16) this->PDF_unc = 0.16/4.21; //"Tbar_s_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==17) this->PDF_unc = 0.8/10.6;  //"Tbar_tW_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==18) this->PDF_unc = 3.0/64.6;  //"Tbar_t_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==19) this->PDF_unc = 0; //"TTWJets_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==20) this->PDF_unc = 0; //"TTZJets_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==21) this->PDF_unc = 14.7/157.5; //"TTJetsLF_HadronicMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==22) this->PDF_unc = 14.7/157.5; //"TTJetsLF_SemiLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==23) this->PDF_unc = 14.7/157.5; //"TTJetsLF_FullLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==24) this->PDF_unc = 14.7/157.5; //"TTJetsCC_HadronicMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==25) this->PDF_unc = 14.7/157.5; //"TTJetsCC_SemiLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==26) this->PDF_unc = 14.7/157.5; //"TTJetsCC_FullLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==27) this->PDF_unc = 14.7/157.5; //"TTJetsBB_HadronicMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==28) this->PDF_unc = 14.7/157.5; //"TTJetsBB_SemiLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==29) this->PDF_unc = 14.7/157.5; //"TTJetsBB_FullLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==30) this->PDF_unc = 0.; //"TTH_Inclusive_M_125_8TeV_53xOn53x";

  if(index==0) this->sampleColor  = kBlack; //"SingleMu_2012ABCD_BEAN_53xOn53x";
  if(index==1) this->sampleColor  = kAzure+2; //"DYJetsToLL_M10To50_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==2) this->sampleColor  = kAzure+2; //"DY1JetsToLL_M50_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==3) this->sampleColor  = kAzure+2; //"DY2JetsToLL_M50_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==4) this->sampleColor  = kAzure+2; //"DY3JetsToLL_M50_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==5) this->sampleColor  = kAzure+2; //"DY4JetsToLL_M50_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==6) this->sampleColor  = kAzure+2; //"W1JetsToLNu_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==7) this->sampleColor  = kAzure+2; //"W2JetsToLNu_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==8) this->sampleColor  = kAzure+2; //"W3JetsToLNu_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==9) this->sampleColor  = kAzure+2; //"W4JetsToLNu_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==10) this->sampleColor = kAzure+2;  //"WW_TuneZ2Star_8TeV_pythia6Tauola_Summer12_53xOn53x";
  if(index==11) this->sampleColor = kAzure+2; //"WZ_TuneZ2Star_8TeV_pythia6Tauola_Summer12_53xOn53x";
  if(index==12) this->sampleColor = kAzure+2; //"ZZ_TuneZ2Star_8TeV_pythia6Tauola_Summer12_53xOn53x";
  if(index==13) this->sampleColor = kMagenta; //"T_s_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==14) this->sampleColor = kMagenta;  //"T_tW_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==15) this->sampleColor = kMagenta;  //"T_t_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==16) this->sampleColor = kMagenta; //"Tbar_s_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==17) this->sampleColor = kMagenta;  //"Tbar_tW_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==18) this->sampleColor = kMagenta;  //"Tbar_t_channel_TuneZ2star_8TeV_powheg_Summer12_53xOn53x";
  if(index==19) this->sampleColor = kBlue-10; //"TTWJets_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==20) this->sampleColor = kBlue-10; //"TTZJets_TuneZ2Star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==21) this->sampleColor = kRed+1; //"TTJetsLF_HadronicMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==22) this->sampleColor = kRed+1; //"TTJetsLF_SemiLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==23) this->sampleColor = kRed+1; //"TTJetsLF_FullLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==24) this->sampleColor = kRed-7; //"TTJetsCC_HadronicMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==25) this->sampleColor = kRed-7; //"TTJetsCC_SemiLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==26) this->sampleColor = kRed-7; //"TTJetsCC_FullLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==27) this->sampleColor = kRed+3; //"TTJetsBB_HadronicMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==28) this->sampleColor = kRed+3; //"TTJetsBB_SemiLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==29) this->sampleColor = kRed+3; //"TTJetsBB_FullLeptMGDecays_TuneZ2star_8TeV_madgraph_Summer12_53xOn53x";
  if(index==30) this->sampleColor = kBlue; //"TTH_Inclusive_M_125_8TeV_53xOn53x";

  if(index==0) this->process = "data";
  if(index==1) this->process = "EWK";
  if(index==2) this->process = "EWK";
  if(index==3) this->process = "EWK";
  if(index==4) this->process = "EWK";
  if(index==5) this->process = "EWK";
  if(index==6) this->process = "EWK";
  if(index==7) this->process = "EWK";
  if(index==8) this->process = "EWK";
  if(index==9) this->process = "EWK";
  if(index==10) this->process = "EWK";
  if(index==11) this->process = "EWK";
  if(index==12) this->process = "EWK";
  if(index==13) this->process = "singleT";
  if(index==14) this->process = "singleT";
  if(index==15) this->process = "singleT";
  if(index==16) this->process = "singleT";
  if(index==17) this->process = "singleT";
  if(index==18) this->process = "singleT";
  if(index==19) this->process = "TTV";
  if(index==20) this->process = "TTV";
  if(index==21) this->process = "TTLF";
  if(index==22) this->process = "TTLF";
  if(index==23) this->process = "TTLF";
  if(index==24) this->process = "TTCC";
  if(index==25) this->process = "TTCC";
  if(index==26) this->process = "TTCC";
  if(index==27) this->process = "TTBB";
  if(index==28) this->process = "TTBB";
  if(index==29) this->process = "TTBB";
  if(index==30) this->process = "signal";

}


void MySample::combine(){


}

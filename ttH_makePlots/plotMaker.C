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
#include "plotMaker.h"


using namespace std;
int main(){

  bool drawAnalysisPlots = true;

  bool drawData = true;
  bool drawLog = true;
  bool drawLegend = true;

  bool noData = false;
  bool noLog = false;
  bool noLegend = false;

  cout << "Set up colours, filenames etc...." << endl;
  for(int i=0; i<Nsample_all;i++) allSamples[i].set_up(lepSel,i);

  std::cout << "Open the files...." << std::endl;
  openTheFiles(lepSel);

  // Choose which categories you would like to make plots for

  std::vector<std::string> cat_labels;
  cat_labels.clear();
  cat_labels.push_back("4j2t");
  cat_labels.push_back("5j2t");
  cat_labels.push_back("6j2t");
  cat_labels.push_back("4j3t");
  cat_labels.push_back("5j3t");
  cat_labels.push_back("6j3t");
  cat_labels.push_back("4j4t");
  cat_labels.push_back("5j4t");
  cat_labels.push_back("6j4t");
  
  int NumCat = int(cat_labels.size());

  // Set up the output directories

  std::string dirprefix = "Images/";
  struct stat st;
  if( stat(dirprefix.c_str(),&st) != 0 )  mkdir(dirprefix.c_str(),0777);

  for( int k=0; k<int(cat_labels.size()); k++ ){
    std::string temp = dirprefix + "/" + cat_labels[k];
    if( stat(temp.c_str(),&st) != 0 )  mkdir(temp.c_str(),0777);
  }
  
  //EXAMPLE: MyPlot *numJets = new MyPlot("h_numJet_ge2tag");
  //EXAMPLE: numJets->setPlotMax(100000.); \\ This is optional
  //EXAMPLE: numJets->drawStack(drawData,setLog,drawLegend);

  if(drawAnalysisPlots){

    MyPlot *numJets = new MyPlot("h_numJet_ge2tag");
    MyPlot *numTags = new MyPlot("h_numTag_full");
    numJets->drawStack(drawData,drawLog,noLegend);
    numTags->drawStack(drawData,drawLog,noLegend);
    
    vector<string> input_vars;
    input_vars.push_back("h_jet_1_pt");
    input_vars.push_back("h_jet_2_pt");
    input_vars.push_back("h_jet_3_pt");
    input_vars.push_back("h_jet_4_pt");
    //input_vars.push_back("h_all_sum_pt_with_met");
    //input_vars.push_back("h_invariant_mass_of_everything");
    //input_vars.push_back("avg_untagged_dijet_mass");
    input_vars.push_back("h_closest_tagged_dijet_mass");
    //input_vars.push_back("h_best_higgs_mass");
    //input_vars.push_back("min_dR_lep_jet");
    //input_vars.push_back("h_min_dR_tag_tag");
    //input_vars.push_back("h_min_dR_lep_jet");
    input_vars.push_back("h_sphericity");
    input_vars.push_back("h_aplanarity");
    input_vars.push_back("h_h0");
    input_vars.push_back("h_h1");
    input_vars.push_back("h_h2");
    input_vars.push_back("h_h3");
    input_vars.push_back("h_avg_btag_disc_btags");
    input_vars.push_back("h_dev_from_avg_disc_btags");
    input_vars.push_back("h_jet_1_csv");
    input_vars.push_back("h_jet_2_csv");
    input_vars.push_back("h_lowest_btag");

    for(int i= 0; i<NumCat;i++){
      for(int k = 0; k<int(input_vars.size()); k++){
	string histname = input_vars[k]+"_"+cat_labels[i];
	MyPlot *plotInCategories = new MyPlot(histname);
	plotInCategories->drawStack(drawData,noLog,drawLegend);
      }
    }

    for(int i= 0; i<NumCat;i++){
      string histname = "h_disc_final10v1_8TeV_CFMlpANN_"+cat_labels[i];
      MyPlot *plotOutCategories = new MyPlot(histname);
      plotOutCategories->drawStack(drawData,noLog,drawLegend);
    } 
  }


  if(drawAnalysisPlots==false){

    MyPlot *numJets = new MyPlot("h_numJet_ge2tag");
    MyPlot *numTags = new MyPlot("h_numTag_full");
   
    numJets->drawStack(drawData,drawLog,noLegend);
    numTags->drawStack(drawData,drawLog,noLegend);

    vector<string> samplesForComparison;
    //samplesForComparison.push_back("data");
    samplesForComparison.push_back("ttbb");
    //samplesForComparison.push_back("ttcc");
    //samplesForComparison.push_back("ttlf");
    //samplesForComparison.push_back("ttV");
    //samplesForComparison.push_back("ewk");
    samplesForComparison.push_back("ttH");
    
    MyPlot *jetcsv_compare = new MyPlot("h_jet_csv_4j2t");
    jetcsv_compare->drawShapeComparison(samplesForComparison);

    vector<string> histoNamesToDraw;
    histoNamesToDraw.push_back("h_jet_csv");
    histoNamesToDraw.push_back("h_tag_csv");
    vector<string> histoNamesToDrawTrue;
    histoNamesToDrawTrue.push_back("h_tag_csvB");
    histoNamesToDrawTrue.push_back("h_jet_csvB");
    histoNamesToDrawTrue.push_back("h_jet_flavour");
    histoNamesToDrawTrue.push_back("h_tag_flavour");
    histoNamesToDrawTrue.push_back("h_untag_flavour");


   vector<MyPlot> histosToDraw;
    
    for(int i= 0; i<NumCat;i++){
      for(int k = 0; k<int(histoNamesToDraw.size()); k++){
	string histname = histoNamesToDraw[k]+"_"+cat_labels[i];
	MyPlot *plotInCategories = new MyPlot(histname);
	plotInCategories->drawStack(drawData,noLog,drawLegend);
      }
    }
  }
    
  closeTheFiles();
 
  PrintMessage();

  return 0;

}


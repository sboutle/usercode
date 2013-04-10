#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH3.h"
#include "TH2F.h"
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TAxis.h"
#include "TMath.h"
#include "TRandom3.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <exception>
#include <cmath> 
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "TGraphAsymmErrors.h"
#include "TVector.h"
#include "TLorentzVector.h"
#include "Math/Interpolator.h"


#ifdef __MAKECINT__
#pragma link C++ class std::vector< TLorentzVector >+; 
#endif


#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


#if !defined(__CINT__) && !defined(__MAKECINT__)

#include "AnalysisCode/LeptonPlusJets/interface/BEANloader.h"
#include "AnalysisCode/LeptonPlusJets/interface/BEANeventVars.h"

#endif


//*****************************************************************************
typedef std::vector< TLorentzVector >          vecTLorentzVector;
typedef std::vector<double>                    vdouble;

double getTopPtWgt( double topPt, int useSys );
int getTopSystem(TLorentzVector leptonV, TLorentzVector metV, vecTLorentzVector jetsV, vdouble btagV,
		 double &minChi, TLorentzVector &hadW, TLorentzVector &lepW, TLorentzVector &hadB, TLorentzVector &lepB, TLorentzVector &hadT, TLorentzVector &lepT);

//*****************************************************************************

void yggdrasil_treeReader_53x_v1( int insample=1, std::string era = "2012_52x", int leptonChoice=0, int maxNentries=-1, int Njobs=1, int jobN=1, bool siteOSU=false, bool isSkim=true, bool longJob=false ) {

  std::cout << "   ===> load the root files! " << std::endl;

  bool isMuon=true;

  sampleInfo sampleLoader = BEAN::loadSamples( insample, era, isMuon, Njobs, jobN, siteOSU, isSkim, false); 
  std::string sampleName = sampleLoader.sampleName;

  std::string sampleType = ( insample>=0 ) ? "mc" : "data";

  std::string str_jobN;
  std::stringstream stream;
  stream << jobN;
  str_jobN = stream.str();

  std::string leptonType;
  if(leptonChoice==0) leptonType = "ele_sel";
  else if(leptonChoice==1) leptonType = "mu_sel";
  else leptonType="lep_sel";

  std::cout << leptonChoice << " == " << leptonType << std::endl;

  //std::string treefilename = "/settebello/users/puigh/samples/53xOn53x/yggdrasil/v7/yggdrasil_treeMaker_" + sampleType + "_" + sampleName + "_job*.root";
  //std::string treefilename = "/data2/ttH/anaTrees_2012_53x_20130323_122853
  std::string treefilename = "/data2/ttH/Trees_withoutCVSreshape4/anaTrees_2012_53x_20130326_120434/yggdrasil_treeMaker_" + sampleType + "_" + sampleName + "_job*.root";
  //std::string histofilename = "/home/gsmith/cms_work/CMSSW_5_2_6/src/Output/HistoFiles/treeReader/yggdrasil_treeReader_53x_newCSV_v14_3rd40_topPtWgt_using_v7_trees_" + sampleType + "_" + sampleName + "_" + leptonType + "_histo_" + str_jobN + ".root";
  std::string histofilename = "/data2/ttH/Trees_withoutCVSreshape4/histos2/yggdrasil_treeReader_53x_newCSV_v14_3rd40_topPtWgt_using_v7_trees_" + sampleType + "_" + sampleName + "_" + leptonType + "_histo_" + str_jobN + ".root";

  std::cout << "  treefilename  = " << treefilename.c_str() << std::endl;
  std::cout << "  histofilename = " << histofilename.c_str() << std::endl;

  TChain *chain = new TChain("worldTree");
  chain->Add(treefilename.c_str());

  //////////////////////////////////////////////////////////////////////////
  ///  Reader information
  //////////////////////////////////////////////////////////////////////////


  std::vector<std::string> weight_dir;
  weight_dir.push_back("42/");
  weight_dir.push_back("52/");
  weight_dir.push_back("62/");
  weight_dir.push_back("43/");
  weight_dir.push_back("53/");
  weight_dir.push_back("63/");
  weight_dir.push_back("44/");
  weight_dir.push_back("54/");
  weight_dir.push_back("64/");

  std::vector<std::string> weight_prefix;
  weight_prefix.push_back("weights/");
  weight_prefix.push_back("weights1/");
  weight_prefix.push_back("weights2/");
  weight_prefix.push_back("weights3/");
  weight_prefix.push_back("weights4/");
  weight_prefix.push_back("weights5/");

  std::vector<std::string> cat_labels;
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


  TMVA::Reader *reader_final10v1_8TeV_CFMLP[NumCat];
  for( int c=0; c<NumCat; c++ ){
    reader_final10v1_8TeV_CFMLP[c] = new TMVA::Reader("!Color:Silent");
  }

  Float_t numJets_float; // due to a bug the reader only takes floats
  Float_t numTags_float;

  Float_t tight_lepton_pt;
  Float_t first_jet_pt;
  Float_t min_dr_tagged_jets;
  Float_t avg_dr_tagged_jets;
  Float_t aplanarity;
  Float_t sphericity;
  Float_t avg_btag_disc_non_btags;
  Float_t MET;
  Float_t second_jet_pt;
  Float_t dr_between_lep_and_closest_jet;
  Float_t h0;
  Float_t avg_btag_disc_btags;
  Float_t dev_from_avg_disc_btags;
  Float_t third_jet_pt;
  Float_t fourth_jet_pt;
  Float_t avg_tagged_dijet_mass;
  Float_t h1;
  Float_t h2;
  Float_t lowest_btag;
  Float_t all_sum_pt_with_met;
  Float_t all_sum_pt_incl_met;
  Float_t all_sum_pt;
  Float_t avg_untagged_dijet_mass;
  Float_t closest_tagged_dijet_mass;
  Float_t h3;
  Float_t h4;
  Float_t first_highest_btag;
  Float_t second_highest_btag;
  Float_t third_highest_btag;
  Float_t fourth_highest_btag;
  Float_t dijet_mass_of_everything;

  Float_t invariant_mass_of_everything;

  Float_t MT_mht;

  Float_t M3_1tag;
  Float_t Mlb;

  Float_t best_higgs_mass;
  Float_t minChi2;
  Float_t dRbb;




  /////////////

  // 4j2t
  reader_final10v1_8TeV_CFMLP[0]->AddVariable( "eve.tight_lepton_pt_[0]", &tight_lepton_pt );
  reader_final10v1_8TeV_CFMLP[0]->AddVariable( "eve.min_dr_tagged_jets_[0]", &min_dr_tagged_jets );
  reader_final10v1_8TeV_CFMLP[0]->AddVariable( "eve.aplanarity_[0]", &aplanarity );
  reader_final10v1_8TeV_CFMLP[0]->AddVariable( "eve.MET_[0]", &MET );
  reader_final10v1_8TeV_CFMLP[0]->AddVariable( "eve.third_jet_pt_[0]", &third_jet_pt );
  reader_final10v1_8TeV_CFMLP[0]->AddVariable( "eve.fourth_jet_pt_[0]", &fourth_jet_pt );
  reader_final10v1_8TeV_CFMLP[0]->AddVariable( "eve.avg_tagged_dijet_mass_[0]", &avg_tagged_dijet_mass );
  reader_final10v1_8TeV_CFMLP[0]->AddVariable( "eve.avg_untagged_dijet_mass_[0]", &avg_untagged_dijet_mass );
  reader_final10v1_8TeV_CFMLP[0]->AddVariable( "eve.invariant_mass_of_everything_[0]", &dijet_mass_of_everything );
  reader_final10v1_8TeV_CFMLP[0]->AddVariable( "sqrt(2*eve.tight_lepton_pt_[0]*eve.MHT_[0]*(1-cos(eve.tight_lepton_phi_[0]-eve.MHT_phi_[0])))", &MT_mht );

  // 5j2t
  reader_final10v1_8TeV_CFMLP[1]->AddVariable( "eve.sphericity_[0]", &sphericity );
  reader_final10v1_8TeV_CFMLP[1]->AddVariable( "eve.MET_[0]", &MET );
  reader_final10v1_8TeV_CFMLP[1]->AddVariable( "eve.dr_between_lep_and_closest_jet_[0]", &dr_between_lep_and_closest_jet );
  reader_final10v1_8TeV_CFMLP[1]->AddVariable( "eve.third_jet_pt_[0]", &third_jet_pt );
  reader_final10v1_8TeV_CFMLP[1]->AddVariable( "eve.fourth_jet_pt_[0]", &fourth_jet_pt );
  reader_final10v1_8TeV_CFMLP[1]->AddVariable( "eve.h1_[0]", &h1 );
  reader_final10v1_8TeV_CFMLP[1]->AddVariable( "eve.all_sum_pt_with_met_[0]", &all_sum_pt_incl_met );
  reader_final10v1_8TeV_CFMLP[1]->AddVariable( "eve.avg_untagged_dijet_mass_[0]", &avg_untagged_dijet_mass );
  reader_final10v1_8TeV_CFMLP[1]->AddVariable( "eve.M3_1tag_[0]", &M3_1tag );
  reader_final10v1_8TeV_CFMLP[1]->AddVariable( "eve.Mlb_[0]", &Mlb );

  // 6j2t
  reader_final10v1_8TeV_CFMLP[2]->AddVariable( "eve.aplanarity_[0]", &aplanarity );
  reader_final10v1_8TeV_CFMLP[2]->AddVariable( "eve.sphericity_[0]", &sphericity );
  reader_final10v1_8TeV_CFMLP[2]->AddVariable( "eve.h0_[0]", &h0 );
  reader_final10v1_8TeV_CFMLP[2]->AddVariable( "eve.avg_btag_disc_btags_[0]", &avg_btag_disc_btags );
  reader_final10v1_8TeV_CFMLP[2]->AddVariable( "eve.third_jet_pt_[0]", &third_jet_pt );
  reader_final10v1_8TeV_CFMLP[2]->AddVariable( "eve.fourth_jet_pt_[0]", &fourth_jet_pt );
  reader_final10v1_8TeV_CFMLP[2]->AddVariable( "eve.h1_[0]", &h1 );
  reader_final10v1_8TeV_CFMLP[2]->AddVariable( "eve.avg_untagged_dijet_mass_[0]", &avg_untagged_dijet_mass );
  reader_final10v1_8TeV_CFMLP[2]->AddVariable( "eve.h3_[0]", &h3 );
  reader_final10v1_8TeV_CFMLP[2]->AddVariable( "eve.invariant_mass_of_everything_[0]", &dijet_mass_of_everything );

  // 4j3t
  reader_final10v1_8TeV_CFMLP[3]->AddVariable( "eve.first_jet_pt_[0]", &first_jet_pt );
  reader_final10v1_8TeV_CFMLP[3]->AddVariable( "eve.second_jet_pt_[0]", &second_jet_pt );
  reader_final10v1_8TeV_CFMLP[3]->AddVariable( "eve.avg_btag_disc_btags_[0]", &avg_btag_disc_btags );
  reader_final10v1_8TeV_CFMLP[3]->AddVariable( "eve.dev_from_avg_disc_btags_[0]", &dev_from_avg_disc_btags );
  reader_final10v1_8TeV_CFMLP[3]->AddVariable( "eve.third_jet_pt_[0]", &third_jet_pt );
  reader_final10v1_8TeV_CFMLP[3]->AddVariable( "eve.fourth_jet_pt_[0]", &fourth_jet_pt );
  reader_final10v1_8TeV_CFMLP[3]->AddVariable( "eve.lowest_btag_[0]", &lowest_btag );
  reader_final10v1_8TeV_CFMLP[3]->AddVariable( "eve.all_sum_pt_with_met_[0]", &all_sum_pt_incl_met );
  reader_final10v1_8TeV_CFMLP[3]->AddVariable( "eve.second_highest_btag_[0]", &second_highest_btag );
  reader_final10v1_8TeV_CFMLP[3]->AddVariable( "eve.invariant_mass_of_everything_[0]", &dijet_mass_of_everything );

  // 5j3t
  reader_final10v1_8TeV_CFMLP[4]->AddVariable( "eve.first_jet_pt_[0]", &first_jet_pt );
  reader_final10v1_8TeV_CFMLP[4]->AddVariable( "eve.min_dr_tagged_jets_[0]", &min_dr_tagged_jets );
  reader_final10v1_8TeV_CFMLP[4]->AddVariable( "eve.second_jet_pt_[0]", &second_jet_pt );
  reader_final10v1_8TeV_CFMLP[4]->AddVariable( "eve.avg_btag_disc_btags_[0]", &avg_btag_disc_btags );
  reader_final10v1_8TeV_CFMLP[4]->AddVariable( "eve.dev_from_avg_disc_btags_[0]", &dev_from_avg_disc_btags );
  reader_final10v1_8TeV_CFMLP[4]->AddVariable( "eve.third_jet_pt_[0]", &third_jet_pt );
  reader_final10v1_8TeV_CFMLP[4]->AddVariable( "eve.fourth_jet_pt_[0]", &fourth_jet_pt );
  reader_final10v1_8TeV_CFMLP[4]->AddVariable( "eve.lowest_btag_[0]", &lowest_btag );
  reader_final10v1_8TeV_CFMLP[4]->AddVariable( "eve.all_sum_pt_with_met_[0]", &all_sum_pt_incl_met );
  reader_final10v1_8TeV_CFMLP[4]->AddVariable( "eve.second_highest_btag_[0]", &second_highest_btag );

  // 6j3t
  reader_final10v1_8TeV_CFMLP[5]->AddVariable( "eve.avg_dr_tagged_jets_[0]", &avg_dr_tagged_jets );
  reader_final10v1_8TeV_CFMLP[5]->AddVariable( "eve.sphericity_[0]", &sphericity );
  reader_final10v1_8TeV_CFMLP[5]->AddVariable( "eve.avg_btag_disc_btags_[0]", &avg_btag_disc_btags );
  reader_final10v1_8TeV_CFMLP[5]->AddVariable( "eve.dev_from_avg_disc_btags_[0]", &dev_from_avg_disc_btags );
  reader_final10v1_8TeV_CFMLP[5]->AddVariable( "eve.h2_[0]", &h2 );
  reader_final10v1_8TeV_CFMLP[5]->AddVariable( "eve.lowest_btag_[0]", &lowest_btag );
  reader_final10v1_8TeV_CFMLP[5]->AddVariable( "eve.avg_untagged_dijet_mass_[0]", &avg_untagged_dijet_mass );
  reader_final10v1_8TeV_CFMLP[5]->AddVariable( "eve.h3_[0]", &h3 );
  reader_final10v1_8TeV_CFMLP[5]->AddVariable( "eve.second_highest_btag_[0]", &second_highest_btag );
  reader_final10v1_8TeV_CFMLP[5]->AddVariable( "eve.invariant_mass_of_everything_[0]", &dijet_mass_of_everything );

  // 4j4t
  reader_final10v1_8TeV_CFMLP[6]->AddVariable( "eve.first_jet_pt_[0]", &first_jet_pt );
  reader_final10v1_8TeV_CFMLP[6]->AddVariable( "eve.h1_[0]", &h1 );
  reader_final10v1_8TeV_CFMLP[6]->AddVariable( "eve.avg_dr_tagged_jets_[0]", &avg_dr_tagged_jets );
  reader_final10v1_8TeV_CFMLP[6]->AddVariable( "eve.aplanarity_[0]", &aplanarity );
  reader_final10v1_8TeV_CFMLP[6]->AddVariable( "eve.avg_btag_disc_btags_[0]", &avg_btag_disc_btags );
  reader_final10v1_8TeV_CFMLP[6]->AddVariable( "eve.dev_from_avg_disc_btags_[0]", &dev_from_avg_disc_btags );
  reader_final10v1_8TeV_CFMLP[6]->AddVariable( "eve.lowest_btag_[0]", &lowest_btag );
  reader_final10v1_8TeV_CFMLP[6]->AddVariable( "eve.all_sum_pt_with_met_[0]", &all_sum_pt_incl_met );
  reader_final10v1_8TeV_CFMLP[6]->AddVariable( "eve.second_highest_btag_[0]", &second_highest_btag );
  reader_final10v1_8TeV_CFMLP[6]->AddVariable( "eve.invariant_mass_of_everything_[0]", &dijet_mass_of_everything );

  // 5j4t
  reader_final10v1_8TeV_CFMLP[7]->AddVariable( "eve.avg_dr_tagged_jets_[0]", &avg_dr_tagged_jets );
  reader_final10v1_8TeV_CFMLP[7]->AddVariable( "eve.dr_between_lep_and_closest_jet_[0]", &dr_between_lep_and_closest_jet );
  reader_final10v1_8TeV_CFMLP[7]->AddVariable( "eve.avg_btag_disc_btags_[0]", &avg_btag_disc_btags );
  reader_final10v1_8TeV_CFMLP[7]->AddVariable( "eve.dev_from_avg_disc_btags_[0]", &dev_from_avg_disc_btags );
  reader_final10v1_8TeV_CFMLP[7]->AddVariable( "eve.third_jet_pt_[0]", &third_jet_pt );
  reader_final10v1_8TeV_CFMLP[7]->AddVariable( "eve.fourth_jet_pt_[0]", &fourth_jet_pt );
  reader_final10v1_8TeV_CFMLP[7]->AddVariable( "eve.lowest_btag_[0]", &lowest_btag );
  reader_final10v1_8TeV_CFMLP[7]->AddVariable( "eve.all_sum_pt_with_met_[0]", &all_sum_pt_incl_met );
  reader_final10v1_8TeV_CFMLP[7]->AddVariable( "eve.first_highest_btag_[0]", &first_highest_btag );
  reader_final10v1_8TeV_CFMLP[7]->AddVariable( "eve.second_highest_btag_[0]", &second_highest_btag );

  // 6j4t
  reader_final10v1_8TeV_CFMLP[8]->AddVariable( "eve.avg_dr_tagged_jets_[0]", &avg_dr_tagged_jets );
  reader_final10v1_8TeV_CFMLP[8]->AddVariable( "eve.sphericity_[0]", &sphericity );
  reader_final10v1_8TeV_CFMLP[8]->AddVariable( "eve.dr_between_lep_and_closest_jet_[0]", &dr_between_lep_and_closest_jet );
  reader_final10v1_8TeV_CFMLP[8]->AddVariable( "eve.avg_btag_disc_btags_[0]", &avg_btag_disc_btags );
  reader_final10v1_8TeV_CFMLP[8]->AddVariable( "eve.second_highest_btag_[0]", &second_highest_btag );
  reader_final10v1_8TeV_CFMLP[8]->AddVariable( "eve.h2_[0]", &h2 );
  reader_final10v1_8TeV_CFMLP[8]->AddVariable( "eve.lowest_btag_[0]", &lowest_btag );
  reader_final10v1_8TeV_CFMLP[8]->AddVariable( "eve.closest_tagged_dijet_mass_[0]", &closest_tagged_dijet_mass );
  reader_final10v1_8TeV_CFMLP[8]->AddVariable( "eve.h3_[0]", &h3 );
  reader_final10v1_8TeV_CFMLP[8]->AddVariable( "eve.invariant_mass_of_everything_[0]", &dijet_mass_of_everything );
  reader_final10v1_8TeV_CFMLP[8]->AddVariable( "eve.best_higgs_mass_[0]", &best_higgs_mass );
  
  

  /////////////



  std::string prefix = "TMVAClassification";

  std::string dir_final10v1_8TeV = "/data2/sboutle/CMSSW_5_3_8_patch1/src/AnalysisCode/LeptonPlusJets/data/tmvaWeights_2012_53Xon53X/nominal_weights_03-15-13/";
  // std::string dir_final10v1_8TeV = "weights/TMVAweights_20121126_8TeV/";
  //std::string dir_final10v1_8TeV = "weights/TMVAweights_20130312_8TeV/";

  for( int c=0; c<NumCat; c++ ){
    reader_final10v1_8TeV_CFMLP[c]->BookMVA("CFMlpANN method", dir_final10v1_8TeV + "weights/" + weight_dir[c] + prefix + "_CFMlpANN.weights.xml" );
  }


  double min_reader_final10v1_8TeV_CFMLP[NumCat];
  double max_reader_final10v1_8TeV_CFMLP[NumCat];


  for( int c=0; c<NumCat; c++ ){

    TFile *f_reader_final10v1_8TeV = new TFile(std::string(dir_final10v1_8TeV + "TMVA_trained_" + weight_dir[c].substr(0,weight_dir[c].size()-1) + ".root").c_str());

    TH1F *rep_hist_reader_final10v1_8TeV_cfmlp = (TH1F*)f_reader_final10v1_8TeV->Get("Method_CFMlpANN/CFMlpANN/MVA_CFMlpANN_S")->Clone();

    min_reader_final10v1_8TeV_CFMLP[c] = rep_hist_reader_final10v1_8TeV_cfmlp->GetBinLowEdge(1);
    max_reader_final10v1_8TeV_CFMLP[c] = rep_hist_reader_final10v1_8TeV_cfmlp->GetXaxis()->GetBinUpEdge(rep_hist_reader_final10v1_8TeV_cfmlp->GetNbinsX());

    delete f_reader_final10v1_8TeV;
  }



  //////////////////////////////////////////////////////////////////////////
  ///  Tree branches/leaves
  //////////////////////////////////////////////////////////////////////////

  eventVars *eve=0;
  chain->SetBranchAddress("eve.", &eve );

  //////////////////////////////////////////////////////////////////////////
  ///  Histogram making
  //////////////////////////////////////////////////////////////////////////


  TFile histofile(histofilename.c_str(),"recreate");
  histofile.cd();


  double NpvMin = 0-0.5;
  double NpvMax = 50+0.5;
  int NpvBins = int( NpvMax - NpvMin + 0.0001 );

  int NjetMin_full = 4;
  int NjetMax_full = 10;
  int NjetBins_full = NjetMax_full - NjetMin_full + 1;

  int NjetMin = 4;
  int NjetMax = 6;
  int NjetBins = NjetMax - NjetMin + 1;
  int NtagMin = 2;
  int NtagMax = 4;
  int NtagBins = NtagMax - NtagMin + 1;

  int NtagMin_full = 2;
  int NtagMax_full = 7;
  int NtagBins_full = NtagMax_full - NtagMin_full + 1;

  double metmax   = 500.;
  double lepPtMax = 350.;
  double jetptmax = 500.;
  double htmax    = 2000.;
  double massmax  = 500.;
  double m3max    = 1000.;
  double massmax_sumTop = 2000.;

  int NmetBins   = int( metmax/10. + 0.0001 );
  int NlepPtBins = int( lepPtMax/10. + 0.0001 );
  int NjetptBins = int( jetptmax/10. + 0.0001 );
  int NhtBins    = int( htmax/50. + 0.0001 );
  int NmassBins  = int( massmax/5. + 0.0001 );
  int Nm3Bins    = int( m3max/5. + 0.0001 );
  int NmassBins_sumTop  = int( massmax/10. + 0.0001 );


  double max_jet_leadCandDistFromPV = 0.2;

  int Nbins_CFMlpANN = 100;

  int NcsvBins = 220;
  int NflavourBins = 51;


  std::vector<std::string> sys_cat_labels;
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


  TH1::SetDefaultSumw2();

  TH1D* h_wgt_lepSF[NumSysCat];
  TH1D* h_wgt_PU[NumSysCat];
  TH1D* h_wgt_bTag[NumSysCat];

  TH1D* h_numPV[NumSysCat];
  TH1D* h_numJet[NumSysCat];
  TH1D* h_numTag[NumSysCat];
  TH1D* h_numTag_full[NumSysCat];
  TH2D* h_numJet_numTag[NumSysCat];
  TH1D* h_numJet_0tag[NumSysCat];
  TH1D* h_numJet_1tag[NumSysCat];
  TH1D* h_numJet_2tag[NumSysCat];
  TH1D* h_numJet_3tag[NumSysCat];
  TH1D* h_numJet_4tag[NumSysCat];
  TH1D* h_numJet_ge2tag[NumSysCat];


  TH1D* h_numPV_cat[NumCat][NumSysCat];

  TH1D* h_mu_pt[NumCat][NumSysCat];
  TH1D* h_mu_phi[NumCat][NumSysCat];
  TH1D* h_mu_eta[NumCat][NumSysCat];

  TH1D* h_ele_pt[NumCat][NumSysCat];
  TH1D* h_ele_phi[NumCat][NumSysCat];
  TH1D* h_ele_eta[NumCat][NumSysCat];

  TH1D* h_lepton_pt[NumCat][NumSysCat];
  TH1D* h_lepton_phi[NumCat][NumSysCat];
  TH1D* h_lepton_eta[NumCat][NumSysCat];


  std::vector<int> numPVbins_lowEdge;
  numPVbins_lowEdge.push_back(0);
  numPVbins_lowEdge.push_back(6);
  numPVbins_lowEdge.push_back(11);
  numPVbins_lowEdge.push_back(16);
  numPVbins_lowEdge.push_back(21);
  numPVbins_lowEdge.push_back(26);

  int NumPV = int(numPVbins_lowEdge.size());

  TH1D* h_jet_pt_NPV[NumPV][NumCat][NumSysCat];
  TH1D* h_jet_1_pt_NPV[NumPV][NumCat][NumSysCat];
  TH1D* h_jet_2_pt_NPV[NumPV][NumCat][NumSysCat];
  TH1D* h_jet_3_pt_NPV[NumPV][NumCat][NumSysCat];
  TH1D* h_jet_4_pt_NPV[NumPV][NumCat][NumSysCat];

  TH1D* h_tag_jet_pt_NPV[NumPV][NumCat][NumSysCat];
  TH1D* h_tag_jet_1_pt_NPV[NumPV][NumCat][NumSysCat];
  TH1D* h_tag_jet_2_pt_NPV[NumPV][NumCat][NumSysCat];
  TH1D* h_tag_jet_3_pt_NPV[NumPV][NumCat][NumSysCat];
  TH1D* h_tag_jet_4_pt_NPV[NumPV][NumCat][NumSysCat];

  TH1D* h_untag_jet_pt_NPV[NumPV][NumCat][NumSysCat];
  TH1D* h_untag_jet_1_pt_NPV[NumPV][NumCat][NumSysCat];
  TH1D* h_untag_jet_2_pt_NPV[NumPV][NumCat][NumSysCat];
  TH1D* h_untag_jet_3_pt_NPV[NumPV][NumCat][NumSysCat];
  TH1D* h_untag_jet_4_pt_NPV[NumPV][NumCat][NumSysCat];


  TH1D* h_jet_eta[NumCat][NumSysCat];
  TH1D* h_jet_csv[NumCat][NumSysCat];
  TH1D* h_jet_flavour[NumCat][NumSysCat];
  TH1D* h_jet_csvB[NumCat][NumSysCat];
  TH1D* h_tag_csv[NumCat][NumSysCat];
  TH1D* h_tag_flavour[NumCat][NumSysCat];
  TH1D* h_tag_csvB[NumCat][NumSysCat];
  TH1D* h_untag_csv[NumCat][NumSysCat];
  TH1D* h_untag_flavour[NumCat][NumSysCat];
  TH1D* h_untag_csvB[NumCat][NumSysCat];

  TH1D* h_jet_1_eta[NumCat][NumSysCat];
  TH1D* h_jet_2_eta[NumCat][NumSysCat];
  TH1D* h_jet_3_eta[NumCat][NumSysCat];
  TH1D* h_jet_4_eta[NumCat][NumSysCat];

  TH1D* h_jet_pt[NumCat][NumSysCat];
  TH1D* h_tag_jet_pt[NumCat][NumSysCat];
  TH1D* h_untag_jet_pt[NumCat][NumSysCat];


  TH1D* h_jet_pt_eta0to1p4[NumCat][NumSysCat];
  TH1D* h_tag_jet_pt_eta0to1p4[NumCat][NumSysCat];
  TH1D* h_untag_jet_pt_eta0to1p4[NumCat][NumSysCat];
  
  TH1D* h_jet_pt_eta1p4to2p4[NumCat][NumSysCat];
  TH1D* h_tag_jet_pt_eta1p4to2p4[NumCat][NumSysCat];
  TH1D* h_untag_jet_pt_eta1p4to2p4[NumCat][NumSysCat];
  


  TH1D* h_jet_pt_CSVL_pass[NumCat][NumSysCat];
  TH1D* h_jet_pt_CSVL_fail[NumCat][NumSysCat];
  TH1D* h_jet_pt_CSVL_all[NumCat][NumSysCat];

  TH1D* h_jet_pt_CSVM_pass[NumCat][NumSysCat];
  TH1D* h_jet_pt_CSVM_fail[NumCat][NumSysCat];
  TH1D* h_jet_pt_CSVM_all[NumCat][NumSysCat];

  TH1D* h_jet_pt_CSVT_pass[NumCat][NumSysCat];
  TH1D* h_jet_pt_CSVT_fail[NumCat][NumSysCat];
  TH1D* h_jet_pt_CSVT_all[NumCat][NumSysCat];



  TH1D* h_jet_1_pt[NumCat][NumSysCat];
  TH1D* h_jet_2_pt[NumCat][NumSysCat];
  TH1D* h_jet_3_pt[NumCat][NumSysCat];
  TH1D* h_jet_4_pt[NumCat][NumSysCat];

  TH1D* h_jet_1_tag_pt[NumCat][NumSysCat];
  TH1D* h_jet_2_tag_pt[NumCat][NumSysCat];
  TH1D* h_jet_3_tag_pt[NumCat][NumSysCat];
  TH1D* h_jet_4_tag_pt[NumCat][NumSysCat];



  TH1D* h_jet_1_pt_eta0to1p4[NumCat][NumSysCat];
  TH1D* h_tag_1_pt_eta0to1p4[NumCat][NumSysCat];
  TH1D* h_untag_1_pt_eta0to1p4[NumCat][NumSysCat];
  
  TH1D* h_jet_1_pt_eta1p4to2p4[NumCat][NumSysCat];
  TH1D* h_tag_1_pt_eta1p4to2p4[NumCat][NumSysCat];
  TH1D* h_untag_1_pt_eta1p4to2p4[NumCat][NumSysCat];
  
  TH1D* h_jet_2_pt_eta0to1p4[NumCat][NumSysCat];
  TH1D* h_tag_2_pt_eta0to1p4[NumCat][NumSysCat];
  TH1D* h_untag_2_pt_eta0to1p4[NumCat][NumSysCat];
  
  TH1D* h_jet_2_pt_eta1p4to2p4[NumCat][NumSysCat];
  TH1D* h_tag_2_pt_eta1p4to2p4[NumCat][NumSysCat];
  TH1D* h_untag_2_pt_eta1p4to2p4[NumCat][NumSysCat];
  
  TH1D* h_jet_3_pt_eta0to1p4[NumCat][NumSysCat];
  TH1D* h_tag_3_pt_eta0to1p4[NumCat][NumSysCat];
  TH1D* h_untag_3_pt_eta0to1p4[NumCat][NumSysCat];
  
  TH1D* h_jet_3_pt_eta1p4to2p4[NumCat][NumSysCat];
  TH1D* h_tag_3_pt_eta1p4to2p4[NumCat][NumSysCat];
  TH1D* h_untag_3_pt_eta1p4to2p4[NumCat][NumSysCat];
  
  
  TH1D* h_jet_4_pt_eta0to1p4[NumCat][NumSysCat];
  TH1D* h_tag_4_pt_eta0to1p4[NumCat][NumSysCat];
  TH1D* h_untag_4_pt_eta0to1p4[NumCat][NumSysCat];
  
  TH1D* h_jet_4_pt_eta1p4to2p4[NumCat][NumSysCat];
  TH1D* h_tag_4_pt_eta1p4to2p4[NumCat][NumSysCat];
  TH1D* h_untag_4_pt_eta1p4to2p4[NumCat][NumSysCat];

  TH1D* h_jet_1_untag_pt[NumCat][NumSysCat];
  TH1D* h_jet_2_untag_pt[NumCat][NumSysCat];
  TH1D* h_jet_3_untag_pt[NumCat][NumSysCat];
  TH1D* h_jet_4_untag_pt[NumCat][NumSysCat];


  TH1D* h_jet_tag_1_pt[NumCat][NumSysCat];
  TH1D* h_jet_tag_2_pt[NumCat][NumSysCat];
  TH1D* h_jet_tag_3_pt[NumCat][NumSysCat];
  TH1D* h_jet_tag_4_pt[NumCat][NumSysCat];

  TH1D* h_jet_untag_1_pt[NumCat][NumSysCat];
  TH1D* h_jet_untag_2_pt[NumCat][NumSysCat];
  TH1D* h_jet_untag_3_pt[NumCat][NumSysCat];
  TH1D* h_jet_untag_4_pt[NumCat][NumSysCat];
  
  TH1D* h_jet_pt30to40_dR2Mean[NumCat][NumSysCat];
  TH1D* h_jet_pt30to40_dRMean[NumCat][NumSysCat];
  TH1D* h_jet_pt30to40_frac01[NumCat][NumSysCat];
  TH1D* h_jet_pt30to40_frac02[NumCat][NumSysCat];
  TH1D* h_jet_pt30to40_frac03[NumCat][NumSysCat];
  TH1D* h_jet_pt30to40_beta[NumCat][NumSysCat];
  TH1D* h_jet_pt30to40_betaStar[NumCat][NumSysCat];
  TH1D* h_jet_pt30to40_leadCandDistFromPV[NumCat][NumSysCat];
  
  
  TH1D* h_jet_1_dR2Mean[NumCat][NumSysCat];
  TH1D* h_jet_2_dR2Mean[NumCat][NumSysCat];
  TH1D* h_jet_3_dR2Mean[NumCat][NumSysCat];
  TH1D* h_jet_4_dR2Mean[NumCat][NumSysCat];
  
  
  TH1D* h_jet_1_dRMean[NumCat][NumSysCat];
  TH1D* h_jet_2_dRMean[NumCat][NumSysCat];
  TH1D* h_jet_3_dRMean[NumCat][NumSysCat];
  TH1D* h_jet_4_dRMean[NumCat][NumSysCat];
  
  
  TH1D* h_jet_1_frac01[NumCat][NumSysCat];
  TH1D* h_jet_2_frac01[NumCat][NumSysCat];
  TH1D* h_jet_3_frac01[NumCat][NumSysCat];
  TH1D* h_jet_4_frac01[NumCat][NumSysCat];
  
  
  TH1D* h_jet_1_frac02[NumCat][NumSysCat];
  TH1D* h_jet_2_frac02[NumCat][NumSysCat];
  TH1D* h_jet_3_frac02[NumCat][NumSysCat];
  TH1D* h_jet_4_frac02[NumCat][NumSysCat];
  
  
  TH1D* h_jet_1_frac03[NumCat][NumSysCat];
  TH1D* h_jet_2_frac03[NumCat][NumSysCat];
  TH1D* h_jet_3_frac03[NumCat][NumSysCat];
  TH1D* h_jet_4_frac03[NumCat][NumSysCat];
  
  
  TH1D* h_jet_1_beta[NumCat][NumSysCat];
  TH1D* h_jet_2_beta[NumCat][NumSysCat];
  TH1D* h_jet_3_beta[NumCat][NumSysCat];
  TH1D* h_jet_4_beta[NumCat][NumSysCat];
  
  
  TH1D* h_jet_1_betaStar[NumCat][NumSysCat];
  TH1D* h_jet_2_betaStar[NumCat][NumSysCat];
  TH1D* h_jet_3_betaStar[NumCat][NumSysCat];
  TH1D* h_jet_4_betaStar[NumCat][NumSysCat];
  
  
  TH1D* h_jet_1_leadCandDistFromPV[NumCat][NumSysCat];
  TH1D* h_jet_2_leadCandDistFromPV[NumCat][NumSysCat];
  TH1D* h_jet_3_leadCandDistFromPV[NumCat][NumSysCat];
  TH1D* h_jet_4_leadCandDistFromPV[NumCat][NumSysCat];
  
  
  TH1D* h_jet_1_csv[NumCat][NumSysCat];
  TH1D* h_jet_2_csv[NumCat][NumSysCat];
  TH1D* h_jet_3_csv[NumCat][NumSysCat];
  TH1D* h_jet_4_csv[NumCat][NumSysCat];

  TH1D* h_maxDEta_jet_aveJetEta[NumCat][NumSysCat];
  TH1D* h_maxDEta_tag_aveJetEta[NumCat][NumSysCat];
  TH1D* h_maxDEta_tag_aveTagEta[NumCat][NumSysCat];
  
  TH1D* h_jet_aveJetEta[NumCat][NumSysCat];
  TH1D* h_jet_aveTagEta[NumCat][NumSysCat];
      
  TH1D* h_met_pt[NumCat][NumSysCat];
  TH1D* h_met_phi[NumCat][NumSysCat];
  TH1D* h_mht_pt[NumCat][NumSysCat];
  TH1D* h_mht_phi[NumCat][NumSysCat];
  TH1D* h_all_sum_pt[NumCat][NumSysCat];
  TH1D* h_all_sum_pt_with_met[NumCat][NumSysCat];
  TH1D* h_invariant_mass_of_everything[NumCat][NumSysCat];

  TH1D* h_HT[NumCat][NumSysCat];

  TH1D* h_transverse_mass_met[NumCat][NumSysCat];
  TH2D* h_transverse_mass_met_vs_met[NumCat][NumSysCat];
  TH1D* h_transverse_mass_mht[NumCat][NumSysCat];
  TH2D* h_transverse_mass_mht_vs_mht[NumCat][NumSysCat];
  
  TH1D* h_HT30[NumCat][NumSysCat];
  
  TH1D* h_min_dR_lep_jet[NumCat][NumSysCat];
  TH1D* h_min_dR_tag_tag[NumCat][NumSysCat];
  TH1D* h_ave_dR_tag_tag[NumCat][NumSysCat];
  
  TH1D* h_ave_mass_tag_tag[NumCat][NumSysCat];
  TH1D* h_ave_mass_untag_untag[NumCat][NumSysCat];
  TH1D* h_ave_mass_untag_untag_reduced[NumCat][NumSysCat];
  
  TH1D* h_M3[NumCat][NumSysCat];
  TH1D* h_M3_1tag[NumCat][NumSysCat];
  TH1D* h_Mlb[NumCat][NumSysCat];
  
  TH1D* h_best_higgs_mass[NumCat][NumSysCat];
  TH1D* h_minChi2[NumCat][NumSysCat];
  TH1D* h_dRbb[NumCat][NumSysCat];
  
  
  TH1D* h_leptonicW_mass[NumCat][NumSysCat];
  TH1D* h_hadronicW_mass[NumCat][NumSysCat];
  TH1D* h_leptonicT_mass[NumCat][NumSysCat];
  TH1D* h_hadronicT_mass[NumCat][NumSysCat];
  
  TH1D* h_hadronicW_pT[NumCat][NumSysCat];
  TH2D* h_hadronicW_pT_mass[NumCat][NumSysCat];
  
  TH1D* h_leptonicT_pT[NumCat][NumSysCat];
  TH2D* h_leptonicT_pT_mass[NumCat][NumSysCat];
  TH1D* h_hadronicT_pT[NumCat][NumSysCat];
  TH2D* h_hadronicT_pT_mass[NumCat][NumSysCat];
  
  TH1D* h_lepT_hadT_DeltaR[NumCat][NumSysCat];
  TH1D* h_lepT_hadT_Angle[NumCat][NumSysCat];
  
  
  TH1D* h_sumTop_pT[NumCat][NumSysCat];
  TH1D* h_sumTop_mass[NumCat][NumSysCat];
  TH2D* h_sumTop_pT_mass[NumCat][NumSysCat];
  
  TH1D* h_minChi2_getTopSystem[NumCat][NumSysCat];
  
  TH2D* h_leptonicW_mass_minChi2_getTopSystem[NumCat][NumSysCat];
  TH2D* h_hadronicW_mass_minChi2_getTopSystem[NumCat][NumSysCat];
  TH2D* h_leptonicT_mass_minChi2_getTopSystem[NumCat][NumSysCat];
  TH2D* h_hadronicT_mass_minChi2_getTopSystem[NumCat][NumSysCat];
  
  
  
  // minChi2 > 30
  TH1D* h_jet_pt_minChi2GT30[NumCat][NumSysCat];
  TH1D* h_jet_1_pt_minChi2GT30[NumCat][NumSysCat];
  
  
  // minChi2 < 30
  TH1D* h_leptonicW_mass_minChi2LT30[NumCat][NumSysCat];
  TH1D* h_hadronicW_mass_minChi2LT30[NumCat][NumSysCat];
  TH1D* h_leptonicT_mass_minChi2LT30[NumCat][NumSysCat];
  TH1D* h_hadronicT_mass_minChi2LT30[NumCat][NumSysCat];
    
  TH1D* h_leptonicT_pT_minChi2LT30[NumCat][NumSysCat];
  TH1D* h_hadronicT_pT_minChi2LT30[NumCat][NumSysCat];
    
  TH1D* h_sumTop_pT_minChi2LT30[NumCat][NumSysCat];
  TH1D* h_sumTop_mass_minChi2LT30[NumCat][NumSysCat];
    
  TH1D* h_jet_pt_minChi2LT30[NumCat][NumSysCat];
  TH1D* h_jet_1_pt_minChi2LT30[NumCat][NumSysCat];
    
  TH1D* h_jet_pt_minChi2LT30_CSVL_pass[NumCat][NumSysCat];
  TH1D* h_jet_pt_minChi2LT30_CSVL_fail[NumCat][NumSysCat];
  TH1D* h_jet_pt_minChi2LT30_CSVL_all[NumCat][NumSysCat];
    
  TH1D* h_jet_pt_minChi2LT30_CSVM_pass[NumCat][NumSysCat];
  TH1D* h_jet_pt_minChi2LT30_CSVM_fail[NumCat][NumSysCat];
  TH1D* h_jet_pt_minChi2LT30_CSVM_all[NumCat][NumSysCat];
    
  TH1D* h_jet_pt_minChi2LT30_CSVT_pass[NumCat][NumSysCat];
  TH1D* h_jet_pt_minChi2LT30_CSVT_fail[NumCat][NumSysCat];
  TH1D* h_jet_pt_minChi2LT30_CSVT_all[NumCat][NumSysCat];
    
    
    // minChi2 < 5
  TH1D* h_leptonicW_mass_minChi2LT5[NumCat][NumSysCat];
  TH1D* h_hadronicW_mass_minChi2LT5[NumCat][NumSysCat];
  TH1D* h_leptonicT_mass_minChi2LT5[NumCat][NumSysCat];
  TH1D* h_hadronicT_mass_minChi2LT5[NumCat][NumSysCat];
    
  TH1D* h_leptonicT_pT_minChi2LT5[NumCat][NumSysCat];
  TH1D* h_hadronicT_pT_minChi2LT5[NumCat][NumSysCat];
    
  TH1D* h_sumTop_pT_minChi2LT5[NumCat][NumSysCat];
  TH1D* h_sumTop_mass_minChi2LT5[NumCat][NumSysCat];
    
  TH1D* h_jet_pt_minChi2LT5[NumCat][NumSysCat];
  TH1D* h_jet_1_pt_minChi2LT5[NumCat][NumSysCat];
    
  TH1D* h_jet_pt_minChi2LT5_CSVL_pass[NumCat][NumSysCat];
  TH1D* h_jet_pt_minChi2LT5_CSVL_fail[NumCat][NumSysCat];
  TH1D* h_jet_pt_minChi2LT5_CSVL_all[NumCat][NumSysCat];
    
  TH1D* h_jet_pt_minChi2LT5_CSVM_pass[NumCat][NumSysCat];
  TH1D* h_jet_pt_minChi2LT5_CSVM_fail[NumCat][NumSysCat];
  TH1D* h_jet_pt_minChi2LT5_CSVM_all[NumCat][NumSysCat];
    
  TH1D* h_jet_pt_minChi2LT5_CSVT_pass[NumCat][NumSysCat];
  TH1D* h_jet_pt_minChi2LT5_CSVT_fail[NumCat][NumSysCat];
  TH1D* h_jet_pt_minChi2LT5_CSVT_all[NumCat][NumSysCat];

  
  
  TH1D* h_aplanarity[NumCat][NumSysCat];
  TH1D* h_sphericity[NumCat][NumSysCat];
  TH1D* h_avg_btag_disc_non_btags[NumCat][NumSysCat];
  TH1D* h_avg_btag_disc_btags[NumCat][NumSysCat];
  TH1D* h_dev_from_avg_disc_btags[NumCat][NumSysCat];
  TH1D* h_lowest_btag[NumCat][NumSysCat];
  TH1D* h_closest_tagged_dijet_mass[NumCat][NumSysCat];

  TH1D* h_h0[NumCat][NumSysCat];
  TH1D* h_h1[NumCat][NumSysCat];
  TH1D* h_h2[NumCat][NumSysCat];
  TH1D* h_h3[NumCat][NumSysCat];
  TH1D* h_h4[NumCat][NumSysCat];


  TH1D* h_disc_final10v1_8TeV_CFMlpANN[NumCat][NumSysCat];
  TH1D* h_disc_final10v1_8TeV_CFMlpANN_minNPV[NumCat][NumSysCat];
  TH1D* h_disc_final10v1_8TeV_CFMlpANN_midNPV[NumCat][NumSysCat];
  TH1D* h_disc_final10v1_8TeV_CFMlpANN_maxNPV[NumCat][NumSysCat];


  TH1D* h_jet_1_chargedHadronEnergyFraction[NumCat][NumSysCat];
  TH1D* h_jet_1_neutralHadronEnergyFraction[NumCat][NumSysCat];
  TH1D* h_jet_1_chargedEmEnergyFraction[NumCat][NumSysCat];
  TH1D* h_jet_1_neutralEmEnergyFraction[NumCat][NumSysCat];
  TH1D* h_jet_1_chargedMultiplicity[NumCat][NumSysCat];
  TH1D* h_jet_1_neutralMultiplicity[NumCat][NumSysCat];
  TH1D* h_jet_1_nconstituents[NumCat][NumSysCat];
    
  TH1D* h_jet_2_chargedHadronEnergyFraction[NumCat][NumSysCat];
  TH1D* h_jet_2_neutralHadronEnergyFraction[NumCat][NumSysCat];
  TH1D* h_jet_2_chargedEmEnergyFraction[NumCat][NumSysCat];
  TH1D* h_jet_2_neutralEmEnergyFraction[NumCat][NumSysCat];
  TH1D* h_jet_2_chargedMultiplicity[NumCat][NumSysCat];
  TH1D* h_jet_2_neutralMultiplicity[NumCat][NumSysCat];
  TH1D* h_jet_2_nconstituents[NumCat][NumSysCat];
    
  TH1D* h_jet_3_chargedHadronEnergyFraction[NumCat][NumSysCat];
  TH1D* h_jet_3_neutralHadronEnergyFraction[NumCat][NumSysCat];
  TH1D* h_jet_3_chargedEmEnergyFraction[NumCat][NumSysCat];
  TH1D* h_jet_3_neutralEmEnergyFraction[NumCat][NumSysCat];
  TH1D* h_jet_3_chargedMultiplicity[NumCat][NumSysCat];
  TH1D* h_jet_3_neutralMultiplicity[NumCat][NumSysCat];
  TH1D* h_jet_3_nconstituents[NumCat][NumSysCat];
    
  TH1D* h_jet_4_chargedHadronEnergyFraction[NumCat][NumSysCat];
  TH1D* h_jet_4_neutralHadronEnergyFraction[NumCat][NumSysCat];
  TH1D* h_jet_4_chargedEmEnergyFraction[NumCat][NumSysCat];
  TH1D* h_jet_4_neutralEmEnergyFraction[NumCat][NumSysCat];
  TH1D* h_jet_4_chargedMultiplicity[NumCat][NumSysCat];
  TH1D* h_jet_4_neutralMultiplicity[NumCat][NumSysCat];
  TH1D* h_jet_4_nconstituents[NumCat][NumSysCat];


  for( int c=0; c<NumCat; c++ ){ 
    double width_final10v1_8TeV_CFMLP = max_reader_final10v1_8TeV_CFMLP[c] - min_reader_final10v1_8TeV_CFMLP[c];

    max_reader_final10v1_8TeV_CFMLP[c] += 0.1*width_final10v1_8TeV_CFMLP;
    min_reader_final10v1_8TeV_CFMLP[c] -= 0.1*width_final10v1_8TeV_CFMLP;
  }



  for( int b=0; b<NumSysCat; b++ ){

    std::string short_suffix = sys_cat_labels[b];

    h_wgt_lepSF[b] = new TH1D((std::string("h_wgt_lepSF" + short_suffix)).c_str(), ";lepton scale factor", 100, 0.5, 1.5 );
    h_wgt_PU[b]    = new TH1D((std::string("h_wgt_PU" + short_suffix)).c_str(), ";PU weight", 200, 0., 2.0 );
    h_wgt_bTag[b]    = new TH1D((std::string("h_wgt_bTag" + short_suffix)).c_str(), ";b-tag weight", 200, 0., 2.0 );

    h_numPV[b]  = new TH1D((std::string("h_numPV" + short_suffix)).c_str(), ";Number of Vertices", NpvBins, NpvMin, NpvMax );
    h_numJet[b] = new TH1D((std::string("h_numJet" + short_suffix)).c_str(), ";Number of Jets", NjetBins_full, NjetMin_full-0.5, NjetMax_full+0.5 );
    h_numTag[b] = new TH1D((std::string("h_numTag" + short_suffix)).c_str(), ";Number of Tags", NtagBins, NtagMin-0.5, NtagMax+0.5 );
    h_numTag_full[b] = new TH1D((std::string("h_numTag_full" + short_suffix)).c_str(), ";Number of Tags", NtagBins_full, NtagMin_full-0.5, NtagMax_full+0.5 );
    h_numJet_numTag[b] = new TH2D((std::string("h_numJet_numTag" + short_suffix)).c_str(), ";Number of Jets;Number of Tags", NjetBins, NjetMin-0.5, NjetMax+0.5, NtagBins, NtagMin-0.5, NtagMax+0.5 );
    /*
    h_numJet_0tag[b] = new TH1D((std::string("h_numJet_0tag" + short_suffix)).c_str(), ";Number of Jets", NjetBins_full, NjetMin_full-0.5, NjetMax_full+0.5 );
    h_numJet_1tag[b] = new TH1D((std::string("h_numJet_1tag" + short_suffix)).c_str(), ";Number of Jets", NjetBins_full, NjetMin_full-0.5, NjetMax_full+0.5 );
    */
    h_numJet_2tag[b] = new TH1D((std::string("h_numJet_2tag" + short_suffix)).c_str(), ";Number of Jets", NjetBins_full, NjetMin_full-0.5, NjetMax_full+0.5 );
    h_numJet_3tag[b] = new TH1D((std::string("h_numJet_3tag" + short_suffix)).c_str(), ";Number of Jets", NjetBins_full, NjetMin_full-0.5, NjetMax_full+0.5 );
    h_numJet_4tag[b] = new TH1D((std::string("h_numJet_4tag" + short_suffix)).c_str(), ";Number of Jets", NjetBins_full, NjetMin_full-0.5, NjetMax_full+0.5 );
    h_numJet_ge2tag[b] = new TH1D((std::string("h_numJet_ge2tag" + short_suffix)).c_str(), ";Number of Jets", NjetBins_full, NjetMin_full-0.5, NjetMax_full+0.5 );



    for( int c=0; c<NumCat; c++ ){ 
      std::string suffix = "_" + cat_labels[c];
      std::string cat_suffix = "_" + cat_labels[c];
      suffix += sys_cat_labels[b];

      if(longJob){
      h_numPV_cat[c][b]  = new TH1D((std::string("h_numPV_cat" + suffix)).c_str(), ";Number of Vertices", NpvBins, NpvMin, NpvMax );
      }

      h_mu_pt[c][b] = new TH1D((std::string("h_mu_pt" + suffix)).c_str(),";#mu p_{T}", NlepPtBins, 0, lepPtMax );
      h_mu_phi[c][b] = new TH1D((std::string("h_mu_phi" + suffix)).c_str(),";#mu #phi", 34, -3.4, 3.4 );
      h_mu_eta[c][b] = new TH1D((std::string("h_mu_eta" + suffix)).c_str(),";#mu #eta", 25, -2.5, 2.5 );

      h_ele_pt[c][b] = new TH1D((std::string("h_ele_pt" + suffix)).c_str(),";ele p_{T}", NlepPtBins, 0, lepPtMax );
      h_ele_phi[c][b] = new TH1D((std::string("h_ele_phi" + suffix)).c_str(),";ele #phi", 34, -3.4, 3.4 );
      h_ele_eta[c][b] = new TH1D((std::string("h_ele_eta" + suffix)).c_str(),";ele #eta", 25, -2.5, 2.5 );

      h_lepton_pt[c][b] = new TH1D((std::string("h_lepton_pt" + suffix)).c_str(),";lepton p_{T}", NlepPtBins, 0, lepPtMax );
      h_lepton_phi[c][b] = new TH1D((std::string("h_lepton_phi" + suffix)).c_str(),";lepton #phi", 34, -3.4, 3.4 );
      h_lepton_eta[c][b] = new TH1D((std::string("h_lepton_eta" + suffix)).c_str(),";lepton #eta", 25, -2.5, 2.5 );

      h_jet_eta[c][b] = new TH1D((std::string("h_jet_eta" + suffix)).c_str(),";jet #eta", 25, -2.5, 2.5 );
      h_jet_csv[c][b] = new TH1D((std::string("h_jet_csv" + suffix)).c_str(),";jet CSV", NcsvBins, -1.1, 1.1 );
      h_jet_flavour[c][b] = new TH1D((std::string("h_jet_flavour" + suffix)).c_str(),";jet Flavour", NflavourBins, -25.5, 25.5 );
      h_jet_csvB[c][b] = new TH1D((std::string("h_jet_csvB" + suffix)).c_str(),";jet CSV B", NcsvBins, -1.1, 1.1 );
      h_tag_flavour[c][b] = new TH1D((std::string("h_tag_flavour" + suffix)).c_str(),";tag Flavour", NflavourBins, -25.5, 25.5 );
      h_tag_csv[c][b] = new TH1D((std::string("h_tag_csv" + suffix)).c_str(),";tag CSV", NcsvBins, -1.1, 1.1 );
      h_tag_csvB[c][b] = new TH1D((std::string("h_tag_csvB" + suffix)).c_str(),";tag CSV B", NcsvBins, -1.1, 1.1 );
      h_untag_flavour[c][b] = new TH1D((std::string("h_untag_flavour" + suffix)).c_str(),";untag Flavour", NflavourBins, -25.5, 25.5 );
      h_untag_csv[c][b] = new TH1D((std::string("h_untag_csv" + suffix)).c_str(),";untag CSV", NcsvBins, -1.1, 1.1 );
      h_untag_csvB[c][b] = new TH1D((std::string("h_untag_csvB" + suffix)).c_str(),";untag CSV B", NcsvBins, -1.1, 1.1 );

      h_jet_1_eta[c][b] = new TH1D((std::string("h_jet_1_eta" + suffix)).c_str(),";jet 1 #eta", 25, -2.5, 2.5 );
      h_jet_2_eta[c][b] = new TH1D((std::string("h_jet_2_eta" + suffix)).c_str(),";jet 2 #eta", 25, -2.5, 2.5 );
      h_jet_3_eta[c][b] = new TH1D((std::string("h_jet_3_eta" + suffix)).c_str(),";jet 3 #eta", 25, -2.5, 2.5 );
      h_jet_4_eta[c][b] = new TH1D((std::string("h_jet_4_eta" + suffix)).c_str(),";jet 4 #eta", 25, -2.5, 2.5 );

      h_jet_pt[c][b] = new TH1D((std::string("h_jet_pt" + suffix)).c_str(),";jet p_{T}", NjetptBins, 0, jetptmax );
      h_tag_jet_pt[c][b] = new TH1D((std::string("h_tag_jet_pt" + suffix)).c_str(),";tagged jet p_{T}", NjetptBins, 0, jetptmax );
      h_untag_jet_pt[c][b] = new TH1D((std::string("h_untag_jet_pt" + suffix)).c_str(),";untagged jet p_{T}", NjetptBins, 0, jetptmax );




      if(longJob){
	h_jet_pt_eta0to1p4[c][b] = new TH1D((std::string("h_jet_pt_eta0to1p4" + suffix)).c_str(),";jet p_{T}", NjetptBins, 0, jetptmax );
	h_tag_jet_pt_eta0to1p4[c][b] = new TH1D((std::string("h_tag_jet_pt_eta0to1p4" + suffix)).c_str(),";tagged jet p_{T}", NjetptBins, 0, jetptmax );
	h_untag_jet_pt_eta0to1p4[c][b] = new TH1D((std::string("h_untag_jet_pt_eta0to1p4" + suffix)).c_str(),";untagged jet p_{T}", NjetptBins, 0, jetptmax );
	
	h_jet_pt_eta1p4to2p4[c][b] = new TH1D((std::string("h_jet_pt_eta1p4to2p4" + suffix)).c_str(),";jet p_{T}", NjetptBins, 0, jetptmax );
	h_tag_jet_pt_eta1p4to2p4[c][b] = new TH1D((std::string("h_tag_jet_pt_eta1p4to2p4" + suffix)).c_str(),";tagged jet p_{T}", NjetptBins, 0, jetptmax );
	h_untag_jet_pt_eta1p4to2p4[c][b] = new TH1D((std::string("h_untag_jet_pt_eta1p4to2p4" + suffix)).c_str(),";untagged jet p_{T}", NjetptBins, 0, jetptmax );
      


	h_jet_pt_CSVL_pass[c][b] = new TH1D((std::string("h_jet_pt_CSVL_pass" + suffix)).c_str(),";jet p_{T}", NjetptBins, 0, jetptmax );
	h_jet_pt_CSVL_fail[c][b] = new TH1D((std::string("h_jet_pt_CSVL_fail" + suffix)).c_str(),";jet p_{T}", NjetptBins, 0, jetptmax );
	h_jet_pt_CSVL_all[c][b]  = new TH1D((std::string("h_jet_pt_CSVL_all" + suffix)).c_str(),";jet p_{T}", NjetptBins, 0, jetptmax );
	
	h_jet_pt_CSVM_pass[c][b] = new TH1D((std::string("h_jet_pt_CSVM_pass" + suffix)).c_str(),";jet p_{T}", NjetptBins, 0, jetptmax );
	h_jet_pt_CSVM_fail[c][b] = new TH1D((std::string("h_jet_pt_CSVM_fail" + suffix)).c_str(),";jet p_{T}", NjetptBins, 0, jetptmax );
	h_jet_pt_CSVM_all[c][b]  = new TH1D((std::string("h_jet_pt_CSVM_all" + suffix)).c_str(),";jet p_{T}", NjetptBins, 0, jetptmax );
	
	h_jet_pt_CSVT_pass[c][b] = new TH1D((std::string("h_jet_pt_CSVT_pass" + suffix)).c_str(),";jet p_{T}", NjetptBins, 0, jetptmax );
	h_jet_pt_CSVT_fail[c][b] = new TH1D((std::string("h_jet_pt_CSVT_fail" + suffix)).c_str(),";jet p_{T}", NjetptBins, 0, jetptmax );
	h_jet_pt_CSVT_all[c][b]  = new TH1D((std::string("h_jet_pt_CSVT_all" + suffix)).c_str(),";jet p_{T}", NjetptBins, 0, jetptmax );
      }


      h_jet_1_pt[c][b] = new TH1D((std::string("h_jet_1_pt" + suffix)).c_str(),";jet 1 p_{T}", NjetptBins, 0, jetptmax );
      h_jet_2_pt[c][b] = new TH1D((std::string("h_jet_2_pt" + suffix)).c_str(),";jet 2 p_{T}", NjetptBins, 0, jetptmax );
      h_jet_3_pt[c][b] = new TH1D((std::string("h_jet_3_pt" + suffix)).c_str(),";jet 3 p_{T}", NjetptBins, 0, jetptmax );
      h_jet_4_pt[c][b] = new TH1D((std::string("h_jet_4_pt" + suffix)).c_str(),";jet 4 p_{T}", NjetptBins, 0, jetptmax );


      h_jet_1_tag_pt[c][b] = new TH1D((std::string("h_jet_1_tag_pt" + suffix)).c_str(),";jet 1 p_{T}", NjetptBins, 0, jetptmax );
      h_jet_2_tag_pt[c][b] = new TH1D((std::string("h_jet_2_tag_pt" + suffix)).c_str(),";jet 2 p_{T}", NjetptBins, 0, jetptmax );
      h_jet_3_tag_pt[c][b] = new TH1D((std::string("h_jet_3_tag_pt" + suffix)).c_str(),";jet 3 p_{T}", NjetptBins, 0, jetptmax );
      h_jet_4_tag_pt[c][b] = new TH1D((std::string("h_jet_4_tag_pt" + suffix)).c_str(),";jet 4 p_{T}", NjetptBins, 0, jetptmax );


      h_jet_1_untag_pt[c][b] = new TH1D((std::string("h_jet_1_untag_pt" + suffix)).c_str(),";jet 1 p_{T}", NjetptBins, 0, jetptmax );
      h_jet_2_untag_pt[c][b] = new TH1D((std::string("h_jet_2_untag_pt" + suffix)).c_str(),";jet 2 p_{T}", NjetptBins, 0, jetptmax );
      h_jet_3_untag_pt[c][b] = new TH1D((std::string("h_jet_3_untag_pt" + suffix)).c_str(),";jet 3 p_{T}", NjetptBins, 0, jetptmax );
      h_jet_4_untag_pt[c][b] = new TH1D((std::string("h_jet_4_untag_pt" + suffix)).c_str(),";jet 4 p_{T}", NjetptBins, 0, jetptmax );


      h_jet_tag_1_pt[c][b] = new TH1D((std::string("h_jet_tag_1_pt" + suffix)).c_str(),";tag jet 1 p_{T}", NjetptBins, 0, jetptmax );
      h_jet_tag_2_pt[c][b] = new TH1D((std::string("h_jet_tag_2_pt" + suffix)).c_str(),";tag jet 2 p_{T}", NjetptBins, 0, jetptmax );
      h_jet_tag_3_pt[c][b] = new TH1D((std::string("h_jet_tag_3_pt" + suffix)).c_str(),";tag jet 3 p_{T}", NjetptBins, 0, jetptmax );
      h_jet_tag_4_pt[c][b] = new TH1D((std::string("h_jet_tag_4_pt" + suffix)).c_str(),";tag jet 4 p_{T}", NjetptBins, 0, jetptmax );


      h_jet_untag_1_pt[c][b] = new TH1D((std::string("h_jet_untag_1_pt" + suffix)).c_str(),";untag jet 1 p_{T}", NjetptBins, 0, jetptmax );
      h_jet_untag_2_pt[c][b] = new TH1D((std::string("h_jet_untag_2_pt" + suffix)).c_str(),";untag jet 2 p_{T}", NjetptBins, 0, jetptmax );
      h_jet_untag_3_pt[c][b] = new TH1D((std::string("h_jet_untag_3_pt" + suffix)).c_str(),";untag jet 3 p_{T}", NjetptBins, 0, jetptmax );
      h_jet_untag_4_pt[c][b] = new TH1D((std::string("h_jet_untag_4_pt" + suffix)).c_str(),";untag jet 4 p_{T}", NjetptBins, 0, jetptmax );

      if(longJob){
	h_jet_1_pt_eta0to1p4[c][b] = new TH1D((std::string("h_jet_1_pt_eta0to1p4" + suffix)).c_str(),";jet 1 p_{T}", NjetptBins, 0, jetptmax );
	h_tag_1_pt_eta0to1p4[c][b] = new TH1D((std::string("h_tag_1_pt_eta0to1p4" + suffix)).c_str(),";jet 1 p_{T}", NjetptBins, 0, jetptmax );
	h_untag_1_pt_eta0to1p4[c][b] = new TH1D((std::string("h_untag_1_pt_eta0to1p4" + suffix)).c_str(),";jet 1 p_{T}", NjetptBins, 0, jetptmax );
	
	h_jet_1_pt_eta1p4to2p4[c][b] = new TH1D((std::string("h_jet_1_pt_eta1p4to2p4" + suffix)).c_str(),";jet 1 p_{T}", NjetptBins, 0, jetptmax );
	h_tag_1_pt_eta1p4to2p4[c][b] = new TH1D((std::string("h_tag_1_pt_eta1p4to2p4" + suffix)).c_str(),";jet 1 p_{T}", NjetptBins, 0, jetptmax );
	h_untag_1_pt_eta1p4to2p4[c][b] = new TH1D((std::string("h_untag_1_pt_eta1p4to2p4" + suffix)).c_str(),";jet 1 p_{T}", NjetptBins, 0, jetptmax );
	
	h_jet_2_pt_eta0to1p4[c][b] = new TH1D((std::string("h_jet_2_pt_eta0to1p4" + suffix)).c_str(),";jet 2 p_{T}", NjetptBins, 0, jetptmax );
	h_tag_2_pt_eta0to1p4[c][b] = new TH1D((std::string("h_tag_2_pt_eta0to1p4" + suffix)).c_str(),";jet 2 p_{T}", NjetptBins, 0, jetptmax );
	h_untag_2_pt_eta0to1p4[c][b] = new TH1D((std::string("h_untag_2_pt_eta0to1p4" + suffix)).c_str(),";jet 2 p_{T}", NjetptBins, 0, jetptmax );
	
	h_jet_2_pt_eta1p4to2p4[c][b] = new TH1D((std::string("h_jet_2_pt_eta1p4to2p4" + suffix)).c_str(),";jet 2 p_{T}", NjetptBins, 0, jetptmax );
	h_tag_2_pt_eta1p4to2p4[c][b] = new TH1D((std::string("h_tag_2_pt_eta1p4to2p4" + suffix)).c_str(),";jet 2 p_{T}", NjetptBins, 0, jetptmax );
	h_untag_2_pt_eta1p4to2p4[c][b] = new TH1D((std::string("h_untag_2_pt_eta1p4to2p4" + suffix)).c_str(),";jet 2 p_{T}", NjetptBins, 0, jetptmax );
	
	h_jet_3_pt_eta0to1p4[c][b] = new TH1D((std::string("h_jet_3_pt_eta0to1p4" + suffix)).c_str(),";jet 3 p_{T}", NjetptBins, 0, jetptmax );
	h_tag_3_pt_eta0to1p4[c][b] = new TH1D((std::string("h_tag_3_pt_eta0to1p4" + suffix)).c_str(),";jet 3 p_{T}", NjetptBins, 0, jetptmax );
	h_untag_3_pt_eta0to1p4[c][b] = new TH1D((std::string("h_untag_3_pt_eta0to1p4" + suffix)).c_str(),";jet 3 p_{T}", NjetptBins, 0, jetptmax );
	
	h_jet_3_pt_eta1p4to2p4[c][b] = new TH1D((std::string("h_jet_3_pt_eta1p4to2p4" + suffix)).c_str(),";jet 3 p_{T}", NjetptBins, 0, jetptmax );
	h_tag_3_pt_eta1p4to2p4[c][b] = new TH1D((std::string("h_tag_3_pt_eta1p4to2p4" + suffix)).c_str(),";jet 3 p_{T}", NjetptBins, 0, jetptmax );
	h_untag_3_pt_eta1p4to2p4[c][b] = new TH1D((std::string("h_untag_3_pt_eta1p4to2p4" + suffix)).c_str(),";jet 3 p_{T}", NjetptBins, 0, jetptmax );
	
	h_jet_4_pt_eta0to1p4[c][b] = new TH1D((std::string("h_jet_4_pt_eta0to1p4" + suffix)).c_str(),";jet 4 p_{T}", NjetptBins, 0, jetptmax );
	h_tag_4_pt_eta0to1p4[c][b] = new TH1D((std::string("h_tag_4_pt_eta0to1p4" + suffix)).c_str(),";jet 4 p_{T}", NjetptBins, 0, jetptmax );
	h_untag_4_pt_eta0to1p4[c][b] = new TH1D((std::string("h_untag_4_pt_eta0to1p4" + suffix)).c_str(),";jet 4 p_{T}", NjetptBins, 0, jetptmax );
	
	h_jet_4_pt_eta1p4to2p4[c][b] = new TH1D((std::string("h_jet_4_pt_eta1p4to2p4" + suffix)).c_str(),";jet 4 p_{T}", NjetptBins, 0, jetptmax );
	h_tag_4_pt_eta1p4to2p4[c][b] = new TH1D((std::string("h_tag_4_pt_eta1p4to2p4" + suffix)).c_str(),";jet 4 p_{T}", NjetptBins, 0, jetptmax );
	h_untag_4_pt_eta1p4to2p4[c][b] = new TH1D((std::string("h_untag_4_pt_eta1p4to2p4" + suffix)).c_str(),";jet 4 p_{T}", NjetptBins, 0, jetptmax );
	
	
	h_jet_pt30to40_dR2Mean[c][b] = new TH1D((std::string("h_jet_pt30to40_dR2Mean" + suffix)).c_str(),";jet dR2Mean", 100, 0, 0.13 );
	h_jet_pt30to40_dRMean[c][b] = new TH1D((std::string("h_jet_pt30to40_dRMean" + suffix)).c_str(),";jet dRMean", 100, 0, 0.35 );
	h_jet_pt30to40_frac01[c][b] = new TH1D((std::string("h_jet_pt30to40_frac01" + suffix)).c_str(),";jet frac01", 100, 0, 1.01 );
	h_jet_pt30to40_frac02[c][b] = new TH1D((std::string("h_jet_pt30to40_frac02" + suffix)).c_str(),";jet frac02", 100, 0, 1.01 );
	h_jet_pt30to40_frac03[c][b] = new TH1D((std::string("h_jet_pt30to40_frac03" + suffix)).c_str(),";jet frac03", 100, 0, 1.01 );
	h_jet_pt30to40_beta[c][b] = new TH1D((std::string("h_jet_pt30to40_beta" + suffix)).c_str(),";jet beta", 100, 0, 1.01 );
	h_jet_pt30to40_betaStar[c][b] = new TH1D((std::string("h_jet_pt30to40_betaStar" + suffix)).c_str(),";jet betaStar", 100, 0, 1.01 );
	h_jet_pt30to40_leadCandDistFromPV[c][b] = new TH1D((std::string("h_jet_pt30to40_leadCandDistFromPV" + suffix)).c_str(),";jet leadCandDistFromPV", 100, 0, max_jet_leadCandDistFromPV );
	
	
	h_jet_1_dR2Mean[c][b] = new TH1D((std::string("h_jet_1_dR2Mean" + suffix)).c_str(),";jet 1 dR2Mean", 100, 0, 0.13 );
	h_jet_2_dR2Mean[c][b] = new TH1D((std::string("h_jet_2_dR2Mean" + suffix)).c_str(),";jet 2 dR2Mean", 100, 0, 0.13 );
	h_jet_3_dR2Mean[c][b] = new TH1D((std::string("h_jet_3_dR2Mean" + suffix)).c_str(),";jet 3 dR2Mean", 100, 0, 0.13 );
	h_jet_4_dR2Mean[c][b] = new TH1D((std::string("h_jet_4_dR2Mean" + suffix)).c_str(),";jet 4 dR2Mean", 100, 0, 0.13 );
	
	
	h_jet_1_dRMean[c][b] = new TH1D((std::string("h_jet_1_dRMean" + suffix)).c_str(),";jet 1 dRMean", 100, 0, 0.35 );
	h_jet_2_dRMean[c][b] = new TH1D((std::string("h_jet_2_dRMean" + suffix)).c_str(),";jet 2 dRMean", 100, 0, 0.35 );
	h_jet_3_dRMean[c][b] = new TH1D((std::string("h_jet_3_dRMean" + suffix)).c_str(),";jet 3 dRMean", 100, 0, 0.35 );
	h_jet_4_dRMean[c][b] = new TH1D((std::string("h_jet_4_dRMean" + suffix)).c_str(),";jet 4 dRMean", 100, 0, 0.35 );
	
	
	h_jet_1_frac01[c][b] = new TH1D((std::string("h_jet_1_frac01" + suffix)).c_str(),";jet 1 frac01", 100, 0, 1.01 );
	h_jet_2_frac01[c][b] = new TH1D((std::string("h_jet_2_frac01" + suffix)).c_str(),";jet 2 frac01", 100, 0, 1.01 );
	h_jet_3_frac01[c][b] = new TH1D((std::string("h_jet_3_frac01" + suffix)).c_str(),";jet 3 frac01", 100, 0, 1.01 );
	h_jet_4_frac01[c][b] = new TH1D((std::string("h_jet_4_frac01" + suffix)).c_str(),";jet 4 frac01", 100, 0, 1.01 );
	
	
	h_jet_1_frac02[c][b] = new TH1D((std::string("h_jet_1_frac02" + suffix)).c_str(),";jet 1 frac02", 100, 0, 1.01 );
	h_jet_2_frac02[c][b] = new TH1D((std::string("h_jet_2_frac02" + suffix)).c_str(),";jet 2 frac02", 100, 0, 1.01 );
	h_jet_3_frac02[c][b] = new TH1D((std::string("h_jet_3_frac02" + suffix)).c_str(),";jet 3 frac02", 100, 0, 1.01 );
	h_jet_4_frac02[c][b] = new TH1D((std::string("h_jet_4_frac02" + suffix)).c_str(),";jet 4 frac02", 100, 0, 1.01 );
	
	
	h_jet_1_frac03[c][b] = new TH1D((std::string("h_jet_1_frac03" + suffix)).c_str(),";jet 1 frac03", 100, 0, 1.01 );
	h_jet_2_frac03[c][b] = new TH1D((std::string("h_jet_2_frac03" + suffix)).c_str(),";jet 2 frac03", 100, 0, 1.01 );
	h_jet_3_frac03[c][b] = new TH1D((std::string("h_jet_3_frac03" + suffix)).c_str(),";jet 3 frac03", 100, 0, 1.01 );
	h_jet_4_frac03[c][b] = new TH1D((std::string("h_jet_4_frac03" + suffix)).c_str(),";jet 4 frac03", 100, 0, 1.01 );
	
	
	h_jet_1_beta[c][b] = new TH1D((std::string("h_jet_1_beta" + suffix)).c_str(),";jet 1 beta", 100, 0, 1.01 );
	h_jet_2_beta[c][b] = new TH1D((std::string("h_jet_2_beta" + suffix)).c_str(),";jet 2 beta", 100, 0, 1.01 );
	h_jet_3_beta[c][b] = new TH1D((std::string("h_jet_3_beta" + suffix)).c_str(),";jet 3 beta", 100, 0, 1.01 );
	h_jet_4_beta[c][b] = new TH1D((std::string("h_jet_4_beta" + suffix)).c_str(),";jet 4 beta", 100, 0, 1.01 );
	
	
	h_jet_1_betaStar[c][b] = new TH1D((std::string("h_jet_1_betaStar" + suffix)).c_str(),";jet 1 betaStar", 100, 0, 1.01 );
	h_jet_2_betaStar[c][b] = new TH1D((std::string("h_jet_2_betaStar" + suffix)).c_str(),";jet 2 betaStar", 100, 0, 1.01 );
	h_jet_3_betaStar[c][b] = new TH1D((std::string("h_jet_3_betaStar" + suffix)).c_str(),";jet 3 betaStar", 100, 0, 1.01 );
	h_jet_4_betaStar[c][b] = new TH1D((std::string("h_jet_4_betaStar" + suffix)).c_str(),";jet 4 betaStar", 100, 0, 1.01 );
	
	
	h_jet_1_leadCandDistFromPV[c][b] = new TH1D((std::string("h_jet_1_leadCandDistFromPV" + suffix)).c_str(),";jet 1 leadCandDistFromPV", 100, 0, max_jet_leadCandDistFromPV );
	h_jet_2_leadCandDistFromPV[c][b] = new TH1D((std::string("h_jet_2_leadCandDistFromPV" + suffix)).c_str(),";jet 2 leadCandDistFromPV", 100, 0, max_jet_leadCandDistFromPV );
	h_jet_3_leadCandDistFromPV[c][b] = new TH1D((std::string("h_jet_3_leadCandDistFromPV" + suffix)).c_str(),";jet 3 leadCandDistFromPV", 100, 0, max_jet_leadCandDistFromPV );
	h_jet_4_leadCandDistFromPV[c][b] = new TH1D((std::string("h_jet_4_leadCandDistFromPV" + suffix)).c_str(),";jet 4 leadCandDistFromPV", 100, 0, max_jet_leadCandDistFromPV );
      }

      h_jet_1_csv[c][b] = new TH1D((std::string("h_jet_1_csv" + suffix)).c_str(),";jet 1 CSV", NcsvBins, -1.1, 1.1 );
      h_jet_2_csv[c][b] = new TH1D((std::string("h_jet_2_csv" + suffix)).c_str(),";jet 2 CSV", NcsvBins, -1.1, 1.1 );
      h_jet_3_csv[c][b] = new TH1D((std::string("h_jet_3_csv" + suffix)).c_str(),";jet 3 CSV", NcsvBins, -1.1, 1.1 );
      h_jet_4_csv[c][b] = new TH1D((std::string("h_jet_4_csv" + suffix)).c_str(),";jet 4 CSV", NcsvBins, -1.1, 1.1 );


      if(longJob){
	h_maxDEta_jet_aveJetEta[c][b] = new TH1D((std::string("h_maxDEta_jet_aveJetEta" + suffix)).c_str(),";max #Delta#eta(jet,ave jet #eta)", 100, 0., 6. );
	h_maxDEta_tag_aveJetEta[c][b] = new TH1D((std::string("h_maxDEta_tag_aveJetEta" + suffix)).c_str(),";max #Delta#eta(tag,ave jet #eta)", 100, 0., 6. );
	h_maxDEta_tag_aveTagEta[c][b] = new TH1D((std::string("h_maxDEta_tag_aveTagEta" + suffix)).c_str(),";max #Delta#eta(tag,ave tag #eta)", 100, 0., 6. );
      	
	
	h_jet_aveJetEta[c][b] = new TH1D((std::string("h_jet_aveJetEta" + suffix)).c_str(),";ave jet #eta", 100, -2.5, 2.5 );
	h_jet_aveTagEta[c][b] = new TH1D((std::string("h_jet_aveTagEta" + suffix)).c_str(),";ave tag #eta", 100, -2.5, 2.5 );
	
      }
      
      h_met_pt[c][b]  = new TH1D((std::string("h_met_pt" + suffix)).c_str(),";MET p_{T}", NmetBins, 0, metmax );
      h_met_phi[c][b] = new TH1D((std::string("h_met_phi" + suffix)).c_str(),";MET #phi", 34, -3.4, 3.4 );
      h_mht_pt[c][b]  = new TH1D((std::string("h_mht_pt" + suffix)).c_str(),";MHT p_{T}", NmetBins, 0, metmax );
      h_mht_phi[c][b] = new TH1D((std::string("h_mht_phi" + suffix)).c_str(),";MHT #phi", 34, -3.4, 3.4 );
      h_all_sum_pt[c][b] = new TH1D((std::string("h_all_sum_pt" + suffix)).c_str(),";sum p_{T} (jets,#mu,ele)", NhtBins, 0, htmax );
      h_all_sum_pt_with_met[c][b] = new TH1D((std::string("h_all_sum_pt_with_met" + suffix)).c_str(),";sum p_{T} (jets,lepton,MET)", NhtBins, 0, htmax );
      h_invariant_mass_of_everything[c][b] = new TH1D((std::string("h_invariant_mass_of_everything" + suffix)).c_str(),";M_{inv} (jets,lepton,MET)", NhtBins, 0, htmax );
      h_HT[c][b] = new TH1D((std::string("h_HT" + suffix)).c_str(),";H_{T} (jets)", NhtBins, 0, htmax );



      if(longJob){

	h_transverse_mass_met[c][b] = new TH1D((std::string("h_transverse_mass_met" + suffix)).c_str(),";M_{T} (lepton,MET)", NmetBins, 0, metmax );
	h_transverse_mass_met_vs_met[c][b] = new TH2D((std::string("h_transverse_mass_met_vs_met" + suffix)).c_str(),";M_{T} (lepton,MET);MET", NmetBins, 0, metmax, NmetBins, 0, metmax );
	h_transverse_mass_mht[c][b] = new TH1D((std::string("h_transverse_mass_mht" + suffix)).c_str(),";M_{T} (lepton,MHT)", NmetBins, 0, metmax );
	h_transverse_mass_mht_vs_mht[c][b] = new TH2D((std::string("h_transverse_mass_mht_vs_mht" + suffix)).c_str(),";M_{T} (lepton,MHT);MHT", NmetBins, 0, metmax, NmetBins, 0, metmax );
	
     	h_HT30[c][b] = new TH1D((std::string("h_HT30" + suffix)).c_str(),";H_{T} (jets)", NhtBins, 0, htmax );
	
	h_min_dR_lep_jet[c][b] = new TH1D((std::string("h_min_dR_lep_jet" + suffix)).c_str(),";min #DeltaR(lep,jet)", 30, 0, 6. );
	h_min_dR_tag_tag[c][b] = new TH1D((std::string("h_min_dR_tag_tag" + suffix)).c_str(),";min #DeltaR(tag,tag)", 30, 0, 6. );
	h_ave_dR_tag_tag[c][b] = new TH1D((std::string("h_ave_dR_tag_tag" + suffix)).c_str(),";ave #DeltaR(tag,tag)", 30, 0, 6. );
	
	h_ave_mass_tag_tag[c][b] = new TH1D((std::string("h_ave_mass_tag_tag" + suffix)).c_str(),";ave mass(tag,tag)", NmassBins, 0, massmax );
	h_ave_mass_untag_untag[c][b] = new TH1D((std::string("h_ave_mass_untag_untag" + suffix)).c_str(),";ave mass(untag,untag)", NmassBins, 0, massmax );
	h_ave_mass_untag_untag_reduced[c][b] = new TH1D((std::string("h_ave_mass_untag_untag_reduced" + suffix)).c_str(),";ave mass(untag,untag)", 200, 0, 200. );
	h_M3[c][b] = new TH1D((std::string("h_M3" + suffix)).c_str(),";M3", Nm3Bins, 0, m3max );
	h_M3_1tag[c][b] = new TH1D((std::string("h_M3_1tag" + suffix)).c_str(),";M3 (1 tag)", Nm3Bins, 0, m3max );
	h_Mlb[c][b] = new TH1D((std::string("h_Mlb" + suffix)).c_str(),";M(lep,tag)", NmassBins, 0, massmax );
	
	h_best_higgs_mass[c][b] = new TH1D((std::string("h_best_higgs_mass" + suffix)).c_str(),";best Higgs mass", NmassBins, 0, massmax );
	h_minChi2[c][b] = new TH1D((std::string("h_minChi2" + suffix)).c_str(),";minChi2", 500, 0, 500 );
	h_dRbb[c][b] = new TH1D((std::string("h_dRbb" + suffix)).c_str(),";best #DeltaR(b,b)", 30, 0, 6. );
	
	
	h_leptonicW_mass[c][b] = new TH1D((std::string("h_leptonicW_mass" + suffix)).c_str(),";lep W mass (minChi2)", NmassBins, 0, massmax );
	h_hadronicW_mass[c][b] = new TH1D((std::string("h_hadronicW_mass" + suffix)).c_str(),";had W mass (minChi2)", NmassBins, 0, massmax );
	h_leptonicT_mass[c][b] = new TH1D((std::string("h_leptonicT_mass" + suffix)).c_str(),";lep top mass (minChi2)", NmassBins, 0, massmax );
	h_hadronicT_mass[c][b] = new TH1D((std::string("h_hadronicT_mass" + suffix)).c_str(),";had top mass (minChi2)", NmassBins, 0, massmax );
	
	h_hadronicW_pT[c][b] = new TH1D((std::string("h_hadronicW_pT" + suffix)).c_str(),";had W mass (minChi2)", NjetptBins, 0, jetptmax );
	h_hadronicW_pT_mass[c][b] = new TH2D((std::string("h_hadronicW_pT_mass" + suffix)).c_str(),";had W mass (minChi2)", NjetptBins, 0, jetptmax, NmassBins, 0, massmax );
	
	h_leptonicT_pT[c][b] = new TH1D((std::string("h_leptonicT_pT" + suffix)).c_str(),";lep top p_{T} (minChi2)", NjetptBins, 0, jetptmax );
	h_leptonicT_pT_mass[c][b] = new TH2D((std::string("h_leptonicT_pT_mass" + suffix)).c_str(),";lep top p_{T} (minChi2);lep top mass (minChi2)", NjetptBins, 0, jetptmax, NmassBins, 0, massmax );
	h_hadronicT_pT[c][b] = new TH1D((std::string("h_hadronicT_pT" + suffix)).c_str(),";had top p_{T} (minChi2)", NjetptBins, 0, jetptmax );
	h_hadronicT_pT_mass[c][b] = new TH2D((std::string("h_hadronicT_pT_mass" + suffix)).c_str(),";had top p_{T} (minChi2);had top mass (minChi2)", NjetptBins, 0, jetptmax, NmassBins, 0, massmax );
	
	h_lepT_hadT_DeltaR[c][b] = new TH1D((std::string("h_lepT_hadT_DeltaR" + suffix)).c_str(),";DeltaR(lepT,hadT)", 70, 0, 7.0 );
	h_lepT_hadT_Angle[c][b]  = new TH1D((std::string("h_lepT_hadT_Angle" + suffix)).c_str(), ";Angle(lepT,hadT)",  30, 0, 3.5 );
	

	h_sumTop_pT[c][b] = new TH1D((std::string("h_sumTop_pT" + suffix)).c_str(),";had + lep top p_{T} (minChi2)", NjetptBins, 0, jetptmax );
	h_sumTop_mass[c][b] = new TH1D((std::string("h_sumTop_mass" + suffix)).c_str(),";had + lep top mass (minChi2)", NmassBins_sumTop, 0, massmax_sumTop );
	h_sumTop_pT_mass[c][b] = new TH2D((std::string("h_sumTop_pT_mass" + suffix)).c_str(),";had + lep top mass (minChi2)", NjetptBins, 0, jetptmax, NmassBins_sumTop, 0, massmax_sumTop );
	
	h_minChi2_getTopSystem[c][b] = new TH1D((std::string("h_minChi2_getTopSystem" + suffix)).c_str(),";minChi2 (getTopSystem)", 500, 0, 500 );


	h_leptonicW_mass_minChi2_getTopSystem[c][b] = new TH2D((std::string("h_leptonicW_mass_minChi2_getTopSystem" + suffix)).c_str(),";lep W mass (minChi2)", NmassBins, 0, massmax, 500, 0, 500 );
	h_hadronicW_mass_minChi2_getTopSystem[c][b] = new TH2D((std::string("h_hadronicW_mass_minChi2_getTopSystem" + suffix)).c_str(),";had W mass (minChi2)", NmassBins, 0, massmax, 500, 0, 500 );
	h_leptonicT_mass_minChi2_getTopSystem[c][b] = new TH2D((std::string("h_leptonicT_mass_minChi2_getTopSystem" + suffix)).c_str(),";lep top mass (minChi2)", NmassBins, 0, massmax, 500, 0, 500 );
	h_hadronicT_mass_minChi2_getTopSystem[c][b] = new TH2D((std::string("h_hadronicT_mass_minChi2_getTopSystem" + suffix)).c_str(),";had top mass (minChi2)", NmassBins, 0, massmax, 500, 0, 500 );
	
	
	// minChi2 > 30
	h_jet_pt_minChi2GT30[c][b] = new TH1D((std::string("h_jet_pt_minChi2GT30" + suffix)).c_str(),";jet p_{T} (minChi2 < 30)", NjetptBins, 0, jetptmax );
	h_jet_1_pt_minChi2GT30[c][b] = new TH1D((std::string("h_jet_1_pt_minChi2GT30" + suffix)).c_str(),";jet 1 p_{T} (minChi2 < 30)", NjetptBins, 0, jetptmax );
	
	// minChi2 < 30
	h_leptonicW_mass_minChi2LT30[c][b] = new TH1D((std::string("h_leptonicW_mass_minChi2LT30" + suffix)).c_str(),";lep W mass (minChi2 < 30)", NmassBins, 0, massmax );
	h_hadronicW_mass_minChi2LT30[c][b] = new TH1D((std::string("h_hadronicW_mass_minChi2LT30" + suffix)).c_str(),";had W mass (minChi2 < 30)", NmassBins, 0, massmax );
	h_leptonicT_mass_minChi2LT30[c][b] = new TH1D((std::string("h_leptonicT_mass_minChi2LT30" + suffix)).c_str(),";lep top mass (minChi2 < 30)", NmassBins, 0, massmax );
	h_hadronicT_mass_minChi2LT30[c][b] = new TH1D((std::string("h_hadronicT_mass_minChi2LT30" + suffix)).c_str(),";had top mass (minChi2 < 30)", NmassBins, 0, massmax );
	
	h_leptonicT_pT_minChi2LT30[c][b] = new TH1D((std::string("h_leptonicT_pT_minChi2LT30" + suffix)).c_str(),";lep top p_{T} (minChi2 < 30)", NjetptBins, 0, jetptmax );
	h_hadronicT_pT_minChi2LT30[c][b] = new TH1D((std::string("h_hadronicT_pT_minChi2LT30" + suffix)).c_str(),";had top p_{T} (minChi2 < 30)", NjetptBins, 0, jetptmax );
	
	h_sumTop_pT_minChi2LT30[c][b] = new TH1D((std::string("h_sumTop_pT_minChi2LT30" + suffix)).c_str(),";had + lep top p_{T} (minChi2 < 30)", NjetptBins, 0, jetptmax );
	h_sumTop_mass_minChi2LT30[c][b] = new TH1D((std::string("h_sumTop_mass_minChi2LT30" + suffix)).c_str(),";had + lep top mass (minChi2 < 30)", NmassBins_sumTop, 0, massmax_sumTop );
	
	
	h_jet_pt_minChi2LT30[c][b] = new TH1D((std::string("h_jet_pt_minChi2LT30" + suffix)).c_str(),";jet p_{T} (minChi2 < 30)", NjetptBins, 0, jetptmax );
	h_jet_1_pt_minChi2LT30[c][b] = new TH1D((std::string("h_jet_1_pt_minChi2LT30" + suffix)).c_str(),";jet 1 p_{T} (minChi2 < 30)", NjetptBins, 0, jetptmax );
	
	h_jet_pt_minChi2LT30_CSVL_pass[c][b] = new TH1D((std::string("h_jet_pt_minChi2LT30_CSVL_pass" + suffix)).c_str(),";jet p_{T} (minChi2 < 30)", NjetptBins, 0, jetptmax );
	h_jet_pt_minChi2LT30_CSVL_fail[c][b] = new TH1D((std::string("h_jet_pt_minChi2LT30_CSVL_fail" + suffix)).c_str(),";jet p_{T} (minChi2 < 30)", NjetptBins, 0, jetptmax );
	h_jet_pt_minChi2LT30_CSVL_all[c][b]  = new TH1D((std::string("h_jet_pt_minChi2LT30_CSVL_all" + suffix)).c_str(),";jet p_{T} (minChi2 < 30)", NjetptBins, 0, jetptmax );
	
	h_jet_pt_minChi2LT30_CSVM_pass[c][b] = new TH1D((std::string("h_jet_pt_minChi2LT30_CSVM_pass" + suffix)).c_str(),";jet p_{T} (minChi2 < 30)", NjetptBins, 0, jetptmax );
	h_jet_pt_minChi2LT30_CSVM_fail[c][b] = new TH1D((std::string("h_jet_pt_minChi2LT30_CSVM_fail" + suffix)).c_str(),";jet p_{T} (minChi2 < 30)", NjetptBins, 0, jetptmax );
	h_jet_pt_minChi2LT30_CSVM_all[c][b]  = new TH1D((std::string("h_jet_pt_minChi2LT30_CSVM_all" + suffix)).c_str(),";jet p_{T} (minChi2 < 30)", NjetptBins, 0, jetptmax );
	
	h_jet_pt_minChi2LT30_CSVT_pass[c][b] = new TH1D((std::string("h_jet_pt_minChi2LT30_CSVT_pass" + suffix)).c_str(),";jet p_{T} (minChi2 < 30)", NjetptBins, 0, jetptmax );
	h_jet_pt_minChi2LT30_CSVT_fail[c][b] = new TH1D((std::string("h_jet_pt_minChi2LT30_CSVT_fail" + suffix)).c_str(),";jet p_{T} (minChi2 < 30)", NjetptBins, 0, jetptmax );
	h_jet_pt_minChi2LT30_CSVT_all[c][b]  = new TH1D((std::string("h_jet_pt_minChi2LT30_CSVT_all" + suffix)).c_str(),";jet p_{T} (minChi2 < 30)", NjetptBins, 0, jetptmax );
	
	
	// minChi2 < 5
	h_leptonicW_mass_minChi2LT5[c][b] = new TH1D((std::string("h_leptonicW_mass_minChi2LT5" + suffix)).c_str(),";lep W mass (minChi2 < 5)", NmassBins, 0, massmax );
	h_hadronicW_mass_minChi2LT5[c][b] = new TH1D((std::string("h_hadronicW_mass_minChi2LT5" + suffix)).c_str(),";had W mass (minChi2 < 5)", NmassBins, 0, massmax );
	h_leptonicT_mass_minChi2LT5[c][b] = new TH1D((std::string("h_leptonicT_mass_minChi2LT5" + suffix)).c_str(),";lep top mass (minChi2 < 5)", NmassBins, 0, massmax );
	h_hadronicT_mass_minChi2LT5[c][b] = new TH1D((std::string("h_hadronicT_mass_minChi2LT5" + suffix)).c_str(),";had top mass (minChi2 < 5)", NmassBins, 0, massmax );
	
	h_leptonicT_pT_minChi2LT5[c][b] = new TH1D((std::string("h_leptonicT_pT_minChi2LT5" + suffix)).c_str(),";lep top p_{T} (minChi2 < 5)", NjetptBins, 0, jetptmax );
	h_hadronicT_pT_minChi2LT5[c][b] = new TH1D((std::string("h_hadronicT_pT_minChi2LT5" + suffix)).c_str(),";had top p_{T} (minChi2 < 5)", NjetptBins, 0, jetptmax );
	
	h_sumTop_pT_minChi2LT5[c][b] = new TH1D((std::string("h_sumTop_pT_minChi2LT5" + suffix)).c_str(),";had + lep top p_{T} (minChi2 < 5)", NjetptBins, 0, jetptmax );
	h_sumTop_mass_minChi2LT5[c][b] = new TH1D((std::string("h_sumTop_mass_minChi2LT5" + suffix)).c_str(),";had + lep top mass (minChi2 < 5)", NmassBins_sumTop, 0, massmax_sumTop );
	
	h_jet_pt_minChi2LT5[c][b] = new TH1D((std::string("h_jet_pt_minChi2LT5" + suffix)).c_str(),";jet p_{T} (minChi2 < 5)", NjetptBins, 0, jetptmax );
	h_jet_1_pt_minChi2LT5[c][b] = new TH1D((std::string("h_jet_1_pt_minChi2LT5" + suffix)).c_str(),";jet 1 p_{T} (minChi2 < 5)", NjetptBins, 0, jetptmax );
	
	h_jet_pt_minChi2LT5_CSVL_pass[c][b] = new TH1D((std::string("h_jet_pt_minChi2LT5_CSVL_pass" + suffix)).c_str(),";jet p_{T} (minChi2 < 5)", NjetptBins, 0, jetptmax );
	h_jet_pt_minChi2LT5_CSVL_fail[c][b] = new TH1D((std::string("h_jet_pt_minChi2LT5_CSVL_fail" + suffix)).c_str(),";jet p_{T} (minChi2 < 5)", NjetptBins, 0, jetptmax );
	h_jet_pt_minChi2LT5_CSVL_all[c][b]  = new TH1D((std::string("h_jet_pt_minChi2LT5_CSVL_all" + suffix)).c_str(),";jet p_{T} (minChi2 < 5)", NjetptBins, 0, jetptmax );
	
	h_jet_pt_minChi2LT5_CSVM_pass[c][b] = new TH1D((std::string("h_jet_pt_minChi2LT5_CSVM_pass" + suffix)).c_str(),";jet p_{T} (minChi2 < 5)", NjetptBins, 0, jetptmax );
	h_jet_pt_minChi2LT5_CSVM_fail[c][b] = new TH1D((std::string("h_jet_pt_minChi2LT5_CSVM_fail" + suffix)).c_str(),";jet p_{T} (minChi2 < 5)", NjetptBins, 0, jetptmax );
	h_jet_pt_minChi2LT5_CSVM_all[c][b]  = new TH1D((std::string("h_jet_pt_minChi2LT5_CSVM_all" + suffix)).c_str(),";jet p_{T} (minChi2 < 5)", NjetptBins, 0, jetptmax );
	
	h_jet_pt_minChi2LT5_CSVT_pass[c][b] = new TH1D((std::string("h_jet_pt_minChi2LT5_CSVT_pass" + suffix)).c_str(),";jet p_{T} (minChi2 < 5)", NjetptBins, 0, jetptmax );
	h_jet_pt_minChi2LT5_CSVT_fail[c][b] = new TH1D((std::string("h_jet_pt_minChi2LT5_CSVT_fail" + suffix)).c_str(),";jet p_{T} (minChi2 < 5)", NjetptBins, 0, jetptmax );
	h_jet_pt_minChi2LT5_CSVT_all[c][b]  = new TH1D((std::string("h_jet_pt_minChi2LT5_CSVT_all" + suffix)).c_str(),";jet p_{T} (minChi2 < 5)", NjetptBins, 0, jetptmax );
      }
	
	
      h_aplanarity[c][b] = new TH1D((std::string("h_aplanarity" + suffix)).c_str(),";aplanarity", 50, 0, 0.5 );
      h_sphericity[c][b] = new TH1D((std::string("h_sphericity" + suffix)).c_str(),";sphericity", 100, 0, 1.0 );
      h_avg_btag_disc_non_btags[c][b] = new TH1D((std::string("h_avg_btag_disc_non_btags" + suffix)).c_str(),";ave btag disc non btags", 1150, -10.5, 1.0 );
      h_avg_btag_disc_btags[c][b] = new TH1D((std::string("h_avg_btag_disc_btags" + suffix)).c_str(),";ave btag disc btags", 1000, -0.01, 1.01 );
      h_dev_from_avg_disc_btags[c][b] = new TH1D((std::string("h_dev_from_avg_disc_btags" + suffix)).c_str(),";dev from avg disc btags", 100, -0.0001, 0.04 );
      h_lowest_btag[c][b] = new TH1D((std::string("h_lowest_btag" + suffix)).c_str(),";lowest btag", 350, 0.65, 1.001 );
      h_closest_tagged_dijet_mass[c][b] = new TH1D((std::string("h_closest_tagged_dijet_mass" + suffix)).c_str(),";mass of closest tagged jets", NmassBins, 0, massmax );

      h_h0[c][b] = new TH1D((std::string("h_h0" + suffix)).c_str(),";h0", 100, 0, 0.5 );
      h_h1[c][b] = new TH1D((std::string("h_h1" + suffix)).c_str(),";h1", 100, -0.3, 0.5 );
      h_h2[c][b] = new TH1D((std::string("h_h2" + suffix)).c_str(),";h2", 100, -0.2, 0.5 );
      h_h3[c][b] = new TH1D((std::string("h_h3" + suffix)).c_str(),";h3", 100, -0.2, 1.1 );
      h_h4[c][b] = new TH1D((std::string("h_h4" + suffix)).c_str(),";h4", 100, -0.2, 0.3 );


      h_disc_final10v1_8TeV_CFMlpANN[c][b] = new TH1D((std::string("h_disc_final10v1_8TeV_CFMlpANN" + suffix)).c_str(),";CFMlpANN response", Nbins_CFMlpANN, min_reader_final10v1_8TeV_CFMLP[c], max_reader_final10v1_8TeV_CFMLP[c] );
      h_disc_final10v1_8TeV_CFMlpANN_minNPV[c][b] = new TH1D((std::string("h_disc_final10v1_8TeV_CFMlpANN_minNPV" + suffix)).c_str(),";CFMlpANN response", Nbins_CFMlpANN, min_reader_final10v1_8TeV_CFMLP[c], max_reader_final10v1_8TeV_CFMLP[c] );
      h_disc_final10v1_8TeV_CFMlpANN_midNPV[c][b] = new TH1D((std::string("h_disc_final10v1_8TeV_CFMlpANN_midNPV" + suffix)).c_str(),";CFMlpANN response", Nbins_CFMlpANN, min_reader_final10v1_8TeV_CFMLP[c], max_reader_final10v1_8TeV_CFMLP[c] );
      h_disc_final10v1_8TeV_CFMlpANN_maxNPV[c][b] = new TH1D((std::string("h_disc_final10v1_8TeV_CFMlpANN_maxNPV" + suffix)).c_str(),";CFMlpANN response", Nbins_CFMlpANN, min_reader_final10v1_8TeV_CFMLP[c], max_reader_final10v1_8TeV_CFMLP[c] );

      if(longJob){
	// jet 1
	h_jet_1_chargedHadronEnergyFraction[c][b] = new TH1D((std::string("h_jet_1_chargedHadronEnergyFraction" + suffix)).c_str(),";jet 1 chargedHadronEnergyFraction", 100, -0.001, 1.001 );
	h_jet_1_neutralHadronEnergyFraction[c][b] = new TH1D((std::string("h_jet_1_neutralHadronEnergyFraction" + suffix)).c_str(),";jet 1 neutralHadronEnergyFraction", 100, -0.001, 1.001 );
	h_jet_1_chargedEmEnergyFraction[c][b] = new TH1D((std::string("h_jet_1_chargedEmEnergyFraction" + suffix)).c_str(),";jet 1 chargedEmEnergyFraction", 100, -0.001, 1.001 );
	h_jet_1_neutralEmEnergyFraction[c][b] = new TH1D((std::string("h_jet_1_neutralEmEnergyFraction" + suffix)).c_str(),";jet 1 neutralEmEnergyFraction", 100, -0.001, 1.001 );
	h_jet_1_chargedMultiplicity[c][b] = new TH1D((std::string("h_jet_1_chargedMultiplicity" + suffix)).c_str(),";jet 1 chargedMultiplicity", 150, 0, 150 );
	h_jet_1_neutralMultiplicity[c][b] = new TH1D((std::string("h_jet_1_neutralMultiplicity" + suffix)).c_str(),";jet 1 neutralMultiplicity", 150, 0, 150 );
	h_jet_1_nconstituents[c][b] = new TH1D((std::string("h_jet_1_nconstituents" + suffix)).c_str(),";jet 1 nconstituents", 150, 0, 150 );
	
	// jet 2
	h_jet_2_chargedHadronEnergyFraction[c][b] = new TH1D((std::string("h_jet_2_chargedHadronEnergyFraction" + suffix)).c_str(),";jet 2 chargedHadronEnergyFraction", 100, -0.001, 1.001 );
	h_jet_2_neutralHadronEnergyFraction[c][b] = new TH1D((std::string("h_jet_2_neutralHadronEnergyFraction" + suffix)).c_str(),";jet 2 neutralHadronEnergyFraction", 100, -0.001, 1.001 );
	h_jet_2_chargedEmEnergyFraction[c][b] = new TH1D((std::string("h_jet_2_chargedEmEnergyFraction" + suffix)).c_str(),";jet 2 chargedEmEnergyFraction", 100, -0.001, 1.001 );
	h_jet_2_neutralEmEnergyFraction[c][b] = new TH1D((std::string("h_jet_2_neutralEmEnergyFraction" + suffix)).c_str(),";jet 2 neutralEmEnergyFraction", 100, -0.001, 1.001 );
	h_jet_2_chargedMultiplicity[c][b] = new TH1D((std::string("h_jet_2_chargedMultiplicity" + suffix)).c_str(),";jet 2 chargedMultiplicity", 150, 0, 150 );
	h_jet_2_neutralMultiplicity[c][b] = new TH1D((std::string("h_jet_2_neutralMultiplicity" + suffix)).c_str(),";jet 2 neutralMultiplicity", 150, 0, 150 );
	h_jet_2_nconstituents[c][b] = new TH1D((std::string("h_jet_2_nconstituents" + suffix)).c_str(),";jet 2 nconstituents", 150, 0, 150 );
	
	// jet 3
	h_jet_3_chargedHadronEnergyFraction[c][b] = new TH1D((std::string("h_jet_3_chargedHadronEnergyFraction" + suffix)).c_str(),";jet 3 chargedHadronEnergyFraction", 100, -0.001, 1.001 );
	h_jet_3_neutralHadronEnergyFraction[c][b] = new TH1D((std::string("h_jet_3_neutralHadronEnergyFraction" + suffix)).c_str(),";jet 3 neutralHadronEnergyFraction", 100, -0.001, 1.001 );
	h_jet_3_chargedEmEnergyFraction[c][b] = new TH1D((std::string("h_jet_3_chargedEmEnergyFraction" + suffix)).c_str(),";jet 3 chargedEmEnergyFraction", 100, -0.001, 1.001 );
	h_jet_3_neutralEmEnergyFraction[c][b] = new TH1D((std::string("h_jet_3_neutralEmEnergyFraction" + suffix)).c_str(),";jet 3 neutralEmEnergyFraction", 100, -0.001, 1.001 );
	h_jet_3_chargedMultiplicity[c][b] = new TH1D((std::string("h_jet_3_chargedMultiplicity" + suffix)).c_str(),";jet 3 chargedMultiplicity", 150, 0, 150 );
	h_jet_3_neutralMultiplicity[c][b] = new TH1D((std::string("h_jet_3_neutralMultiplicity" + suffix)).c_str(),";jet 3 neutralMultiplicity", 150, 0, 150 );
	h_jet_3_nconstituents[c][b] = new TH1D((std::string("h_jet_3_nconstituents" + suffix)).c_str(),";jet 3 nconstituents", 150, 0, 150 );
	
	// jet 4
	h_jet_4_chargedHadronEnergyFraction[c][b] = new TH1D((std::string("h_jet_4_chargedHadronEnergyFraction" + suffix)).c_str(),";jet 4 chargedHadronEnergyFraction", 100, -0.001, 1.001 );
	h_jet_4_neutralHadronEnergyFraction[c][b] = new TH1D((std::string("h_jet_4_neutralHadronEnergyFraction" + suffix)).c_str(),";jet 4 neutralHadronEnergyFraction", 100, -0.001, 1.001 );
	h_jet_4_chargedEmEnergyFraction[c][b] = new TH1D((std::string("h_jet_4_chargedEmEnergyFraction" + suffix)).c_str(),";jet 4 chargedEmEnergyFraction", 100, -0.001, 1.001 );
	h_jet_4_neutralEmEnergyFraction[c][b] = new TH1D((std::string("h_jet_4_neutralEmEnergyFraction" + suffix)).c_str(),";jet 4 neutralEmEnergyFraction", 100, -0.001, 1.001 );
	h_jet_4_chargedMultiplicity[c][b] = new TH1D((std::string("h_jet_4_chargedMultiplicity" + suffix)).c_str(),";jet 4 chargedMultiplicity", 150, 0, 150 );
	h_jet_4_neutralMultiplicity[c][b] = new TH1D((std::string("h_jet_4_neutralMultiplicity" + suffix)).c_str(),";jet 4 neutralMultiplicity", 150, 0, 150 );
	h_jet_4_nconstituents[c][b] = new TH1D((std::string("h_jet_4_nconstituents" + suffix)).c_str(),";jet 4 nconstituents", 150, 0, 150 );
	
	
	for( int iPV=0; iPV<NumPV; iPV++ ){
	  std::string npv_suffix = "_";
	  TString PVbin_start = Form("%d",numPVbins_lowEdge[iPV]);
	  TString PVbin_end   = ( iPV<(NumPV-1) ) ? Form("%d",numPVbins_lowEdge[iPV+1]) : "Inf";
	  
	  npv_suffix += std::string( "NPV" + PVbin_start + "to" + PVbin_end );
	  npv_suffix += "_";
	  npv_suffix += cat_labels[c];
	  npv_suffix += sys_cat_labels[b];
	  
	  h_jet_pt_NPV[iPV][c][b] = new TH1D((std::string("h_jet_pt" + npv_suffix)).c_str(),";jet p_{T}", NjetptBins, 0, jetptmax );
	  h_jet_1_pt_NPV[iPV][c][b] = new TH1D((std::string("h_jet_1_pt" + npv_suffix)).c_str(),";highest jet p_{T}", NjetptBins, 0, jetptmax );
	  h_jet_2_pt_NPV[iPV][c][b] = new TH1D((std::string("h_jet_2_pt" + npv_suffix)).c_str(),";second-highest jet p_{T}", NjetptBins, 0, jetptmax );
	  h_jet_3_pt_NPV[iPV][c][b] = new TH1D((std::string("h_jet_3_pt" + npv_suffix)).c_str(),";third-highest jet p_{T}", NjetptBins, 0, jetptmax );
	  h_jet_4_pt_NPV[iPV][c][b] = new TH1D((std::string("h_jet_4_pt" + npv_suffix)).c_str(),";fourth-highest jet p_{T}", NjetptBins, 0, jetptmax );
	  
	  h_tag_jet_pt_NPV[iPV][c][b] = new TH1D((std::string("h_tag_jet_pt" + npv_suffix)).c_str(),";tagged jet p_{T}", NjetptBins, 0, jetptmax );
	  h_tag_jet_1_pt_NPV[iPV][c][b] = new TH1D((std::string("h_tag_jet_1_pt" + npv_suffix)).c_str(),";highest tagged jet p_{T}", NjetptBins, 0, jetptmax );
	  h_tag_jet_2_pt_NPV[iPV][c][b] = new TH1D((std::string("h_tag_jet_2_pt" + npv_suffix)).c_str(),";second-highest tagged jet p_{T}", NjetptBins, 0, jetptmax );
	  h_tag_jet_3_pt_NPV[iPV][c][b] = new TH1D((std::string("h_tag_jet_3_pt" + npv_suffix)).c_str(),";third-highest tagged jet p_{T}", NjetptBins, 0, jetptmax );
	  h_tag_jet_4_pt_NPV[iPV][c][b] = new TH1D((std::string("h_tag_jet_4_pt" + npv_suffix)).c_str(),";fourth-highest tagged jet p_{T}", NjetptBins, 0, jetptmax );
	  
	  h_untag_jet_pt_NPV[iPV][c][b] = new TH1D((std::string("h_untag_jet_pt" + npv_suffix)).c_str(),";untagged jet p_{T}", NjetptBins, 0, jetptmax );
	  h_untag_jet_1_pt_NPV[iPV][c][b] = new TH1D((std::string("h_untag_jet_1_pt" + npv_suffix)).c_str(),";highest untagged jet p_{T}", NjetptBins, 0, jetptmax );
	  h_untag_jet_2_pt_NPV[iPV][c][b] = new TH1D((std::string("h_untag_jet_2_pt" + npv_suffix)).c_str(),";second-highest untagged jet p_{T}", NjetptBins, 0, jetptmax );
	  h_untag_jet_3_pt_NPV[iPV][c][b] = new TH1D((std::string("h_untag_jet_3_pt" + npv_suffix)).c_str(),";third-highest untagged jet p_{T}", NjetptBins, 0, jetptmax );
	  h_untag_jet_4_pt_NPV[iPV][c][b] = new TH1D((std::string("h_untag_jet_4_pt" + npv_suffix)).c_str(),";fourth-highest untagged jet p_{T}", NjetptBins, 0, jetptmax );
	}
      }

    } // end loop over categories

  } // end loop over systematics



  //////////////////////////////////////////////////////////////////////////
  ///  Weights
  //////////////////////////////////////////////////////////////////////////
  
  ////////////////////////

  /*
  int numhists[NumCat];
  for( int c=0; c<NumCat-1; c++ ){
    numhists[c] = 10;
  }
  numhists[NumCat-1] = 10;

  double histxmin[NumCat][11];
  double histxmax[NumCat][11];
  int histbins[NumCat][11];
  std::string histXaxis[NumCat][11];


  // 4j2t
  histbins[0][0] = NlepPtBins;  histxmin[0][0] = 0; histxmax[0][0] = lepPtMax; histXaxis[0][0] = "tight_lepton_pt";
  histbins[0][1] = 30;          histxmin[0][1] = 0; histxmax[0][1] = 6.;       histXaxis[0][1] = "min_dr_tagged_jets";
  histbins[0][2] = 25;          histxmin[0][2] = 0; histxmax[0][2] = 0.5;      histXaxis[0][2] = "aplanarity";
  histbins[0][3] = NmetBins/2;  histxmin[0][3] = 0; histxmax[0][3] = metmax;   histXaxis[0][3] = "MET";
  histbins[0][4] = NjetptBins;  histxmin[0][4] = 0; histxmax[0][4] = jetptmax; histXaxis[0][4] = "third_jet_pt";
  histbins[0][5] = NjetptBins;  histxmin[0][5] = 0; histxmax[0][5] = jetptmax; histXaxis[0][5] = "fourth_jet_pt";
  histbins[0][6] = NmassBins/2; histxmin[0][6] = 0; histxmax[0][6] = massmax;  histXaxis[0][6] = "avg_tagged_dijet_mass";
  histbins[0][7] = NmassBins/2; histxmin[0][7] = 0; histxmax[0][7] = massmax;  histXaxis[0][7] = "avg_untagged_dijet_mass";
  histbins[0][8] = NhtBins;     histxmin[0][8] = 0; histxmax[0][8] = htmax;    histXaxis[0][8] = "invariant_mass_of_everything";
  histbins[0][9] = NmetBins;    histxmin[0][9] = 0; histxmax[0][9] = metmax;   histXaxis[0][9] = "MT_mht";

  // 5j2t
  histbins[1][0] = 20;          histxmin[1][0] = 0;    histxmax[1][0] = 1.0;      histXaxis[1][0] = "sphericity";
  histbins[1][1] = NmetBins/2;  histxmin[1][1] = 0;    histxmax[1][1] = metmax;   histXaxis[1][1] = "MET";
  histbins[1][2] = 30;          histxmin[1][2] = 0;    histxmax[1][2] = 6.;       histXaxis[1][2] = "dr_between_lep_and_closest_jet";
  histbins[1][3] = NjetptBins;  histxmin[1][3] = 0;    histxmax[1][3] = jetptmax; histXaxis[1][3] = "third_jet_pt";
  histbins[1][4] = NjetptBins;  histxmin[1][4] = 0;    histxmax[1][4] = jetptmax; histXaxis[1][4] = "fourth_jet_pt";
  histbins[1][5] = 20;          histxmin[1][5] = -0.3; histxmax[1][5] = 0.5;      histXaxis[1][5] = "h1";
  histbins[1][6] = NhtBins/2;   histxmin[1][6] = 0;    histxmax[1][6] = htmax;    histXaxis[1][6] = "all_sum_pt_with_met";
  histbins[1][7] = NmassBins/2; histxmin[1][7] = 0;    histxmax[1][7] = massmax;  histXaxis[1][7] = "avg_untagged_dijet_mass";
  histbins[1][8] = Nm3Bins/5;   histxmin[1][8] = 0;    histxmax[1][8] = m3max;    histXaxis[1][8] = "M3_1tag";
  histbins[1][9] = NmassBins/2; histxmin[1][9] = 0;    histxmax[1][9] = massmax;  histXaxis[1][9] = "Mlb";

  // 6j2t
  histbins[2][0] = 25;          histxmin[2][0] = 0;     histxmax[2][0] = 0.5;      histXaxis[2][0] = "aplanarity";
  histbins[2][1] = 20;          histxmin[2][1] = 0;     histxmax[2][1] = 1.0;      histXaxis[2][1] = "sphericity";
  histbins[2][2] = 20;          histxmin[2][2] = 0.;    histxmax[2][2] = 0.5;      histXaxis[2][2] = "h0";
  histbins[2][3] = 40;          histxmin[2][3] = 0.601; histxmax[2][3] = 1.01;     histXaxis[2][3] = "avg_btag_disc_btags";
  histbins[2][4] = NjetptBins;  histxmin[2][4] = 0;     histxmax[2][4] = jetptmax; histXaxis[2][4] = "third_jet_pt";
  histbins[2][5] = NjetptBins;  histxmin[2][5] = 0;     histxmax[2][5] = jetptmax; histXaxis[2][5] = "fourth_jet_pt";
  histbins[2][6] = 20;          histxmin[2][6] = -0.3;  histxmax[2][6] = 0.5;      histXaxis[2][6] = "h1";
  histbins[2][7] = NmassBins/2; histxmin[2][7] = 0;     histxmax[2][7] = massmax;  histXaxis[2][7] = "avg_untagged_dijet_mass";
  histbins[2][8] = 20;          histxmin[2][8] = -0.2;  histxmax[2][8] = 1.1;      histXaxis[2][8] = "h3";
  histbins[2][9] = NhtBins;     histxmin[2][9] = 0;     histxmax[2][9] = htmax;    histXaxis[2][9] = "invariant_mass_of_everything";

  // 4j3t
  histbins[3][0] = NjetptBins;  histxmin[3][0] = 0;     histxmax[3][0] = jetptmax; histXaxis[3][0] = "first_jet_pt";
  histbins[3][1] = NjetptBins;  histxmin[3][1] = 0;     histxmax[3][1] = jetptmax; histXaxis[3][1] = "second_jet_pt";
  histbins[3][2] = 40;          histxmin[3][2] = 0.601; histxmax[3][2] = 1.01;     histXaxis[3][2] = "avg_btag_disc_btags";
  histbins[3][3] = 12;          histxmin[3][3] = 0;     histxmax[3][3] = 0.04;     histXaxis[3][3] = "dev_from_avg_disc_btags";
  histbins[3][4] = NjetptBins;  histxmin[3][4] = 0;     histxmax[3][4] = jetptmax; histXaxis[3][4] = "third_jet_pt";
  histbins[3][5] = NjetptBins;  histxmin[3][5] = 0;     histxmax[3][5] = jetptmax; histXaxis[3][5] = "fourth_jet_pt";
  histbins[3][6] = 40;          histxmin[3][6] = 0.601; histxmax[3][6] = 1.01;     histXaxis[3][6] = "lowest_btag";
  histbins[3][7] = NhtBins/2;   histxmin[3][7] = 0;     histxmax[3][7] = htmax;    histXaxis[3][7] = "all_sum_pt_with_met";
  histbins[3][8] = 40;          histxmin[3][8] = 0.601; histxmax[3][8] = 1.01;     histXaxis[3][8] = "second_highest_btag";
  histbins[3][9] = NhtBins;     histxmin[3][9] = 0;     histxmax[3][9] = htmax;    histXaxis[3][9] = "invariant_mass_of_everything";

  // 5j3t
  histbins[4][0] = NjetptBins;  histxmin[4][0] = 0;     histxmax[4][0] = jetptmax; histXaxis[4][0] = "first_jet_pt";
  histbins[4][1] = 30;          histxmin[4][1] = 0;     histxmax[4][1] = 6.;       histXaxis[4][1] = "min_dr_tagged_jets";
  histbins[4][2] = NjetptBins;  histxmin[4][2] = 0;     histxmax[4][2] = jetptmax; histXaxis[4][2] = "second_jet_pt";
  histbins[4][3] = 40;          histxmin[4][3] = 0.601; histxmax[4][3] = 1.01;     histXaxis[4][3] = "avg_btag_disc_btags";
  histbins[4][4] = 12;          histxmin[4][4] = 0;     histxmax[4][4] = 0.04;     histXaxis[4][4] = "dev_from_avg_disc_btags";
  histbins[4][5] = NjetptBins;  histxmin[4][5] = 0;     histxmax[4][5] = jetptmax; histXaxis[4][5] = "third_jet_pt";
  histbins[4][6] = NjetptBins;  histxmin[4][6] = 0;     histxmax[4][6] = jetptmax; histXaxis[4][6] = "fourth_jet_pt";
  histbins[4][7] = 40;          histxmin[4][7] = 0.601; histxmax[4][7] = 1.01;     histXaxis[4][7] = "lowest_btag";
  histbins[4][8] = NhtBins/2;   histxmin[4][8] = 0;     histxmax[4][8] = htmax;    histXaxis[4][8] = "all_sum_pt_with_met";
  histbins[4][9] = 40;          histxmin[4][9] = 0.601; histxmax[4][9] = 1.01;     histXaxis[4][9] = "second_highest_btag";

  // 6j3t
  histbins[5][0] = 30;          histxmin[5][0] = 0;     histxmax[5][0] = 6.;       histXaxis[5][0] = "avg_dr_tagged_jets";
  histbins[5][1] = 20;          histxmin[5][1] = 0;     histxmax[5][1] = 1.0;      histXaxis[5][1] = "sphericity";
  histbins[5][2] = 40;          histxmin[5][2] = 0.601; histxmax[5][2] = 1.01;     histXaxis[5][2] = "avg_btag_disc_btags";
  histbins[5][3] = 12;          histxmin[5][3] = 0;     histxmax[5][3] = 0.04;     histXaxis[5][3] = "dev_from_avg_disc_btags";
  histbins[5][4] = 20;          histxmin[5][4] = -0.2;  histxmax[5][4] = 0.5;      histXaxis[5][4] = "h2";
  histbins[5][5] = 40;          histxmin[5][5] = 0.601; histxmax[5][5] = 1.01;     histXaxis[5][5] = "lowest_btag";
  histbins[5][6] = NmassBins/2; histxmin[5][6] = 0;     histxmax[5][6] = massmax;  histXaxis[5][6] = "avg_untagged_dijet_mass";
  histbins[5][7] = 20;          histxmin[5][7] = -0.2;  histxmax[5][7] = 1.1;      histXaxis[5][7] = "h3";
  histbins[5][8] = 40;          histxmin[5][8] = 0.601; histxmax[5][8] = 1.01;     histXaxis[5][8] = "second_highest_btag";
  histbins[5][9] = NhtBins;     histxmin[5][9] = 0;     histxmax[5][9] = htmax;    histXaxis[5][9] = "invariant_mass_of_everything";

  // 4j4t
  histbins[6][0] = NjetptBins;  histxmin[6][0] = 0;     histxmax[6][0] = jetptmax; histXaxis[6][0] = "first_jet_pt";
  histbins[6][1] = 20;          histxmin[6][1] = -0.3;  histxmax[6][1] = 0.5;      histXaxis[6][1] = "h1";
  histbins[6][2] = 30;          histxmin[6][2] = 0;     histxmax[6][2] = 6.;       histXaxis[6][2] = "avg_dr_tagged_jets";
  histbins[6][3] = 25;          histxmin[6][3] = 0;     histxmax[6][3] = 0.5;      histXaxis[6][3] = "aplanarity";
  histbins[6][4] = 40;          histxmin[6][4] = 0.601; histxmax[6][4] = 1.01;     histXaxis[6][4] = "avg_btag_disc_btags";
  histbins[6][5] = 12;          histxmin[6][5] = 0;     histxmax[6][5] = 0.04;     histXaxis[6][5] = "dev_from_avg_disc_btags";
  histbins[6][6] = 40;          histxmin[6][6] = 0.601; histxmax[6][6] = 1.01;     histXaxis[6][6] = "lowest_btag";
  histbins[6][7] = NhtBins/2;   histxmin[6][7] = 0;     histxmax[6][7] = htmax;    histXaxis[6][7] = "all_sum_pt_with_met";
  histbins[6][8] = 40;          histxmin[6][8] = 0.601; histxmax[6][8] = 1.01;     histXaxis[6][8] = "second_highest_btag";
  histbins[6][9] = NhtBins;     histxmin[6][9] = 0;     histxmax[6][9] = htmax;    histXaxis[6][9] = "invariant_mass_of_everything";

  // 5j4t
  histbins[7][0] = 30;          histxmin[7][0] = 0;     histxmax[7][0] = 6.;       histXaxis[7][0] = "avg_dr_tagged_jets";
  histbins[7][1] = 30;          histxmin[7][1] = 0;     histxmax[7][1] = 6.;       histXaxis[7][1] = "dr_between_lep_and_closest_jet";
  histbins[7][2] = 40;          histxmin[7][2] = 0.601; histxmax[7][2] = 1.01;     histXaxis[7][2] = "avg_btag_disc_btags";
  histbins[7][3] = 12;          histxmin[7][3] = 0;     histxmax[7][3] = 0.04;     histXaxis[7][3] = "dev_from_avg_disc_btags";
  histbins[7][4] = NjetptBins;  histxmin[7][4] = 0;     histxmax[7][4] = jetptmax; histXaxis[7][4] = "third_jet_pt";
  histbins[7][5] = NjetptBins;  histxmin[7][5] = 0;     histxmax[7][5] = jetptmax; histXaxis[7][5] = "fourth_jet_pt";
  histbins[7][6] = 40;          histxmin[7][6] = 0.601; histxmax[7][6] = 1.01;     histXaxis[7][6] = "lowest_btag";
  histbins[7][7] = NhtBins/2;   histxmin[7][7] = 0;     histxmax[7][7] = htmax;    histXaxis[7][7] = "all_sum_pt_with_met";
  histbins[7][8] = 40;          histxmin[7][8] = 0.601; histxmax[7][8] = 1.01;     histXaxis[7][8] = "first_highest_btag";
  histbins[7][9] = 40;          histxmin[7][9] = 0.601; histxmax[7][9] = 1.01;     histXaxis[7][9] = "second_highest_btag";

  // 6j4t
  histbins[8][0] = 30;          histxmin[8][0] = 0;     histxmax[8][0] = 6.;       histXaxis[8][0] = "avg_dr_tagged_jets";
  histbins[8][1] = 20;          histxmin[8][1] = 0;     histxmax[8][1] = 1.0;      histXaxis[8][1] = "sphericity";
  histbins[8][2] = 30;          histxmin[8][2] = 0;     histxmax[8][2] = 6.;       histXaxis[8][2] = "dr_between_lep_and_closest_jet";
  histbins[8][3] = 40;          histxmin[8][3] = 0.601; histxmax[8][3] = 1.01;     histXaxis[8][3] = "avg_btag_disc_btags";  
  histbins[8][4] = 40;          histxmin[8][4] = 0.601; histxmax[8][4] = 1.01;     histXaxis[8][4] = "second_highest_btag";
  histbins[8][5] = 20;          histxmin[8][5] = -0.2;  histxmax[8][5] = 0.5;      histXaxis[8][5] = "h2";
  histbins[8][6] = 40;          histxmin[8][6] = 0.601; histxmax[8][6] = 1.01;     histXaxis[8][6] = "lowest_btag";
  histbins[8][7] = NmassBins/2; histxmin[8][7] = 0;     histxmax[8][7] = massmax;  histXaxis[8][7] = "closest_tagged_dijet_mass";
  histbins[8][8] = 20;          histxmin[8][8] = -0.2;  histxmax[8][8] = 1.1;      histXaxis[8][8] = "h3";
  histbins[8][9] = NhtBins;     histxmin[8][9] = 0;     histxmax[8][9] = htmax;    histXaxis[8][9] = "invariant_mass_of_everything";

  /////////////


  TH2D* h_profiles[NumCat][55];
  for( int c=0; c<NumCat; c++ ){
    //if( c!=0 ) continue;
    int p=0;
    std::string cat_suffix = "_" + cat_labels[c];

    for( int i=0; i<numhists[c]; i++ ){
      for( int j=i+1; j<numhists[c]; j++ ){
  	std::string title = std::string("h_" + histXaxis[c][i] + "_" + histXaxis[c][j] + cat_suffix);
  	std::string axistitle = std::string(";" + histXaxis[c][i] + ";" + histXaxis[c][j]);
  	h_profiles[c][p] = new TH2D(title.c_str(), axistitle.c_str(), histbins[c][i], histxmin[c][i], histxmax[c][i], histbins[c][j], histxmin[c][j], histxmax[c][j] );
  	p++;
      }
    }
  }
  */



  int nentries = chain->GetEntries();
  std::cout << "\n\t Number of entries = " << nentries << std::endl;
  std::cout << "\t Max number of entries = " << maxNentries << std::endl;
  std::cout << "\n" << std::endl;

  int use_nentries = std::max( maxNentries, nentries);

  int NeventsPerJob = int( double(use_nentries)/double(Njobs) + 0.000001 ) + 1;

  int firstEvent = (jobN-1)*NeventsPerJob + 1;
  int lastEvent  = firstEvent + NeventsPerJob;
  if( jobN==Njobs ) lastEvent = -1;


  std::cout << "========  Starting Event Loop  ========" << std::endl;
  for (Long64_t ievt=0; ievt<chain->GetEntries();ievt++) {    //Long64_t

    if( ievt<firstEvent ) continue;
    if( ievt==lastEvent ) break;

    if( ievt==1 )        std::cout << "     Event " << ievt << std::endl;
    if( ievt%10000==0 && ievt!=1 ) std::cout << "           " << ievt << "\t" 
					     << int(double(ievt-firstEvent)/double(NeventsPerJob)*100) << "% done" << std::endl;

    if( ievt==(maxNentries+1) ) break;

    chain->GetEntry(ievt);

    bool passEvent = false;
    //if( isMuon && eve->leptonType_==1 )       passEvent = eve->passMuonTrigger_;
    //else if( !isMuon && eve->leptonType_==0 ) passEvent = eve->passElectronTrigger_;

    if(leptonChoice==0){
      if(eve->leptonType_==0 && eve->passElectronTrigger_==1) passEvent=true;
    }
    else if(leptonChoice==1){
      if(eve->leptonType_==1 && eve->passMuonTrigger_==1) passEvent=true;
    }
    else {
      if((eve->leptonType_==1 && eve->passMuonTrigger_==1) || (eve->leptonType_==0 && eve->passElectronTrigger_==1)) passEvent=true;
    }


    if( !passEvent ) continue;


    int numPVs = eve->numPVs_;
    bool isMinNPV = ( numPVs>=0 && numPVs<=12 );
    bool isMidNPV = ( numPVs>12 && numPVs<18 );


    int run = eve->run_;
    int lumi = eve->lumi_;
    long evt = eve->evt_;


    double wgt_top_pt = 1;
    double wgt_top_pt_up = 1;
    double wgt_top_pt_down = 1;
    if( sampleLoader.isTTJets ){
      wgt_top_pt = getTopPtWgt( eve->top_pt_, 0 );
      wgt_top_pt_up = getTopPtWgt( eve->top_pt_, 1 );
      wgt_top_pt_down = getTopPtWgt( eve->top_pt_, 2 );
    }


    for(int iSys=0; iSys<kNumSys; iSys++){

      numJets_float = eve->numJets_float_[iSys];
      numTags_float = eve->numTags_float_[iSys];

      // convert floats to ints for nJets, nTags
      int numJets = int( numJets_float + 0.001 );
      int numTags = int( numTags_float + 0.001 );

      int njet = ( numJets>NjetMax ) ? NjetMax : numJets;
      int ntag = ( numTags>NtagMax ) ? NtagMax : numTags;
      int njet_full = ( numJets>NjetMax_full ) ? NjetMax_full : numJets;
      int ntag_full = ( numTags>NtagMax_full ) ? NtagMax_full : numTags;



      // 
      // Associate the variables with the MVA reader with the vars in the tree
      //
      tight_lepton_pt           = eve->tight_lepton_pt_;
      first_jet_pt              = eve->first_jet_pt_[iSys];
      min_dr_tagged_jets        = eve->min_dr_tagged_jets_[iSys];
      avg_dr_tagged_jets        = eve->avg_dr_tagged_jets_[iSys];
      aplanarity                = eve->aplanarity_[iSys];
      sphericity                = eve->sphericity_[iSys];
      MET                       = eve->MET_[iSys];
      second_jet_pt             = eve->second_jet_pt_[iSys];
      avg_btag_disc_btags       = eve->avg_btag_disc_btags_[iSys]; 
      dev_from_avg_disc_btags   = eve->dev_from_avg_disc_btags_[iSys];
      third_jet_pt              = eve->third_jet_pt_[iSys];
      fourth_jet_pt             = eve->fourth_jet_pt_[iSys];
      avg_tagged_dijet_mass     = eve->avg_tagged_dijet_mass_[iSys];
      lowest_btag               = eve->lowest_btag_[iSys];
      all_sum_pt_incl_met       = eve->all_sum_pt_with_met_[iSys];
      avg_untagged_dijet_mass   = eve->avg_untagged_dijet_mass_[iSys];
      closest_tagged_dijet_mass = eve->closest_tagged_dijet_mass_[iSys];
      first_highest_btag        = eve->first_highest_btag_[iSys];
      second_highest_btag       = eve->second_highest_btag_[iSys];
      third_highest_btag        = eve->third_highest_btag_[iSys];
      fourth_highest_btag       = eve->fourth_highest_btag_[iSys];
      best_higgs_mass           = eve->best_higgs_mass_[iSys];
      invariant_mass_of_everything   = eve->invariant_mass_of_everything_[iSys];
      dr_between_lep_and_closest_jet = eve->dr_between_lep_and_closest_jet_[iSys];
      
      dijet_mass_of_everything = invariant_mass_of_everything;

      h0 = eve->h0_[iSys];
      h1 = eve->h1_[iSys];
      h2 = eve->h2_[iSys];
      h3 = eve->h3_[iSys];
      h4 = eve->h4_[iSys];

      minChi2 = eve->minChi2_[iSys];
      dRbb = eve->dRbb_[iSys];
      avg_btag_disc_non_btags = eve->avg_btag_disc_non_btags_[iSys];
      all_sum_pt = eve->all_sum_pt_[iSys];

      all_sum_pt_with_met = all_sum_pt_incl_met;

      double Mlb     = eve->Mlb_[iSys];
      double M3      = eve->M3_[iSys];
      double M3_1tag = eve->M3_1tag_[iSys];

      double MHT     = eve->MHT_[iSys];
      double MHT_phi = eve->MHT_phi_[iSys];

      double MET_phi = eve->MET_phi_[iSys];
      double tight_lepton_phi = eve->tight_lepton_phi_;
      double tight_lepton_eta = eve->tight_lepton_eta_;

      double jet_tag_1_pt = eve->jet_tag_1_pt_[iSys];
      double jet_tag_2_pt = eve->jet_tag_2_pt_[iSys];
      double jet_tag_3_pt = eve->jet_tag_3_pt_[iSys];
      double jet_tag_4_pt = eve->jet_tag_4_pt_[iSys];

      double jet_1_chargedHadronEnergyFraction = eve->jet_1_chargedHadronEnergyFraction_[iSys];
      double jet_1_neutralHadronEnergyFraction = eve->jet_1_neutralHadronEnergyFraction_[iSys];
      double jet_1_chargedEmEnergyFraction = eve->jet_1_chargedEmEnergyFraction_[iSys];
      double jet_1_neutralEmEnergyFraction = eve->jet_1_neutralEmEnergyFraction_[iSys];
      double jet_1_chargedMultiplicity = eve->jet_1_chargedMultiplicity_[iSys];
      double jet_1_neutralMultiplicity = eve->jet_1_neutralMultiplicity_[iSys];
      double jet_1_nconstituents = eve->jet_1_nconstituents_[iSys];

      double jet_2_chargedHadronEnergyFraction = eve->jet_2_chargedHadronEnergyFraction_[iSys];
      double jet_2_neutralHadronEnergyFraction = eve->jet_2_neutralHadronEnergyFraction_[iSys];
      double jet_2_chargedEmEnergyFraction = eve->jet_2_chargedEmEnergyFraction_[iSys];
      double jet_2_neutralEmEnergyFraction = eve->jet_2_neutralEmEnergyFraction_[iSys];
      double jet_2_chargedMultiplicity = eve->jet_2_chargedMultiplicity_[iSys];
      double jet_2_neutralMultiplicity = eve->jet_2_neutralMultiplicity_[iSys];
      double jet_2_nconstituents = eve->jet_2_nconstituents_[iSys];

      double jet_3_chargedHadronEnergyFraction = eve->jet_3_chargedHadronEnergyFraction_[iSys];
      double jet_3_neutralHadronEnergyFraction = eve->jet_3_neutralHadronEnergyFraction_[iSys];
      double jet_3_chargedEmEnergyFraction = eve->jet_3_chargedEmEnergyFraction_[iSys];
      double jet_3_neutralEmEnergyFraction = eve->jet_3_neutralEmEnergyFraction_[iSys];
      double jet_3_chargedMultiplicity = eve->jet_3_chargedMultiplicity_[iSys];
      double jet_3_neutralMultiplicity = eve->jet_3_neutralMultiplicity_[iSys];
      double jet_3_nconstituents = eve->jet_3_nconstituents_[iSys];

      double jet_4_chargedHadronEnergyFraction = eve->jet_4_chargedHadronEnergyFraction_[iSys];
      double jet_4_neutralHadronEnergyFraction = eve->jet_4_neutralHadronEnergyFraction_[iSys];
      double jet_4_chargedEmEnergyFraction = eve->jet_4_chargedEmEnergyFraction_[iSys];
      double jet_4_neutralEmEnergyFraction = eve->jet_4_neutralEmEnergyFraction_[iSys];
      double jet_4_chargedMultiplicity = eve->jet_4_chargedMultiplicity_[iSys];
      double jet_4_neutralMultiplicity = eve->jet_4_neutralMultiplicity_[iSys];
      double jet_4_nconstituents = eve->jet_4_nconstituents_[iSys];


      if( third_jet_pt<40. ) continue;

      double wgt = eve->wgt_[iSys]*eve->wgt_btagSF_[iSys];

      int this_category = -1;
      if( numJets==4 && numTags==2) this_category=0;
      if( numJets==5 && numTags==2) this_category=1;
      if( numJets>=6 && numTags==2) this_category=2;	
      if( numJets==4 && numTags==3) this_category=3;
      if( numJets==5 && numTags==3) this_category=4;
      if( numJets>=6 && numTags==3) this_category=5;
      if( numJets==4 && numTags>=4) this_category=6;
      if( numJets==5 && numTags>=4) this_category=7;
      if( numJets>=6 && numTags>=4) this_category=8;


      double MT_met = sqrt( 2 * tight_lepton_pt * MET * (1-cos(tight_lepton_phi - MET_phi)) );
      double MT_mht = sqrt( 2 * tight_lepton_pt * MHT * (1-cos(tight_lepton_phi - MHT_phi)) );

      /*
      double variable[NumCat][11];
      // 4j2t
      variable[0][0] = tight_lepton_pt;
      variable[0][1] = min_dr_tagged_jets;
      variable[0][2] = aplanarity;
      variable[0][3] = MET;
      variable[0][4] = third_jet_pt;
      variable[0][5] = fourth_jet_pt;
      variable[0][6] = avg_tagged_dijet_mass;
      variable[0][7] = avg_untagged_dijet_mass;
      variable[0][8] = invariant_mass_of_everything;
      variable[0][9] = MT_mht;

      // 5j2t
      variable[1][0] = sphericity;
      variable[1][1] = MET;
      variable[1][2] = dr_between_lep_and_closest_jet;
      variable[1][3] = third_jet_pt;
      variable[1][4] = fourth_jet_pt;
      variable[1][5] = h1;
      variable[1][6] = all_sum_pt_with_met;
      variable[1][7] = avg_untagged_dijet_mass;
      variable[1][8] = M3_1tag;
      variable[1][9] = Mlb;

      // 6j2t
      variable[2][0] = aplanarity;
      variable[2][1] = sphericity;
      variable[2][2] = h0;
      variable[2][3] = avg_btag_disc_btags;
      variable[2][4] = third_jet_pt;
      variable[2][5] = fourth_jet_pt;
      variable[2][6] = h1;
      variable[2][7] = avg_untagged_dijet_mass;
      variable[2][8] = h3;
      variable[2][9] = invariant_mass_of_everything;

      // 4j3t
      variable[3][0] = first_jet_pt;
      variable[3][1] = second_jet_pt;
      variable[3][2] = avg_btag_disc_btags;
      variable[3][3] = dev_from_avg_disc_btags;
      variable[3][4] = third_jet_pt;
      variable[3][5] = fourth_jet_pt;
      variable[3][6] = lowest_btag;
      variable[3][7] = all_sum_pt_with_met;
      variable[3][8] = second_highest_btag;
      variable[3][9] = invariant_mass_of_everything;

      // 5j3t
      variable[4][0] = first_jet_pt;
      variable[4][1] = min_dr_tagged_jets;
      variable[4][2] = second_jet_pt;
      variable[4][3] = avg_btag_disc_btags;
      variable[4][4] = dev_from_avg_disc_btags;
      variable[4][5] = third_jet_pt;
      variable[4][6] = fourth_jet_pt;
      variable[4][7] = lowest_btag;
      variable[4][8] = all_sum_pt_with_met;
      variable[4][9] = second_highest_btag;

      // 6j3t
      variable[5][0] = avg_dr_tagged_jets;
      variable[5][1] = sphericity;
      variable[5][2] = avg_btag_disc_btags;
      variable[5][3] = dev_from_avg_disc_btags;
      variable[5][4] = h2;
      variable[5][5] = lowest_btag;
      variable[5][6] = avg_untagged_dijet_mass;
      variable[5][7] = h3;
      variable[5][8] = second_highest_btag;
      variable[5][9] = invariant_mass_of_everything;

      // 4j4t
      variable[6][0] = first_jet_pt;
      variable[6][1] = h1;
      variable[6][2] = avg_dr_tagged_jets;
      variable[6][3] = aplanarity;
      variable[6][4] = avg_btag_disc_btags;
      variable[6][5] = dev_from_avg_disc_btags;
      variable[6][6] = lowest_btag;
      variable[6][7] = all_sum_pt_with_met;
      variable[6][8] = second_highest_btag;
      variable[6][9] = invariant_mass_of_everything;

      // 5j4t
      variable[7][0] = avg_dr_tagged_jets;
      variable[7][1] = dr_between_lep_and_closest_jet;
      variable[7][2] = avg_btag_disc_btags;
      variable[7][3] = dev_from_avg_disc_btags;
      variable[7][4] = third_jet_pt;
      variable[7][5] = fourth_jet_pt;
      variable[7][6] = lowest_btag;
      variable[7][7] = all_sum_pt_with_met;
      variable[7][8] = first_highest_btag;
      variable[7][9] = second_highest_btag;

      // 6j4t
      variable[8][0] = avg_dr_tagged_jets;
      variable[8][1] = sphericity;
      variable[8][2] = dr_between_lep_and_closest_jet;
      variable[8][3] = avg_btag_disc_btags;
      variable[8][4] = second_highest_btag;
      variable[8][5] = h2;
      variable[8][6] = lowest_btag;
      variable[8][7] = closest_tagged_dijet_mass;
      variable[8][8] = h3;
      variable[8][9] = invariant_mass_of_everything;
      */

      /////////////

      h_wgt_lepSF[iSys]->Fill(eve->wgt_lepSF_);
      h_wgt_PU[iSys]->Fill(eve->wgt_pu_);
      h_wgt_bTag[iSys]->Fill(eve->wgt_btagSF_[0]);

      h_numJet[iSys]->Fill(njet_full,wgt);


      if( numJets>=4 ){
	h_numPV[iSys]->Fill(numPVs,wgt);
	h_numTag[iSys]->Fill(ntag,wgt);
	h_numTag_full[iSys]->Fill(ntag_full,wgt);
      }

      //if( numTags==0 ) h_numJet_0tag[iSys]->Fill(njet_full,wgt);
      //if( numTags==1 ) h_numJet_1tag[iSys]->Fill(njet_full,wgt);
      if( numTags==2 ) h_numJet_2tag[iSys]->Fill(njet_full,wgt);
      if( numTags==3 ) h_numJet_3tag[iSys]->Fill(njet_full,wgt);
      if( numTags>=4 ) h_numJet_4tag[iSys]->Fill(njet_full,wgt);

      if( numTags>=2 ){
	h_numJet_ge2tag[iSys]->Fill(njet_full,wgt);
      }

      h_numJet_numTag[iSys]->Fill(njet,ntag,wgt);


      double HT = eve->HT_[iSys];



      vvdouble jet_vect_TLV = eve->jet_vect_TLV_[iSys];
      vdouble jet_CSV = eve->jet_CSV_[iSys];
      vint jet_flavour = eve->jet_flavour_[iSys];

      vdouble jet_dR2Mean = eve->jet_dR2Mean_[iSys];
      vdouble jet_dRMean  = eve->jet_dRMean_[iSys];
      vdouble jet_frac01  = eve->jet_frac01_[iSys];
      vdouble jet_frac02  = eve->jet_frac02_[iSys];
      vdouble jet_frac03  = eve->jet_frac03_[iSys];
      vdouble jet_beta    = eve->jet_beta_[iSys];
      vdouble jet_betaStar = eve->jet_betaStar_[iSys];
      vdouble jet_leadCandDistFromPV = eve->jet_leadCandDistFromPV_[iSys];



      for( int category=0; category<NumCat; category++ ){
	if( this_category!=category ) continue;

	h_met_pt[category][iSys]->Fill(std::min(double(MET),metmax)-0.00001,wgt);
	h_met_phi[category][iSys]->Fill(MET_phi,wgt);

	h_mht_pt[category][iSys]->Fill(std::min(double(MHT),metmax)-0.00001,wgt);
	h_mht_phi[category][iSys]->Fill(MHT_phi,wgt);

	if(longJob){
	  h_transverse_mass_met[category][iSys]->Fill(std::min(double(MT_met),metmax)-0.00001,wgt);
	  h_transverse_mass_met_vs_met[category][iSys]->Fill(std::min(double(MT_met),metmax)-0.00001,std::min(double(MET),metmax)-0.00001,wgt);
	  
	  h_transverse_mass_mht[category][iSys]->Fill(std::min(double(MT_mht),metmax)-0.00001,wgt);
	  h_transverse_mass_mht_vs_mht[category][iSys]->Fill(std::min(double(MT_mht),metmax)-0.00001,std::min(double(MHT),metmax)-0.00001,wgt);
	}
	
	if(longJob) h_numPV_cat[category][iSys]->Fill(numPVs,wgt);

	//if( leptonChoice==1 ){
	if( leptonChoice==1 ){
	  h_mu_pt[category][iSys]->Fill(std::min(double(tight_lepton_pt),lepPtMax)-0.00001,wgt);
	  h_mu_eta[category][iSys]->Fill(tight_lepton_eta,wgt);
	  h_mu_phi[category][iSys]->Fill(tight_lepton_phi,wgt);
	}
	//else if ( leptonChoice==0) {
	else if ( leptonChoice==0) {
	  h_ele_pt[category][iSys]->Fill(std::min(double(tight_lepton_pt),lepPtMax)-0.00001,wgt);
	  h_ele_eta[category][iSys]->Fill(tight_lepton_eta,wgt);
	  h_ele_phi[category][iSys]->Fill(tight_lepton_phi,wgt);
	}

	h_lepton_pt[category][iSys]->Fill(std::min(double(tight_lepton_pt),lepPtMax)-0.00001,wgt);
	h_lepton_eta[category][iSys]->Fill(tight_lepton_eta,wgt);
	h_lepton_phi[category][iSys]->Fill(tight_lepton_phi,wgt);

	h_HT[category][iSys]->Fill(std::min(double(HT),htmax)-0.00001,wgt);

	if(longJob){
	  h_min_dR_lep_jet[category][iSys]->Fill(dr_between_lep_and_closest_jet,wgt);
	  h_all_sum_pt[category][iSys]->Fill(std::min(double(all_sum_pt),htmax)-0.00001,wgt);
	  h_all_sum_pt_with_met[category][iSys]->Fill(std::min(double(all_sum_pt_with_met),htmax)-0.00001,wgt);
	  h_invariant_mass_of_everything[category][iSys]->Fill(std::min(double(invariant_mass_of_everything+MET),htmax)-0.00001,wgt);

	  double HT30 = HT;
	  h_HT30[category][iSys]->Fill(std::min(HT30,htmax)-0.00001,wgt);
	}

	double sumJetEta = 0;
	double sumTagEta = 0;
	double cntJetEta = 0;
	double cntTagEta = 0;

	int myTags=0, myUnTags=0;
	vecTLorentzVector jetsV;
	for( int iJet=0; iJet<int(jet_vect_TLV.size()); iJet++ ){
	  TLorentzVector myJet;
	  myJet.SetPxPyPzE( jet_vect_TLV[iJet][0], jet_vect_TLV[iJet][1], jet_vect_TLV[iJet][2], jet_vect_TLV[iJet][3] );

	  jetsV.push_back(myJet);
	  
	  double myCSV = jet_CSV[iJet];
	  int myflavour = jet_flavour[iJet];
	  double myJetPT = myJet.Pt();
	  double myJetEta = myJet.Eta();
	  double myjet_leadCandDistFromPV = std::min(jet_leadCandDistFromPV[iJet],max_jet_leadCandDistFromPV)-0.000001;

	  h_jet_eta[category][iSys]->Fill(myJet.Eta(),wgt);
	  h_jet_csv[category][iSys]->Fill(myCSV,wgt);
	  h_jet_flavour[category][iSys]->Fill(myflavour,wgt);
	  if(myflavour==5 || myflavour==-5){
	    h_jet_csvB[category][iSys]->Fill(myCSV,wgt);
	    if(myCSV>0.679) h_tag_csvB[category][iSys]->Fill(myCSV,wgt);
	    else h_untag_csvB[category][iSys]->Fill(myCSV,wgt);
	  }

	  if(myCSV>0.679){
	    h_tag_csv[category][iSys]->Fill(myCSV,wgt);
	    h_tag_flavour[category][iSys]->Fill(myflavour,wgt);
	  }	    
	  else{
	    h_untag_csv[category][iSys]->Fill(myCSV,wgt);
	    h_untag_flavour[category][iSys]->Fill(myflavour,wgt);

	  }

	  if(longJob){

	    h_jet_pt[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    if( myCSV>0.679 ) h_tag_jet_pt[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    else              h_untag_jet_pt[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    
	    
	    h_jet_pt_CSVL_all[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    if( myCSV>0.244 ) h_jet_pt_CSVL_pass[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    else              h_jet_pt_CSVL_fail[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    
	    h_jet_pt_CSVM_all[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    if( myCSV>0.679 ) h_jet_pt_CSVM_pass[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    else              h_jet_pt_CSVM_fail[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    
	    h_jet_pt_CSVT_all[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    if( myCSV>0.898 ) h_jet_pt_CSVT_pass[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    else              h_jet_pt_CSVT_fail[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    
	    if( fabs( myJetEta )<1.4 ){
	      h_jet_pt_eta0to1p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      if( myCSV>0.679 ) h_tag_jet_pt_eta0to1p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      else              h_untag_jet_pt_eta0to1p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    }
	    else{
	      h_jet_pt_eta1p4to2p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      if( myCSV>0.679 ) h_tag_jet_pt_eta1p4to2p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      else              h_untag_jet_pt_eta1p4to2p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    }
	  }

	  if(longJob){
	    
	    for( int iPV=0; iPV<NumPV; iPV++ ){
	      int NPVlow = numPVbins_lowEdge[iPV];
	      int NPVhigh = ( iPV<(NumPV-1) ) ? numPVbins_lowEdge[iPV+1] : 99999;
	      if( numPVs>=NPVlow && numPVs<NPVhigh ){
		h_jet_pt_NPV[iPV][category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
		if( iJet==0 )      h_jet_1_pt_NPV[iPV][category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
		else if( iJet==1 ) h_jet_2_pt_NPV[iPV][category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
		else if( iJet==2 ) h_jet_3_pt_NPV[iPV][category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
		else if( iJet==3 ) h_jet_4_pt_NPV[iPV][category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
		
		if( myCSV>0.679 ){
		  h_tag_jet_pt_NPV[iPV][category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
		  if( myTags==0 )      h_tag_jet_1_pt_NPV[iPV][category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
		  else if( myTags==1 ) h_tag_jet_2_pt_NPV[iPV][category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
		  else if( myTags==2 ) h_tag_jet_3_pt_NPV[iPV][category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
		  else if( myTags==3 ) h_tag_jet_4_pt_NPV[iPV][category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
		}
		else {
		  h_untag_jet_pt_NPV[iPV][category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
		  if( myUnTags==0 )      h_untag_jet_1_pt_NPV[iPV][category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
		  else if( myUnTags==1 ) h_untag_jet_2_pt_NPV[iPV][category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
		  else if( myUnTags==2 ) h_untag_jet_3_pt_NPV[iPV][category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
		  else if( myUnTags==3 ) h_untag_jet_4_pt_NPV[iPV][category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
		}
	      }
	    }
	    
	    if( myJetPT>=30 && myJetPT<40. ){
	      h_jet_pt30to40_dR2Mean[category][iSys]->Fill(jet_dR2Mean[iJet],wgt);
	      h_jet_pt30to40_dRMean[category][iSys]->Fill(jet_dRMean[iJet],wgt);
	      h_jet_pt30to40_frac01[category][iSys]->Fill(jet_frac01[iJet],wgt);
	      h_jet_pt30to40_frac02[category][iSys]->Fill(jet_frac02[iJet],wgt);
	      h_jet_pt30to40_frac03[category][iSys]->Fill(jet_frac03[iJet],wgt);
	      h_jet_pt30to40_beta[category][iSys]->Fill(jet_beta[iJet],wgt);
	      h_jet_pt30to40_betaStar[category][iSys]->Fill(jet_betaStar[iJet],wgt);
	      h_jet_pt30to40_leadCandDistFromPV[category][iSys]->Fill(myjet_leadCandDistFromPV,wgt);
	    }
	  }
	  sumJetEta += myJet.Eta();
	  cntJetEta += 1.;
	  
	  if( myCSV>0.679 ){
	    sumTagEta += myJet.Eta();
	    cntTagEta += 1.;
	  }
	  
	  if( iJet==0 ){
	    h_jet_1_eta[category][iSys]->Fill(myJet.Eta(),wgt);
	    h_jet_1_pt[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    if(longJob){
	      if( fabs( myJet.Eta() ) < 1.4 ) h_jet_1_pt_eta0to1p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      else                            h_jet_1_pt_eta1p4to2p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      h_jet_1_dR2Mean[category][iSys]->Fill(jet_dR2Mean[iJet],wgt);
	      h_jet_1_dRMean[category][iSys]->Fill(jet_dRMean[iJet],wgt);
	      h_jet_1_frac01[category][iSys]->Fill(jet_frac01[iJet],wgt);
	      h_jet_1_frac02[category][iSys]->Fill(jet_frac02[iJet],wgt);
	      h_jet_1_frac03[category][iSys]->Fill(jet_frac03[iJet],wgt);
	      h_jet_1_beta[category][iSys]->Fill(jet_beta[iJet],wgt);
	      h_jet_1_betaStar[category][iSys]->Fill(jet_betaStar[iJet],wgt);
	      h_jet_1_leadCandDistFromPV[category][iSys]->Fill(myjet_leadCandDistFromPV,wgt);
	    }
	    if( myCSV>0.679 ) h_jet_1_tag_pt[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    else              h_jet_1_untag_pt[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	  }
	  else if( iJet==1 ){
	    h_jet_2_eta[category][iSys]->Fill(myJet.Eta(),wgt);
	    h_jet_2_pt[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    if(longJob){
	      if( fabs( myJet.Eta() ) < 1.4 ) h_jet_2_pt_eta0to1p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      else                            h_jet_2_pt_eta1p4to2p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      h_jet_2_dR2Mean[category][iSys]->Fill(jet_dR2Mean[iJet],wgt);
	      h_jet_2_dRMean[category][iSys]->Fill(jet_dRMean[iJet],wgt);
	      h_jet_2_frac01[category][iSys]->Fill(jet_frac01[iJet],wgt);
	      h_jet_2_frac02[category][iSys]->Fill(jet_frac02[iJet],wgt);
	      h_jet_2_frac03[category][iSys]->Fill(jet_frac03[iJet],wgt);
	      h_jet_2_beta[category][iSys]->Fill(jet_beta[iJet],wgt);
	      h_jet_2_betaStar[category][iSys]->Fill(jet_betaStar[iJet],wgt);
	      h_jet_2_leadCandDistFromPV[category][iSys]->Fill(myjet_leadCandDistFromPV,wgt);
	    }
	    if( myCSV>0.679 ) h_jet_2_tag_pt[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    else              h_jet_2_untag_pt[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	  }
	  else if( iJet==2 ){
	    h_jet_3_eta[category][iSys]->Fill(myJet.Eta(),wgt);
	    h_jet_3_pt[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    if(longJob){
	      if( fabs( myJet.Eta() ) < 1.4 ) h_jet_3_pt_eta0to1p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      else                            h_jet_3_pt_eta1p4to2p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      h_jet_3_dR2Mean[category][iSys]->Fill(jet_dR2Mean[iJet],wgt);
	      h_jet_3_dRMean[category][iSys]->Fill(jet_dRMean[iJet],wgt);
	      h_jet_3_frac01[category][iSys]->Fill(jet_frac01[iJet],wgt);
	      h_jet_3_frac02[category][iSys]->Fill(jet_frac02[iJet],wgt);
	      h_jet_3_frac03[category][iSys]->Fill(jet_frac03[iJet],wgt);
	      h_jet_3_beta[category][iSys]->Fill(jet_beta[iJet],wgt);
	      h_jet_3_betaStar[category][iSys]->Fill(jet_betaStar[iJet],wgt);
	      h_jet_3_leadCandDistFromPV[category][iSys]->Fill(myjet_leadCandDistFromPV,wgt);
	    }
	    if( myCSV>0.679 ) h_jet_3_tag_pt[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    else              h_jet_3_untag_pt[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	  }
	  else if( iJet==3 ){
	    h_jet_4_eta[category][iSys]->Fill(myJet.Eta(),wgt);
	    h_jet_4_pt[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    if(longJob){
	      if( fabs( myJet.Eta() ) < 1.4 ) h_jet_4_pt_eta0to1p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      else                            h_jet_4_pt_eta1p4to2p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      h_jet_4_dR2Mean[category][iSys]->Fill(jet_dR2Mean[iJet],wgt);
	      h_jet_4_dRMean[category][iSys]->Fill(jet_dRMean[iJet],wgt);
	      h_jet_4_frac01[category][iSys]->Fill(jet_frac01[iJet],wgt);
	      h_jet_4_frac02[category][iSys]->Fill(jet_frac02[iJet],wgt);
	      h_jet_4_frac03[category][iSys]->Fill(jet_frac03[iJet],wgt);
	      h_jet_4_beta[category][iSys]->Fill(jet_beta[iJet],wgt);
	      h_jet_4_betaStar[category][iSys]->Fill(jet_betaStar[iJet],wgt);
	      h_jet_4_leadCandDistFromPV[category][iSys]->Fill(myjet_leadCandDistFromPV,wgt);
	    }
	    if( myCSV>0.679 ) h_jet_4_tag_pt[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    else              h_jet_4_untag_pt[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	  }

	  if( myCSV>0.679 ){
	    if( myTags==0 ){
	      h_jet_tag_1_pt[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      if(longJob){
		if( fabs( myJet.Eta() ) < 1.4 ) h_tag_1_pt_eta0to1p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
		else                            h_tag_1_pt_eta1p4to2p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      }
	    }
	    else if( myTags==1 ){
	      h_jet_tag_2_pt[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      if(longJob){
		if( fabs( myJet.Eta() ) < 1.4 ) h_tag_2_pt_eta0to1p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
		else                            h_tag_2_pt_eta1p4to2p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      }
	    }
	    else if( myTags==2 ){
	      h_jet_tag_3_pt[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      if(longJob){
		if( fabs( myJet.Eta() ) < 1.4 ) h_tag_3_pt_eta0to1p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
		else                            h_tag_3_pt_eta1p4to2p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      }
	    }
	    else if( myTags==3 ){
	      h_jet_tag_4_pt[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      if(longJob){
		if( fabs( myJet.Eta() ) < 1.4 ) h_tag_4_pt_eta0to1p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
		else                            h_tag_4_pt_eta1p4to2p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      }
	    }
	  }
	  else {
	    if( myUnTags==0 ){
	      h_jet_untag_1_pt[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      if(longJob){
		if( fabs( myJet.Eta() ) < 1.4 ) h_untag_1_pt_eta0to1p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
		else                            h_untag_1_pt_eta1p4to2p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      }
	    }
	    if( myUnTags==1 ){
	      h_jet_untag_2_pt[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      if(longJob){
		if( fabs( myJet.Eta() ) < 1.4 ) h_untag_2_pt_eta0to1p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
		else                            h_untag_2_pt_eta1p4to2p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      }
	    }
	    if( myUnTags==2 ){
	      h_jet_untag_3_pt[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      if(longJob){
		if( fabs( myJet.Eta() ) < 1.4 ) h_untag_3_pt_eta0to1p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
		else                            h_untag_3_pt_eta1p4to2p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      }
	    }
	    if( myUnTags==3 ){
	      h_jet_untag_4_pt[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      if(longJob){
		if( fabs( myJet.Eta() ) < 1.4 ) h_untag_4_pt_eta0to1p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
		else                            h_untag_4_pt_eta1p4to2p4[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      }
	    }
	  }


	  if( myCSV>0.679 ) myTags++;
	  else              myUnTags++;

	} // end loop over jets



	if(longJob){
	  double aveJetEta = ( cntJetEta>0 ) ? sumJetEta/cntJetEta : -999;
	  double aveTagEta = ( cntTagEta>0 ) ? sumTagEta/cntTagEta : -999;
	  
	  double maxDEta_jet_aveJetEta = -1;
	  double maxDEta_tag_aveJetEta = -1;
	  double maxDEta_tag_aveTagEta = -1;
	  
	  for( int iJet=0; iJet<int(jet_vect_TLV.size()); iJet++ ){
	    TLorentzVector myJet;
	    myJet.SetPxPyPzE( jet_vect_TLV[iJet][0], jet_vect_TLV[iJet][1], jet_vect_TLV[iJet][2], jet_vect_TLV[iJet][3] );
	    
	    double myCSV = jet_CSV[iJet];
	    double myJetEta = myJet.Eta();
	    
	    maxDEta_jet_aveJetEta = std::max( maxDEta_jet_aveJetEta, fabs(myJetEta - aveJetEta) );
	    if( myCSV>0.679 ){
	      maxDEta_tag_aveJetEta = std::max( maxDEta_tag_aveJetEta, fabs(myJetEta - aveJetEta) );
	      maxDEta_tag_aveTagEta = std::max( maxDEta_tag_aveTagEta, fabs(myJetEta - aveTagEta) );
	    }
	  }
	  
	  h_maxDEta_jet_aveJetEta[category][iSys]->Fill(maxDEta_jet_aveJetEta,wgt);
	  h_maxDEta_tag_aveJetEta[category][iSys]->Fill(maxDEta_tag_aveJetEta,wgt);
	  h_maxDEta_tag_aveTagEta[category][iSys]->Fill(maxDEta_tag_aveTagEta,wgt);
	  
	  h_jet_aveJetEta[category][iSys]->Fill(aveJetEta,wgt);
	  h_jet_aveTagEta[category][iSys]->Fill(aveTagEta,wgt);
	}
	
	double CFMlpANN_final10v1_8TeV = reader_final10v1_8TeV_CFMLP[category]->EvaluateMVA( "CFMlpANN method" );
	CFMlpANN_final10v1_8TeV = std::min( CFMlpANN_final10v1_8TeV, max_reader_final10v1_8TeV_CFMLP[category]-0.0001 );
	CFMlpANN_final10v1_8TeV = std::max( CFMlpANN_final10v1_8TeV, min_reader_final10v1_8TeV_CFMLP[category]+0.0001 );
	
	h_disc_final10v1_8TeV_CFMlpANN[category][iSys]->Fill(CFMlpANN_final10v1_8TeV,wgt);
	
	
	h_jet_1_csv[category][iSys]->Fill(first_highest_btag,wgt);
	h_jet_2_csv[category][iSys]->Fill(second_highest_btag,wgt);
	h_jet_3_csv[category][iSys]->Fill(third_highest_btag,wgt);
	h_jet_4_csv[category][iSys]->Fill(fourth_highest_btag,wgt);
	

	if(longJob){
	  h_min_dR_tag_tag[category][iSys]->Fill(min_dr_tagged_jets,wgt);
	  h_ave_dR_tag_tag[category][iSys]->Fill(avg_dr_tagged_jets,wgt);
	  h_ave_mass_tag_tag[category][iSys]->Fill(std::min(double(avg_tagged_dijet_mass),massmax)-0.00001,wgt);
	  h_ave_mass_untag_untag[category][iSys]->Fill(std::min(double(avg_untagged_dijet_mass),massmax)-0.00001,wgt);
	  h_ave_mass_untag_untag_reduced[category][iSys]->Fill(double(avg_untagged_dijet_mass),wgt);
	  
	  h_best_higgs_mass[category][iSys]->Fill(std::min(double(best_higgs_mass),massmax)-0.00001,wgt);
	  h_minChi2[category][iSys]->Fill(minChi2,wgt);
	  h_dRbb[category][iSys]->Fill(dRbb,wgt);
	  
	  
	  
	  TLorentzVector mydummyguy;
	  mydummyguy.SetPxPyPzE(0,0,0,0);
	  TLorentzVector lepW = mydummyguy;
	  TLorentzVector hadW = mydummyguy;
	  TLorentzVector lepB = mydummyguy;
	  TLorentzVector hadB = mydummyguy;
	  TLorentzVector lepT = mydummyguy;
	  TLorentzVector hadT = mydummyguy;
	  double minChi2_getTopSystem = 100000000;
	  
	  vdouble lepton_vect = eve->lepton_TLV_[iSys];
	  TLorentzVector leptonV;
	  leptonV.SetPxPyPzE( lepton_vect[0], lepton_vect[1], lepton_vect[2], lepton_vect[3] );
	  
	  TLorentzVector metV;
	  metV.SetPxPyPzE( MET*cos(MET_phi), MET*sin(MET_phi), 0, MET );
	  
	  
	  int result_getTopSystem = getTopSystem( leptonV, metV, jetsV, jet_CSV, minChi2_getTopSystem, hadW, lepW, hadB, lepB, hadT, lepT );
	  
	  TLorentzVector sumT = hadT + lepT;
	  
	  if( result_getTopSystem<0 ) std::cout << " ERROR !! result_getTopSystem = " << result_getTopSystem << std::endl;
	  
	  h_leptonicW_mass[category][iSys]->Fill(std::min(lepW.M(),massmax)-0.00001,wgt);
	  h_hadronicW_mass[category][iSys]->Fill(std::min(hadW.M(),massmax)-0.00001,wgt);
	  h_leptonicT_mass[category][iSys]->Fill(std::min(lepT.M(),massmax)-0.00001,wgt);
	  h_hadronicT_mass[category][iSys]->Fill(std::min(hadT.M(),massmax)-0.00001,wgt);
	  
	  h_hadronicW_pT[category][iSys]->Fill(std::min(hadW.Pt(),massmax)-0.00001,wgt);
	  h_hadronicW_pT_mass[category][iSys]->Fill(std::min(hadW.Pt(),massmax)-0.00001,std::min(hadW.M(),massmax)-0.00001,wgt);
	  
	  h_leptonicT_pT[category][iSys]->Fill(std::min(lepT.Pt(),massmax)-0.00001,wgt);
	  h_leptonicT_pT_mass[category][iSys]->Fill(std::min(lepT.Pt(),massmax)-0.00001,std::min(lepT.M(),massmax)-0.00001,wgt);
	  h_hadronicT_pT[category][iSys]->Fill(std::min(hadT.Pt(),massmax)-0.00001,wgt);
	  h_hadronicT_pT_mass[category][iSys]->Fill(std::min(hadT.Pt(),massmax)-0.00001,std::min(hadT.M(),massmax)-0.00001,wgt);
	  
	  h_lepT_hadT_DeltaR[category][iSys]->Fill(std::min(lepT.DeltaR(hadT),7.0)-0.00001,wgt);
	  h_lepT_hadT_Angle[category][iSys]->Fill(lepT.Angle(hadT.Vect()),wgt);
	  
	  h_sumTop_pT[category][iSys]->Fill(std::min(sumT.Pt(),jetptmax)-0.00001,wgt);
	  h_sumTop_mass[category][iSys]->Fill(std::min(sumT.M(),massmax_sumTop)-0.00001,wgt);
	  h_sumTop_pT_mass[category][iSys]->Fill(std::min(sumT.Pt(),jetptmax)-0.00001,std::min(sumT.M(),massmax_sumTop)-0.00001,wgt);
	  
	  h_minChi2_getTopSystem[category][iSys]->Fill(minChi2_getTopSystem,wgt);
	  
	  h_leptonicW_mass_minChi2_getTopSystem[category][iSys]->Fill(std::min(lepW.M(),massmax)-0.00001,minChi2_getTopSystem,wgt);
	  h_hadronicW_mass_minChi2_getTopSystem[category][iSys]->Fill(std::min(hadW.M(),massmax)-0.00001,minChi2_getTopSystem,wgt);
	  h_leptonicT_mass_minChi2_getTopSystem[category][iSys]->Fill(std::min(lepT.M(),massmax)-0.00001,minChi2_getTopSystem,wgt);
	  h_hadronicT_mass_minChi2_getTopSystem[category][iSys]->Fill(std::min(hadT.M(),massmax)-0.00001,minChi2_getTopSystem,wgt);
	  
	  
	  if( minChi2_getTopSystem>30 ){
	    for( int iJet=0; iJet<int(jetsV.size()); iJet++ ){
	      h_jet_pt_minChi2GT30[category][iSys]->Fill(std::min(jetsV[iJet].Pt(),jetptmax)-0.00001,wgt);
	      if( iJet==0 ) h_jet_1_pt_minChi2GT30[category][iSys]->Fill(std::min(jetsV[iJet].Pt(),jetptmax)-0.00001,wgt);
	    }
	  }

	  if( minChi2_getTopSystem<=30 ){
	    h_leptonicW_mass_minChi2LT30[category][iSys]->Fill(std::min(lepW.M(),massmax)-0.00001,wgt);
	    h_hadronicW_mass_minChi2LT30[category][iSys]->Fill(std::min(hadW.M(),massmax)-0.00001,wgt);
	    h_leptonicT_mass_minChi2LT30[category][iSys]->Fill(std::min(lepT.M(),massmax)-0.00001,wgt);
	    h_hadronicT_mass_minChi2LT30[category][iSys]->Fill(std::min(hadT.M(),massmax)-0.00001,wgt);
	    
	    h_leptonicT_pT_minChi2LT30[category][iSys]->Fill(std::min(lepT.Pt(),massmax)-0.00001,wgt);
	    h_hadronicT_pT_minChi2LT30[category][iSys]->Fill(std::min(hadT.Pt(),massmax)-0.00001,wgt);
	    
	    h_sumTop_pT_minChi2LT30[category][iSys]->Fill(std::min(sumT.Pt(),jetptmax)-0.00001,wgt);
	    h_sumTop_mass_minChi2LT30[category][iSys]->Fill(std::min(sumT.M(),massmax_sumTop)-0.00001,wgt);
	    
	    for( int iJet=0; iJet<int(jetsV.size()); iJet++ ){
	      h_jet_pt_minChi2LT30[category][iSys]->Fill(std::min(jetsV[iJet].Pt(),jetptmax)-0.00001,wgt);
	      if( iJet==0 ) h_jet_1_pt_minChi2LT30[category][iSys]->Fill(std::min(jetsV[iJet].Pt(),jetptmax)-0.00001,wgt);
	      
	      double myCSV = jet_CSV[iJet];
	      double myJetPT = jetsV[iJet].Pt();
	      
	      h_jet_pt_minChi2LT30_CSVL_all[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      if( myCSV>0.244 ) h_jet_pt_minChi2LT30_CSVL_pass[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      else              h_jet_pt_minChi2LT30_CSVL_fail[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      
	      h_jet_pt_minChi2LT30_CSVM_all[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      if( myCSV>0.679 ) h_jet_pt_minChi2LT30_CSVM_pass[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      else              h_jet_pt_minChi2LT30_CSVM_fail[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      
	      h_jet_pt_minChi2LT30_CSVT_all[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      if( myCSV>0.898 ) h_jet_pt_minChi2LT30_CSVT_pass[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      else              h_jet_pt_minChi2LT30_CSVT_fail[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    }
	  }
	  
	  if( minChi2_getTopSystem<=5 ){
	    h_leptonicW_mass_minChi2LT5[category][iSys]->Fill(std::min(lepW.M(),massmax)-0.00001,wgt);
	    h_hadronicW_mass_minChi2LT5[category][iSys]->Fill(std::min(hadW.M(),massmax)-0.00001,wgt);
	    h_leptonicT_mass_minChi2LT5[category][iSys]->Fill(std::min(lepT.M(),massmax)-0.00001,wgt);
	    h_hadronicT_mass_minChi2LT5[category][iSys]->Fill(std::min(hadT.M(),massmax)-0.00001,wgt);
	    
	    h_leptonicT_pT_minChi2LT5[category][iSys]->Fill(std::min(lepT.Pt(),massmax)-0.00001,wgt);
	    h_hadronicT_pT_minChi2LT5[category][iSys]->Fill(std::min(hadT.Pt(),massmax)-0.00001,wgt);
	    
	    h_sumTop_pT_minChi2LT5[category][iSys]->Fill(std::min(sumT.Pt(),jetptmax)-0.00001,wgt);
	    h_sumTop_mass_minChi2LT5[category][iSys]->Fill(std::min(sumT.M(),massmax_sumTop)-0.00001,wgt);
	    
	    for( int iJet=0; iJet<int(jetsV.size()); iJet++ ){
	      h_jet_pt_minChi2LT5[category][iSys]->Fill(std::min(jetsV[iJet].Pt(),jetptmax)-0.00001,wgt);
	      if( iJet==0 ) h_jet_1_pt_minChi2LT5[category][iSys]->Fill(std::min(jetsV[iJet].Pt(),jetptmax)-0.00001,wgt);
	      
	      double myCSV = jet_CSV[iJet];
	      double myJetPT = jetsV[iJet].Pt();
	      
	      h_jet_pt_minChi2LT5_CSVL_all[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      if( myCSV>0.244 ) h_jet_pt_minChi2LT5_CSVL_pass[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      else              h_jet_pt_minChi2LT5_CSVL_fail[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      
	      h_jet_pt_minChi2LT5_CSVM_all[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      if( myCSV>0.679 ) h_jet_pt_minChi2LT5_CSVM_pass[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      else              h_jet_pt_minChi2LT5_CSVM_fail[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      
	      h_jet_pt_minChi2LT5_CSVT_all[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      if( myCSV>0.898 ) h_jet_pt_minChi2LT5_CSVT_pass[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	      else              h_jet_pt_minChi2LT5_CSVT_fail[category][iSys]->Fill(std::min(myJetPT,jetptmax)-0.00001,wgt);
	    }
	  }
	}


	if(longJob){
	  h_M3[category][iSys]->Fill(std::min(double(M3),m3max)-0.00001,wgt);
	  h_M3_1tag[category][iSys]->Fill(std::min(double(M3_1tag),m3max)-0.00001,wgt);
	  h_Mlb[category][iSys]->Fill(std::min(double(Mlb),massmax)-0.00001,wgt);
	}

	h_aplanarity[category][iSys]->Fill(aplanarity,wgt);
	h_sphericity[category][iSys]->Fill(sphericity,wgt);
	h_avg_btag_disc_non_btags[category][iSys]->Fill(avg_btag_disc_non_btags,wgt);
	h_avg_btag_disc_btags[category][iSys]->Fill(avg_btag_disc_btags,wgt);
	h_dev_from_avg_disc_btags[category][iSys]->Fill(dev_from_avg_disc_btags,wgt);
	h_lowest_btag[category][iSys]->Fill(lowest_btag,wgt);
	h_closest_tagged_dijet_mass[category][iSys]->Fill(closest_tagged_dijet_mass,wgt);


	h_h0[category][iSys]->Fill(h0,wgt);
	h_h1[category][iSys]->Fill(h1,wgt);
	h_h2[category][iSys]->Fill(h2,wgt);
	h_h3[category][iSys]->Fill(h3,wgt);
	h_h4[category][iSys]->Fill(h4,wgt);


	if(longJob){
	  h_jet_1_chargedHadronEnergyFraction[category][iSys]->Fill(jet_1_chargedHadronEnergyFraction,wgt);
	  h_jet_1_neutralHadronEnergyFraction[category][iSys]->Fill(jet_1_neutralHadronEnergyFraction,wgt);
	  h_jet_1_chargedEmEnergyFraction[category][iSys]->Fill(jet_1_chargedEmEnergyFraction,wgt);
	  h_jet_1_neutralEmEnergyFraction[category][iSys]->Fill(jet_1_neutralEmEnergyFraction,wgt);
	  h_jet_1_chargedMultiplicity[category][iSys]->Fill(jet_1_chargedMultiplicity,wgt);
	  h_jet_1_neutralMultiplicity[category][iSys]->Fill(jet_1_neutralMultiplicity,wgt);
	  h_jet_1_nconstituents[category][iSys]->Fill(jet_1_nconstituents,wgt);
	  
	  h_jet_2_chargedHadronEnergyFraction[category][iSys]->Fill(jet_2_chargedHadronEnergyFraction,wgt);
	  h_jet_2_neutralHadronEnergyFraction[category][iSys]->Fill(jet_2_neutralHadronEnergyFraction,wgt);
	  h_jet_2_chargedEmEnergyFraction[category][iSys]->Fill(jet_2_chargedEmEnergyFraction,wgt);
	  h_jet_2_neutralEmEnergyFraction[category][iSys]->Fill(jet_2_neutralEmEnergyFraction,wgt);
	  h_jet_2_chargedMultiplicity[category][iSys]->Fill(jet_2_chargedMultiplicity,wgt);
	  h_jet_2_neutralMultiplicity[category][iSys]->Fill(jet_2_neutralMultiplicity,wgt);
	  h_jet_2_nconstituents[category][iSys]->Fill(jet_2_nconstituents,wgt);
	  
	  h_jet_3_chargedHadronEnergyFraction[category][iSys]->Fill(jet_3_chargedHadronEnergyFraction,wgt);
	  h_jet_3_neutralHadronEnergyFraction[category][iSys]->Fill(jet_3_neutralHadronEnergyFraction,wgt);
	  h_jet_3_chargedEmEnergyFraction[category][iSys]->Fill(jet_3_chargedEmEnergyFraction,wgt);
	  h_jet_3_neutralEmEnergyFraction[category][iSys]->Fill(jet_3_neutralEmEnergyFraction,wgt);
	  h_jet_3_chargedMultiplicity[category][iSys]->Fill(jet_3_chargedMultiplicity,wgt);
	  h_jet_3_neutralMultiplicity[category][iSys]->Fill(jet_3_neutralMultiplicity,wgt);
	  h_jet_3_nconstituents[category][iSys]->Fill(jet_3_nconstituents,wgt);
	  
	  h_jet_4_chargedHadronEnergyFraction[category][iSys]->Fill(jet_4_chargedHadronEnergyFraction,wgt);
	  h_jet_4_neutralHadronEnergyFraction[category][iSys]->Fill(jet_4_neutralHadronEnergyFraction,wgt);
	  h_jet_4_chargedEmEnergyFraction[category][iSys]->Fill(jet_4_chargedEmEnergyFraction,wgt);
	  h_jet_4_neutralEmEnergyFraction[category][iSys]->Fill(jet_4_neutralEmEnergyFraction,wgt);
	  h_jet_4_chargedMultiplicity[category][iSys]->Fill(jet_4_chargedMultiplicity,wgt);
	  h_jet_4_neutralMultiplicity[category][iSys]->Fill(jet_4_neutralMultiplicity,wgt);
	  h_jet_4_nconstituents[category][iSys]->Fill(jet_4_nconstituents,wgt);
	}



	if( isMinNPV ){
	  h_disc_final10v1_8TeV_CFMlpANN_minNPV[category][iSys]->Fill(CFMlpANN_final10v1_8TeV,wgt);
	}
	else if( isMidNPV ){
	  h_disc_final10v1_8TeV_CFMlpANN_midNPV[category][iSys]->Fill(CFMlpANN_final10v1_8TeV,wgt);
	}
	else{
	  h_disc_final10v1_8TeV_CFMlpANN_maxNPV[category][iSys]->Fill(CFMlpANN_final10v1_8TeV,wgt);
	}


 

	// /// correlation plots
	// if( iSys==0 ){
	//   int pp=0;
	//   for( int i=0; i<numhists[category]; i++ ){
	//     for( int j=i+1; j<numhists[category]; j++ ){
	//       h_profiles[category][pp++]->Fill(variable[category][i],variable[category][j],wgt);
	//     }
	//   }
	// }

	if( numJets>=6 && numTags>=4 && CFMlpANN_final10v1_8TeV>0.5 && iSys==0 && insample<0 ){
	  std::cout << " *************************************************************** " << std::endl;
	  if( leptonChoice==1 ) std::cout << " \t   muon   + jets event " << std::endl;
	  else if( leptonChoice==0 ) std::cout << " \t electron + jets event " << std::endl;
	  else std::cout << " \t lepton + jets event " << std::endl;
	  std::cout << "\t grab " << run << ":" << lumi << ":" << evt << std::endl;
	  std::cout << "  run = " << run << ",\t lumi = " << lumi << ",\t event = " << evt << std::endl;
	  std::cout << "  number of jets = " << numJets << ",\t number of tags = " << numTags << std::endl;
	  std::cout << "  ANN output (v1 8TeV) = " << CFMlpANN_final10v1_8TeV << std::endl;
	  std::cout << "  lepton: pt = " << tight_lepton_pt << ",\t eta = " << tight_lepton_eta << ",\t phi = " << tight_lepton_phi << std::endl;
	  std::cout << "  min DR (lep,jet) = " << dr_between_lep_and_closest_jet << std::endl;
	  std::cout << "  first:  jet pt = " << first_jet_pt  << ",\t tag jet pt = " << jet_tag_1_pt << ",\t jet csv = " << first_highest_btag   << std::endl;
	  std::cout << "  second: jet pt = " << second_jet_pt << ",\t tag jet pt = " << jet_tag_2_pt << ",\t jet csv = " << second_highest_btag  << std::endl;
	  std::cout << "  third:  jet pt = " << third_jet_pt  << ",\t tag jet pt = " << jet_tag_3_pt << ",\t jet csv = " << third_highest_btag  << std::endl;
	  std::cout << "  fourth: jet pt = " << fourth_jet_pt << ",\t tag jet pt = " << jet_tag_4_pt << ",\t jet csv = " << fourth_highest_btag << std::endl;
	  std::cout << "  lowest csv     = " << lowest_btag << std::endl;
	  std::cout << "  ave CSV (tags) = " << avg_btag_disc_btags << ",\t (untags) = " << avg_btag_disc_non_btags << std::endl;
	  std::cout << "  dev from ave CSV (tags) = " << dev_from_avg_disc_btags << std::endl;
	  std::cout << "  min DR (tag,tag) = " << min_dr_tagged_jets << ",\t ave DR (tag,tag) = " << avg_dr_tagged_jets << std::endl;
	  std::cout << "  MET:    pt = " << MET << ",\t phi = " << MET_phi << std::endl;
	  std::cout << "  MHT:    pt = " << MHT << ",\t phi = " << MHT_phi << std::endl;
	  std::cout << "  Transverse mass (MET) = " << MT_met << ",\t (MHT) = " << MT_mht << std::endl;
	  std::cout << "  Mlb = " << Mlb << ",\t M3 (1 tag) = " << M3_1tag << ",\t M3 = " << M3 << std::endl;
	  std::cout << "  aplanarity = " << aplanarity << ",\t sphericity = " << sphericity << std::endl;
	  std::cout << "  closest mass (tag,tag) = " << closest_tagged_dijet_mass << ",\t ave mass (tag,tag) = " << avg_tagged_dijet_mass << ",\t (untag,untag) = " << avg_untagged_dijet_mass << std::endl;
	  std::cout << "  sum pt (lepton,jets,met) = " << all_sum_pt_with_met << ",\t sum pt (lepton,jets) = " << all_sum_pt << std::endl;
	  std::cout << "  mass (lepton,jets,met) = " << invariant_mass_of_everything << std::endl;
	  std::cout << "  best_higgs_mass = " << best_higgs_mass << std::endl;
	  std::cout << " *************************************************************** " << std::endl;
	}
      } // end loop over categories


    } // end loop over systematics
  } // end loop over events

  std::cout << " Done! " << std::endl;

  histofile.Write();
  histofile.Close();

}



double getTopPtWgt( double topPt, int useSys ){
  if( topPt<0 ) return 1;

  if( topPt>714 ) topPt = 714;

  double result = 1.4e-6 * topPt * topPt - 2.0e-3 * topPt + 1.2;

  double newresult = result;
  if( useSys==1 )      newresult = 2*(result-1) + 1;
  else if( useSys==2 ) newresult = 1.;
  return newresult;
}



int getTopSystem(TLorentzVector lepton, TLorentzVector met, vecTLorentzVector jets, vdouble btag,
		 double &minChi, TLorentzVector &hadW, TLorentzVector &lepW, TLorentzVector &hadB, TLorentzVector &lepB, TLorentzVector &hadT, TLorentzVector &lepT){

  
  int nJets = int(jets.size());

  double chi_top_lep=10000;
  double chi_top_had=10000;
  //double chi_W_lep=10000; //isn't really used
  double chi_W_had=10000;

  minChi = 1000000;
  double btagCut = 0.679;
  double W_mass = 80.4;
  double top_mass = 172.5;

  // from Darren 2/22/2013
  double sigma_hadW   = 10.51;
  double sigma_hadTop = 17.97;
  double sigma_lepTop = 10.02;

  // updated 8/22/2012 from J. Timcheck
  //sigma's from >=6j >=4t, muon, no imaginary neutrino pz ttH
  // double sigma_hadW   = 12.77;
  // double sigma_hadTop = 18.9;
  // double sigma_lepTop = 32.91;

  // //sigma's from >=6j >=4t, muon, no imaginary neutrino pz ttH
  // double sigma_hadW   = 12.59;
  // double sigma_hadTop = 19.9;
  // double sigma_lepTop = 39.05;

  //sigma's from >=6j >=4t, muon, no imaginary neutrino pz ttJets
  //double sigma_hadW		= 12.72,
  //sigma_hadTop	= 18.12,
  //sigma_lepTop	= 38.72;

  double metPz[2];
  double chi=999999;


  int nBtags = 0;
  for(int i=0;i<nJets;i++){
    if(btag[i]>btagCut) nBtags++;
  }

  if( !(nJets>=4 && nBtags>=2) ) return -1;

  int nUntags = nJets-nBtags;

  double lowest_btag = 99.;
  double second_lowest_btag = 999.;
  int ind_lowest_btag = 999;
  int ind_second_lowest_btag = 999;

  //std::cout << "\t nJets = " << nJets << ",\t nbtags = " << btag.size() << ",\t nBtags = " << nBtags << ",\t nUntags = " << nUntags << std::endl;

  if( nUntags<2 ){
    for(int i=0;i<nJets;i++){
      if( btag[i]<lowest_btag ){
	second_lowest_btag = lowest_btag;
	ind_second_lowest_btag = ind_lowest_btag;
	
	lowest_btag = btag[i];
	ind_lowest_btag = i;
      }
      else if( btag[i]<second_lowest_btag ){
	second_lowest_btag = btag[i];
	ind_second_lowest_btag = i;
      }
    }
  }


  // First get the neutrino z
  double energyLep = lepton.E();
  double a = (W_mass*W_mass/(2.0*energyLep)) + (lepton.Px()*met.Px() + lepton.Py()*met.Py())/energyLep;
  double radical = (2.0*lepton.Pz()*a/energyLep)*(2.0*lepton.Pz()*a/energyLep);
  radical = radical - 4.0*(1.0 - (lepton.Pz()/energyLep)*(lepton.Pz()/energyLep))*(met.Px()*met.Px() + met.Py()*met.Py()- a*a);
  if (radical < 0.0) radical = 0.0;
  metPz[0] = (lepton.Pz()*a/energyLep) + 0.5*sqrt(radical);
  metPz[0] = metPz[0] / (1.0 - (lepton.Pz()/energyLep)*(lepton.Pz()/energyLep));
  metPz[1] = (lepton.Pz()*a/energyLep) - 0.5*sqrt(radical);
  metPz[1] = metPz[1] / (1.0 - (lepton.Pz()/energyLep)*(lepton.Pz()/energyLep));


  // Loop over all jets, both Pz, calcaulte chi-square
  TLorentzVector metNew;
  for( int ipznu=0; ipznu<2; ipznu++ ){
    metNew.SetXYZM(met.Px(),met.Py(),metPz[ipznu],0.0); //neutrino has mass 0
    //with b-tag info
    if( nJets>=4 && nBtags>=2 ){
      vecTLorentzVector not_b_tagged,b_tagged;
      //fill not_b_tagged and b_tagged
      for( int i=0;i<nJets;i++ ){
	if( btag[i]>btagCut && i!=ind_second_lowest_btag && i!=ind_lowest_btag) b_tagged.push_back(jets[i]);
	else not_b_tagged.push_back(jets[i]);
      }
      //first make possible t_lep's with b-tagged jets (includes making W_lep)
      for( int i=0; i<int(b_tagged.size()); i++ ){
	TLorentzVector W_lep=metNew+lepton; //used for histogram drawing only
	TLorentzVector top_lep=metNew+lepton+b_tagged.at(i);
	chi_top_lep=pow((top_lep.M()-top_mass)/sigma_lepTop,2);
	//next make possible W_had's with not b-tagged jets
	for( int j=0; j<int(not_b_tagged.size()); j++ ){
	  for( int k=0; k<int(not_b_tagged.size()); k++ ){
	    if( j!=k ){
	      TLorentzVector W_had=not_b_tagged.at(j)+not_b_tagged.at(k);
	      chi_W_had=pow((W_had.M()-W_mass)/sigma_hadW,2);
	      //now make possible top_had's (using the W_had + some b-tagged jet)
	      for( int l=0; l<int(b_tagged.size()); l++ ){
		if( l!=i ){
		  TLorentzVector top_had=W_had+b_tagged.at(l);
		  chi_top_had=pow((top_had.M()-top_mass)/sigma_hadTop,2);
		  chi=chi_top_lep+chi_W_had+chi_top_had;
		  //accept the lowest chi
		  if( chi<minChi ){
		    minChi=chi;
		    hadW = W_had;
		    lepW = W_lep;
		    hadB = b_tagged.at(l);
		    lepB = b_tagged.at(i);
		    hadT = top_had;
		    lepT = top_lep;
		  } // end if chi<minChi
		}
	      }
	    }
	  }
	}
      }
    }
  }

  return 1;
}




/*

  TFile *f_old = new TFile("HistoFiles/results/farm_yggdrasil_treeReader_v52_3rd40_mc_TTJets_TuneZ2_7TeV_madgraph_Fall11_BEANV05_prune_lep_sel_histo_ttbarPlusOther.root");
  TFile *f_new = new TFile("HistoFiles/results/farm_yggdrasil_treeReader_v54_3rd40_mc_TTJets_TuneZ2_7TeV_madgraph_Fall11_BEANV05_wQ2_prune_lep_sel_histo_ttbarPlusOther.root");

  TH1D* h_old = (TH1D*)f_old->Get("h_numJet_numTag");
  TH1D* h_new = (TH1D*)f_new->Get("h_numJet_numTag");

  TH1D* h_rat = (TH1D*)h_new->Clone("h_rat");

  h_rat->Divide(h_new,h_old);

  h_rat->Draw("textcolz");




  TFile *f_old = new TFile("lepJets_result_2012_53Xon52X_jet3Pt40_26Oct2012_v1_topPagCSV.root");
  TFile *f_new = new TFile("ttH_8TeV_split_byCat_manyQ2.root");

  TH1D* h_old = (TH1D*)f_old->Get("ttbar_CFMlpANN_ljets_jge6_tge4");
  TH1D* h_new = (TH1D*)f_new->Get("ttbar_CFMlpANN_ljets_jge6_tge4");

  TH1D* h_rat = (TH1D*)h_new->Clone("h_rat");

  h_rat->Divide(h_new,h_old);

  h_rat->Draw();


 */

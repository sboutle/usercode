#ifndef AnalysisCode_LeptonPlusJets_BEANeventVarsAH_h
#define AnalysisCode_LeptonPlusJets_BEANeventVarsAH_h

//
// Dependencies (#includes)
//
#include <iostream>
#include <vector>
#include "TLorentzVector.h"

#ifdef __MAKECINT__
#pragma link C++ class std::vector< TLorentzVector >+; 
#endif

using namespace std;



typedef std::vector<std::vector<double> > vvdouble;
typedef std::vector<std::vector<std::string> > vvstring;
typedef std::vector<std::vector<int> > vvint;
typedef std::vector<double> vdouble;
typedef std::vector<string> vstring;
typedef std::vector<bool> vbool;
typedef std::vector<int> vint;
typedef std::vector< TLorentzVector > vecTLorentzVector;

//
// Utility Class for Handling Event Variables
//

const Int_t kNumSys = 17;

struct eventVars{


  //////////////////////////////////////////////////////////////////////////
  ///  Tree branches/leaves
  //////////////////////////////////////////////////////////////////////////

  //explicit eventVars(): { }

  int runNum_;
  int lumiNum_;
  long evtNum_;

  double top_pt_;
  double antitop_pt_;

  double Q2ScaleUpWgt_;
  double Q2ScaleDownWgt_;

  Int_t   numPVs_;
  Int_t   numSys_;

  Float_t wgt_lumi_;
  Float_t wgt_xs_;
  Float_t wgt_nGen_;
  Float_t wgt_lepSF_;
  Float_t wgt_pu_;
  Float_t wgt_puUp_;
  Float_t wgt_puDown_;
  Float_t wgt_pu_2012AB_;
  Float_t wgt_puUp_2012AB_;
  Float_t wgt_puDown_2012AB_;
  Float_t wgt_pu_2012ABCD_;
  Float_t wgt_puUp_2012ABCD_;
  Float_t wgt_puDown_2012ABCD_;
  Float_t wgt_topPt_;
  Float_t wgt_topPtUp_;
  Float_t wgt_topPtDown_;
  double  wgt_[kNumSys];  
  double  wgt_btagSF_[kNumSys];  
  
  vvdouble thejets_[kNumSys];
  vdouble CSVdisc_[kNumSys];
  vdouble met_[kNumSys];
  vint jetflavourmatch_[kNumSys];
  vint jet_genId_[kNumSys];
  vint jet_genParentId_[kNumSys];
  vint jet_genGrandParentId_[kNumSys];
  // vstring jet_assignment_from_bhm_[kNumSys]; // leave out for now

  int isSignal_;
  vint topWdecay_; 
  vint antitopWdecay_;
  vint higgsdecay_; 
  vint higgsWplusdecay_;
  vint higgsWminusdecay_;
  vvdouble higgsdecayb_; 

  vvdouble mcparticles_;
  vint mcparticleID_;
  vint mcparticleINDEX_;
  vint mcparticleMOTHER_;
  vvint mcparticleDAUGHTERS_;


  int IsTauTauLeptonEvent_[kNumSys];

  void initialize();

};


typedef std::vector<eventVars> veventVars;


void eventVars::initialize(){

  runNum_  = -99;
  lumiNum_ = -99;
  evtNum_ = -99;

  top_pt_     = -99.9;
  antitop_pt_ = -99.9;

  Q2ScaleUpWgt_   = -99.9;
  Q2ScaleDownWgt_ = -99.9;

  numPVs_ = -99;
  numSys_ = -99;

  wgt_lumi_             = -99.9;
  wgt_xs_               = -99.9;
  wgt_nGen_             = -99.9;
  wgt_lepSF_            = -99.9;
  wgt_pu_               = -99.9;
  wgt_puUp_             = -99.9;
  wgt_puDown_           = -99.9;
  wgt_pu_2012AB_        = -99.9;
  wgt_puUp_2012AB_      = -99.9;
  wgt_puDown_2012AB_    = -99.9;
  wgt_pu_2012ABCD_      = -99.9;
  wgt_puUp_2012ABCD_    = -99.9;
  wgt_puDown_2012ABCD_  = -99.9;;
  wgt_topPt_            = -99.9;
  wgt_topPtUp_          = -99.9;
  wgt_topPtDown_        = -99.9;

  isSignal_ = 0;

  higgsdecay_.clear();
  higgsdecayb_.clear();
  higgsWplusdecay_.clear();
  higgsWminusdecay_.clear();
  topWdecay_.clear();
  antitopWdecay_.clear();

  mcparticles_.clear();
  mcparticleID_.clear();
  mcparticleINDEX_.clear();
  mcparticleMOTHER_.clear();
  mcparticleDAUGHTERS_.clear();
 
  for(int iSys=0; iSys<kNumSys; iSys++){

    wgt_[iSys] = -99.9;  
    wgt_btagSF_[iSys] = -99.9;  

    thejets_[iSys].clear();
    CSVdisc_[iSys].clear();
    jetflavourmatch_[iSys].clear();
    jet_genId_[iSys].clear();
    jet_genParentId_[iSys].clear();
    jet_genGrandParentId_[iSys].clear();
    // jet_assignment_from_bhm_[iSys].clear(); // leave out for now

    IsTauTauLeptonEvent_[iSys] = -99;
  }


  return;
}

  

#endif

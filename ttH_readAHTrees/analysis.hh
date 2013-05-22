#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TFile.h"
#include "TRandom3.h"
#include <fstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <TROOT.h>
#include "THStack.h"
#include "TLegend.h"
#include "eventVars.h"

using namespace std;

class SimpleAnalysis {

public:
  SimpleAnalysis();
  //~SimpleAnalysis();
  void BeginJob(TTree *t1) ;
  void EndJob() ;
  void Loop(TTree *t1) ;
  void Event() ;
  void BranchDefine();

  eventVars *eve_;

private:
  // ->>>>
  TH1F*              h_jets_pt;
  TH1F*              h_njets;
  TH1F*              h_njets_ge6;
  TH1F*              h_ntags_ge6;
  TH1F*              h_avecsv_tagged;
  TH1F*              h_avecsv_untagged;
  TH1F*              h_avecsv_all;
  TH1F*              h_csv_alljets;
  TH1F*              h_csv_Hbjets;
  TH1F*              h_csv_Hbjets_true;
  TH1F*              h_csv_Tbjets;
  TH1F*              h_csv_Tbjets_true;
  TH1F*              h_csv_Wjets;
  TH1F*              h_csv_Wjets_true;
  TH1F*              h_csv_Wcjets_true;

  TH1F*              h_Wdecay_true_8j4t;
  TH1F*              h_Wdecay_true;

  TH1F*              h_besthadWmass;
  TH1F*              h_secondbesthadWmass;
  TH1F*              h_besthadWTopmass;
  TH1F*              h_secondbesthadWTopmass;
  TH1F*              h_finalJetsmass;
  TH1F*              h_finalJetsmass_0Hb;
  TH1F*              h_finalJetsmass_1Hb;
  TH1F*              h_finalJetsmass_2Hb;
  TH1F*              h_finalJetsmass_0Hb99;
  TH1F*              h_finalJetsmass_1Hb99;
  TH1F*              h_finalJetsmass_2Hb99;
  TH1F*              h_finalJetsmass_2Hb99_nomerge;
  TH1F*              h_finalJetsmass_2Hb_nomerge;

  TH1F*              h_HbtoOther_dr_all;
  TH1F*              h_HbtoOther_dr_2Hb99;
  TH1F*              h_HbtoOther_dr_2Hb;

  TH1F*              h_finalJetsmass_0Hbanywhere;
  TH1F*              h_finalJetsmass_1Hbanywhere;
  TH1F*              h_finalJetsmass_2Hbanywhere;
  TH1F*              h_finalJetsmass_2Hbanywhere_allevents;
  TH1F*              h_finalJetsmass_2Hbanywhere_allevents_nomerge;

  TH1F*              h_W1_chiSq;
  TH1F*              h_W2_chiSq;
  TH1F*              h_T1_chiSq;
  TH1F*              h_T2_chiSq;
  TH1F*              h_TOT_chiSq;
  TH1F*              h_allcombo_chiSq;


  TH1F*              h_jet_genID;
  TH1F*              h_jet_genMotherID;

  TH2F*              h_b_destination;
  TH1F*              h_singleb_destination;
  TH2F*              h_higgsdecay;
  TH1F*              h_higgsdecayb_eta;
  TH1F*              h_higgsdecayb_dr;
  TH1F*              h_higgsdecayb_eta_missing;
  TH1F*              h_higgsdecayb_pt_missing;

  // ->>>>


};

SimpleAnalysis::SimpleAnalysis(){
}

vector<TLorentzVector> getAnaObjectVector(vector<vector<double> > myTreeObject)
{
  vector<TLorentzVector> tempVLV;
  TLorentzVector tempLV;
  tempVLV.clear();

  //loop over the vector of vectors and make TLorentzVectors

  vector<vector<double> >::iterator it =myTreeObject.begin();
  for(; it!=myTreeObject.end(); ++it){
  
    vector<double> myPhysOb=*it;
    if(myPhysOb.size()!=4){
      cout << "ERROR: vector not filled correctly. " << endl;
      return tempVLV;
    }
    
    double myX=myPhysOb[0];     
    double myY=myPhysOb[1];     
    double myZ=myPhysOb[2];     
    double myE=myPhysOb[3];     

    if(myE < 8.0){
      //cout<<"WARNING: Object energy is "<<myE<<endl;
      if (myE<0.){
	cout<<"WARNING: Object has negative energy"<<endl;
      }
    }
    tempLV.SetPxPyPzE(myX,myY,myZ,myE);
    tempVLV.push_back(tempLV);

  }

  return tempVLV;
}

TLorentzVector getAnaObject(vector<double> myTreeObject){

  TLorentzVector tempLV;

  if(myTreeObject.size()!=4){
    cout << "ERROR: vector not filled correctly. " << endl;
    return tempLV;
  }
    
  double myX=myTreeObject[0];     
  double myY=myTreeObject[1];     
  double myZ=myTreeObject[2];     
  double myE=myTreeObject[3];     
  
  if (myE<0.){
    cout<<"WARNING: Object has negative energy"<<endl;
  }

  tempLV.SetPxPyPzE(myX,myY,myZ,myE);
  return tempLV;


}

void SimpleAnalysis::BeginJob(TTree *t1){

  cout << "Beginning Job......." << endl;

  eve_ = 0;
  t1->SetBranchAddress("eve.",&eve_);

}
    
void SimpleAnalysis::EndJob(){

  cout << "Ending Job......." << endl;

}

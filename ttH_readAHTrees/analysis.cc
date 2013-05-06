//#ifdef __CINT__
#include <iostream>
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <fstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <TROOT.h>
#include <TChain.h>
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <math.h>
#include "TStyle.h"
#include "analysis.hh"

using namespace std;

int main()
{

  //TFile *f = new TFile("yggdrasil_treeMaker_v3_mc_TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_Summer12_53xOn52x_mu_sel_1.root");
  //TFile *f = new TFile("yggdrasil_treeMaker_mc_TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_53xOn53x_job1of1.root");
  //TString f1 = "/data2/ttH/allHad_trees/yggdrasil_treeMaker_mc_TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_53xOn53x_job1of5.root";
  TString f1 = "/data2/ttH/allHad_trees/yggdrasil_treeMaker_mc_TTH_Inclusive_M_125_8TeV_53xOn53x_job1of1.root";
  TString tree2 = "worldTree";
  TChain *c2 = new TChain(tree2);
  c2->Add(f1);
  
  SimpleAnalysis* myana = new typename SimpleAnalysis::SimpleAnalysis();

  myana->BeginJob(c2);
  myana->Loop(c2);
  myana->EndJob();

}

void SimpleAnalysis::Loop(TTree *t1){

  // Initialise histograms

  h_jets_pt = new TH1F("Jet PT (GeV)", "Jet PT (GeV)", 100, 0.0, 1000.0);
  h_njets = new TH1F("N jets", "N jets", 7, 3.5, 10.5);
  h_jet_genID = new TH1F("Jet genID", "Jet genID", 101, -50.5, 50.5);
  h_jet_genMotherID = new TH1F("Jet genMotherID", "Jet genMotherID", 101, -50.5, 50.5);

  vector<TH1F*> jetpthistos;
  TH1F *h1 = new TH1F("jet1pT", "Jet 1 pT(GeV) for all events",100,0.0,500.0);
  TH1F *h2 = new TH1F("jet2pT", "Jet 2 pT(GeV) for all events",100,0.0,400.0);
  TH1F *h3 = new TH1F("jet3pT", "Jet 3 pT(GeV) for all events",100,0.0,300.0);
  TH1F *h4 = new TH1F("jet4pT", "Jet 4 pT(GeV) for all events",100,0.0,250.0);
  TH1F *h5 = new TH1F("jet5pT", "Jet 5 pT(GeV) for all events",100,0.0,200.0);
  TH1F *h6 = new TH1F("jet6pT", "Jet 6 pT(GeV) for all events",100,0.0,150.0);
  TH1F *h7 = new TH1F("jet7pT", "Jet 7 pT(GeV) for all events",100,0.0,100.0);
  TH1F *h8 = new TH1F("jet8pT", "Jet 8 pT(GeV) for all events",100,0.0,100.0);
  TH1F *h9 = new TH1F("jet9pT", "Jet 9 pT(GeV) for all events",100,0.0,100.0);
  TH1F *h10 = new TH1F("jet10pT", "Jet 10 pT(GeV) for all events",100,0.0,100.0);
  jetpthistos.push_back(h1);
  jetpthistos.push_back(h2);
  jetpthistos.push_back(h3);
  jetpthistos.push_back(h4);
  jetpthistos.push_back(h5);
  jetpthistos.push_back(h6);
  jetpthistos.push_back(h7);
  jetpthistos.push_back(h8);
  jetpthistos.push_back(h9);
  jetpthistos.push_back(h10);

  // Loop over Events

  cout << " in the loop " << endl;

  t1->GetEntry(0);   
  cout << "First event number " << eve_->evtNum_ << endl;
    
  Int_t nentries = (Int_t)t1->GetEntries();

  cout << "Number of events: " << nentries << endl;

  for (Int_t i=0; i<nentries; i++) {    
    t1->GetEntry(i);   

    int njets = 0;    
    int ntags = 0;

    // Get the jets, muons, electrons etc.


    vector<vector<double> > thejets = eve_->thejets_[0];
    vector<double> myDetBTagDiscr   = eve_->CSVdisc_[0];          // B-tag discriminant value with jet index the same as "thejets"
    vector<int> jet_partonID        = eve_->jet_genId_[0];          // B-tag discriminant value with jet index the same as "thejets"
    vector<int> jet_motherID        = eve_->jet_genParentId_[0];          // B-tag discriminant value with jet index the same as "thejets"
    vector<int> flavours            = eve_->jetflavourmatch_[0];  // Flavour of generator parton matched to jet with jet index the same as "thejets"
    vector<double> met              = eve_->met_[0];              // Missing transverse energy 
    int isSignal                    = eve_->isSignal_;            // ==1 if ttH ==0 if ttbar
    int NPV                         = eve_->numPVs_;              // Number of reconstructed primary vertices
    double PUwgt                    = eve_->wgt_pu_;              // Monte Carlo PU weight
    vector<int> higgsdecayproducts  = eve_->higgsdecay_;          // Size 2 vector with decay products of the Higgs: 0:b, 1:c 2:tau, 3:W, 4:gamma, 5:Z, 6:gluon, 7:other
    vector<int> topWdecayproducts   = eve_->topWdecay_;           // Size 2 vector with decay products of the W from top: 0:e, 1:mu, 2:tau, 3:d, 4:u, 5:s, 6:c, 7:b, 8:t, 9: neutrino
    vector<int> antitopWdecayproducts   = eve_->antitopWdecay_;       // Size 2 vector with decay products of the W from top: 0:e, 1:mu, 2:tau, 3:d, 4:u, 5:s, 6:c, 7:b, 8:t, 9: neutrino
    vector<vector<double> > higgsb  = eve_->higgsdecayb_;         // Generator level b's from the Higgs 
    int eventNum                    = eve_->evtNum_;              // Event number
    int runNum                      = eve_->runNum_;              // Run number

    
    vector<vector<double> > mcTrue  = eve_->mcparticles_;
    vector<int> mcTrueID  = eve_->mcparticleID_;
    vector<int> mcTrueINDEX  = eve_->mcparticleINDEX_;
    vector<int> mcTrueMOTHER  = eve_->mcparticleMOTHER_;

    vector<TLorentzVector> mcParticles;
    mcParticles.clear();
    mcParticles = getAnaObjectVector(mcTrue);

    vector<TLorentzVector>::iterator mcit = mcParticles.begin();    
    for(; mcit!=mcParticles.end(); ++mcit){
      int mcindex = mcit - mcParticles.begin(); // Get index
      TLorentzVector myMCparticle=*mcit; 
      int myMCparticleID = mcTrueID[mcindex];
      int myMCparticleINDEX = mcTrueINDEX[mcindex];
      //cout << "index: " << myMCparticleINDEX << " ID: " << myMCparticleID << " PT: " << myMCparticle.Pt() << endl;
      for(int i=0; i<mcTrueINDEX.size();i++){
	if(mcTrueINDEX[i]==mcTrueMOTHER[mcindex]){
	  TLorentzVector myMCparticleMother=mcParticles[i];
	  //cout << "Mother index: " << mcTrueINDEX[i] << " MOTHER ID: " << mcTrueID[i] << " MOTHER PT: " << myMCparticleMother.Pt() << endl;
	}
      } 
      for(int j=0; j < mcTrueINDEX.size(); j++){
	if(mcTrueMOTHER[j]==mcTrueINDEX[mcindex]){	
	  TLorentzVector myMCparticleDaughter=mcParticles[j];
	  //cout << "Daughter index: " << mcTrueINDEX[j] << " DAUGHTER ID: " << mcTrueID[j] << " DAUGHTER PT: " << myMCparticleDaughter.Pt() << endl;
	}
      }
    }
    
    vector<vector<double> > nomJets = eve_->thejets_[0];
    vector<TLorentzVector> detJets;
    detJets.clear();
    detJets = getAnaObjectVector(thejets);
    vector<TLorentzVector> detBtags;
    detBtags.clear();
    detBtags = getAnaObjectVector(higgsb);
    TLorentzVector myDetMET = getAnaObject(met);

    vector<TLorentzVector>::iterator it =detJets.begin();    
    for(; it!=detJets.end(); ++it){
      int jetindex = it - detJets.begin(); // Get index for current jet
      TLorentzVector myDetJet=*it; 

      if ((myDetJet.Pt()>55) || (myDetJet.Pt()>30 && njets>4)){
	njets++; 
	h_jets_pt->Fill(myDetJet.Pt());
	if(jetindex<10) jetpthistos[jetindex]->Fill(myDetJet.Pt());
	h_jet_genID->Fill(jet_partonID[jetindex]);
	h_jet_genMotherID->Fill(jet_motherID[jetindex]);
      }
    }

  
    // Fill the event-level histograms
    h_njets->Fill(njets);
 }   

  // Write the histograms to a file
  
  TFile* outfile = new TFile("histograms_ttH.root","RECREATE");
  h_jets_pt->Write();
  h_njets->Write();
  h_jet_genID->Write();
  h_jet_genMotherID->Write();
  for (unsigned int i = 0; i < jetpthistos.size(); ++i) jetpthistos[i]->Write();
 
  return;
}

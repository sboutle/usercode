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
  TString f1 = "/data2/ttH/allHad_trees/yggdrasil_treeMaker_mc_TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_53xOn53x_job*of5.root";
  //TString f1 = "/data2/ttH/allHad_trees/yggdrasil_treeMaker_mc_TTH_Inclusive_M_125_8TeV_53xOn53x_job1of1.root";
  //TString f1 = "/data2/ttH/allHad_trees/trees_25GeV/yggdrasil_treeMaker_mc_TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_53xOn53x_job*of5.root";

  //TString f1 = "/data2/ttH/allHad_trees/trees_25GeV/yggdrasil_treeMaker_mc_TTH_Inclusive_M_125_8TeV_53xOn53x_job1of1.root";
  //TString f1 = "/data2/ttH/allHad_trees/TTJets_100eventtest.root";
  TString tree2 = "worldTree";
  TChain *c2 = new TChain(tree2);
  c2->Add(f1);
  
  SimpleAnalysis* myana = new typename SimpleAnalysis::SimpleAnalysis();

  myana->BeginJob(c2);
  myana->Loop(c2);
  myana->EndJob();

}

void SimpleAnalysis::Loop(TTree *t1){

  bool isLinear = false;
  bool kinfit = true;


  // Initialise histograms

  h_jets_pt = new TH1F("Jet PT (GeV)", "Jet PT (GeV)", 100, 0.0, 1000.0);
  h_njets = new TH1F("N jets", "N jets", 9, 3.5, 12.5);
  h_njets_ge6 = new TH1F("N jets (>=6)", "N jets (>=6)", 7, 5.5, 12.5);
  h_ntags_ge6 = new TH1F("N tags (>=6)", "N tags (>=6)", 6, 1.5, 7.5);

  h_avecsv_tagged = new TH1F("Ave. CSV discr. (tagged jets)", "Ave. CSV discr. (tagged jets)", 40, 0.6, 1.0);
  h_avecsv_untagged = new TH1F("Ave. CSV discr. (untagged jets)", "Ave. CSV discr. (untagged jets)", 70, 0.0, 0.7);
  h_avecsv_all = new TH1F("Ave. CSV discr. (all jets)", "Ave. CSV discr. (all jets)", 100, 0.0, 1.0);

  h_csv_alljets = new TH1F("CSV discr. (all jets)", "CSV discr. (all jets)", 100, 0.0, 1.0);
  h_csv_Hbjets = new TH1F("CSV discr. (Jets assigned to Higgs)", "CSV discr. (Jets assigned to Higgs)", 100, 0.0, 1.0);
  h_csv_Hbjets_true = new TH1F("CSV discr. (Higgs b jets)", "CSV discr. (Higgs b jets)", 100, 0.0, 1.0);
  h_csv_Tbjets = new TH1F("CSV discr. (Jets assigned to Top b)", "CSV discr. (Jets assigned to Top b)", 100, 0.0, 1.0);
  h_csv_Tbjets_true = new TH1F("CSV discr. (Top b jets)", "CSV discr. (Top b jets)", 100, 0.0, 1.0);
  h_csv_Wjets = new TH1F("CSV discr. (Jets assigned to W)", "CSV discr. (Jets assigned to W)", 100, 0.0, 1.0);
  h_csv_Wjets_true = new TH1F("CSV discr. (W decay jets)", "CSV discr. (W decay jets)", 100, 0.0, 1.0);
  h_csv_Wcjets_true = new TH1F("CSV discr. (W decay c jets)", "CSV discr. (W decay c jets)", 100, 0.0, 1.0);

  h_besthadWmass = new TH1F("Best had W mass", "Best had W mass", 200, 0.0, 200.0);
  h_secondbesthadWmass = new TH1F("Second Best had W mass", "Second best had W mass", 200, 0.0, 200.0);

  h_besthadWTopmass = new TH1F("Best had W Top mass", "Best had W Top mass", 400, 0.0, 400.0);
  h_secondbesthadWTopmass = new TH1F("Second Best had W Top mass", "Second best had W Top mass", 400, 0.0, 400.0);

  h_finalJetsmass = new TH1F("Mass of remaining 2 jets", "Mass of remaining 2 jets", 400, 0.0, 400.0);
  h_finalJetsmass_1Hb = new TH1F("Mass of remaining 2 jets for ==1 matched to b", "Mass of remaining 2 jets for ==1 matched to b", 400, 0.0, 400.0);
  h_finalJetsmass_2Hb = new TH1F("Mass of remaining 2 jets for ==2 matched to b", "Mass of remaining 2 jets for ==2 matched to b", 400, 0.0, 400.0);
  h_finalJetsmass_0Hb = new TH1F("Mass of remaining 2 jets for ==0 matched to b", "Mass of remaining 2 jets for ==0 matched to b", 400, 0.0, 400.0);


  h_finalJetsmass_0Hbanywhere = new TH1F("Mass of remaining 2 jets for ==0 bs from Higgs in the event", "Mass of remaining 2 jets for ==0 bs from Higgs in the event", 400, 0.0, 400.0);
  h_finalJetsmass_1Hbanywhere = new TH1F("Mass of remaining 2 jets for ==1 bs from Higgs in the event", "Mass of remaining 2 jets for ==1 bs from Higgs in the event", 400, 0.0, 400.0);

  h_finalJetsmass_2Hbanywhere = new TH1F("Mass of remaining 2 jets for ==2 bs from Higgs in the event", "Mass of remaining 2 jets for ==2 bs from Higgs in the event", 400, 0.0, 400.0);

  h_b_destination = new TH2F("Destination of the Higgs b", "Destination of the Higgs b", 4, 0.5, 4.5, 4, 0.5, 4.5);
  h_singleb_destination = new TH1F("Destination of the single Higgs b", "Destination of the single Higgs b", 4, 0.5, 4.5);

  h_higgsdecay = new TH2F("Higgs_decay_products", "Higgs_decay_products", 8, -0.5, 7.5, 8, -0.5, 7.5);
  h_higgsdecayb_eta = new TH1F("Higgs decay b eta", "Higgs decay b eta", 48, -3.0, 3.0);
  h_higgsdecayb_dr = new TH1F("Higgs decay b DeltaR", "Higgs decay b DeltaR", 48, 0.0, 6.0);
  h_higgsdecayb_eta_missing = new TH1F("Higgs decay b eta for events with missing Higgs jets", "Higgs decay b eta for events with missing Higgs jets", 48, -3.0, 3.0);
  h_higgsdecayb_pt_missing = new TH1F("Higgs decay b pt for events with missing Higgs jets", "Higgs decay b pt for events with missing Higgs jets", 100, 0.0, 1000.0);

  h_W1_chiSq = new TH1F("W1 chisq", "W1 chisq", 100, 0.0, 50.0);
  h_W2_chiSq = new TH1F("W2 chisq", "W2 chisq", 100, 0.0, 50.0);
  h_T1_chiSq = new TH1F("T1 chisq", "T1 chisq", 100, 0.0, 50.0);
  h_T2_chiSq = new TH1F("T2 chisq", "T2 chisq", 100, 0.0, 50.0);
  h_TOT_chiSq = new TH1F("TOT chisq", "TOT chisq", 100, 0.0, 200.0);
  h_allcombo_chiSq = new TH1F("all possible chisq", "all possible chisq", 100, 0.0, 2000.0);

  h_jet_genID = new TH1F("Jet genID", "Jet genID", 101, -50.5, 50.5);
  h_jet_genMotherID = new TH1F("Jet genMotherID", "Jet genMotherID", 101, -50.5, 50.5);

  h_Wdecay_true_8j4t = new TH1F("True W decay 8j4t", "True W decay 8j4t", 41, -20.5, 20.5);
  h_Wdecay_true = new TH1F("True W decay", "True W decay", 41, -20.5, 20.5);


  vector<TH1F*> jetpthistos;
  TH1F *pT_h1 = new TH1F("jet1pT", "Jet 1 pT(GeV) for all events",100,0.0,500.0);
  TH1F *pT_h2 = new TH1F("jet2pT", "Jet 2 pT(GeV) for all events",100,0.0,400.0);
  TH1F *pT_h3 = new TH1F("jet3pT", "Jet 3 pT(GeV) for all events",100,0.0,300.0);
  TH1F *pT_h4 = new TH1F("jet4pT", "Jet 4 pT(GeV) for all events",100,0.0,250.0);
  TH1F *pT_h5 = new TH1F("jet5pT", "Jet 5 pT(GeV) for all events",100,0.0,200.0);
  TH1F *pT_h6 = new TH1F("jet6pT", "Jet 6 pT(GeV) for all events",100,0.0,150.0);
  TH1F *pT_h7 = new TH1F("jet7pT", "Jet 7 pT(GeV) for all events",100,0.0,100.0);
  TH1F *pT_h8 = new TH1F("jet8pT", "Jet 8 pT(GeV) for all events",100,0.0,100.0);
  TH1F *pT_h9 = new TH1F("jet9pT", "Jet 9 pT(GeV) for all events",100,0.0,100.0);
  TH1F *pT_h10 = new TH1F("jet10pT", "Jet 10 pT(GeV) for all events",100,0.0,100.0);
  jetpthistos.push_back(pT_h1);
  jetpthistos.push_back(pT_h2);
  jetpthistos.push_back(pT_h3);
  jetpthistos.push_back(pT_h4);
  jetpthistos.push_back(pT_h5);
  jetpthistos.push_back(pT_h6);
  jetpthistos.push_back(pT_h7);
  jetpthistos.push_back(pT_h8);
  jetpthistos.push_back(pT_h9);
  jetpthistos.push_back(pT_h10);

  vector<TH1F*> jetetahistos;
  TH1F *eta_h1 = new TH1F("jet1eta", "Jet 1 eta(GeV) for all events",30,-3.0,3.0);
  TH1F *eta_h2 = new TH1F("jet2eta", "Jet 2 eta(GeV) for all events",30,-3.0,3.0);
  TH1F *eta_h3 = new TH1F("jet3eta", "Jet 3 eta(GeV) for all events",30,-3.0,3.0);
  TH1F *eta_h4 = new TH1F("jet4eta", "Jet 4 eta(GeV) for all events",30,-3.0,3.0);
  TH1F *eta_h5 = new TH1F("jet5eta", "Jet 5 eta(GeV) for all events",30,-3.0,3.0);
  TH1F *eta_h6 = new TH1F("jet6eta", "Jet 6 eta(GeV) for all events",30,-3.0,3.0);
  TH1F *eta_h7 = new TH1F("jet7eta", "Jet 7 eta(GeV) for all events",30,-3.0,3.0);
  TH1F *eta_h8 = new TH1F("jet8eta", "Jet 8 eta(GeV) for all events",30,-3.0,3.0);
  TH1F *eta_h9 = new TH1F("jet9eta", "Jet 9 eta(GeV) for all events",30,-3.0,3.0);
  TH1F *eta_h10 = new TH1F("jet10eta", "Jet 10 eta(GeV) for all events",30,-3.0,3.0);
  jetetahistos.push_back(eta_h1);
  jetetahistos.push_back(eta_h2);
  jetetahistos.push_back(eta_h3);
  jetetahistos.push_back(eta_h4);
  jetetahistos.push_back(eta_h5);
  jetetahistos.push_back(eta_h6);
  jetetahistos.push_back(eta_h7);
  jetetahistos.push_back(eta_h8);
  jetetahistos.push_back(eta_h9);
  jetetahistos.push_back(eta_h10);

  vector<TH1F*> jetcsvhistos;
  TH1F *csv_h1 = new TH1F("jet1csv", "Jet 1 csv for all events",100,0.0,1.0);
  TH1F *csv_h2 = new TH1F("jet2csv", "Jet 2 csv for all events",100,0.0,1.0);
  TH1F *csv_h3 = new TH1F("jet3csv", "Jet 3 csv for all events",100,0.0,1.0);
  TH1F *csv_h4 = new TH1F("jet4csv", "Jet 4 csv for all events",100,0.0,1.0);
  TH1F *csv_h5 = new TH1F("jet5csv", "Jet 5 csv for all events",100,0.0,1.0);
  TH1F *csv_h6 = new TH1F("jet6csv", "Jet 6 csv for all events",100,0.0,1.0);
  TH1F *csv_h7 = new TH1F("jet7csv", "Jet 7 csv for all events",100,0.0,1.0);
  TH1F *csv_h8 = new TH1F("jet8csv", "Jet 8 csv for all events",100,0.0,1.0);
  TH1F *csv_h9 = new TH1F("jet9csv", "Jet 9 csv for all events",100,0.0,1.0);
  TH1F *csv_h10 = new TH1F("jet10csv", "Jet 10 csv for all events",100,0.0,1.0);
  jetcsvhistos.push_back(csv_h1);
  jetcsvhistos.push_back(csv_h2);
  jetcsvhistos.push_back(csv_h3);
  jetcsvhistos.push_back(csv_h4);
  jetcsvhistos.push_back(csv_h5);
  jetcsvhistos.push_back(csv_h6);
  jetcsvhistos.push_back(csv_h7);
  jetcsvhistos.push_back(csv_h8);
  jetcsvhistos.push_back(csv_h9);
  jetcsvhistos.push_back(csv_h10);

  vector<TH1F*> jetcsvorderedhistos;
  TH1F *csvordered_h1 = new TH1F("jet1csvordered", "Jet 1 csvordered for all events",100,0.0,1.0);
  TH1F *csvordered_h2 = new TH1F("jet2csvordered", "Jet 2 csvordered for all events",100,0.0,1.0);
  TH1F *csvordered_h3 = new TH1F("jet3csvordered", "Jet 3 csvordered for all events",100,0.0,1.0);
  TH1F *csvordered_h4 = new TH1F("jet4csvordered", "Jet 4 csvordered for all events",100,0.0,1.0);
  TH1F *csvordered_h5 = new TH1F("jet5csvordered", "Jet 5 csvordered for all events",100,0.0,1.0);
  TH1F *csvordered_h6 = new TH1F("jet6csvordered", "Jet 6 csvordered for all events",100,0.0,1.0);
  TH1F *csvordered_h7 = new TH1F("jet7csvordered", "Jet 7 csvordered for all events",100,0.0,1.0);
  TH1F *csvordered_h8 = new TH1F("jet8csvordered", "Jet 8 csvordered for all events",100,0.0,1.0);
  TH1F *csvordered_h9 = new TH1F("jet9csvordered", "Jet 9 csvordered for all events",100,0.0,1.0);
  TH1F *csvordered_h10 = new TH1F("jet10csvordered", "Jet 10 csvordered for all events",100,0.0,1.0);
  jetcsvorderedhistos.push_back(csvordered_h1);
  jetcsvorderedhistos.push_back(csvordered_h2);
  jetcsvorderedhistos.push_back(csvordered_h3);
  jetcsvorderedhistos.push_back(csvordered_h4);
  jetcsvorderedhistos.push_back(csvordered_h5);
  jetcsvorderedhistos.push_back(csvordered_h6);
  jetcsvorderedhistos.push_back(csvordered_h7);
  jetcsvorderedhistos.push_back(csvordered_h8);
  jetcsvorderedhistos.push_back(csvordered_h9);
  jetcsvorderedhistos.push_back(csvordered_h10);

  // Loop over Events

  cout << " in the loop " << endl;

  t1->GetEntry(0);   
  cout << "First event number " << eve_->evtNum_ << endl;
    
  Int_t nentries = (Int_t)t1->GetEntries();

  cout << "Number of events: " << nentries << endl;

  int WProb_counter = 0;

  for (Int_t i=0; i<nentries; i++) {    
    //for (Int_t i=0; i<100; i++) {    
    t1->GetEntry(i);   

    int njets = 0;    
    int ntags = 0;

    // Get the jets, muons, electrons etc.


    vector<vector<double> > thejets = eve_->thejets_[0];
    vector<vector<double> > jesUpJets = eve_->thejets_[11];
    vector<vector<double> > jesDownJets = eve_->thejets_[12];

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

    
  

    vector<TLorentzVector> detJets;
    detJets.clear();
    detJets = getAnaObjectVector(thejets);

    vector<TLorentzVector> JESUpJets;
    JESUpJets.clear();
    JESUpJets = getAnaObjectVector(jesUpJets);

    vector<TLorentzVector> JESDownJets;
    JESDownJets.clear();
    JESDownJets = getAnaObjectVector(jesDownJets);

    vector<TLorentzVector> detBs;
    detBs.clear();
    detBs = getAnaObjectVector(higgsb);

    TLorentzVector myDetMET = getAnaObject(met);

    double sum_discr_tagged = 0.0;
    double sum_discr_untagged = 0.0;

    vector<TLorentzVector>::iterator it =detJets.begin();    

    for(; it!=detJets.end(); ++it){
      int jetindex = it - detJets.begin(); // Get index for current jet
      TLorentzVector myDetJet=*it; 

      if ((myDetJet.Pt()>55) || (myDetJet.Pt()>30 && njets>3)){
	njets++; 
	h_jets_pt->Fill(myDetJet.Pt());

	if(jetindex<10){
	  jetpthistos[jetindex]->Fill(myDetJet.Pt());
	  jetetahistos[jetindex]->Fill(myDetJet.Eta());
	  jetcsvhistos[jetindex]->Fill(myDetBTagDiscr[jetindex]);
	}

	h_jet_genID->Fill(jet_partonID[jetindex]);
	h_jet_genMotherID->Fill(jet_motherID[jetindex]);

	if(myDetBTagDiscr[jetindex]>=0.679){
	  ntags++;
	  sum_discr_tagged = sum_discr_tagged + myDetBTagDiscr[jetindex];
	}
	else{
	  sum_discr_untagged = sum_discr_untagged + myDetBTagDiscr[jetindex];
	}
      }
    }


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
	    if(fabs(mcTrueID[i]) == 24) h_Wdecay_true->Fill(myMCparticleID);
	    if(fabs(mcTrueID[i]) == 24 && njets==8 && ntags>=4) h_Wdecay_true_8j4t->Fill(myMCparticleID);

	  }
	} 
	/*
	for(int j=0; j < mcTrueINDEX.size(); j++){
	  if(mcTrueMOTHER[j]==mcTrueINDEX[mcindex]){	
	  TLorentzVector myMCparticleDaughter=mcParticles[j];
	  //cout << "Daughter index: " << mcTrueINDEX[j] << " DAUGHTER ID: " << mcTrueID[j] << " DAUGHTER PT: " << myMCparticleDaughter.Pt() << endl;
	  }
	}
	*/
      }

    vector<int> btag_order;
    btag_order.clear();
    vector<TLorentzVector>::iterator it2 =detJets.begin();    
    vector<TLorentzVector>::iterator it22 =detJets.begin();    
    vector<TLorentzVector>::iterator it3 =detJets.begin();    
    vector<TLorentzVector>::iterator it32 =detJets.begin();    
    vector<TLorentzVector>::iterator it4 =detJets.begin();    
    vector<TLorentzVector>::iterator it5 =detJets.begin();    
    vector<TLorentzVector>::iterator it6 =detJets.begin();    

    //if((njets==8 || njets==9) && ntags>=4){
    if(njets==8 && ntags>=4){
      h_avecsv_tagged->Fill(sum_discr_tagged/ntags);
      h_avecsv_untagged->Fill(sum_discr_untagged/(njets-ntags));
      h_avecsv_all->Fill((sum_discr_tagged+sum_discr_untagged)/njets);

      if(isSignal){
	h_higgsdecay->Fill(higgsdecayproducts[0],higgsdecayproducts[1]);
	//cout << "NEW EVENT =======================================================" << endl;
	
	
	if(higgsdecayproducts[0]==0 &&higgsdecayproducts[1]==0 ){
	  h_higgsdecayb_eta->Fill(detBs[0].Eta());
	  h_higgsdecayb_eta->Fill(detBs[1].Eta());
	  
	  h_higgsdecayb_dr->Fill(detBs[0].DeltaR(detBs[1]));
	}
      }


      if(isLinear==false){



	for(; it2!=detJets.end(); ++it2){
	  int jetindex2 = it2 - detJets.begin(); // Get index for current jet 
	  TLorentzVector myDetJet2=*it2;
	  for(it3 = it2+1; it3!=detJets.end(); ++it3){
	    int jetindex3 = it3 - detJets.begin(); // Get index for current jet
	    
	  }
	}
	
	int combo_counter = 0;
	double chisq = 9998.;
	TLorentzVector hadW1 = 9999.;
	TLorentzVector hadW2 = 9999.;
	TLorentzVector Top1 = 9999.;
	TLorentzVector Top2 = 9999.;
	TLorentzVector leftover1 = 9999.;
	TLorentzVector leftover2 = 9999.;
	TLorentzVector leftover3 = 9999.;
	TLorentzVector hadW1_JESUp = 9999.;
	TLorentzVector hadW2_JESUp = 9999.;
	TLorentzVector Top1_JESUp = 9999.;
	TLorentzVector Top2_JESUp = 9999.;
	TLorentzVector leftover1_JESUp = 9999.;
	TLorentzVector leftover2_JESUp = 9999.;
	TLorentzVector leftover3_JESUp = 9999.;
	TLorentzVector hadW1_JESDown = 9999.;
	TLorentzVector hadW2_JESDown = 9999.;
	TLorentzVector Top1_JESDown = 9999.;
	TLorentzVector Top2_JESDown = 9999.;
	TLorentzVector leftover1_JESDown = 9999.;
	TLorentzVector leftover2_JESDown = 9999.;
	TLorentzVector leftover3_JESDown = 9999.;
	int chosen_jet1 = 99;
	int chosen_jet2 = 99;
	int chosen_jet3 = 99;
	int W_jet_11 = 99;
	int W_jet_12 = 99;
	int W_jet_21 = 99;
	int W_jet_22 = 99;
	int T_jet_1 = 99;
	int T_jet_2 = 99;

	
	if(njets==8 && kinfit){
	  for(int i = 0; i < njets; i++){
	    for(int j = i+1; j<njets; j++){
	      if(i != j){
		//cout << "i " << i << " j " << j << endl;
		bool foundk = false;
	      
		for(int k = 0; k<njets; k++){
		  if(k!=i && k!=j && foundk==false){
		    foundk = true;

		    for(int l = k+1; l<njets; l++ ){
		      if(l!=i && l!=j){
			for(int m = l+1; m<njets; m++){
			  if(m!=i && m!=j){
			    for(int n=0; n<njets; n++){
			      if(n!=i && n!=j && n!=k && n!=l && n!=m){
				for(int o=n+1; o<njets; o++){
				  if(o!=i && o!=j && o!=k && o!=l && o!=m){
				    for(int p=o+1; p<njets; p++){
				      if(p!=i && p!=j && p!=k && p!=l && p!=m){
					for(int s = 0; s<9; s++){
				     
					  int j1 = 99;
					  int j2 = 99;
					  int j3 = 99;
					  int j4 = 99;
					  int j5 = 99;
					  int j6 = 99;

					  if(s==0) { j1 = k; j2 = l; j3 = m; j4 = n; j5 = o; j6 = p;}
					  if(s==1) { j1 = k; j2 = l; j3 = m; j4 = o; j5 = p; j6 = n;}
					  if(s==2) { j1 = k; j2 = l; j3 = m; j4 = p; j5 = n; j6 = o;}
					  if(s==3) { j1 = l; j2 = m; j3 = k; j4 = n; j5 = o; j6 = p;}
					  if(s==4) { j1 = l; j2 = m; j3 = k; j4 = o; j5 = p; j6 = n;}
					  if(s==5) { j1 = l; j2 = m; j3 = k; j4 = p; j5 = n; j6 = o;}
					  if(s==6) { j1 = m; j2 = k; j3 = l; j4 = n; j5 = o; j6 = p;}
					  if(s==7) { j1 = m; j2 = k; j3 = l; j4 = o; j5 = p; j6 = n;}
					  if(s==8) { j1 = m; j2 = k; j3 = l; j4 = p; j5 = n; j6 = o;}
				      
					  double W1unc = TMath::Max(fabs((detJets[j1]+detJets[j2]).M()-(JESUpJets[j1]+JESUpJets[j2]).M()),fabs((detJets[j1]+detJets[j2]).M()-(JESDownJets[j1]+JESDownJets[j2]).M()));
					  double chisqW1 = ((detJets[j1]+detJets[j2]).M()-80.1)*((detJets[j1]+detJets[j2]).M()-80.1)/(W1unc*W1unc);
					  double W2unc = TMath::Max(fabs((detJets[j4]+detJets[j5]).M()-(JESUpJets[j4]+JESUpJets[j5]).M()),fabs((detJets[j4]+detJets[j5]).M()-(JESDownJets[j4]+JESDownJets[j5]).M()));
					  double chisqW2 = ((detJets[j4]+detJets[j5]).M()-80.1)*((detJets[j4]+detJets[j5]).M()-80.1)/(W2unc*W2unc);

					  double T1unc = TMath::Max(fabs((detJets[j1]+detJets[j2]+detJets[j3]).M()-(JESUpJets[j1]+JESUpJets[j2]+JESUpJets[j3]).M()),fabs((detJets[j1]+detJets[j2]+detJets[j3]).M()-(JESDownJets[j1]+JESDownJets[j2]+JESDownJets[j3]).M()));
					  double chisqT1 = ((detJets[j1]+detJets[j2]+detJets[j3]).M()-173.)*((detJets[j1]+detJets[j2]+detJets[j3]).M()-173.)/(T1unc*T1unc);
					  double T2unc = TMath::Max(fabs((detJets[j4]+detJets[j5]+detJets[j6]).M()-(JESUpJets[j4]+JESUpJets[j5]+JESUpJets[j6]).M()),fabs((detJets[j4]+detJets[j5]+detJets[j6]).M()-(JESDownJets[j4]+JESDownJets[j5]+JESDownJets[j6]).M()));
					  double chisqT2 = ((detJets[j4]+detJets[j5]+detJets[j6]).M()-173.)*((detJets[j4]+detJets[j5]+detJets[j6]).M()-173.)/(T2unc*T2unc);
					
					  double chisq_temp = chisqW1 + chisqW2 + chisqT1 + chisqT2;
					  h_allcombo_chiSq->Fill(chisq_temp);
					  if(myDetBTagDiscr[j1]<0.95 && myDetBTagDiscr[j2]<0.95 && myDetBTagDiscr[j4]<0.95 && myDetBTagDiscr[j5]<0.95 && myDetBTagDiscr[j3]>0.679 && myDetBTagDiscr[j4]>0.679){
					  //if(myDetBTagDiscr[j3]>0.679 && myDetBTagDiscr[j6]>0.679){
					    if(chisq_temp < chisq ){
					      chisq = chisq_temp;
					      if(chisqW1<chisqW2){
						hadW1 = detJets[j1]+detJets[j2];
						hadW2 = detJets[j4]+detJets[j5];
						W_jet_11 = j1;
						W_jet_12 = j2;
						W_jet_21 = j4;
						W_jet_22 = j5;
					      }
					      else{
						hadW1 = detJets[j4]+detJets[j5];
						hadW2 = detJets[j1]+detJets[j2];
						W_jet_11 = j4;
						W_jet_12 = j5;
						W_jet_21 = j1;
						W_jet_22 = j2;
					      }
					      if(chisqT1<chisqT1){
						Top1 = detJets[j1]+detJets[j2]+detJets[j3];
						Top2 = detJets[j4]+detJets[j5]+detJets[j6];
						T_jet_1 = j3;
						T_jet_2 = j6;
					      }
					      else{
						Top1 = detJets[j4]+detJets[j5]+detJets[j6];
						Top2 = detJets[j1]+detJets[j2]+detJets[j3];
						T_jet_1 = j6;
						T_jet_2 = j3;
					      }
					      leftover1 = detJets[i];
					      leftover2 = detJets[j];
					      chosen_jet1 = i;
					      chosen_jet2 = j;
					    }
					  }
					}
				      }
				    }
				  }
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}


	if(njets==9 && kinfit){
	  for(int h = 0; h < njets; h++){
	    for(int i = h+1; i < njets; i++){
	      for(int j = i+1; j<njets; j++){
		if(i != j && j!=h && i!=h){
		  //cout << "i " << i << " j " << j << endl;
		  bool foundk = false;
	      
		  for(int k = 0; k<njets; k++){
		    if(k!=i && k!=j && k!=h && foundk==false){
		      foundk = true;

		      for(int l = k+1; l<njets; l++ ){
			if(l!=i && l!=j && l!=h){
			  for(int m = l+1; m<njets; m++){
			    if(m!=i && m!=j && m!=h){
			      for(int n=0; n<njets; n++){
				if(n!=h && n!=i && n!=j && n!=k && n!=l && n!=m){
				  for(int o=n+1; o<njets; o++){
				    if(o!=h && o!=i && o!=j && o!=k && o!=l && o!=m){
				      for(int p=o+1; p<njets; p++){
					if(p!=h && p!=i && p!=j && p!=k && p!=l && p!=m){
					  for(int s = 0; s<9; s++){
				     
					    int j1 = 99;
					    int j2 = 99;
					    int j3 = 99;
					    int j4 = 99;
					    int j5 = 99;
					    int j6 = 99;

					    if(s==0) { j1 = k; j2 = l; j3 = m; j4 = n; j5 = o; j6 = p;}
					    if(s==1) { j1 = k; j2 = l; j3 = m; j4 = o; j5 = p; j6 = n;}
					    if(s==2) { j1 = k; j2 = l; j3 = m; j4 = p; j5 = n; j6 = o;}
					    if(s==3) { j1 = l; j2 = m; j3 = k; j4 = n; j5 = o; j6 = p;}
					    if(s==4) { j1 = l; j2 = m; j3 = k; j4 = o; j5 = p; j6 = n;}
					    if(s==5) { j1 = l; j2 = m; j3 = k; j4 = p; j5 = n; j6 = o;}
					    if(s==6) { j1 = m; j2 = k; j3 = l; j4 = n; j5 = o; j6 = p;}
					    if(s==7) { j1 = m; j2 = k; j3 = l; j4 = o; j5 = p; j6 = n;}
					    if(s==8) { j1 = m; j2 = k; j3 = l; j4 = p; j5 = n; j6 = o;}
					  
				      
					    double W1unc = TMath::Max(fabs((detJets[j1]+detJets[j2]).M()-(JESUpJets[j1]+JESUpJets[j2]).M()),fabs((detJets[j1]+detJets[j2]).M()-(JESDownJets[j1]+JESDownJets[j2]).M()));
					    double chisqW1 = ((detJets[j1]+detJets[j2]).M()-80.1)*((detJets[j1]+detJets[j2]).M()-80.1)/W1unc;
					    double W2unc = TMath::Max(fabs((detJets[j4]+detJets[j5]).M()-(JESUpJets[j4]+JESUpJets[j5]).M()),fabs((detJets[j4]+detJets[j5]).M()-(JESDownJets[j4]+JESDownJets[j5]).M()));
					    double chisqW2 = ((detJets[j4]+detJets[j5]).M()-80.1)*((detJets[j4]+detJets[j5]).M()-80.1)/W2unc;
					    
					    double T1unc = TMath::Max(fabs((detJets[j1]+detJets[j2]+detJets[j3]).M()-(JESUpJets[j1]+JESUpJets[j2]+JESUpJets[j3]).M()),fabs((detJets[j1]+detJets[j2]+detJets[j3]).M()-(JESDownJets[j1]+JESDownJets[j2]+JESDownJets[j3]).M()));
					    double chisqT1 = ((detJets[j1]+detJets[j2]+detJets[j3]).M()-173.)*((detJets[j1]+detJets[j2]+detJets[j3]).M()-173.)/T1unc;
					    double T2unc = TMath::Max(fabs((detJets[j4]+detJets[j5]+detJets[j6]).M()-(JESUpJets[j4]+JESUpJets[j5]+JESUpJets[j6]).M()),fabs((detJets[j4]+detJets[j5]+detJets[j6]).M()-(JESDownJets[j4]+JESDownJets[j5]+JESDownJets[j6]).M()));
					    double chisqT2 = ((detJets[j4]+detJets[j5]+detJets[j6]).M()-173.)*((detJets[j4]+detJets[j5]+detJets[j6]).M()-173.)/T2unc;
					    double chisq_temp = chisqW1 + chisqW2 + chisqT1 + chisqT2;
					    combo_counter++;
					    h_allcombo_chiSq->Fill(chisq_temp);
					  
					    if(chisq_temp < chisq){
					      chisq = chisq_temp;
					      if(chisqW1<chisqW2){
						hadW1 = detJets[j1]+detJets[j2];
						hadW2 = detJets[j4]+detJets[j5];
						W_jet_11 = j1;
						W_jet_12 = j2;
						W_jet_21 = j4;
						W_jet_22 = j5;
					      }
					      else{
						hadW1 = detJets[j4]+detJets[j5];
						hadW2 = detJets[j1]+detJets[j2];
						W_jet_11 = j4;
						W_jet_12 = j5;
						W_jet_21 = j1;
						W_jet_22 = j2;
					      }
					      if(chisqT1<chisqT1){
						Top1 = detJets[j1]+detJets[j2]+detJets[j3];
						Top2 = detJets[j4]+detJets[j5]+detJets[j6];
						T_jet_1 = j3;
						T_jet_2 = j6;
					      }
					      else{
						Top1 = detJets[j4]+detJets[j5]+detJets[j6];
						Top2 = detJets[j1]+detJets[j2]+detJets[j3];
						T_jet_1 = j6;
						T_jet_2 = j3;
					      }

					      leftover1 = detJets[h];
					      leftover2 = detJets[i];
					      leftover3 = detJets[j];
					      chosen_jet1 = h;
					      chosen_jet2 = i;
					      chosen_jet3 = j;
					    }
					  }
					}
				      }
				    }
				  }
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}














	
	//cout << "number of combinations = " << combo_counter << endl;
	h_besthadWmass->Fill(hadW1.M());
	h_secondbesthadWmass->Fill(hadW2.M());
	h_besthadWTopmass->Fill(Top1.M());
	h_secondbesthadWTopmass->Fill(Top2.M());
	h_TOT_chiSq->Fill(chisq);
	
	if(isSignal){
	  if(higgsdecayproducts[0]==0 &&higgsdecayproducts[1]==0){
	    h_finalJetsmass->Fill((leftover1+leftover2).M());
	    if(jet_motherID[chosen_jet1]!=25 && jet_motherID[chosen_jet1]!=-99 && jet_motherID[chosen_jet2]!=25 && jet_motherID[chosen_jet2]!=-99){
	      //if(jet_motherID[chosen_jet1]!=25 && jet_motherID[chosen_jet2]!=25){
	      h_finalJetsmass_0Hb->Fill((leftover1+leftover2).M());
	    }
	    if((jet_motherID[chosen_jet1]==25 || jet_motherID[chosen_jet1]==-99) ^ (jet_motherID[chosen_jet2]==25 || jet_motherID[chosen_jet2]==-99)){
	      //if(jet_motherID[chosen_jet1]==25 ^ jet_motherID[chosen_jet2]==25){
	      h_finalJetsmass_1Hb->Fill((leftover1+leftover2).M());
	      if(jet_motherID[W_jet_11]==25 || jet_motherID[W_jet_12]==25 ) h_singleb_destination->Fill(1);
	      if(jet_motherID[W_jet_21]==25 || jet_motherID[W_jet_22]==25 ) h_singleb_destination->Fill(2);
	      if(jet_motherID[T_jet_1]==25) h_singleb_destination->Fill(3);
	      if(jet_motherID[T_jet_2]==25) h_singleb_destination->Fill(4);
	    }
	    if((jet_motherID[chosen_jet1]==25 || jet_motherID[chosen_jet1]==-99) && (jet_motherID[chosen_jet2]==25 || jet_motherID[chosen_jet2]==-99)){
	      //if(jet_motherID[chosen_jet1]==25 && jet_motherID[chosen_jet2]==25){
	      h_finalJetsmass_2Hb->Fill((leftover1+leftover2).M());
	    }
	    
	    
	    if(jet_motherID[W_jet_11]==25 && jet_motherID[W_jet_12]==25) h_b_destination->Fill(1,1);
	    
	    if(jet_motherID[W_jet_11]==25 && jet_motherID[W_jet_21]==25) h_b_destination->Fill(1,2);
	    if(jet_motherID[W_jet_12]==25 && jet_motherID[W_jet_21]==25) h_b_destination->Fill(1,2);
	    if(jet_motherID[W_jet_12]==25 && jet_motherID[W_jet_22]==25) h_b_destination->Fill(1,2);
	    if(jet_motherID[W_jet_11]==25 && jet_motherID[W_jet_22]==25) h_b_destination->Fill(1,2);
	    
	    if(jet_motherID[W_jet_21]==25 && jet_motherID[W_jet_22]==25) h_b_destination->Fill(2,2);
	    
	    if(jet_motherID[W_jet_11]==25 && jet_motherID[T_jet_1]==25) h_b_destination->Fill(1,3);
	    if(jet_motherID[W_jet_12]==25 && jet_motherID[T_jet_1]==25) h_b_destination->Fill(1,3);
	    
	    if(jet_motherID[W_jet_11]==25 && jet_motherID[T_jet_2]==25) h_b_destination->Fill(1,4);
	    if(jet_motherID[W_jet_12]==25 && jet_motherID[T_jet_2]==25) h_b_destination->Fill(1,4);
	    
	    if(jet_motherID[W_jet_21]==25 && jet_motherID[T_jet_1]==25) h_b_destination->Fill(2,3);
	    if(jet_motherID[W_jet_22]==25 && jet_motherID[T_jet_1]==25) h_b_destination->Fill(2,3);
	    
	    if(jet_motherID[W_jet_21]==25 && jet_motherID[T_jet_2]==25) h_b_destination->Fill(2,4);
	    if(jet_motherID[W_jet_22]==25 && jet_motherID[T_jet_2]==25) h_b_destination->Fill(2,4);
	    
	    if(jet_motherID[T_jet_1]==25 && jet_motherID[T_jet_2]==25) h_b_destination->Fill(3,4);
	    
	    int Hbcounter = 0;
	    
	    for(int ijet = 0; ijet<njets; ijet++){
	      h_csv_alljets->Fill(myDetBTagDiscr[ijet]);
	      if(jet_motherID[ijet]==25){
		Hbcounter++;
		h_csv_Hbjets_true->Fill(myDetBTagDiscr[ijet]);
	      }
	      
	      if(ijet==chosen_jet1 || ijet==chosen_jet2) h_csv_Hbjets->Fill(myDetBTagDiscr[ijet]);
	      //if(ijet==chosen_jet1 || ijet==chosen_jet2) cout << "ijet " << ijet << " jet 1 " << chosen_jet1 << " jet 2 " << chosen_jet2 << endl;
	      if(fabs(jet_motherID[ijet])==24) h_csv_Wjets_true->Fill(myDetBTagDiscr[ijet]);
	      if(fabs(jet_motherID[ijet])==24 && abs(jet_partonID[ijet])==4) h_csv_Wcjets_true->Fill(myDetBTagDiscr[ijet]);
	      if(ijet==W_jet_11 || ijet==W_jet_12 || ijet==W_jet_21 || ijet==W_jet_22) h_csv_Wjets->Fill(myDetBTagDiscr[ijet]);
	      
	      if(fabs(jet_motherID[ijet])==6 && fabs(jet_partonID[ijet])==5) h_csv_Tbjets_true->Fill(myDetBTagDiscr[ijet]);
	      if(ijet==T_jet_1 || ijet==T_jet_2 ) h_csv_Tbjets->Fill(myDetBTagDiscr[ijet]);
	      
	    }
	    
	    if(Hbcounter==0){
	      h_finalJetsmass_0Hbanywhere->Fill((leftover1+leftover2).M());
	      h_higgsdecayb_eta_missing->Fill(detBs[0].Eta());
	      h_higgsdecayb_eta_missing->Fill(detBs[1].Eta());
	      h_higgsdecayb_pt_missing->Fill(detBs[0].Pt());
	      h_higgsdecayb_pt_missing->Fill(detBs[1].Pt());
	      
	    }
	    if(Hbcounter==1){
	      h_finalJetsmass_1Hbanywhere->Fill((leftover1+leftover2).M());
	      
	    /*
	    cout << "EVENT: " << eventNum << " ===============================" << endl;
	    cout << "b quark 1: pt = " << detBs[0].Pt() << " eta = " << detBs[0].Eta() << endl;
	    cout << "b quark 2: pt = " << detBs[1].Pt() << " eta = " << detBs[1].Eta() << endl;

	  
	    for(int ijet = 0; ijet<8; ijet++){
	      cout << "jet " << ijet << ": pt = " << detJets[ijet].Pt() << " eta = " << detJets[ijet].Eta() << " partonID = " << jet_partonID[ijet] << " mother ID = " << jet_motherID[ijet] << " deltaR to b1 = " << detJets[ijet].DeltaR(detBs[0]) << " deltaR to b2 = " << detJets[ijet].DeltaR(detBs[1]) << endl;
	    }

	    cout << "Jets used in the mass calculation: " << chosen_jet1 << " " << chosen_jet2 << endl;
	    cout << endl;
	    */
	    }
	    if(Hbcounter==2){
	      h_finalJetsmass_2Hbanywhere->Fill((leftover1+leftover2).M());
	    }
	  }
	  
	}
	else{
	  h_finalJetsmass->Fill((leftover1+leftover2).M());
	  

	  for(int ijet = 0; ijet<njets; ijet++){
	    h_csv_alljets->Fill(myDetBTagDiscr[ijet]);
	    if(jet_motherID[ijet]==25){
	      h_csv_Hbjets_true->Fill(myDetBTagDiscr[ijet]);
	    }
	    if(ijet==chosen_jet1 || ijet==chosen_jet2) h_csv_Hbjets->Fill(myDetBTagDiscr[ijet]);
	    //if(ijet==chosen_jet1 || ijet==chosen_jet2) 
	    //cout << "ijet " << ijet << " jet 1 " << chosen_jet1 << " jet 2 " << chosen_jet2 << endl;
	    if(fabs(jet_motherID[ijet])==24) h_csv_Wjets_true->Fill(myDetBTagDiscr[ijet]);
	    if(fabs(jet_motherID[ijet])==24 && abs(jet_partonID[ijet])==4) h_csv_Wcjets_true->Fill(myDetBTagDiscr[ijet]);
	    if(ijet==W_jet_11 || ijet==W_jet_12 || ijet==W_jet_21 || ijet==W_jet_22) h_csv_Wjets->Fill(myDetBTagDiscr[ijet]);
	    if(fabs(jet_motherID[ijet])==6 && fabs(jet_partonID[ijet])==5) h_csv_Tbjets_true->Fill(myDetBTagDiscr[ijet]);
	    if(ijet==T_jet_1 || ijet==T_jet_2 ) h_csv_Tbjets->Fill(myDetBTagDiscr[ijet]);
	    
	    
	  }


	}

	  
    

      }

      if(isLinear){
	double hadW1_chi2 = 9998.;
	double hadW2_chi2 = 9999.;
	TLorentzVector hadW1 = 9999.;
	TLorentzVector hadW2 = 9999.;
	TLorentzVector hadW1_JESUp = 9999.;
	TLorentzVector hadW2_JESUp = 9999.;
	TLorentzVector hadW1_JESDown = 9999.;
	TLorentzVector hadW2_JESDown = 9999.;
	
	int W1J1 = 99.;
	int W1J2 = 99.;
	int W2J1 = 99.;
	int W2J2 = 99.;
	
	
	
	for(; it2!=detJets.end(); ++it2){
	  int jetindex2 = it2 - detJets.begin(); // Get index for current jet
	  
	  TLorentzVector myDetJet2=*it2;
	  TLorentzVector JESUpJet2= JESUpJets[jetindex2];
	  TLorentzVector JESDownJet2= JESDownJets[jetindex2];
	  
	  if ((myDetJet2.Pt()>55) || (myDetJet2.Pt()>30 && njets>3)){
	    int nhigher = 0;
	    for(int jtag=0; jtag<njets && jtag<10;jtag++){
	      if(myDetBTagDiscr[jetindex2]<myDetBTagDiscr[jtag])
		nhigher++;
	    }
	    btag_order.push_back(nhigher);
	    
	    for(it3 = it2+1; it3!=detJets.end(); ++it3){
	      int jetindex3 = it3 - detJets.begin(); // Get index for current jet
	      
	      TLorentzVector myDetJet3=*it3;
	      TLorentzVector JESUpJet3= JESUpJets[jetindex3];
	      TLorentzVector JESDownJet3= JESDownJets[jetindex3];
	      
	      
	      if ((myDetJet3.Pt()>55) || (myDetJet3.Pt()>30 && njets>3)){
		
		if(jetindex2>njets || jetindex3>njets) cout << "TOO MANY JET ITERATIONS" << endl;
		
		TLorentzVector hadW1_temp = myDetJet2+myDetJet3;
		TLorentzVector hadW1_JESUp_temp = JESUpJet2+JESUpJet3;
		TLorentzVector hadW1_JESDown_temp = JESDownJet2+JESDownJet3;
		
		double dijetM_UpDiff = fabs(hadW1_JESUp_temp.M() - hadW1_temp.M());
		double dijetM_DownDiff = fabs(hadW1_JESDown_temp.M() - hadW1_temp.M());
		double dijetM_unc = TMath::Max(dijetM_UpDiff,dijetM_DownDiff);
		
		double hadW1_chi2_temp = fabs(hadW1_temp.M()-80.1)/dijetM_unc; 
		
		//cout << "i = " << jetindex2 << " j = " << jetindex3 << " mass = " <<  dijetM << " chisq = " << chi2_hadW << endl;

		//cout << "dijetM_secondlowestchi2 " << dijetM_secondlowestchi2 << " dijetM_lowestchi2 " << dijetM_lowestchi2 << endl; 	      
		//cout << "secondlowestchi2 " << secondlowestchi2 << " lowestchi2 " << lowestchi2 << endl;
	      
		if(hadW1_chi2_temp < hadW1_chi2){
		  hadW1 = hadW1_temp;
		  hadW1_JESUp = hadW1_JESUp_temp;
		  hadW1_JESDown = hadW1_JESDown_temp;
		  hadW1_chi2 = hadW1_chi2_temp;
		  W1J1 = jetindex2;
		  W1J2 = jetindex3;
		}
	      }
	    }// end third loop on jets
	  }
	}//end second loop on jets
   
	for(; it22!=detJets.end(); ++it22){
	  int jetindex22 = it22 - detJets.begin(); // Get index for current jet
	  if(jetindex22 != W1J1 && jetindex22 != W1J2 ){

	    TLorentzVector myDetJet22=*it22;
	    TLorentzVector JESUpJet22= JESUpJets[jetindex22];
	    TLorentzVector JESDownJet22= JESDownJets[jetindex22];
	  
	    if ((myDetJet22.Pt()>55) || (myDetJet22.Pt()>30 && njets>3)){
	   	    
	      for(it32 = it22+1; it32!=detJets.end(); ++it32){
		int jetindex32 = it32 - detJets.begin(); // Get index for current jet
		if(jetindex32 != W1J1 && jetindex32 != W1J2 ){
	    
		  TLorentzVector myDetJet32=*it32;
		  TLorentzVector JESUpJet32= JESUpJets[jetindex32];
		  TLorentzVector JESDownJet32= JESDownJets[jetindex32];

		  if ((myDetJet32.Pt()>55) || (myDetJet32.Pt()>30 && njets>3)){
		 		  
		    if(jetindex22>njets || jetindex32>njets) cout << "TOO MANY JET ITERATIONS" << endl;
		  
		    TLorentzVector hadW2_temp = myDetJet22+myDetJet32;
		    TLorentzVector hadW2_JESUp_temp = JESUpJet22+JESUpJet32;
		    TLorentzVector hadW2_JESDown_temp = JESDownJet22+JESDownJet32;
		  
		    double dijetM_UpDiff = fabs(hadW2_JESUp_temp.M() - hadW2_temp.M());
		    double dijetM_DownDiff = fabs(hadW2_JESDown_temp.M() - hadW2_temp.M());
		    double dijetM_unc = TMath::Max(dijetM_UpDiff,dijetM_DownDiff);
		  
		    double hadW2_chi2_temp = fabs(hadW2_temp.M()-80.1)/dijetM_unc; 
		  
		    //cout << "i = " << jetindex2 << " j = " << jetindex3 << " mass = " <<  dijetM << " chisq = " << chi2_hadW << endl;
		  
		    //cout << "dijetM_secondlowestchi2 " << dijetM_secondlowestchi2 << " dijetM_lowestchi2 " << dijetM_lowestchi2 << endl; 	      
		    //cout << "secondlowestchi2 " << secondlowestchi2 << " lowestchi2 " << lowestchi2 << endl;
		  
		    if(hadW2_chi2_temp < hadW2_chi2){
		      hadW2 = hadW2_temp;
		      hadW2_JESUp = hadW2_JESUp_temp;
		      hadW2_JESDown = hadW2_JESDown_temp;
		      hadW2_chi2 = hadW2_chi2_temp;
		      W2J1 = jetindex22;
		      W2J2 = jetindex32;
		    }
		  }
		}
	      }// end third loop on jets
	    }
	  }
	}//end second loop on jets
   
	if(W1J1==W2J1 || W1J1==W2J2 || W1J2==W2J1 || W1J2==W2J2 ) WProb_counter++;


	h_besthadWmass->Fill(hadW1.M());
	h_secondbesthadWmass->Fill(hadW2.M());

	//cout << "WINNER: " << dijetM_lowestchi2 << " runner-up: " <<  dijetM_secondlowestchi2 << endl;

	double Top1_chi2 = 9999.0;
	double Top2_chi2 = 9999.0;
	TLorentzVector Top1 = 9999.;
	TLorentzVector Top2 = 9999.;
	TLorentzVector Top1_JESUp = 9999.;
	TLorentzVector Top2_JESUp = 9999.;
	TLorentzVector Top1_JESDown = 9999.;
	TLorentzVector Top2_JESDown = 9999.;
 
	int W1T1 = 99;
	int W2T2 = 99;

	for(; it4!=detJets.end(); ++it4){
	  int jetindex4 = it4 - detJets.begin(); // Get index for current jet
	  TLorentzVector myDetJet4=*it4;
	  TLorentzVector JESUpJet4= JESUpJets[jetindex4];
	  TLorentzVector JESDownJet4= JESDownJets[jetindex4];

	  if ((myDetJet4.Pt()>55) || (myDetJet4.Pt()>30 && njets>3)){
	    if(jetindex4!=W1J1 && jetindex4!=W1J2 && jetindex4!=W2J1 && jetindex4!=W2J2){

	      TLorentzVector Top1_temp = myDetJet4 + hadW1;
	      TLorentzVector Top1_JESUp_temp = JESUpJet4 + hadW1_JESUp;
	      TLorentzVector Top1_JESDown_temp = JESDownJet4 + hadW1_JESDown;
	    
	      double topmass_UpDiff = fabs(Top1_JESUp_temp.M() - Top1_temp.M());
	      double topmass_DownDiff = fabs(Top1_JESDown_temp.M() - Top1_temp.M());
	      double topmass_unc = TMath::Max(topmass_UpDiff,topmass_DownDiff);
		  
	      double Top1_chi2_temp = fabs(Top1_temp.M()-173.1)/topmass_unc; // SB Change me 

	      if(Top1_chi2_temp < Top1_chi2){
		Top1 = Top1_temp;
		Top1_JESUp = Top1_JESUp_temp;
		Top1_JESDown = Top1_JESDown_temp;
		Top1_chi2 = Top1_chi2_temp;
		W1T1 = jetindex4;
	      }
	    }
	  }
	}
	for(; it5!=detJets.end(); ++it5){
	  int jetindex5 = it5 - detJets.begin(); // Get index for current jet
	  TLorentzVector myDetJet5=*it5;
	  TLorentzVector JESUpJet5= JESUpJets[jetindex5];
	  TLorentzVector JESDownJet5= JESDownJets[jetindex5];
	  if ((myDetJet5.Pt()>55) || (myDetJet5.Pt()>30 && njets>3)){
	    if(jetindex5!=W1J1 && jetindex5!=W1J2 && jetindex5!=W2J1 && jetindex5!=W2J2 && jetindex5!=W1T1){

	      TLorentzVector Top2_temp = myDetJet5 + hadW2;
	      TLorentzVector Top2_JESUp_temp = JESUpJet5 + hadW2_JESUp;
	      TLorentzVector Top2_JESDown_temp = JESDownJet5 + hadW2_JESDown;
	    
	      double topmass_UpDiff = fabs(Top2_JESUp_temp.M() - Top2_temp.M());
	      double topmass_DownDiff = fabs(Top2_JESDown_temp.M() - Top2_temp.M());
	      double topmass_unc = TMath::Max(topmass_UpDiff,topmass_DownDiff);
		  
	      double Top2_chi2_temp = fabs(Top2_temp.M()-173.1)/topmass_unc; // SB Change me 

	      if(Top2_chi2_temp < Top2_chi2){
		Top2 = Top2_temp;
		Top2_JESUp = Top2_JESUp_temp;
		Top2_JESDown = Top2_JESDown_temp;
		Top2_chi2 = Top2_chi2_temp;
		W2T2 = jetindex5;
	      }
	    }
	  }
	}
	h_besthadWTopmass->Fill(Top1.M());
	h_secondbesthadWTopmass->Fill(Top2.M());

	int num_left = 0;
	TLorentzVector finalJet1 = 9999.;
	TLorentzVector finalJet2 = 9999.;
	TLorentzVector finalJets = 9999.;

	for(; it6!=detJets.end(); ++it6){
	  int jetindex6 = it6 - detJets.begin(); // Get index for current jet
	  TLorentzVector myDetJet6=*it6;

	  if(jetindex6!=W1J1 && jetindex6!=W1J2 && jetindex6!=W2J1 && jetindex6!=W2J2 && jetindex6!=W1T1 && jetindex6!=W2T2)num_left++;
	  if(num_left>2) cout << "TOO MANY JETS LEFT AFTER MAKING TOPS" << endl;

	  if(num_left==1) finalJet1 = myDetJet6;
	  if(num_left==2) finalJet2 = myDetJet6;
	}

	finalJets = finalJet1 + finalJet2;

	if((hadW1_chi2+hadW2_chi2+Top1_chi2+Top2_chi2)<15.){
	  h_finalJetsmass->Fill(finalJets.M());
	}
	h_W1_chiSq->Fill(hadW1_chi2);
	h_W2_chiSq->Fill(hadW2_chi2);
	h_T1_chiSq->Fill(Top1_chi2);
	h_T2_chiSq->Fill(Top2_chi2);
	h_TOT_chiSq->Fill(hadW1_chi2+hadW2_chi2+Top1_chi2+Top2_chi2);
      } 

    }//end 8j4t selection

    for(int itag=0; btag_order.size()>0 && itag<btag_order.size();itag++){
      int my_position = btag_order[itag];
      jetcsvorderedhistos[my_position]->Fill(myDetBTagDiscr[itag]);


      //cout << "For i=" << itag << " rank: " << btag_order[itag] << " discr. " << myDetBTagDiscr[itag] << endl;
    }
    
    // Fill the event-level histograms
    h_njets->Fill(njets);
    if(njets>=6) h_njets_ge6->Fill(njets);
    if(njets>=6) h_ntags_ge6->Fill(ntags);

 }   

  // Write the histograms to a file
  

  cout << "number of jets used in W1 AND W2 = " << WProb_counter << endl;

  TFile* outfile = new TFile("histograms_ttJets_massivefitter_8j4t_oldfile_tagVetoW_tagTopB.root","RECREATE");
  //TFile* outfile = new TFile("histograms_ttH_massivefitter_8j4t_oldfile.root","RECREATE");
  h_jets_pt->Write();
  h_njets->Write();
  h_njets_ge6->Write();
  h_ntags_ge6->Write();
  h_avecsv_tagged->Write();
  h_avecsv_untagged->Write();
  h_avecsv_all->Write();
  h_csv_alljets->Write();
  h_csv_Hbjets->Write();
  h_csv_Hbjets_true->Write();
  h_csv_Wjets->Write();
  h_csv_Wjets_true->Write();
  h_csv_Wcjets_true->Write();
  h_csv_Tbjets->Write();
  h_csv_Tbjets_true->Write();

  h_besthadWmass->Write();
  h_secondbesthadWmass->Write();
  h_besthadWTopmass->Write();
  h_secondbesthadWTopmass->Write();
  h_finalJetsmass->Write();
  h_finalJetsmass_0Hbanywhere->Write();
  h_finalJetsmass_1Hbanywhere->Write();
  h_finalJetsmass_2Hbanywhere->Write();
  h_finalJetsmass_0Hb->Write();
  h_finalJetsmass_1Hb->Write();
  h_finalJetsmass_2Hb->Write();
  h_W1_chiSq->Write();
  h_W2_chiSq->Write();
  h_T1_chiSq->Write();
  h_T2_chiSq->Write();
  h_TOT_chiSq->Write();
  h_allcombo_chiSq->Write();
  h_higgsdecay->Write();
  h_higgsdecayb_eta->Write();
  h_higgsdecayb_dr->Write();
  h_higgsdecayb_pt_missing->Write();
  h_higgsdecayb_eta_missing->Write();
  h_b_destination->Write();
  h_singleb_destination->Write();

  h_Wdecay_true_8j4t->Write();
  h_Wdecay_true->Write();

  h_Wdecay_true_8j4t->Print("all");
  h_Wdecay_true->Print("all");


  //h_higgsdecay->Print("all");


  h_jet_genID->Write();
  h_jet_genMotherID->Write();
  for (unsigned int i = 0; i < jetpthistos.size(); ++i) jetpthistos[i]->Write();
  for (unsigned int i = 0; i < jetetahistos.size(); ++i) jetetahistos[i]->Write();
  for (unsigned int i = 0; i < jetcsvhistos.size(); ++i) jetcsvhistos[i]->Write();
  for (unsigned int i = 0; i < jetcsvorderedhistos.size(); ++i) jetcsvorderedhistos[i]->Write();

  return;
}

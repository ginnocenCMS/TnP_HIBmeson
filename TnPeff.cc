#include "format.h"
#include "utilities.h"

using namespace std;

//const int nMuPtBin = 1;
const int nMuPtBin = 7;
//const double MuPtBin[nMuPtBin+1] = {0,30};
const double MuPtBin[nMuPtBin+1] = {0.0,1.5,3.0,4.5,6.0,9.0,20.0,30.0};
////TH1 setting of mumu
int nBin = 140;
double blow = 2.6;
double bhigh = 3.5;
double measure = (bhigh-blow)/nBin;
//width 2.9~3.3
double sideband_width = 0.2;
bool quiet = false;
int numCPU = 6;

//TTree variable
double massTrg, massTrk, massID;
double passTrg, passTrk, passID;
double ptTrg, ptTrk, ptID;
double etaTrg, etaTrk, etaID;

void doTnP(TH1D* hmumu_trg_all[], TH1D* hmumu_trg_pass[], TH1D* hmumu_trk_all[], TH1D* hmumu_trk_pass[], TH1D* hmumu_id_all[], TH1D* hmumu_id_pass[], 
TTree* TnPtreeTrg, TTree* TnPtreeTrk, TTree* TnPtreeID){
  TChain *nt = new TChain("demo/root");
  //TChain *hlt = new TChain("hltanalysis/HltTree");
  TChain *hlt = new TChain("demo/HltTree");
  string inf = "../BfinderBoostedMC_20140628_inclBtoPsiMuMu_pa_STARTHI53_V27-v1_00000_reduce.root";
  nt->Add(inf.c_str());
  hlt->Add(inf.c_str());

  nt->SetBranchStatus("*" ,0);
  nt->SetBranchStatus("MuonInfo*" ,1);
  hlt->SetBranchStatus("*" ,0);
  hlt->SetBranchStatus("HLT_PAMu3_v1" ,1);              
  nt->AddFriend(hlt);

  MuonInfoBranches MuonInfo;
  MuonInfo.setbranchadd(nt);

  Int_t HLT_PAMu3_v1; 
  nt->SetBranchAddress("HLT_PAMu3_v1",&HLT_PAMu3_v1);                         

  TLorentzVector mu1Vec;
  TLorentzVector mu2Vec;
  TLorentzVector jpsiVec;
  int nevents_total = nt->GetEntries();                                                  
//  nevents_total = 100000;
  for(int entry=0; entry<nevents_total; entry++){
    if ((entry%10000) == 0) printf("Loading event #%d of %d.\n",entry,nevents_total);
    nt->GetEntry(entry);
    if(MuonInfo.size<2) continue;
    if(!HLT_PAMu3_v1) continue;
	for(int mu1 = 0; mu1 < MuonInfo.size; mu1++){
	  for(int mu2 = mu1+1; mu2 < MuonInfo.size; mu2++){
		if(MuonInfo.charge[mu1]*MuonInfo.charge[mu2] > 0) continue;
        mu1Vec.SetPtEtaPhiM(MuonInfo.pt[mu1], MuonInfo.eta[mu1], MuonInfo.phi[mu1], MUON_MASS);
        mu2Vec.SetPtEtaPhiM(MuonInfo.pt[mu2], MuonInfo.eta[mu2], MuonInfo.phi[mu2], MUON_MASS);
		jpsiVec = mu1Vec+mu2Vec;
		////Trigger
		bool mu1Trg = false;
		bool mu2Trg = false;
		//std::cout<<MuTrgMatchTrgObjE->at(mu1).at(0)<<std::endl;
		if(MuonInfo.MuTrgMatchTrgObjE->at(mu1).at(0) > 0) mu1Trg = true;
		if(MuonInfo.MuTrgMatchTrgObjE->at(mu2).at(0) > 0) mu2Trg = true;
		////
		////Trigger efficiency
		//bool pickmu1 = false;
		//if(int(MuonInfo.pt[mu1]*10)%2==1) pickmu1 = true;
		//if(pickmu1){
		if(MuonInfo.isTrackerMuon[mu1] && KisooTrackSel(MuonInfo, mu1) && KisooGlobalSel(MuonInfo, mu1) && Acc(MuonInfo, mu1) && mu1Trg &&
		   MuonInfo.isTrackerMuon[mu2] && KisooTrackSel(MuonInfo, mu2) && KisooGlobalSel(MuonInfo, mu2) && Acc(MuonInfo, mu2) 
		){
		  //Fill Tree
		  massTrg = jpsiVec.Mag(); ptTrg = MuonInfo.pt[mu2]; etaTrg = MuonInfo.eta[mu2]; passTrg = 0;
          if(mu2Trg) {passTrg = 1;}
		  TnPtreeTrg->Fill();
		  //Fill Histogram
		  for(int i = 0; i < nMuPtBin; i++){
            if(MuonInfo.pt[mu2] > MuPtBin[i] && MuonInfo.pt[mu2] < MuPtBin[i+1]){
		      hmumu_trg_all[i]->Fill(jpsiVec.Mag());
			    if(mu2Trg){
                hmumu_trg_pass[i]->Fill(jpsiVec.Mag());
			    }
              break;
            }
          }
        }
		//if(!pickmu1){
		if(MuonInfo.isTrackerMuon[mu2] && KisooTrackSel(MuonInfo, mu2) && KisooGlobalSel(MuonInfo, mu2) && Acc(MuonInfo, mu2) && mu2Trg &&
		   MuonInfo.isTrackerMuon[mu1] && KisooTrackSel(MuonInfo, mu1) && KisooGlobalSel(MuonInfo, mu1) && Acc(MuonInfo, mu1) 
		){
		  //Fill Tree
		  massTrg = jpsiVec.Mag(); ptTrg = MuonInfo.pt[mu1]; etaTrg = MuonInfo.eta[mu1]; passTrg = 0;
          if(mu1Trg) {passTrg = 1;}
		  TnPtreeTrg->Fill();
		  //Fill Histogram
		  for(int i = 0; i < nMuPtBin; i++){
            if(MuonInfo.pt[mu1] > MuPtBin[i] && MuonInfo.pt[mu1] < MuPtBin[i+1]){
		      hmumu_trg_all[i]->Fill(jpsiVec.Mag());
			    if(mu1Trg){
                hmumu_trg_pass[i]->Fill(jpsiVec.Mag());
			    }
              break;
            }
          }
        }
        ////
		////Tracking efficiency
		if(MuonInfo.isTrackerMuon[mu1] && KisooTrackSel(MuonInfo, mu1) && KisooGlobalSel(MuonInfo, mu1) && Acc(MuonInfo, mu1) && mu1Trg	&&
           MuonInfo.outerTrackisNonnull[mu2]
		){
		  //Fill Tree
		  massTrk = jpsiVec.Mag(); ptTrk = MuonInfo.pt[mu2]; etaTrk = MuonInfo.eta[mu2]; passTrk = 0;
          if(MuonInfo.isTrackerMuon[mu2] && KisooTrackSel(MuonInfo, mu2)) {passTrk = 1;}
		  TnPtreeTrk->Fill();
		  //Fill Histogram
		  for(int i = 0; i < nMuPtBin; i++){
            if(MuonInfo.pt[mu2] > MuPtBin[i] && MuonInfo.pt[mu2] < MuPtBin[i+1]){
		      hmumu_trk_all[i]->Fill(jpsiVec.Mag());
			    if(MuonInfo.isTrackerMuon[mu2] && KisooTrackSel(MuonInfo, mu2)){
                hmumu_trk_pass[i]->Fill(jpsiVec.Mag());
			    }
              break;
            }
          }
        }
		if(MuonInfo.isTrackerMuon[mu2] && KisooTrackSel(MuonInfo, mu2) && KisooGlobalSel(MuonInfo, mu2) && Acc(MuonInfo, mu2) && mu2Trg	&&
           MuonInfo.outerTrackisNonnull[mu1]  
		){
		  //Fill Tree
		  massTrk = jpsiVec.Mag(); ptTrk = MuonInfo.pt[mu1]; etaTrk = MuonInfo.eta[mu1]; passTrk = 0;
          if(MuonInfo.isTrackerMuon[mu1] && KisooTrackSel(MuonInfo, mu1)) {passTrk = 1;}
		  TnPtreeTrk->Fill();
		  //Fill Histogram
		  for(int i = 0; i < nMuPtBin; i++){
            if(MuonInfo.pt[mu1] > MuPtBin[i] && MuonInfo.pt[mu1] < MuPtBin[i+1]){
		      hmumu_trk_all[i]->Fill(jpsiVec.Mag());
			    if(MuonInfo.isTrackerMuon[mu1] && KisooTrackSel(MuonInfo, mu1)){
                hmumu_trk_pass[i]->Fill(jpsiVec.Mag());
			    }
              break;
            }
          }
        }
		////
		////ID efficiency
        if(MuonInfo.isTrackerMuon[mu1] && KisooTrackSel(MuonInfo, mu1) && KisooGlobalSel(MuonInfo, mu1) && Acc(MuonInfo, mu1) && mu1Trg &&
           (MuonInfo.type[mu2] & (1<<4)) && KisooTrackSel(MuonInfo, mu2) && Acc(MuonInfo, mu2) 
        ){  
		  //Fill Tree
		  massID = jpsiVec.Mag(); ptID = MuonInfo.pt[mu2]; etaID = MuonInfo.eta[mu2]; passID = 0;
          if(MuonInfo.isTrackerMuon[mu2] && KisooTrackSel(MuonInfo, mu2) && KisooGlobalSel(MuonInfo, mu2) && Acc(MuonInfo, mu2)) {passID = 1;}
		  TnPtreeID->Fill();
		  //Fill Histogram
		  for(int i = 0; i < nMuPtBin; i++){
            if(MuonInfo.pt[mu2] > MuPtBin[i] && MuonInfo.pt[mu2] < MuPtBin[i+1]){
		      hmumu_id_all[i]->Fill(jpsiVec.Mag());
			    if(MuonInfo.isTrackerMuon[mu2] && KisooTrackSel(MuonInfo, mu2) && KisooGlobalSel(MuonInfo, mu2) && Acc(MuonInfo, mu2)){
                hmumu_id_pass[i]->Fill(jpsiVec.Mag());
			    }
              break;
            }
          }
        }
        if(MuonInfo.isTrackerMuon[mu2] && KisooTrackSel(MuonInfo, mu2) && KisooGlobalSel(MuonInfo, mu2) && Acc(MuonInfo, mu2) && mu2Trg &&
           (MuonInfo.type[mu1] & (1<<4)) && KisooTrackSel(MuonInfo, mu1) && Acc(MuonInfo, mu1) 
        ){  
		  //Fill Tree
		  massID = jpsiVec.Mag(); ptID = MuonInfo.pt[mu1]; etaID = MuonInfo.eta[mu1]; passID = 0;
          if(MuonInfo.isTrackerMuon[mu1] && KisooTrackSel(MuonInfo, mu1) && KisooGlobalSel(MuonInfo, mu1) && Acc(MuonInfo, mu1)) {passID = 1;}
		  TnPtreeID->Fill();
		  //Fill Histogram
		  for(int i = 0; i < nMuPtBin; i++){
            if(MuonInfo.pt[mu1] > MuPtBin[i] && MuonInfo.pt[mu1] < MuPtBin[i+1]){
		      hmumu_id_all[i]->Fill(jpsiVec.Mag());
			    if(MuonInfo.isTrackerMuon[mu1] && KisooTrackSel(MuonInfo, mu1) && KisooGlobalSel(MuonInfo, mu1) && Acc(MuonInfo, mu1)){
                hmumu_id_pass[i]->Fill(jpsiVec.Mag());
			    }
              break;
            }
          }
        }
	    ////
      }
    }
  }
}

int TnPeff(){
  gStyle->SetOptStat("0");
  TTree* TnPtreeTrg = new TTree("TnPtreeTrg","");
  TTree* TnPtreeTrk = new TTree("TnPtreeTrk","");
  TTree* TnPtreeID = new TTree("TnPtreeID","");
  TnPtreeTrg->Branch("mass",&massTrg);
  TnPtreeTrg->Branch("pt",&ptTrg);
  TnPtreeTrg->Branch("eta",&etaTrg);
  TnPtreeTrg->Branch("pass",&passTrg);

  TnPtreeTrk->Branch("mass",&massTrk);
  TnPtreeTrk->Branch("pt",&ptTrk);
  TnPtreeTrk->Branch("eta",&etaTrk);
  TnPtreeTrk->Branch("pass",&passTrk);

  TnPtreeID->Branch("mass",&massID);
  TnPtreeID->Branch("pt",&ptID);
  TnPtreeID->Branch("eta",&etaID);
  TnPtreeID->Branch("pass",&passID);
  
  TH1D* hmumu_trg_all[nMuPtBin];
  TH1D* hmumu_trg_pass[nMuPtBin];
  TH1D* hmumu_trk_all[nMuPtBin];
  TH1D* hmumu_trk_pass[nMuPtBin];
  TH1D* hmumu_id_all[nMuPtBin];
  TH1D* hmumu_id_pass[nMuPtBin];
  TF1* tf1mumu_trg_all[nMuPtBin];
  TF1* tf1mumu_trg_pass[nMuPtBin];
  TF1* tf1mumu_trk_all[nMuPtBin];
  TF1* tf1mumu_trk_pass[nMuPtBin];
  TF1* tf1mumu_id_all[nMuPtBin];
  TF1* tf1mumu_id_pass[nMuPtBin];
  TH1D* eff_trg = myTH1D("eff_trg", "Mu Pt (GeV)", "efficiency", 2, nMuPtBin, MuPtBin);
  TH1D* eff_trk = myTH1D("eff_trk", "Mu Pt (GeV)", "efficiency", 3, nMuPtBin, MuPtBin);
  TH1D* eff_id = myTH1D("eff_id", "Mu Pt (GeV)", "efficiency", 4, nMuPtBin, MuPtBin);
  for(int i = 0; i < nMuPtBin; i++){
    hmumu_trg_all[i] = new TH1D(Form("hmumu_trg_all%d",i),"",nBin,blow,bhigh);
    hmumu_trg_pass[i] = new TH1D(Form("hmumu_trg_pass%d",i),"",nBin,blow,bhigh);
    hmumu_trk_all[i] = new TH1D(Form("hmumu_trk_all%d",i),"",nBin,blow,bhigh);
    hmumu_trk_pass[i] = new TH1D(Form("hmumu_trk_pass%d",i),"",nBin,blow,bhigh);
    hmumu_id_all[i] = new TH1D(Form("hmumu_id_all%d",i),"",nBin,blow,bhigh);
    hmumu_id_pass[i] = new TH1D(Form("hmumu_id_pass%d",i),"",nBin,blow,bhigh);

	tf1mumu_trg_all[i] = myTF1(Form("tf1_trg_all%d",i), blow, bhigh);
	tf1mumu_trg_pass[i] = myTF1(Form("tf1_trg_pass%d",i), blow, bhigh);
	tf1mumu_trk_all[i] = myTF1(Form("tf1_trk_all%d",i), blow, bhigh);
	tf1mumu_trk_pass[i] = myTF1(Form("tf1_trk_pass%d",i), blow, bhigh);
	tf1mumu_id_all[i] = myTF1(Form("tf1_id_all%d",i), blow, bhigh);
	tf1mumu_id_pass[i] = myTF1(Form("tf1_id_pass%d",i), blow, bhigh);
  }
  TCanvas *c= new TCanvas("c","",600,600);
  c->SetTopMargin(0.07504363);
  c->cd();
  doTnP(hmumu_trg_all, hmumu_trg_pass, hmumu_trk_all, hmumu_trk_pass, hmumu_id_all, hmumu_id_pass, TnPtreeTrg, TnPtreeTrk, TnPtreeID);
  TF1 *signal = new TF1("signal","[0]*([4]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[4])*Gaus(x,[1],[3])/(sqrt(2*3.14159)*[3]))");
  signal->SetRange(blow, bhigh);
  for(int i = 0; i < nMuPtBin; i++){
	cout<<"Mu Pt Bin(GeV): " <<MuPtBin[i]<<endl;
	////side band method
	int siglow = -1;
	int sighigh = -1;
	for(int b = 0; b < nBin; b++){
	 if(siglow == -1 && hmumu_trg_all[i]->GetBinCenter(b+1) > JPSI_MASS-sideband_width) siglow = b+1; 
	 if(hmumu_trg_all[i]->GetBinCenter(b+1) < JPSI_MASS+sideband_width) sighigh = b+1; 
	}
	int sigwidth = int((sighigh-siglow+1)/2);
	int sblow = siglow-sigwidth;
	int sbhigh = sighigh+sigwidth;
	if(sblow < 0) sblow = 0;
	if(sbhigh > nBin) sbhigh = nBin;
	//cout<<siglow-sigwidth<<endl;cout<<siglow-1<<endl;cout<<sighigh+1<<endl;cout<<sighigh+sigwidth<<endl;

	double nall, npass;
	//Using rootfit
//	tf1mumu_trg_all[i]->SetParameter(0,hmumu_trg_all[i]->GetEntries());
//	hmumu_trg_all[i]->Fit(Form("tf1_trg_all%d",i), "L q", "", blow, bhigh);
//	hmumu_trg_all[i]->Fit(Form("tf1_trg_all%d",i), "L m q", "", blow, bhigh);
//	for(int k = 0; k < 5; k++){
//	  signal->SetParameter(k, tf1mumu_trg_all[i]->GetParameter(k));}
//	nall = signal->Integral(blow, bhigh)/measure;
//    tf1mumu_trg_pass[i]->SetParLimits(0,0,nall);
//	hmumu_trg_pass[i]->Fit(Form("tf1_trg_pass%d",i), "L q", "", blow, bhigh);
//	hmumu_trg_pass[i]->Fit(Form("tf1_trg_pass%d",i), "L m q", "", blow, bhigh);
//	for(int k = 0; k < 5; k++){
//	  signal->SetParameter(k, tf1mumu_trg_pass[i]->GetParameter(k));}
//	npass = signal->Integral(blow, bhigh)/measure;

	//Using sideband subtraction
	nall = hmumu_trg_all[i]->Integral(siglow,sighigh)-(hmumu_trg_all[i]->Integral(siglow-sigwidth, siglow-1)+hmumu_trg_all[i]->Integral(sighigh+1,sighigh+sigwidth));
	npass = hmumu_trg_pass[i]->Integral(siglow,sighigh)-(hmumu_trg_pass[i]->Integral(siglow-sigwidth, siglow-1)+hmumu_trg_pass[i]->Integral(sighigh+1,sighigh+sigwidth));
	eff_trg->SetBinContent(i+1, npass/nall);
	//eff_trg->SetBinError(i+1, tf1mumu_trg_all[i]->GetParError(0)/measure);
	eff_trg->SetBinError(i+1, 0.00001);
	cout<<"his trg all : "<<hmumu_trg_all[i]->GetEntries()<<endl;
	cout<<"his trg pass: "<<hmumu_trg_pass[i]->GetEntries()<<endl;
    cout<<"trg all : "<<nall<<endl;
    cout<<"trg pass: "<<npass<<endl;
	cout<<"eff trg: "<<npass/nall<<endl;

	//Using sideband subtraction
	nall = hmumu_trk_all[i]->Integral(siglow,sighigh)-(hmumu_trk_all[i]->Integral(siglow-sigwidth, siglow-1)+hmumu_trk_all[i]->Integral(sighigh+1,sighigh+sigwidth));
	npass = hmumu_trk_pass[i]->Integral(siglow,sighigh)-(hmumu_trk_pass[i]->Integral(siglow-sigwidth, siglow-1)+hmumu_trk_pass[i]->Integral(sighigh+1,sighigh+sigwidth));
	eff_trk->SetBinContent(i+1, npass/nall);
	eff_trk->SetBinError(i+1, 0.00001);
	cout<<"his trk all : "<<hmumu_trk_all[i]->GetEntries()<<endl;
	cout<<"his trk pass: "<<hmumu_trk_pass[i]->GetEntries()<<endl;
    cout<<"trk all : "<<nall<<endl;
    cout<<"trk pass: "<<npass<<endl;
	cout<<"eff trk: "<<npass/nall<<endl;

	//Using sideband subtraction
	nall = hmumu_id_all[i]->Integral(siglow,sighigh)-(hmumu_id_all[i]->Integral(siglow-sigwidth, siglow-1)+hmumu_id_all[i]->Integral(sighigh+1,sighigh+sigwidth));
	npass = hmumu_id_pass[i]->Integral(siglow,sighigh)-(hmumu_id_pass[i]->Integral(siglow-sigwidth, siglow-1)+hmumu_id_pass[i]->Integral(sighigh+1,sighigh+sigwidth));
	eff_id->SetBinContent(i+1, npass/nall);
	eff_id->SetBinError(i+1, 0.00001);
    cout<<"his id all : "<<hmumu_id_all[i]->GetEntries()<<endl;
    cout<<"his id pass: "<<hmumu_id_pass[i]->GetEntries()<<endl;
    cout<<"id all : "<<nall<<endl;
    cout<<"id pass: "<<npass<<endl;
    cout<<"eff id: "<<npass/nall<<endl;
  }
  eff_trg->Draw("pe");
  eff_trk->Draw("pe same");
  eff_id->Draw("pe same");
  TLegend *leg = myLegend(0.6778523,0.399651,0.8775168,0.6073298);
  leg->AddEntry(eff_trg, "Trg", "p");
  leg->AddEntry(eff_trk, "Trk", "p");
  leg->AddEntry(eff_id, "ID", "p");
  leg->Draw("same");
  TFile *outf = new TFile("../TnPOutMC.root","recreate");
  outf->cd();
  TnPtreeTrg->Write();
  TnPtreeTrk->Write();
  TnPtreeID->Write();
  eff_trg->Write();
  eff_trk->Write();
  eff_id->Write();
  c->Write();

  for(int i = 0; i < nMuPtBin; i++){
    hmumu_trg_all[i]->Write();
    hmumu_trg_pass[i]->Write();
    hmumu_trk_all[i]->Write();
    hmumu_trk_pass[i]->Write();
    hmumu_id_all[i]->Write();
    hmumu_id_pass[i]->Write();
  }
  outf->Close();
  return 0;
}


#include <RAT/DU/DSReader.hh>
#include <RAT/DS/FitVertex.hh>
#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DU/LightPathCalculator.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DS/FitResult.hh>
#include <TGraph2D.h>
#include <string>
#include <TMultiGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFrame.h>
#include <TGraph.h>
#include <TPaveText.h>
#include <TCanvas.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLine.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TStyle.h>
#include <TPad.h>


#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <iostream>
#include <fstream>
/*
 * The point of this is that you want to use root instead of echidna to look at the energy recostruction over z dimension
 * the way to do this is to give a function a couple of histrograms.
 */
/*Function signitures*/


vector<string> glob( const string& path, const string& start )
{
  vector<string> result;
  TSystemDirectory dir(path.c_str(), path.c_str());
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(".root") && fname.BeginsWith( start ) ) {
        stringstream fullPath; fullPath << path << "/" << fname.Data();
        result.push_back(fullPath.str());
      }
    }
  }
  return result;
}






void findTres(string filename,TH1D* hist2){

	RAT::DU::LightPathCalculator fLightPath = RAT::DU::Utility::Get()->GetLightPathCalculator();
	RAT::DU::DSReader dsreader2( filename);
	double time_res;
	int entrycount2= dsreader2.GetEntryCount();
	for( int ent=0; ent</*entrycount*/500;ent++){

		const RAT::DS::Entry& ds=dsreader2.GetEntry(ent);

		for(int iEv=0; iEv<ds.GetEVCount(); iEv++){
			const RAT::DS::EV& rev = ds.GetEV(iEv);

			const RAT::DS::FitVertex& vertex= rev.GetFitResult("scintFitter").GetVertex(0);

			TVector3 fEventPosition= vertex.GetPosition();
			double fEventTime= vertex.GetTime();
			const RAT::DS::CalPMTs& calpmts = rev.GetCalPMTs();

			for(int iPMT=0;iPMT<calpmts.GetCount();iPMT++){
				const RAT::DS::PMTCal& pmt = calpmts.GetPMT(iPMT);

				const TVector3 pmtPosition = RAT::DU::Utility::Get()->GetPMTInfo().GetPosition(pmt.GetID());
				const float pmtTime = pmt.GetTime();
				fLightPath.CalcByPosition(fEventPosition, pmtPosition);
				double distInInnerAV = fLightPath.GetDistInInnerAV();
				double distInAV = fLightPath.GetDistInAV();
				double distInWater = fLightPath.GetDistInWater();

				const double flightTime = RAT::DU::Utility::Get()->GetEffectiveVelocity().CalcByDistance(distInInnerAV, distInAV, distInWater);

				time_res=pmtTime-flightTime-fEventTime;
				cout<< "Time "<< time_res<< endl;
				hist2->Fill(time_res);
			}	

		}
	}
}


void go(){

	int i=0;

	vector<string> biFileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Bi210/root/","SolarBi210");
	TH1D* hist= new TH1D("hist","",100,-100,300);
	hist->SetLineColor(kBlue);

	findTres(biFileList[0],hist);

	vector<string> poFileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po210/root/","SolarPo210");
	TH1D* hist2= new TH1D("hist2","",100,-100,300);
	hist2->SetLineColor(kBlack);

	findTres(poFileList[0],hist2);
	{
	TAxis* xaxis = hist->GetXaxis();
	double biTotal=hist->Intergral(xaxis->GetBin(-100),xaxis->GetBin(300));
	TAxis* xaxis = hist2->GetXaxis();
	double poTotal=hist2->Intergral(xaxis->GetBin(-100),xaxis->GetBin(300));
	cout<<"biTotal "<<biTotal<<" poTotal "<<poTotal<<endl;
	TCanvas *c1 = new TCanvas();
	//gStyle->SetOptStat(0);
	hist->Draw();
	hist2->Draw("same");
	TLegend * leg= new TLegend(0.5,0.7,0.9,0.9);
	leg->AddEntry(hist,"Bi210","l");
	leg->AddEntry(hist2,"Po210","l");
	leg->Draw();
	c1->Print("TimingResiudals.png");
	}

	

}





















void FillHists( TFile* f,TH2D* RecoVsZ,TH2D* McVsZ,TH2D* NhitsVsR)
{

//cout << fileName << endl;

    TNtuple *T = (TNtuple *) f->Get("output");
    
    Int_t runID; // Run number
    TString ratversion; // version of rat used
    Int_t mcIndex; //MC event number
    // identifies a trigger within a MC event
    // counts from 0 to GetEvCount() -1 within each MC events
    Int_t evIndex;
    Bool_t mc; // tells whether there is MC data associated to a trigger
    Int_t eventID; // Event number (event with trigger within a MC event)
    // MC variables
    Int_t pdg1; // PDG code of the most energetic particle in an event
    Int_t pdg2; // PDG code of the second most energetic particle in an event
    Int_t parentpdg1; // PDG code of the parent of the most energetic particle in an event
    Int_t parentpdg2; // PDG code of the parent of the second most energetic particle in an event
    Double_t mcke1;
    Double_t mcke2;
    Double_t mcPosx,mcPosy,mcPosz, mcPosr;
    Int_t nhits;
    Double_t q;
    ULong64_t clockCount10,clockCount50;
    Int_t uTDays,uTSecs,uTNSecs;
    Bool_t scintFit,waterFit,partialFit,fitValid;
    Double_t posx,posy,posz, posr;
    Double_t energy;
    Double_t time;
    Double_t timingPeaks;
    Double_t meanTime;
    Double_t biPoCumul;
    Double_t biPoLikelihood212;
    Double_t biPoLikelihood214;
    Double_t mcEdepQuenched;
    Double_t alphaBeta212;
    Double_t alphaBeta214;

    // Number of ntuples entries
    Int_t nentries = T->GetEntries();
    cout << "Total number of entries " << nentries << endl;

    // Set branch addresses
    T->SetBranchAddress("runID",&runID);
    T->SetBranchAddress("mcIndex",&mcIndex);
    T->SetBranchAddress("evIndex",&evIndex);
    T->SetBranchAddress("mc",&mc);
    // MC
    T->SetBranchAddress("pdg1",&pdg1);
    T->SetBranchAddress("pdg2",&pdg2);
    T->SetBranchAddress("parentpdg1",&parentpdg1);
    T->SetBranchAddress("parentpdg2",&parentpdg2);
    T->SetBranchAddress("mcke1",&mcke1);
    T->SetBranchAddress("mcke2",&mcke2);
    T->SetBranchAddress("mcPosx",&mcPosx);
    T->SetBranchAddress("mcPosy",&mcPosy);
    T->SetBranchAddress("mcPosz",&mcPosz);
    // Exists if trigger
    T->SetBranchAddress("eventID",&eventID); // 0 if no trigger
    T->SetBranchAddress("uTDays",&uTDays);
    T->SetBranchAddress("uTSecs",&uTSecs);
    T->SetBranchAddress("uTNSecs",&uTNSecs);
    T->SetBranchAddress("scintFit",&scintFit);
    T->SetBranchAddress("waterFit",&waterFit);
    T->SetBranchAddress("partialFit",&partialFit);
    T->SetBranchAddress("fitValid",&fitValid);
    T->SetBranchAddress("posx",&posx);
    T->SetBranchAddress("posy",&posy);
    T->SetBranchAddress("posz",&posz);
    T->SetBranchAddress("posr",&posr);
    T->SetBranchAddress("energy",&energy);
    T->SetBranchAddress("time",&time);
    T->SetBranchAddress("clockCount10",&clockCount10);
    T->SetBranchAddress("clockCount50",&clockCount50);
    T->SetBranchAddress("nhits",&nhits);
    T->SetBranchAddress("q",&q);
    T->SetBranchAddress("mcEdepQuenched",&mcEdepQuenched);
    // Classifier valuse
    T->SetBranchAddress("meanTime",&meanTime);
    T->SetBranchAddress("timingPeaks",&timingPeaks);
    T->SetBranchAddress("biPoCumul",&biPoCumul);
    T->SetBranchAddress("biPoLikelihood212",&biPoLikelihood212);
    T->SetBranchAddress("biPoLikelihood214",&biPoLikelihood214);
    T->SetBranchAddress("alphaBeta212",&alphaBeta212);
    T->SetBranchAddress("alphaBeta214",&alphaBeta214);
    int totalEvents = 0;
    
    int i=0;
    while(i<nentries){
        T->GetEntry(i);
                RecoVsZ->Fill(energy, posz);
                McVsZ->Fill( mcEdepQuenched, posz);
                NhitsVsR->Fill( nhits, posr);
   i++;
    }// loop over entries
   f->Close();
}



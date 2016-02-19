#include <RAT/DU/DSReader.hh>
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
void all(string path,string name);





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


void FillHist(TFile* file,TH1D * hist, string entry){

	TTree* Tree = (TTree*) file->Get("output");
	Int_t para;

	Tree->SetBranchAddress(entry.c_str(),&para);
	Int_t n = (Int_t)Tree->GetEntries();

	for( Int_t i =0;i<n;i++){
	Tree->GetEntry(i);	
	hist->Fill(para);
	}

}

void FillHist(TFile* file,TH1D * hist){

	TTree* Tree = (TTree*) file->Get("output");
	Double_t para;

	Tree->SetBranchAddress("berkeleyAlphaBeta",&para);
	Int_t n = (Int_t)Tree->GetEntries();

	for( Int_t i =0;i<n;i++){
	Tree->GetEntry(i);	
	hist->Fill(para);
	}

}
void FindHits(string filename,TH1D * scaleHist ){
	vector<double> scaleList;
	double scale;
	cout<<"filename : "<<filename<<endl;
	RAT::DU::DSReader dsReader( filename );
	// Loop through entrys in rootfile
	for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ ){
		int firstWindowHits;	
		int unbuiltHits;	
		int allHits=0;	
		const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );

		size_t count=rDS.GetMCEVCount();
		for( int iMCEV=0; iMCEV<count; iMCEV++){
			const RAT::DS::MCHits& hits=rDS.GetMCEV(iMCEV).GetMCHits();
			cout <<"count : "<<iMCEV<< " hits : "<< hits.GetCount() <<endl;	
			if (iMCEV==0){
				firstWindowHits=hits.GetCount();
			}
		allHits+=hits.GetCount();
		}
		
		const RAT::DS::MCHits& unhits=rDS.GetMC().GetUnbuiltMCHits();
		unbuiltHits=unhits.GetCount();	
		scale=(double)firstWindowHits/(allHits+unbuiltHits);
		cout<<"ratio between 1st window/all hits : "<< scale<<endl;
		//scaleList.push_back(scale);
		scaleHist->Fill(scale);
	}
scaleHist->Draw();
}

void alphas(){
	all("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po210/root","SolarPo210");
	all("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po212/root","SolarPo212");
	all("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po214/root","SolarPo214");
}

void all(string path,string name){

	TH1D * scaleHist= new TH1D("scaleHist","",100,0,1);
	//vector<string> poFileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po210/root/","SolarPo210");
	vector<string> poFileList= glob(path,name);
	if (poFileList.size()>20){

	for( int i=0; i<20/*poFileList.size()*/; i++ ){
		TFile * file= TFile::Open(poFileList[i].c_str());	
		FindHits( poFileList[i],  scaleHist);
	  }
	}else{

	for( int i=0; i<poFileList.size(); i++ ){
		TFile * file= TFile::Open(poFileList[i].c_str());	
		FindHits( poFileList[i],  scaleHist);
	  }
	}

// Plotting things
	{
	TCanvas* c1= new TCanvas();
	c1->cd();
	TAxis * xAxis=   scaleHist->GetXaxis();
	TAxis * yAxis=    scaleHist->GetYaxis();
	//McVsZ->Draw("surf1");
	scaleHist->SetTitle(name.c_str());
	xAxis->SetTitle("Fraction in 1st window");
	//scaleHist->SetXaxis("Fraction in 1st window");
	scaleHist->Draw();
	c1->Print((name+".png").c_str());
	c1->Print((name+".eps").c_str());
	
	TFile fileout("scales.root","UPDATE");
	scaleHist->Write();
	fileout.Close();
	}
}
void all(){

	TH1D * scaleHist= new TH1D("scaleHist","",100,0,1);
	vector<string> poFileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po210/root/","SolarPo210");

	for( int i=0; i<20/*poFileList.size()*/; i++ ){
		TFile * file= TFile::Open(poFileList[i].c_str());	
		FindHits( poFileList[i],  scaleHist);
	  }
// Plotting things
	{
	TCanvas* c1= new TCanvas();
	c1->cd();
	TAxis * xAxis=   scaleHist->GetXaxis();
	TAxis * yAxis=    scaleHist->GetYaxis();
	//McVsZ->Draw("surf1");
	scaleHist->Draw();
	c1->Print("test.png");
	c1->Print("test.eps");
	
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


int func(){
	
	gStyle->SetOptStat(0);
	TH2D * RecoVsZ   = new TH2D("RecoVsZ","Reconstrution energy vs z coordinate",20,0,1,6000,-6000,6000);
	TH2D * McVsZ   = new TH2D("McVsZ","MC energy vs z coordinate",100,0,1,100,-7000,7000);
	TH2D * NhitsVsR   = new TH2D("NhitsVsR","Nhits vs r",50,0, 1000,6000,0,6000);
	
	TLegend* t1 = new TLegend( 0.6, 0.7, 0.89, 0.88 );
	vector<string> biFileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Bi210/root/","SolarBi210");
	vector<string> poFileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po210/root/","SolarPo210");
//	for( int i=0; i<biFileList.size(); i++ ){
//		TFile * file= TFile::Open(biFileList[i].c_str());	
//		FillHist( file, hBi210,"nhits");
//	  }

	cout<<"Got Here"<<endl;
	for( int i=0; i<poFileList.size(); i++ ){
		TFile * file= TFile::Open(poFileList[i].c_str());	
		//FillHist( file, hPo210,"nhits");

		FillHists( file, RecoVsZ,McVsZ,NhitsVsR);
	  }
// Plotting things
	{
	TCanvas* c1= new TCanvas();
	c1->cd();
	TAxis * xAxis= NhitsVsR->GetXaxis();
	TAxis * yAxis= NhitsVsR->GetYaxis();
	TH1D * projection= NhitsVsR->ProjectionX("projection",xAxis->FindBin(0.0), xAxis->FindBin(1.0),"posr<1000");
	NhitsVsR->SetXTitle("Reco Energy");
	NhitsVsR->SetYTitle("posz");
	NhitsVsR->Draw("contz0");
	//McVsZ->Draw("surf1");
	c1->Print("test.png");
	c1->Print("test.eps");
	
	TCanvas* c2= new TCanvas();
	c2->cd();
	projection->Draw();
	
	}

	{
	TCanvas* c1= new TCanvas();
	c1->cd();
	TAxis * xAxis= McVsZ->GetXaxis();
	TH1D * projection= McVsZ->ProjectionX("projection",xAxis->FindBin(0.0), xAxis->FindBin(1.0),"posz<6000");
	McVsZ->SetYTitle("posz");
	McVsZ->SetXTitle("Reco Energy");
	McVsZ->Draw("contz0");
	//McVsZ->Draw("surf1");
	c1->Print("test.png");
	c1->Print("test.eps");
	
	projection->Draw();
	
	}
//	McVsZ->Draw();
//	NhitsVsR->Draw();
//	hPo210->Draw("same");
//
//	t1->AddEntry( hBi210, "Bi 210","f");
//	t1->AddEntry( hPo210, "Po 210","f");
//	t1->Draw();
//
//	TFile fileout("plot_nhits.root","RECREATE");
//	fileout.cd();
//	hBi210->Write();
//	hPo210->Write();
//	fileout.Close();
	return 0;

}


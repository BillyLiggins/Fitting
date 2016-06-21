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
#include <TGraphErrors.h>
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

#define SSTR( x ) static_cast< std::ostringstream & >( \
		( std::ostringstream() << std::dec << x ) ).str()
/*
 * The point of this is that you want to use root instead of echidna to look at the energy recostruction over z dimension
 * the way to do this is to give a function a couple of histrograms.
 */
/*Function signitures*/

double mean( vector<double> v ){
	double sum=0;
	for( int i=0 ; i<v.size();i++){
		sum=sum +v[i];
	}

	double mean = sum / v.size();
	return mean;
}

double variance ( vector<double> v , double mean )
{
	double sum = 0.0;
	double temp =0.0;
	double var =0.0;

	for ( int j =0; j <v.size(); j++)
	{
		temp = pow(v[j] - mean,2);
		sum += temp;
	}

	return var = sum/(v.size() -2);
}




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


void FillHist(TFile* file,TH1D * hist){

	TTree* Tree = (TTree*) file->Get("output");

	Double_t MC, RE;
	Int_t evIndex;
	Bool_t fitValid;

	Tree->SetBranchAddress("mcEdepQuenched",&MC);
	Tree->SetBranchAddress("energy",&RE);
	Tree->SetBranchAddress("evIndex",&evIndex);
	Tree->SetBranchAddress("fitValid",&fitValid);
	Int_t n = (Int_t)Tree->GetEntries();
	for( Int_t i =0;i<n;i++){
		Tree->GetEntry(i);	
		if (evIndex==0 && fitValid) {
		// std::cout << "evIndex = "<<evIndex << std::endl;
		hist->Fill(MC-RE);
		}
	}

}
void FillHist(TFile* file,TH1D * hist, const string& entry){

	TTree* Tree = (TTree*) file->Get("output");

	Double_t para;
	Int_t evIndex;

	Tree->SetBranchAddress(entry.c_str(),&para);
	Tree->SetBranchAddress("evIndex",&evIndex);
	Int_t n = (Int_t)Tree->GetEntries();
	for( Int_t i =0;i<n;i++){
		Tree->GetEntry(i);	
		if (evIndex==0) {
		// std::cout << "evIndex = "<<evIndex << std::endl;
		hist->Fill(para);
		}
	}

}




void AlphaEnergyReconstructionDiscrepancy(){
	//This function will plot the discrepancy between reconstruct energy and MC energy. 

	TH1D* MC = new TH1D("MC","MC data",100,0,2);
	TH1D* RE = new TH1D("RE","RE data",100,0,2);
	TH1D* DIF = new TH1D("DIF","",100,1,1);


	vector<string> fileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po212","SolarPo");

	for (int i = 0; i < fileList.size(); i++) {
		TFile* tfile = new TFile(fileList[i].c_str());
		FillHist(tfile,MC,"mcEdepQuenched");
		FillHist(tfile,RE,"energy");
		FillHist(tfile,DIF);
		tfile->Close();
	}

	{
		TCanvas * c1 = new TCanvas();
		c1->SetGrid();
		gStyle->SetOptStat(false);
		MC->SetFillColorAlpha(kGreen,0.1);
		RE->SetFillColorAlpha(kRed,0.1);
		MC->SetTitle("MC and Reconstruction energy of Po212");
		MC->GetXaxis()->SetTitle("Energy (MeV)");
		MC->GetYaxis()->SetTitle("Counts/0.02 MeV bin");
		MC->Draw();
		RE->Draw("same");
		TLegend * leg =new TLegend(0.7,0.7,0.9,0.9);
		leg->AddEntry(MC,"MC","f");
		leg->AddEntry(RE,"RECON","f");
		leg->Draw();
		c1->Print("EnergyComparsion.png");
		c1->Print("EnergyComparsion.tex");
	}
	{
		TCanvas * c1 = new TCanvas();
		c1->SetGrid();
		
		gStyle->SetOptStat(false);
		// DIF->SetFillColorAlpha(kGreen,0.1);
		DIF->SetLineWidth(2);
		DIF->SetTitle("MC-RECON Energy");
		DIF->GetXaxis()->SetTitle("Energy (MeV)");
		DIF->GetYaxis()->SetTitle("Counts/0.02 MeV bin");
		DIF->Draw();
		// TLegend * leg =new TLegend(0.7,0.7,0.9,0.9);
		// leg->AddEntry(MC,"MC","f");
		// leg->AddEntry(RE,"RECON","f");
		// leg->Draw();
		c1->Print("EnergyDifference.png");
		c1->Print("EnergyDifference.tex");
	}

}


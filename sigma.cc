/*
 *
 * TODO:
 * 	-You need to find the rejection and efficacies.
 *
 */
#include "util.h"
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
#include <TAxis.h>


#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <iostream>
#include <fstream>



#include "util.h"

/*
 * The point of this is that you want to use root instead of echidna to look at the energy recostruction over z dimension
 * the way to do this is to give a function a couple of histrograms.
 */
/*Function signitures*/
/********************/

#define SSTR( x ) static_cast< std::ostringstream & >( \
		        ( std::ostringstream() << std::dec << x ) ).str()

string partFlag;



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


void FillHist(string filename,TH2D* hDeltaE, TH2D * hNormE){
	//cout<<"filename : "<<filename<<endl;

	TFile * file = TFile::Open(filename.c_str());
        TTree* Tree = (TTree*) file->Get("output");
	Double_t energy;
	Double_t energyReco;
	Double_t mcPosz, mcPosr;
        Int_t nhits, evIndex;
        Bool_t fitValid;
        Int_t flag=0;
        Int_t pdg1, pdg2;
        Int_t parentpdg1,parentpdg2;

        Tree->SetBranchAddress("pdg1",&pdg1);
        Tree->SetBranchAddress("pdg2",&pdg2);
        Tree->SetBranchAddress("parentpdg1",&parentpdg1);
        Tree->SetBranchAddress("parentpdg2",&parentpdg2);
        Tree->SetBranchAddress("nhits",&nhits);
	Tree->SetBranchAddress("mcEdepQuenched",&energy);
	Tree->SetBranchAddress("energy",&energyReco);
        Tree->SetBranchAddress("fitValid",&fitValid);
        Tree->SetBranchAddress("evIndex",&evIndex);
        Tree->SetBranchAddress("mcPosz",&mcPosz);
        Tree->SetBranchAddress("mcPosr",&mcPosr);
        Int_t n = (Int_t)Tree->GetEntries();
        //Int_t n = 9;

	if ( partFlag=="Bi") flag=11;
	if ( partFlag=="Po" ) flag=1000020040;

	for( Int_t iEntry = 0; iEntry < n; iEntry++ ){
		Tree->GetEntry(iEntry);
		if (fitValid !=0 && energyReco>0.2 && evIndex==0 && pdg1==flag && mcPosr<5000){// || pdg2==flag)){
				hDeltaE->Fill(energy,energy-energyReco);	
				hNormE->Fill(energy,(energy-energyReco)/energy);
		}
	}
	file->Close();
	
}

/* void FillHistrograms(string filename,double E_low,double E_high,TH1D* hNhits,TH1D* hRes,TH1D* hRMS,TH1D* hDelta, TH1D * normal){ */
/* 	//cout<<"filename : "<<filename<<endl; */

/* 	TFile * file = TFile::Open(filename.c_str()); */
/*         TTree* Tree = (TTree*) file->Get("output"); */
/* 	Double_t energy; */
/* 	Double_t energyReco; */
/* 	Double_t mcPosz, mcPosr; */
/*         Int_t nhits, evIndex; */
/*         Bool_t fitValid; */
/*         Int_t flag=0; */
/*         Int_t pdg1, pdg2; */
/*         Int_t parentpdg1,parentpdg2; */

/*         Tree->SetBranchAddress("pdg1",&pdg1); */
/*         Tree->SetBranchAddress("pdg2",&pdg2); */
/*         Tree->SetBranchAddress("parentpdg1",&parentpdg1); */
/*         Tree->SetBranchAddress("parentpdg2",&parentpdg2); */
/*         Tree->SetBranchAddress("nhits",&nhits); */
/* 	Tree->SetBranchAddress("mcEdepQuenched",&energy); */
/* 	Tree->SetBranchAddress("energy",&energyReco); */
/*         Tree->SetBranchAddress("fitValid",&fitValid); */
/*         Tree->SetBranchAddress("evIndex",&evIndex); */
/*         Tree->SetBranchAddress("mcPosz",&mcPosz); */
/*         Tree->SetBranchAddress("mcPosr",&mcPosr); */
/*         Int_t n = (Int_t)Tree->GetEntries(); */
/*         //Int_t n = 9; */

/* 	if ( partFlag=="Bi") flag=11; */
/* 	if ( partFlag=="Po" ) flag=1000020040; */

/* 	for( Int_t iEntry = 0; iEntry < n; iEntry++ ){ */
/* 		Tree->GetEntry(iEntry); */
/* 		if (fitValid !=0 && energyReco>0.2 && evIndex==0 && pdg1==flag && mcPosr<5000){// || pdg2==flag)){ */
/* 			if(E_low<energy && energy< E_high){ */
/* 				hNhits->Fill(nhits); */
/* 				hRes->Fill(1/sqrt(nhits)); */
/* 				hDelta->Fill(energy-energyReco); */	
/* 				normal->Fill((energy-energyReco)/energy); */
/* 			} */
/* 		} */
/* 	} */
/* 	file->Close(); */
	
/* } */

int main(){
	//gStyle->SetOptStat(0);

	int n=25;

	vector<string> pureAlpha= glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output/ntuple","alpha");
	vector<string> pureElectron= glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output_electron/ntuple","electron");
	double E_bins=100;
	double E_low=0;
	double E_high=10;

	//YOU SHOULD MAKE IT LIKE THE abB STUFF AND TAKE THE FOR LOOP LIMITS FROM THE TH1D LIMTS.
	TH2D* Histrogram_deltaE_bi  = new TH2D("Histrogram_deltaE_bi","",E_bins,E_low,E_high,100,-0.4,0.4);
	Histrogram_deltaE_bi->SetLineColor(kBlue);Histrogram_deltaE_bi->SetLineWidth(2);

	TH2D* Histrogram_deltaE_po  = new TH2D("Histrogram_deltaE_po","",E_bins,E_low,E_high,100,-0.4,0.4);
	Histrogram_deltaE_po->SetLineColor(kRed);Histrogram_deltaE_po->SetLineWidth(2);

	TH2D* Histrogram_NormE_bi  = new TH2D("Histrogram_NormE_bi","",E_bins,E_low,E_high,100,-0.4,0.4);
	Histrogram_NormE_bi->SetLineColor(kBlue);Histrogram_NormE_bi->SetLineWidth(2);

	TH2D* Histrogram_NormE_po  = new TH2D("Histrogram_NormE_po","",E_bins,E_low,E_high,100,-0.4,0.4);
	Histrogram_NormE_po->SetLineColor(kRed);Histrogram_NormE_po->SetLineWidth(2);

	for( int i=0; i<pureElectron.size(); i++ ){
		partFlag="Bi";
		FillHist(pureElectron[i],Histrogram_deltaE_bi,Histrogram_NormE_po);
	  }
	
	for( int i=0; i<pureAlpha.size(); i++ ){
		partFlag="Po";
		FillHist(pureAlpha[i],Histrogram_deltaE_po,Histrogram_NormE_bi);
	  }

	double E[n],delta_bi_RMS[n],delta_po_RMS[n],NormE_bi_RMS[n],NormE_po_RMS[n];
	double error_E[n],delta_bi_error[n],delta_po_error[n],NormE_RMS_bi_error[n],NormE_RMS_po_error[n];

	/* TCanvas* TotalCanvas = new TCanvas(); */
	/* TCanvas* TotalCanvasPo = new TCanvas(); */
	/* cout<<"before TCanvas "<<endl; */
	TCanvas* TotalCanvas_deltaE = new TCanvas();
	TCanvas* TotalCanvasPo_deltaE = new TCanvas();
	/* cout<<"After TCanvas "<<endl; */

	for (double energy=0;energy<25;energy++){

		double energy_real= energy/10;
		E[(int)energy]=(energy/10)+0.05;
		error_E[(int)energy]=0;

		TH1D * slicedBi_deltaE = Histrogram_deltaE_bi->ProjectionY("slicedBi_deltaE", Histrogram_deltaE_bi->GetXaxis()->FindBin(energy_real),Histrogram_deltaE_bi->GetXaxis()->FindBin(energy_real+0.1));
		TH1D * slicedPo_deltaE = Histrogram_deltaE_po->ProjectionY("slicedPo_deltaE", Histrogram_deltaE_po->GetXaxis()->FindBin(energy_real),Histrogram_deltaE_po->GetXaxis()->FindBin(energy_real+0.1));

		TH1D * slicedBi_normE = Histrogram_NormE_bi->ProjectionY("slicedBi_normE", Histrogram_NormE_bi->GetXaxis()->FindBin(energy_real),Histrogram_NormE_bi->GetXaxis()->FindBin(energy_real+0.1));
		TH1D * slicedPo_normE = Histrogram_NormE_po->ProjectionY("slicedPo_normE", Histrogram_NormE_po->GetXaxis()->FindBin(energy_real),Histrogram_NormE_po->GetXaxis()->FindBin(energy_real+0.1));

		bool first=true;
		if(slicedBi_normE->GetEntries()>200 && slicedPo_normE->GetEntries()>200  ){	


			// slicedBi_deltaE ->Sumw2();
			// slicedPo_deltaE->Sumw2();
			// Double_t norm =slicedBi_deltaE ->GetEntries();
			// slicedBi_deltaE ->Scale(1/norm);
			// norm = slicedPo_deltaE->GetEntries();
			// slicedPo_deltaE->Scale(1/norm);

			slicedBi_deltaE ->Fit("gaus");
			slicedPo_deltaE->Fit("gaus");

			TF1 *fit1 = (TF1*)slicedBi_deltaE ->GetFunction("gaus");
			TF1 *fit2 = (TF1*)slicedPo_deltaE->GetFunction("gaus");
			fit1->SetLineColor(kBlack);
			fit2->SetLineColor(kBlack);

			if (true){
				slicedBi_deltaE ->GetXaxis()->SetTitle("#Delta E");
				slicedPo_deltaE->GetXaxis()->SetTitle("#Delta E");
				slicedBi_deltaE ->SetTitle(("deltaE_bi_"+SSTR(energy_real)+" < E <"+SSTR(energy_real+0.1)+" MeV").c_str());
				slicedPo_deltaE->SetTitle(("deltaE_po_"+SSTR(energy_real)+" < E <"+SSTR(energy_real+0.1)+" MeV").c_str());
				TCanvas* checkerCanvas = new TCanvas();
				slicedBi_deltaE ->Draw();
				fit1->Draw("same");
				checkerCanvas->Print(("plotChecker/bi/DeltaE_between_"+SSTR(energy_real)+"0_to"+SSTR(energy_real+0.1)+".png").c_str());

				TCanvas* checkerCanvas_2= new TCanvas();
				slicedPo_deltaE->Draw();
				fit2->Draw("same");
				checkerCanvas_2->Print(("plotChecker/po/DeltaE_between_"+SSTR(energy_real)+"0_to"+SSTR(energy_real+0.1)+".png").c_str());
			}
			int index = (int)energy;
			delta_bi_RMS[index]=fit1->GetParameter(2);
			delta_po_RMS[index]=fit2->GetParameter(2);
			delta_bi_error[index]=fit1->GetParError(2);
			delta_po_error[index]=fit2->GetParError(2);

			if(false){
				//This is to compute the all on one plots.
				TotalCanvas_deltaE->cd();
				if (first=true)slicedBi_deltaE ->Draw();
				slicedBi_deltaE ->Draw("same");

				TotalCanvasPo_deltaE->cd();
				if (first=true) slicedPo_deltaE->Draw();
				slicedPo_deltaE->Draw("same");
				first=false;
			}
		
//+++++++++++++++++++++++++NormE++++++++++++++++++++++++++++



			// slicedBi_normE->Sumw2();
			// slicedPo_normE->Sumw2();
			// norm = slicedBi_normE->GetEntries();
			// slicedBi_normE->Scale(1/norm);
			// norm = slicedPo_normE->GetEntries();
			// slicedPo_normE->Scale(1/norm);

			slicedBi_normE->Fit("gaus");
			slicedPo_normE->Fit("gaus");

			fit1 = (TF1*)slicedBi_normE->GetFunction("gaus");
			fit2 = (TF1*)slicedPo_normE->GetFunction("gaus");
			fit1->SetLineColor(kBlack);
			fit2->SetLineColor(kBlack);

			if (true){
				slicedBi_normE->GetXaxis()->SetTitle("#Delta E");
				slicedPo_normE->GetXaxis()->SetTitle("#Delta E/E_{MC}");
				slicedBi_normE->SetTitle(("NormE_bi_"+SSTR(energy_real)+" < E <"+SSTR(energy_real+0.1)+" MeV").c_str());
				slicedPo_normE->SetTitle(("NormE_po_"+SSTR(energy_real)+" < E <"+SSTR(energy_real+0.1)+" MeV").c_str());
				TCanvas* checkerCanvas = new TCanvas();
				slicedBi_normE->Draw();
				fit1->Draw("same");
				checkerCanvas->Print(("plotChecker/bi/NormE_between_"+SSTR(energy_real)+"0_to"+SSTR(energy_real+0.1)+".png").c_str());

				TCanvas* checkerCanvas_2= new TCanvas();
				slicedPo_normE->Draw();
				fit2->Draw("same");
				checkerCanvas_2->Print(("plotChecker/po/NormE_between_"+SSTR(energy_real)+"0_to"+SSTR(energy_real+0.1)+".png").c_str());
			}

			NormE_bi_RMS[(int)energy]=fit1->GetParameter(2);
			NormE_po_RMS[(int)energy]=fit2->GetParameter(2);
			NormE_RMS_bi_error[(int)energy]=fit1->GetParError(2);
			NormE_RMS_po_error[(int)energy]=fit2->GetParError(2);

		}else{
			int index = (int)energy;
			delta_bi_RMS[index]=0;
			delta_po_RMS[index]=0;
			delta_bi_error[index]=0;
			delta_po_error[index]=0;

			NormE_bi_RMS[(int)energy]=0;
			NormE_po_RMS[(int)energy]=0;
			NormE_RMS_bi_error[(int)energy]=0;
			NormE_RMS_po_error[(int)energy]=0;
		}
	
	}

	/* TotalCanvas->Print("All_electrons_nhits.png"); */
	/* TotalCanvasPo->Print("All_alpha_nhits.png"); */
	/* TotalCanvas_deltaE->Print("All_electrons_deltaE.png"); */
	/* TotalCanvasPo_deltaE->Print("All_alpha_deltaE.png"); */


//=========================NhitVsEnergy==============================
	//        /* TCanvas* error= new TCanvas(); */
	/* TMultiGraph *mg = new TMultiGraph(); */
	/* error->cd(); */
        /* TGraphErrors * gr_nhits_bi =new TGraphErrors(n,E,nhits_bi,error_E,error_nhits_bi); */
	/* gr_nhits_bi ->SetMarkerColor(kBlue); */
	/* gr_nhits_bi->SetLineColor(kBlue); */
	/* gr_nhits_bi->SetFillColor(kBlue); */
        /* mg->Add(gr_nhits_bi); */

        /* TGraphErrors * gr_nhits_po =new TGraphErrors(n,E,nhits_po,error_E,error_nhits_po); */
	/* gr_nhits_po ->SetMarkerColor(kRed); */
	/* gr_nhits_po->SetLineColor(kRed); */
	/* gr_nhits_po->SetFillColor(kRed); */
        /* mg->Add(gr_nhits_po); */

	/* mg->Draw("ap"); */
/* //        mg->GetXaxis()->SetRangeUser(0.0,3.5); */
/* //        mg->GetYaxis()->SetRangeUser(0.0,800); */
	/* mg->SetTitle("Mean nhits across energy"); */
        /* mg->GetXaxis()->SetTitle("Energy MC (MeV)"); */
       	/* mg->GetYaxis()->SetTitle("nhits"); */

	/* TLegend * errorLeg = new TLegend(0.7,0.7,0.9,0.9); */
	/* errorLeg->AddEntry(gr_nhits_bi,"Beta","f"); */
	/* errorLeg->AddEntry(gr_nhits_po,"Alpha","f"); */
	/* //errorLeg->AddEntry(gr2,"Scaled Po","f"); */
	/* errorLeg->Draw(); */

	/* error->Print("NhitsVsEnergy.png"); */

//==============RMS===================================================================================

        TCanvas* RMSCan= new TCanvas();
	TMultiGraph *mg_RMS = new TMultiGraph();
	RMSCan->cd();
        TGraphErrors * gr_RMS_bi =new TGraphErrors(n,E,delta_bi_RMS,error_E,delta_bi_error);
	gr_RMS_bi ->SetMarkerColor(kBlue);
	gr_RMS_bi->SetLineColor(kBlue);
	gr_RMS_bi->SetFillColor(kBlue);
        mg_RMS->Add(gr_RMS_bi);

        TGraphErrors * gr_RMS_po =new TGraphErrors(n,E,delta_po_RMS,error_E,delta_po_error);
	gr_RMS_po ->SetMarkerColor(kRed);
	gr_RMS_po->SetLineColor(kRed);
	gr_RMS_po->SetFillColor(kRed);
        mg_RMS->Add(gr_RMS_po);


	mg_RMS->Draw("a.");
//        mg_res->GetXaxis()->SetRangeUser(0.0,3.5);
//        mg_res->GetYaxis()->SetRangeUser(0.0,0.12);
	mg_RMS->SetTitle("#Delta E RMS across energy");
        mg_RMS->GetXaxis()->SetTitle("Energy MC (MeV)");
       	mg_RMS->GetYaxis()->SetTitle("#Delta E RMS");

	TLegend * RMSLeg = new TLegend(0.7,0.7,0.9,0.9);
	RMSLeg->AddEntry(gr_RMS_bi,"Beta","f");
	RMSLeg->AddEntry(gr_RMS_po,"Alpha","f");
	//errorLeg->AddEntry(gr2,"Scaled Po","f");
	RMSLeg->Draw();

	RMSCan->Print("deltaE_RMSVsEnergy.png");
	
//==============NormRMS===================================================================================

        TCanvas* NormRMSCan= new TCanvas();
	TMultiGraph *mg_NormE_RMS = new TMultiGraph();
	NormRMSCan->cd();
        TGraphErrors * gr_NormE_RMS_bi =new TGraphErrors(n,E,NormE_bi_RMS,error_E,NormE_RMS_bi_error);
	gr_NormE_RMS_bi ->SetMarkerColor(kBlue);
	gr_NormE_RMS_bi->SetLineColor(kBlue);
	gr_NormE_RMS_bi->SetFillColor(kBlue);
        mg_NormE_RMS->Add(gr_NormE_RMS_bi);

        TGraphErrors * gr_NormE_RMS_po =new TGraphErrors(n,E,NormE_po_RMS,error_E,NormE_RMS_po_error);
	gr_NormE_RMS_po ->SetMarkerColor(kRed);
	gr_NormE_RMS_po->SetLineColor(kRed);
	gr_NormE_RMS_po->SetFillColor(kRed);
        mg_NormE_RMS->Add(gr_NormE_RMS_po);


	mg_NormE_RMS->Draw("a.");
//        mg_res->GetXaxis()->SetRangeUser(0.0,3.5);
//        mg_res->GetYaxis()->SetRangeUser(0.0,0.12);
	mg_NormE_RMS->SetTitle("#Delta E/E_{MC} RMS across energy");
        mg_NormE_RMS->GetXaxis()->SetTitle("Energy MC (MeV)");
       	mg_NormE_RMS->GetYaxis()->SetTitle("#Delta E/E_{MC} RMS");

	TLegend * NormE_RMSLeg = new TLegend(0.7,0.7,0.9,0.9);
	NormE_RMSLeg->AddEntry(gr_NormE_RMS_bi,"Beta","f");
	NormE_RMSLeg->AddEntry(gr_NormE_RMS_po,"Alpha","f");
	//errorLeg->AddEntry(gr2,"Scaled Po","f");
	NormE_RMSLeg->Draw();

	NormRMSCan->Print("NormE_RMSVsEnergy.png");
	return 0;
}

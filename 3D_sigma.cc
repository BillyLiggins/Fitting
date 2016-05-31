/*
 *
 * TODO:
 * 	-You need to find the rejection and efficacies.
 *
 */
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
#include <TH3D.h>
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
#include <TF1.h>

#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <iostream>
#include <fstream>

#define SSTR( x ) static_cast< std::ostringstream & >( \
		        ( std::ostringstream() << std::dec << x ) ).str()
/*
 * the point of this is that you want to use root instead of echidna to look at the energy recostruction over z dimension
 * the way to do this is to give a function a couple of histrograms.
 */
/*function signitures*/
/********************/

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

//================================================================================================================

void FillHistrograms(string filename,TH3D *hist){
	//cout<<"filename : "<<filename<<endl;

	TFile * file = TFile::Open(filename.c_str());
        TTree* Tree = (TTree*) file->Get("output");
	Double_t energy;
	Double_t energyReco;
	Double_t mcPosz, mcPosr,posr;
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
        Tree->SetBranchAddress("posr",&posr);
        Int_t n = (Int_t)Tree->GetEntries();
        //Int_t n = 9;

	if ( partFlag=="Bi") flag=11;
	if ( partFlag=="Po" ) flag=1000020040;

	for( Int_t iEntry = 0; iEntry < n; iEntry++ ){
		Tree->GetEntry(iEntry);
		if (fitValid !=0 && energyReco>0.2 && evIndex==0 && pdg1==flag && mcPosr<5000){// || pdg2==flag)){
			double deltaE=(energy-energyReco)/energy;
				hist->Fill(energy,posr,deltaE);
		}
	}
	file->Close();
	
}

void second(){
	//gStyle->SetOptStat(0);


	/* vector<string> bi210FileList= glob("/data/snoplus/OfficialProcessing/production_5_3_0/Bi210","SolarBi210_r"); */
	/* vector<string> po210FileList= glob("/data/snoplus/OfficialProcessing/production_5_3_0/Po210","SolarPo210_r"); */

	/* vector<string> bi212FileList= glob("/data/snoplus/OfficialProcessing/production_5_3_0/Bi212","SolarBi212_r"); */
	/* vector<string> po212FileList= glob("/data/snoplus/OfficialProcessing/production_5_3_0/Po212","SolarPo212_r"); */
//	vector<string> bi212FileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Bi212","SolarBi212_r");
//	vector<string> po212FileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po212","SolarPo212_r");

	/* vector<string> bi214FileList= glob("/data/snoplus/OfficialProcessing/production_5_3_0/Bi214","SolarBi214_r"); */
	/* vector<string> po214FileList= glob("/data/snoplus/OfficialProcessing/production_5_3_0/Po214","SolarPo214_r"); */
//	vector<string> bi214FileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Bi214","SolarBi214_r");
//	vector<string> po214FileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po214","SolarPo214_r");

	vector<string> pureAlpha= glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output/ntuple","alpha");
	vector<string> pureElectron= glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output_electron/ntuple","electron");


	TH3D* hist_complete_ele= new TH3D("hist_complete_ele","",25,0,2.5,6,0,6000,100,-0.4,0.4);
	hist_complete_ele->SetLineColor(kBlue);
	hist_complete_ele->SetLineWidth(2);

	TH3D* hist_complete_alp= new TH3D("hist_complete_alp","",25,0,2.5,6,0,6000,100,-0.4,0.4);
	hist_complete_alp->SetLineColor(kBlue);
	hist_complete_alp->SetLineWidth(2);



	for( int i=0; i<pureElectron.size(); i++ ){
		partFlag="Bi";
		FillHistrograms(pureElectron[i],hist_complete_ele);
	  }
	
	for( int i=0; i<pureAlpha.size(); i++ ){
		partFlag="Po";
		FillHistrograms(pureAlpha[i],hist_complete_alp);
	  }
	double radialStep=1000;
	double energyStep=0.1;

		cout<<"hist_complete_ele = "<<hist_complete_ele->GetEntries()<<endl;
		cout<<"hist_complete_alp = "<<hist_complete_alp->GetEntries()<<endl;

	for(double radius=0;radius<6000;radius+=radialStep){

		double lower_radius=radius;
		double upper_radius=lower_radius+radialStep;
	
		//How many energy slices.
		int n=25;
		double E[n],NormE_ele[n],NormE_alp[n];
		double E_error[n],NormE_ele_error[n],NormE_alp_error[n];

		for (double energy=0;energy<n;energy++){
			/* TH1::SetDefaultSumw2(); */
			double lower_energy=energy/10;
			double upper_energy=lower_energy+energyStep;

			cout<<"Energy = "<<lower_energy<<" < E <"<<upper_energy<<endl;
			cout<<"Radius = "<<lower_radius<<" < R <"<<upper_radius<<endl;
			E[(int)energy]=(energy/10)+0.05;
			E_error[(int)energy]=0;

			double energy_real= energy/10;

			TAxis* x_axis = hist_complete_ele->GetXaxis();
			TAxis* y_axis = hist_complete_ele->GetYaxis();
			TH1D* hist_sample_ele;
			hist_sample_ele=hist_complete_ele->ProjectionZ("hist_sample_ele",x_axis->FindBin(lower_energy),x_axis->FindBin(upper_energy), y_axis->FindBin(lower_radius), y_axis->FindBin(upper_radius));
			hist_sample_ele->Sumw2();
			
			cout<<"Printed the sample plot"<<endl;
			x_axis =hist_complete_alp->GetXaxis();
			y_axis =hist_complete_alp->GetYaxis();
			TH1D* hist_sample_alp=hist_complete_alp->ProjectionZ("hist_sample_alp",x_axis->FindBin(lower_energy),x_axis->FindBin(upper_energy), y_axis->FindBin(lower_radius), y_axis->FindBin(upper_radius));
			hist_sample_alp->Sumw2();

			TCanvas *c1=new TCanvas();
			c1->cd();
			hist_sample_ele->Draw();
			hist_sample_alp->Draw("same");
			c1->Print("sample.png");

			bool first=true;


			cout<<"hist_sample_ele = "<<hist_sample_ele->GetEntries()<<endl;
			cout<<"hist_sample_alp = "<<hist_sample_alp->GetEntries()<<endl;
			if(hist_sample_ele->GetEntries()>100 && hist_sample_alp->GetEntries()>100){	

				Double_t norm =hist_sample_ele->GetEntries();
				hist_sample_ele->Scale(1/norm);
				norm =hist_sample_alp->GetEntries();
				hist_sample_alp->Scale(1/norm);

				TCanvas *c1=new TCanvas();
				c1->cd();
				hist_sample_ele->Draw();
				hist_sample_alp->Draw("same");
				c1->Print("sample_inside.png");

				cout<<"from above fit hist_sample_ele = "<<hist_sample_ele->GetEntries()<<endl;
				cout<<"from above fit hist_sample_alp = "<<hist_sample_alp->GetEntries()<<endl;

				hist_sample_ele->Fit("gaus");
				TF1* fit1 = (TF1*)hist_sample_ele->GetFunction("gaus");
				NormE_ele[(int)energy]=fit1->GetParameter(2);
				cout<<"After beta fit"<<endl;

				hist_sample_alp->Fit("gaus");
				TF1* fit2 = (TF1*)hist_sample_alp->GetFunction("gaus");
				NormE_alp[(int)energy]=fit2->GetParameter(2);
				cout<<"After alpaha fit"<<endl;

				/* fit1->SetLineColor(kBlack); */
				/* fit2->SetLineColor(kBlack); */

				if (true){
					cout<<"Printing plot checkers"<<endl;
					TCanvas* checkerCanvas = new TCanvas();
					hist_sample_ele->Draw();
				cout<<"GotHere 2"<<endl;
					fit1->Draw("same");
				cout<<"GotHere 2"<<endl;
					checkerCanvas->Print(("shellPlots/bi/NormE_between_"+SSTR(energy_real)+"0_to"+SSTR(energy_real+0.1)+".png").c_str());

					TCanvas* checkerCanvas_2= new TCanvas();
					hist_sample_alp->Draw();
					fit2->Draw("same");
					checkerCanvas_2->Print(("shellPlots/po/NormE_between_"+SSTR(energy_real)+"0_to"+SSTR(energy_real+0.1)+".png").c_str());
				}

				cout<<"GotHere above getParameter"<<endl;

				NormE_ele_error[(int)energy]=fit1->GetParError(2);
				NormE_alp_error[(int)energy]=fit2->GetParError(2);

				cout<<"GotHere below getParameter"<<endl;
				/* if(false){ */
				/* 	//This is to compute the all on one plots. */
				/* 	TotalCanvas_deltaE->cd(); */
				/* 	if (first=true) Histrogram_deltaE_bi->Draw(); */
				/* 	Histrogram_deltaE_bi->Draw("same"); */

				/* 	TotalCanvasPo_deltaE->cd(); */
				/* 	if (first=true) Histrogram_deltaE_po->Draw(); */
				/* 	Histrogram_deltaE_po->Draw("same"); */
				/* 	first=false; */
				/* } */
			
	//+++++++++++++++++++++++++NormE++++++++++++++++++++++++++++

			}else{
				NormE_ele[(int)energy]=0;
				NormE_alp[(int)energy]=0;
				NormE_ele_error[(int)energy]=0;
				NormE_alp_error[(int)energy]=0;
			}
		
		}

		/* TotalCanvas->Print("All_electrons_nhits.png"); */
		/* TotalCanvasPo->Print("All_alpha_nhits.png"); */
		/* TotalCanvas_deltaE->Print("All_electrons_deltaE.png"); */
		/* TotalCanvasPo_deltaE->Print("All_alpha_deltaE.png"); */


		TCanvas* NormRMSCan= new TCanvas();
		TMultiGraph *mg_NormE_RMS = new TMultiGraph();
		NormRMSCan->cd();
		TGraphErrors * gr_NormE_RMS_bi =new TGraphErrors(n,E,NormE_ele,E_error,NormE_ele_error);
		gr_NormE_RMS_bi ->SetMarkerColor(kBlue);
		gr_NormE_RMS_bi->SetLineColor(kBlue);
		gr_NormE_RMS_bi->SetFillColor(kBlue);
		mg_NormE_RMS->Add(gr_NormE_RMS_bi);

		TGraphErrors * gr_NormE_RMS_po =new TGraphErrors(n,E,NormE_alp,E_error,NormE_alp_error);
		gr_NormE_RMS_po ->SetMarkerColor(kRed);
		gr_NormE_RMS_po->SetLineColor(kRed);
		gr_NormE_RMS_po->SetFillColor(kRed);
		mg_NormE_RMS->Add(gr_NormE_RMS_po);


		mg_NormE_RMS->Draw("a*");
	//        mg_res->GetXaxis()->SetRangeUser(0.0,3.5);
	//        mg_res->GetYaxis()->SetRangeUser(0.0,0.12);
		mg_NormE_RMS->SetTitle(("NormE {"+SSTR(lower_radius)+"< posr <"+SSTR(upper_radius)+"}").c_str());
		mg_NormE_RMS->GetXaxis()->SetTitle("Energy MC (MeV)");
		mg_NormE_RMS->GetYaxis()->SetTitle("#Delta E/E_{MC} RMS");

		TLegend * NormE_RMSLeg = new TLegend(0.7,0.7,0.9,0.9);
		NormE_RMSLeg->AddEntry(gr_NormE_RMS_bi,"Beta","f");
		NormE_RMSLeg->AddEntry(gr_NormE_RMS_po,"Alpha","f");
		//errorLeg->AddEntry(gr2,"Scaled Po","f");
		NormE_RMSLeg->Draw();

		/* NormRMSCan->Print(("shellPlots/NormE_E_ "+SSTR(lower_energy)+"_to_"+SSTR(upper_energy)+"_r_"+SSTR(lower_radius)+"_to_"+SSTR(upper_radius)+".png").c_str()); */
		NormRMSCan->Print(("shellPlots/NormE_r_"+SSTR(lower_radius)+"_to_"+SSTR(upper_radius)+".png").c_str());
	}




}




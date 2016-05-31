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

#define SSTR( x ) static_cast< std::ostringstream & >( \
		        ( std::ostringstream() << std::dec << x ) ).str()
/*
 * The point of this is that you want to use root instead of echidna to look at the energy recostruction over z dimension
 * the way to do this is to give a function a couple of histrograms.
 */
/*Function signitures*/
void nhits(string filename, TH1D* sigma,TH2D* plot);
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


void nhitsMod(string filename,TH2D* plot,Double_t scale){
	cout<<"filename : "<<filename<<endl;

	TFile * file = TFile::Open(filename.c_str());
        TTree* Tree = (TTree*) file->Get("output");
	Double_t energy;
	Double_t energyReco;
	Double_t mcPosz;
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
        Int_t n = (Int_t)Tree->GetEntries();

	if ( partFlag=="Bi") flag=11;
	if ( partFlag=="Po" ) flag=1000020040;

	for( Int_t iEntry = 0; iEntry < n; iEntry++ ){
			Tree->GetEntry(iEntry);
			if (fitValid !=0 && energyReco>0.2 && evIndex==0 && pdg1==flag && mcPosz<5000){
				plot->Fill(energy,scale/sqrt(nhits));
			}
	}

}

void nhitsMod(string filename,TH2D* plot){
	nhitsMod(filename, plot,1.0);
}

void FindNhits(string filename,TH2D* plot){
	cout<<"filename : "<<filename<<endl;

	TFile * file = TFile::Open(filename.c_str());
        TTree* Tree = (TTree*) file->Get("output");
	Double_t energy;
	Double_t energyReco;
	Double_t mcPosz;
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
        Int_t n = (Int_t)Tree->GetEntries();
        //Int_t n = 9;

	if ( partFlag=="Bi") flag=11;
	if ( partFlag=="Po" ) flag=1000020040;

	for( Int_t iEntry = 0; iEntry < n; iEntry++ ){
			Tree->GetEntry(iEntry);
			if (fitValid !=0 && energyReco>0.2 && evIndex==0 && pdg1==flag && mcPosz<5000){// || pdg2==flag)){
				plot->Fill(energy, nhits);

			}
	}

}
void nhits(string filename, TH1D* sigma,TH2D* plot){
	cout<<"filename : "<<filename<<endl;

	TFile * file = TFile::Open(filename.c_str());
        TTree* Tree = (TTree*) file->Get("output");
	Double_t energy;
	Double_t energyReco;
	Double_t mcPosz;
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
        Int_t n = (Int_t)Tree->GetEntries();

	if ( partFlag=="Bi") flag=11;
	if ( partFlag=="Po" ) flag=1000020040;

	for( Int_t iEntry = 0; iEntry < n; iEntry++ ){
			Tree->GetEntry(iEntry);
			if (fitValid !=0 && energyReco>0.2 && evIndex==0 && pdg1==flag && mcPosz<5000){// || pdg2==flag)){
				//cout<<"The flag is "<<flag<<endl;
				sigma->Fill(1/sqrt(nhits) );
				plot->Fill(energy, 1/sqrt(nhits));
			}
	}
}

void spray(string biString,string poString){
	gStyle->SetOptStat(0);

	int n=35;
	TH1D* bi210Hits = new TH1D("bi210Hits","",100,0,0.5);
	TH2D* bi210HitsVsEnergy = new TH2D("bi210HitsVsEnergy","",35,0,3.5,100,0,0.2);
	bi210HitsVsEnergy->SetMarkerColor(kRed);bi210HitsVsEnergy->SetFillColor(kRed);bi210HitsVsEnergy->SetLineColor(kRed);

	TH1D* po210Hits = new TH1D("po210Hits","",100,0,0.5);
	TH2D* po210HitsVsEnergy = new TH2D("po210HitsVsEnergy","",35,0,3.5,100,0,0.2);
	po210HitsVsEnergy->SetMarkerColor(kRed);po210HitsVsEnergy->SetFillColor(kRed);po210HitsVsEnergy->SetLineColor(kRed);

	TH1D* bi212Hits = new TH1D("bi212Hits","",100,0,0.5);
	TH2D* bi212HitsVsEnergy = new TH2D("bi212HitsVsEnergy","",35,0,3.5,100,0,0.2);
	bi212HitsVsEnergy->SetMarkerColor(kBlue);bi212HitsVsEnergy->SetFillColor(kBlue);bi212HitsVsEnergy->SetLineColor(kBlue);

	TH1D* po212Hits = new TH1D("po212Hits","",100,0,0.5);
	TH2D* po212HitsVsEnergy = new TH2D("po212HitsVsEnergy","",35,0,3.5,100,0,0.2);
	po212HitsVsEnergy->SetMarkerColor(kBlue);po212HitsVsEnergy->SetFillColor(kBlue);po212HitsVsEnergy->SetLineColor(kBlue);

	TH1D* bi214Hits = new TH1D("bi214Hits","",100,0,0.5);
	TH2D* bi214HitsVsEnergy = new TH2D("bi214HitsVsEnergy","",35,0,3.5,100,0,0.2);
	bi214HitsVsEnergy->SetMarkerColor(kBlack);bi214HitsVsEnergy->SetFillColor(kBlack);bi214HitsVsEnergy->SetLineColor(kBlack);

	TH1D* po214Hits = new TH1D("po214Hits","",100,0,0.5);
	TH2D* po214HitsVsEnergy = new TH2D("po214HitsVsEnergy","",35,0,3.5,100,0,0.2);
	po214HitsVsEnergy->SetMarkerColor(kGreen+3); po214HitsVsEnergy->SetFillColor(kGreen+3);po214HitsVsEnergy->SetLineColor(kGreen+3);

	TH2D* po210HitsVsEnergyMod = new TH2D("po210HitsVsEnergyMod","",35,0,3.5,40,0,0.2);
	po210HitsVsEnergyMod->SetMarkerColor(kRed);po210HitsVsEnergyMod->SetFillColor(kRed);po210HitsVsEnergyMod->SetLineColor(kRed);

	TH2D* po212HitsVsEnergyMod = new TH2D("po212HitsVsEnergyMod","",35,0,3.5,40,0,0.2);
	po212HitsVsEnergyMod->SetMarkerColor(kBlue);po212HitsVsEnergyMod->SetFillColor(kBlue);po212HitsVsEnergyMod->SetLineColor(kBlue);

	TH2D* po214HitsVsEnergyMod = new TH2D("po214HitsVsEnergyMod","",35,0,3.5,40,0,0.2);
	po214HitsVsEnergyMod ->SetMarkerColor(kGreen+3); po214HitsVsEnergyMod->SetFillColor(kGreen+3);po214HitsVsEnergyMod->SetLineColor(kGreen+3);

	TH2D* bi210NHitsVsEnergy = new TH2D("bi210NHitsVsEnergy","",35,0,3.5,2000,0,2000);
	bi210NHitsVsEnergy->SetMarkerColor(kRed); bi210NHitsVsEnergy->SetFillColor(kRed);bi210NHitsVsEnergy->SetLineColor(kRed);
	bi210NHitsVsEnergy->SetLineWidth(3);

	TH2D* bi212NHitsVsEnergy = new TH2D("bi212NHitsVsEnergy","",35,0,3.5,2000,0,2000);
	bi212NHitsVsEnergy->SetMarkerColor(kBlue); bi212NHitsVsEnergy->SetFillColor(kBlue);bi212NHitsVsEnergy->SetLineColor(kBlue);
	bi212NHitsVsEnergy->SetLineWidth(3);

	TH2D* bi214NHitsVsEnergy = new TH2D("bi214NHitsVsEnergy","",35,0,3.5,2000,0,2000);
	bi214NHitsVsEnergy->SetMarkerColor(kGreen+3); bi214NHitsVsEnergy->SetFillColor(kGreen+3);bi214NHitsVsEnergy->SetLineColor(kGreen+3);
	bi214NHitsVsEnergy->SetLineWidth(3);

	TH2D* po210NHitsVsEnergy = new TH2D("po210NHitsVsEnergy","",35,0,3.5,2000,0,2000);
	po210NHitsVsEnergy->SetMarkerColor(kRed); po210NHitsVsEnergy->SetFillColor(kRed);po210NHitsVsEnergy->SetLineColor(kRed);
	po210NHitsVsEnergy->SetLineWidth(3);

	TH2D* po212NHitsVsEnergy = new TH2D("po212NHitsVsEnergy","",35,0,3.5,2000,0,2000);
	po212NHitsVsEnergy->SetMarkerColor(kBlue); po212NHitsVsEnergy->SetFillColor(kBlue);po212NHitsVsEnergy->SetLineColor(kBlue);
	po212NHitsVsEnergy->SetLineWidth(3);

	TH2D* po214NHitsVsEnergy = new TH2D("po214NHitsVsEnergy","",35,0,3.5,2000,0,2000);
	po214NHitsVsEnergy->SetMarkerColor(kGreen+3); po214NHitsVsEnergy->SetFillColor(kGreen+3);po214NHitsVsEnergy->SetLineColor(kGreen+3);
	po214NHitsVsEnergy->SetLineWidth(3);

	TH2D* allHist= new TH2D("allHist","",35,0,3.5,100,0,0.2);

	vector<string> bi210FileList= glob("/data/snoplus/OfficialProcessing/production_5_3_0/Bi210","SolarBi210_r");
	vector<string> po210FileList= glob("/data/snoplus/OfficialProcessing/production_5_3_0/Po210","SolarPo210_r");

	vector<string> bi212FileList= glob("/data/snoplus/OfficialProcessing/production_5_3_0/Bi212","SolarBi212_r");
	vector<string> po212FileList= glob("/data/snoplus/OfficialProcessing/production_5_3_0/Po212","SolarPo212_r");
//	vector<string> bi212FileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Bi212","SolarBi212_r");
//	vector<string> po212FileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po212","SolarPo212_r");

	vector<string> bi214FileList= glob("/data/snoplus/OfficialProcessing/production_5_3_0/Bi214","SolarBi214_r");
	vector<string> po214FileList= glob("/data/snoplus/OfficialProcessing/production_5_3_0/Po214","SolarPo214_r");
//	vector<string> bi214FileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Bi214","SolarBi214_r");
//	vector<string> po214FileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po214","SolarPo214_r");

	//for( int i=0; i<bi210FileList.size(); i++ ){
	for( int i=0; i<5; i++ ){
		partFlag="Bi";
		nhits( bi210FileList[i], bi210Hits,bi210HitsVsEnergy);
		FindNhits( bi210FileList[i],bi210NHitsVsEnergy);
	  }
	for( int i=0; i<bi212FileList.size(); i++ ){
		partFlag="Bi";
		nhits( bi212FileList[i], bi212Hits,bi212HitsVsEnergy);
		FindNhits( bi212FileList[i],bi212NHitsVsEnergy);
	  }
	for( int i=0; i<bi214FileList.size(); i++ ){
		partFlag="Bi";
		nhits( bi214FileList[i], bi214Hits,bi214HitsVsEnergy);
		FindNhits( bi214FileList[i],bi214NHitsVsEnergy);
	  }

	//for( int i=0; i<po210FileList.size(); i++ ){
	for( int i=0; i<5; i++ ){
		partFlag="Po";
		nhits( po210FileList[i], po210Hits,po210HitsVsEnergy);
		nhitsMod( po210FileList[i],po210HitsVsEnergyMod);
		FindNhits( po210FileList[i],po210NHitsVsEnergy);
	  }
	for( int i=0; i<po212FileList.size(); i++ ){
		partFlag="Po";
		nhits( po212FileList[i], po212Hits,po212HitsVsEnergy);
		nhitsMod( po212FileList[i],po212HitsVsEnergyMod);
		FindNhits( po212FileList[i],po212NHitsVsEnergy);
	  }
	for( int i=0; i<po214FileList.size(); i++ ){
		partFlag="Po";
		nhits( po214FileList[i], po214Hits,po214HitsVsEnergy);
		nhitsMod( po214FileList[i],po214HitsVsEnergyMod);
		FindNhits( po214FileList[i],po214NHitsVsEnergy);
	  }

	TAxis* eAxispo = po210HitsVsEnergy->GetXaxis();
	TAxis* nhitsAxispo =  po210HitsVsEnergy->GetYaxis();

//	TH1D* tester = po210HitsVsEnergy->ProjectionY("tester",eAxispo->FindBin(0.5),eAxispo->FindBin(0.6) );	
//	tester->Draw();


        double E[n],nhits[n],eE[n],enhits[n];
        double normalNhits[n], enormalNhits[n];

	TAxis* eAxis = bi210HitsVsEnergy->GetXaxis();
	TAxis* nhitsAxis = bi210HitsVsEnergy->GetYaxis();
	double temp,temp2;
	TCanvas* ccc=new TCanvas();
	TH1D* blankkk= new TH1D("blankkk","",600,0,600);
	blankkk->Fill(2.5);
	blankkk->Draw();
        for( double i=0; i<n; i++ ){
		//TH1D* projectX  =  bi210HitsVsEnergy->ProjectionX( "projectX", eAxis->FindBin( i/10 ), eAxis->FindBin((i+1)/10) );	
		TH1D* projectY  =  bi210HitsVsEnergy->ProjectionY( "projectY", eAxis->FindBin( i/10 ), eAxis->FindBin((i+1)/10) );	
		projectY->Add(bi212HitsVsEnergy->ProjectionY( "hold1", eAxis->FindBin( i/10 ), eAxis->FindBin((i+1)/10)));
		projectY->Add(bi214HitsVsEnergy->ProjectionY( "holder1", eAxis->FindBin( i/10 ), eAxis->FindBin((i+1)/10)));
		TH1D* projectY2  =  bi210NHitsVsEnergy->ProjectionY( "projectY2", eAxis->FindBin( i/10 ), eAxis->FindBin((i+1)/10) );	
		projectY2->Add(bi212NHitsVsEnergy->ProjectionY( "hold11", eAxis->FindBin( i/10 ), eAxis->FindBin((i+1)/10)));
		projectY2->Add(bi214NHitsVsEnergy->ProjectionY( "holder11", eAxis->FindBin( i/10 ), eAxis->FindBin((i+1)/10)));
		if (i =22) projectY2->Draw();
		//project->GetMean();
//		E[i]=projectX->GetMean();
//		eE[i]=projectX->GetRMS();
		E[(int)i]=(i/10)+0.05;
		eE[(int)i]=0;
		temp=projectY->GetMean();
		temp2=projectY->GetRMS();
		nhits[(int)i]=temp;
		enhits[(int)i]=temp2;
		normalNhits[(int)i]=projectY2->GetMean();
		enormalNhits[(int)i]=projectY2->GetRMS();
	}	



        double Epo[n],nhitspo[n],eEpo[n],enhitspo[n],nhitspo2[n];
        double normalNhitspo[n], normalNhitspo_scale[n],enormalNhitspo[n];

	int counter=0;
        for( double i=0; i<n; i++ ){
		//TH1D* projectX  =  bi210HitsVsEnergy->ProjectionX( "projectX", eAxis->FindBin( i/10 ), eAxis->FindBin((i+1)/10) );	
		//TCanvas* temptemp= new TCanvas();
		TH1D* projectYpo  =  po210HitsVsEnergy->ProjectionY( "projectYpo", eAxispo->FindBin( i/10 ), eAxispo->FindBin((i+1)/10) );	
		projectYpo->Add(po212HitsVsEnergy->ProjectionY( "hold",  eAxispo->FindBin( i/10 ), eAxispo->FindBin((i+1)/10)));
		projectYpo->Add(po214HitsVsEnergy->ProjectionY( "holder",  eAxispo->FindBin( i/10 ), eAxispo->FindBin((i+1)/10)));
		TH1D* projectYpo2  =  po210NHitsVsEnergy->ProjectionY( "projectY2po", eAxispo->FindBin( i/10 ), eAxispo->FindBin((i+1)/10) );	
		projectYpo2->Add(po212NHitsVsEnergy->ProjectionY( "hold11", eAxispo->FindBin( i/10 ), eAxispo->FindBin((i+1)/10)));
		projectYpo2->Add(po214NHitsVsEnergy->ProjectionY( "holder11", eAxispo->FindBin( i/10 ), eAxispo->FindBin((i+1)/10)));
		//E[i]=projectX->GetMean();
		//eE[i]=projectX->GetRMS();
		Epo[(int)i]=(i/10)+0.05;
		eEpo[(int)i]=0;
		temp=projectYpo->GetMean();
		temp2=projectYpo->GetRMS();
		if( temp>0. && projectYpo->GetEntries()>20){
			counter++;
			nhitspo[(int)i]=temp;
			nhitspo2[(int)i]=temp * 1/sqrt(95./86.);
			enhitspo[(int)i]=temp2;
		}

		normalNhitspo[(int)i]=projectYpo2->GetMean();
		normalNhitspo_scale[(int)i]=projectYpo2->GetMean()*95./86.;
		enormalNhitspo[(int)i]=projectYpo2->GetRMS();
		cout<<"this is the "<<i<<"th nhits mean = "<< normalNhitspo[(int)i]<<endl;
		cout<<"this is the "<<i<<"th nhits rms = "<< enormalNhitspo[(int)i]<<endl;
	}	
        
	TMultiGraph *mg = new TMultiGraph();
        TCanvas* error= new TCanvas();
	error->cd();
        TGraphErrors * gr =new TGraphErrors(n,E,nhits,eE,enhits);
	gr->SetMarkerColor(kBlue);
	gr->SetLineColor(kBlue);
	gr->SetFillColor(kBlue);
        mg->Add(gr);

        TGraphErrors * gr1 =new TGraphErrors(n,E,nhitspo,eEpo,enhitspo);
	gr1->SetMarkerColor(kRed);
	gr1->SetLineColor(kRed);
	gr1->SetFillColor(kRed);
        mg->Add(gr1);

        TGraphErrors * gr2 =new TGraphErrors(n,E,nhitspo2,eEpo,enhitspo);
	gr2->SetMarkerColor(kBlack);
	gr2->SetLineColor(kBlack);
	gr2->SetFillColor(kBlack);
        mg->Add(gr2);


	mg->Draw("ap");
        mg->GetXaxis()->SetRangeUser(0.0,3.5);
        mg->GetYaxis()->SetRangeUser(0.0,0.12);
	mg->SetTitle("Resolution across energy");
        mg->GetXaxis()->SetTitle("Energy MC (MeV)");
       	mg->GetYaxis()->SetTitle("1/sqrt(nhits)");



	TLegend * errorLeg = new TLegend(0.7,0.7,0.9,0.9);
	errorLeg->AddEntry(gr,"Bi","f");
	errorLeg->AddEntry(gr1,"Po","f");
	errorLeg->AddEntry(gr2,"Scaled Po","f");
	errorLeg->Draw();
	error->Print("MeanResolution.png");

//=================================================================================================

	TMultiGraph *nhitmg = new TMultiGraph();
        TCanvas* nhitCan= new TCanvas();
 	nhitCan->cd();
        TGraphErrors * nhits_gr_bi =new TGraphErrors(n,E,normalNhits,eE,enormalNhits);
	nhits_gr_bi ->SetMarkerColor(kBlue);
	nhits_gr_bi->SetLineColor(kBlue);
	nhits_gr_bi->SetFillColor(kBlue);
        nhitmg->Add(nhits_gr_bi);

        TGraphErrors * nhits_gr_po =new TGraphErrors(n,E,normalNhitspo,eE,enormalNhitspo);
	nhits_gr_po->SetMarkerColor(kRed);
	nhits_gr_po->SetLineColor(kRed);
	nhits_gr_po->SetFillColor(kRed);
        nhitmg->Add(nhits_gr_po);

        TGraphErrors * nhits_gr_po2 =new TGraphErrors(n,E,normalNhitspo_scale,eE,enormalNhitspo);
	nhits_gr_po2->SetMarkerColor(kBlack);
	nhits_gr_po2->SetLineColor(kBlack);
	nhits_gr_po2->SetFillColor(kBlack);
        nhitmg->Add(nhits_gr_po2);

	nhitmg->Draw("ap");
        //nhitmg->GetXaxis()->SetRangeUser(0.0,3.5);
        //nhitmg->GetYaxis()->SetRangeUser(0.0,0.12);
	nhitmg->SetTitle("Mean nhits across energy");
        nhitmg->GetXaxis()->SetTitle("Energy MC (MeV)");
       	nhitmg->GetYaxis()->SetTitle("Mean nhits");

	TLegend * nhitsLeg = new TLegend(0.2,0.7,0.4,0.9);
	nhitsLeg->AddEntry(nhits_gr_bi,"Bi","f");
	nhitsLeg->AddEntry(nhits_gr_po,"Po","f");
	nhitsLeg->AddEntry(nhits_gr_po2,"Scaled Po","f");
	nhitsLeg->Draw();
	nhitCan->Print("MeanNhits.png");
	
//=================================================================================================

	/*TCanvas* c1 =new TCanvas();
	TH2D* can = new TH2D("can","",100,0,5,100,0,5);
	bi210HitsVsEnergy->SetTitle( "Combined #beta Spectrum " );
	bi210HitsVsEnergy->GetXaxis()->SetTitle( "E_{true}" );
	bi210HitsVsEnergy->GetYaxis()->SetTitle( "1/sqrt(hits)" );
//	bi212HitsVsEnergy->Draw();
	bi210HitsVsEnergy->Draw();
	bi212HitsVsEnergy->Draw("same");
	bi214HitsVsEnergy->Draw("same");
	TLegend * histLeg = new TLegend(0.7,0.7,0.9,0.9);
	histLeg->AddEntry(bi210HitsVsEnergy,"Bi210","f");
	histLeg->AddEntry(bi212HitsVsEnergy,"Bi212","f");
	histLeg->AddEntry(bi214HitsVsEnergy,"Bi214","f");
	histLeg->Draw();
	//c1->Print("212BetaRes.png");
	c1->Print("BetaRes.png");


	TCanvas* c2 =new TCanvas();
	c2->cd();
	po210HitsVsEnergy->SetTitle( "Combined #alpha Spectrum" );
	po210HitsVsEnergy->GetXaxis()->SetTitle( "E_{true}" );
	po210HitsVsEnergy->GetYaxis()->SetTitle( "1/sqrt(hits)" );
	po210HitsVsEnergy->Draw();
	po212HitsVsEnergy->Draw("same");
	po214HitsVsEnergy->Draw("same");
	TLegend * histLegPo = new TLegend(0.7,0.7,0.9,0.9);
	histLegPo->AddEntry(po210HitsVsEnergy,"Po210","f");
	histLegPo->AddEntry(po212HitsVsEnergy,"Po212","f");
	histLegPo->AddEntry(po214HitsVsEnergy,"Po214","f");
	histLegPo->Draw();
	c2->Print("AlphaRes.png");
*/
	TCanvas* c3 =new TCanvas();
	c3->cd();
	bi210HitsVsEnergy->SetMarkerColor(kBlack);bi210HitsVsEnergy->SetFillColor(kBlack);bi210HitsVsEnergy->SetLineColor(kBlack);
	bi212HitsVsEnergy->SetMarkerColor(kBlack);bi212HitsVsEnergy->SetFillColor(kBlack);bi212HitsVsEnergy->SetLineColor(kBlack);
	bi214HitsVsEnergy->SetMarkerColor(kBlack);bi214HitsVsEnergy->SetFillColor(kBlack);bi214HitsVsEnergy->SetLineColor(kBlack);

	allHist->SetTitle( "Alpha Resolution" );
	allHist->GetXaxis()->SetTitle( "E_{true}" );
	allHist->GetYaxis()->SetTitle( "1/sqrt(hits)" );
	allHist->Draw();
	bi210HitsVsEnergy->Draw("same");
	bi212HitsVsEnergy->Draw("same");
	bi214HitsVsEnergy->Draw("same");
	po210HitsVsEnergy->Draw("same");
	po212HitsVsEnergy->Draw("same");
	po214HitsVsEnergy->Draw("same");
	TLegend * histLegPo2 = new TLegend(0.7,0.7,0.9,0.9);
	histLegPo2->AddEntry(po210HitsVsEnergy,"Po210","f");
	histLegPo2->AddEntry(po212HitsVsEnergy,"Po212","f");
	histLegPo2->AddEntry(po214HitsVsEnergy,"Po214","f");
	histLegPo2->AddEntry(bi210HitsVsEnergy,"Combined #beta spectrum","f");
	histLegPo2->Draw();
	c3->Print("Beta-AlphaRes.png");

	TCanvas* c4 =new TCanvas();
	c4->cd();
	bi210HitsVsEnergy->SetMarkerColor(kBlack);bi210HitsVsEnergy->SetFillColor(kBlack);bi210HitsVsEnergy->SetLineColor(kBlack);
	bi212HitsVsEnergy->SetMarkerColor(kBlack);bi212HitsVsEnergy->SetFillColor(kBlack);bi212HitsVsEnergy->SetLineColor(kBlack);
	bi214HitsVsEnergy->SetMarkerColor(kBlack);bi214HitsVsEnergy->SetFillColor(kBlack);bi214HitsVsEnergy->SetLineColor(kBlack);
	allHist->SetTitle( "Modified Alpha Resolution" );
	allHist->GetXaxis()->SetTitle( "E_{true}" );
	allHist->GetYaxis()->SetTitle( "1/sqrt(hits)" );
	allHist->Draw();
	bi210HitsVsEnergy->Draw("same");
	bi212HitsVsEnergy->Draw("same");
	bi214HitsVsEnergy->Draw("same");
	po210HitsVsEnergyMod->Draw("same");
	po212HitsVsEnergyMod->Draw("same");
	po214HitsVsEnergyMod->Draw("same");
	TLegend * histLegPoMod = new TLegend(0.7,0.7,0.9,0.9);
	histLegPoMod->AddEntry(po210HitsVsEnergyMod,"Scaled Po210","f");
	histLegPoMod->AddEntry(po212HitsVsEnergyMod,"Scaled Po212","f");
	histLegPoMod->AddEntry(po214HitsVsEnergyMod,"Scaled Po214","f");
	histLegPoMod->AddEntry(bi210HitsVsEnergy,"Combined #beta spectrum","f");
	histLegPoMod->Draw();
	c4->Print("ModAlphaRes.png");

	TFile fileout1("allHist.root","RECREATE");

	fileout1.cd();
	bi210Hits->Write();
	bi212Hits->Write();
	bi214Hits->Write();
	po210Hits->Write();
	po212Hits->Write();
	po214Hits->Write();
	fileout1.Close();

}


//================================================================================================================


void FillHistrograms(string filename,double E_low,double E_high,TH1D* hNhits,TH1D* hRes,TH1D* hRMS,TH1D* hDelta, TH1D * normal){
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
			if(E_low<energy && energy< E_high){
				hNhits->Fill(nhits);
				hRes->Fill(1/sqrt(nhits));
				hDelta->Fill(energy-energyReco);	
				normal->Fill((energy-energyReco)/energy);
			}
		}
	}
	file->Close();
	
}

void second(){
	//gStyle->SetOptStat(0);

	int n=25;

	vector<string> pureAlpha= glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output/ntuple","alpha");
	vector<string> pureElectron= glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output_electron/ntuple","electron");

	double E[n],nhits_bi[n],nhits_po[n],res_bi[n],res_po[n],RMS_bi[n],RMS_po[n],delta_bi_RMS[n],delta_po_RMS[n],NormE_bi_RMS[n],NormE_po_RMS[n];
	double error_E[n],error_nhits_po[n],error_res_po[n],error_nhits_bi[n],error_res_bi[n],RMS_bi_error[n],RMS_po_error[n],NormE_RMS_bi_error[n],NormE_RMS_po_error[n];

	/* TCanvas* TotalCanvas = new TCanvas(); */
	/* TCanvas* TotalCanvasPo = new TCanvas(); */
	/* cout<<"before TCanvas "<<endl; */
	TCanvas* TotalCanvas_deltaE = new TCanvas();
	TCanvas* TotalCanvasPo_deltaE = new TCanvas();
	/* cout<<"After TCanvas "<<endl; */

	for (double energy=0;energy<25;energy++){
		TH1D* Histrogram_nhits_bi  = new TH1D("Histrogram_nhits_bi","",2000,0,2000);
		Histrogram_nhits_bi->SetLineColor(kBlue);Histrogram_nhits_bi->SetLineWidth(3);

		TH1D* Histrogram_nhits_po  = new TH1D("Histrogram_nhits_po","",2000,0,2000);
		Histrogram_nhits_po->SetLineColor(kRed);Histrogram_nhits_po->SetLineWidth(2);
		
		TH1D* Histrogram_res_bi  = new TH1D("Histrogram_res_bi","",20,0,2.);
		Histrogram_res_bi->SetLineColor(kBlue);Histrogram_res_bi->SetLineWidth(2);

		TH1D* Histrogram_res_po  = new TH1D("Histrogram_res_po","",20,0,2);
		Histrogram_res_po->SetLineColor(kRed);Histrogram_res_po->SetLineWidth(2);
		
		TH1D* Histrogram_RMS_bi  = new TH1D("Histrogram_RMS_bi","",20,0,2.);
		Histrogram_res_bi->SetLineColor(kBlue);Histrogram_res_bi->SetLineWidth(2);

		TH1D* Histrogram_RMS_po  = new TH1D("Histrogram_RMS_po","",20,0,2);
		Histrogram_res_po->SetLineColor(kRed);Histrogram_res_po->SetLineWidth(2);

		TH1D* Histrogram_deltaE_bi  = new TH1D("Histrogram_deltaE_bi","",100,-0.4,0.4);
		Histrogram_deltaE_bi->SetLineColor(kBlue);Histrogram_deltaE_bi->SetLineWidth(2);

		TH1D* Histrogram_deltaE_po  = new TH1D("Histrogram_deltaE_po","",100,-0.4,0.4);
		Histrogram_deltaE_po->SetLineColor(kRed);Histrogram_deltaE_po->SetLineWidth(2);

		TH1D* Histrogram_NormE_bi  = new TH1D("Histrogram_NormE_bi","",100,-0.4,0.4);
		Histrogram_NormE_bi->SetLineColor(kBlue);Histrogram_NormE_bi->SetLineWidth(2);

		TH1D* Histrogram_NormE_po  = new TH1D("Histrogram_NormE_po","",100,-0.4,0.4);
		Histrogram_NormE_po->SetLineColor(kRed);Histrogram_NormE_po->SetLineWidth(2);

		double energy_real= energy/10;

		for( int i=0; i<pureElectron.size(); i++ ){
			partFlag="Bi";
			FillHistrograms(pureElectron[i],energy_real,energy_real+0.1,Histrogram_nhits_bi,Histrogram_res_bi,Histrogram_RMS_bi,Histrogram_deltaE_bi,Histrogram_NormE_po);
		  }
		
		for( int i=0; i<pureAlpha.size(); i++ ){
			partFlag="Po";
			FillHistrograms(pureAlpha[i],energy_real,energy_real+0.1,Histrogram_nhits_po,Histrogram_res_po,Histrogram_RMS_po,Histrogram_deltaE_po,Histrogram_NormE_bi);
		  }

			E[(int)energy]=(energy/10)+0.05;
			error_E[(int)energy]=0;

		/* if(Histrogram_nhits_bi->GetEntries()>200 && Histrogram_nhits_po->GetEntries()>200  ){ */	

		/* 	nhits_bi[(int)energy]=Histrogram_nhits_bi->GetMean(); */
		/* 	error_nhits_bi[(int)energy]=Histrogram_nhits_bi->GetRMS(); */
		/* 	nhits_po[(int)energy]=Histrogram_nhits_po->GetMean(); */
		/* 	error_nhits_po[(int)energy]=Histrogram_nhits_po->GetRMS(); */

		/* 	res_bi[(int)energy]=Histrogram_res_bi->GetMean(); */
		/* 	error_res_bi[(int)energy]=Histrogram_res_bi->GetRMS(); */
		/* 	res_po[(int)energy]=Histrogram_res_po->GetMean(); */
		/* 	error_res_po[(int)energy]=Histrogram_res_po->GetRMS(); */
			

		/* 	TCanvas* nhitsCan= new TCanvas(); */
		/* 	Histrogram_nhits_bi->GetXaxis()->SetTitle("nhits"); */
		/* 	Histrogram_nhits_po->GetXaxis()->SetTitle("nhits"); */
		/* 	Histrogram_nhits_bi->SetTitle(("nhits_bi_"+SSTR(energy_real)+" < E <"+SSTR(energy_real+0.1)+" MeV").c_str()); */
		/* 	Histrogram_nhits_po->SetTitle(("nhits_po_"+SSTR(energy_real)+" < E <"+SSTR(energy_real+0.1)+" MeV").c_str()); */

		/* 	Double_t norm = Histrogram_nhits_bi->GetEntries(); */
		/* 	Histrogram_nhits_bi->Scale(1/norm); */
		/* 	norm = Histrogram_nhits_po->GetEntries(); */
		/* 	Histrogram_nhits_po->Scale(1/norm); */

		/* 	Histrogram_nhits_bi->Fit("gaus"); */
		/* 	Histrogram_nhits_po->Fit("gaus"); */
		/* 	TF1 *fit1 = (TF1*)Histrogram_nhits_bi->GetFunction("gaus"); */
		/* 	TF1 *fit2 = (TF1*)Histrogram_nhits_po->GetFunction("gaus"); */
		/* 	fit1->SetLineColor(kBlack); */
		/* 	fit2->SetLineColor(kBlack); */

		/* 	TCanvas* checkerCanvas = new TCanvas(); */
		/* 	Histrogram_nhits_bi->Draw(); */
		/* 	fit1->Draw("same"); */
		/* 	nhitsCan->Print(("plotChecker/bi/nhits_between_"+SSTR(energy_real)+"to"+SSTR(energy_real+0.1)+".png").c_str()); */
		/* 	Histrogram_nhits_po->Draw(); */
		/* 	fit2->Draw("same"); */
		/* 	nhitsCan->Print(("plotChecker/po/nhits_between_"+SSTR(energy_real)+"to"+SSTR(energy_real+0.1)+".png").c_str()); */

		/* 	TotalCanvas->cd(); */
		/* 	if (energy==0) Histrogram_nhits_bi->Draw(); */
		/* 	Histrogram_nhits_bi->Draw("same"); */

		/* 	TotalCanvasPo->cd(); */
		/* 	if (energy==0) Histrogram_nhits_po->Draw(); */
		/* 	Histrogram_nhits_po->Draw("same"); */
		/* }else{ */
		/* 	nhits_bi[(int)energy]=0; */
		/* 	error_nhits_bi[(int)energy]=0; */
		/* 	nhits_po[(int)energy]=0; */
		/* 	error_nhits_po[(int)energy]=0; */

		/* 	res_bi[(int)energy]=0;; */
		/* 	error_res_bi[(int)energy]=0; */
		/* 	res_po[(int)energy]=0; */
		/* 	error_res_po[(int)energy]=0; */
			
		/* 	RMS_bi[(int)energy]=0; */
		/* 	RMS_po[(int)energy]=0; */
		/* } */

		bool first=true;
		if(Histrogram_nhits_bi->GetEntries()>200 && Histrogram_nhits_po->GetEntries()>200  ){	

			Histrogram_deltaE_bi->GetXaxis()->SetTitle("#Delta E");
			Histrogram_deltaE_po->GetXaxis()->SetTitle("#Delta E");
			Histrogram_deltaE_bi->SetTitle(("deltaE_bi_"+SSTR(energy_real)+" < E <"+SSTR(energy_real+0.1)+" MeV").c_str());
			Histrogram_deltaE_po->SetTitle(("deltaE_po_"+SSTR(energy_real)+" < E <"+SSTR(energy_real+0.1)+" MeV").c_str());

			Histrogram_deltaE_bi->Sumw2();
			Histrogram_deltaE_po->Sumw2();
			Double_t norm = Histrogram_deltaE_bi->GetEntries();
			Histrogram_deltaE_bi->Scale(1/norm);
			norm = Histrogram_deltaE_po->GetEntries();
			Histrogram_deltaE_po->Scale(1/norm);

			Histrogram_deltaE_bi->Fit("gaus");
			Histrogram_deltaE_po->Fit("gaus");

			TF1 *fit1 = (TF1*)Histrogram_deltaE_bi->GetFunction("gaus");
			TF1 *fit2 = (TF1*)Histrogram_deltaE_po->GetFunction("gaus");
			fit1->SetLineColor(kBlack);
			fit2->SetLineColor(kBlack);

			if (false){
				TCanvas* checkerCanvas = new TCanvas();
				Histrogram_deltaE_bi->Draw();
				fit1->Draw("same");
				checkerCanvas->Print(("plotChecker/bi/DeltaE_between_"+SSTR(energy_real)+"0_to"+SSTR(energy_real+0.1)+".png").c_str());

				TCanvas* checkerCanvas_2= new TCanvas();
				Histrogram_deltaE_po->Draw();
				fit2->Draw("same");
				checkerCanvas_2->Print(("plotChecker/po/DeltaE_between_"+SSTR(energy_real)+"0_to"+SSTR(energy_real+0.1)+".png").c_str());
			}

			RMS_bi[(int)energy]=fit1->GetParameter(2);
			RMS_po[(int)energy]=fit2->GetParameter(2);
			RMS_bi_error[(int)energy]=fit1->GetParError(2);
			RMS_po_error[(int)energy]=fit2->GetParError(2);

			if(false){
				//This is to compute the all on one plots.
				TotalCanvas_deltaE->cd();
				if (first=true) Histrogram_deltaE_bi->Draw();
				Histrogram_deltaE_bi->Draw("same");

				TotalCanvasPo_deltaE->cd();
				if (first=true) Histrogram_deltaE_po->Draw();
				Histrogram_deltaE_po->Draw("same");
				first=false;
			}
		
//+++++++++++++++++++++++++NormE++++++++++++++++++++++++++++


			Histrogram_NormE_bi->GetXaxis()->SetTitle("#Delta E");
			Histrogram_NormE_po->GetXaxis()->SetTitle("#Delta E/E_{MC}");
			Histrogram_NormE_bi->SetTitle(("NormE_bi_"+SSTR(energy_real)+" < E <"+SSTR(energy_real+0.1)+" MeV").c_str());
			Histrogram_NormE_po->SetTitle(("NormE_po_"+SSTR(energy_real)+" < E <"+SSTR(energy_real+0.1)+" MeV").c_str());

			Histrogram_NormE_bi->Sumw2();
			Histrogram_NormE_po->Sumw2();
			norm = Histrogram_NormE_bi->GetEntries();
			Histrogram_NormE_bi->Scale(1/norm);
			norm = Histrogram_NormE_po->GetEntries();
			Histrogram_NormE_po->Scale(1/norm);

			Histrogram_NormE_bi->Fit("gaus");
			Histrogram_NormE_po->Fit("gaus");

			fit1 = (TF1*)Histrogram_NormE_bi->GetFunction("gaus");
			fit2 = (TF1*)Histrogram_NormE_po->GetFunction("gaus");
			fit1->SetLineColor(kBlack);
			fit2->SetLineColor(kBlack);

			if (false){
				TCanvas* checkerCanvas = new TCanvas();
				Histrogram_NormE_bi->Draw();
				fit1->Draw("same");
				checkerCanvas->Print(("plotChecker/bi/NormE_between_"+SSTR(energy_real)+"0_to"+SSTR(energy_real+0.1)+".png").c_str());

				TCanvas* checkerCanvas_2= new TCanvas();
				Histrogram_NormE_po->Draw();
				fit2->Draw("same");
				checkerCanvas_2->Print(("plotChecker/po/NormE_between_"+SSTR(energy_real)+"0_to"+SSTR(energy_real+0.1)+".png").c_str());
			}

			NormE_bi_RMS[(int)energy]=fit1->GetParameter(2);
			NormE_po_RMS[(int)energy]=fit2->GetParameter(2);
			NormE_RMS_bi_error[(int)energy]=fit1->GetParError(2);
			NormE_RMS_po_error[(int)energy]=fit2->GetParError(2);

		}else{
			RMS_bi[(int)energy]=0;
			RMS_po[(int)energy]=0;
			RMS_bi_error[(int)energy]=0;
			RMS_po_error[(int)energy]=0;
			
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
        TGraphErrors * gr_RMS_bi =new TGraphErrors(n,E,RMS_bi,error_E,RMS_bi_error);
	gr_RMS_bi ->SetMarkerColor(kBlue);
	gr_RMS_bi->SetLineColor(kBlue);
	gr_RMS_bi->SetFillColor(kBlue);
        mg_RMS->Add(gr_RMS_bi);

        TGraphErrors * gr_RMS_po =new TGraphErrors(n,E,RMS_po,error_E,RMS_po_error);
	gr_RMS_po ->SetMarkerColor(kRed);
	gr_RMS_po->SetLineColor(kRed);
	gr_RMS_po->SetFillColor(kRed);
        mg_RMS->Add(gr_RMS_po);


	mg_RMS->Draw("a*");
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


	mg_NormE_RMS->Draw("a*");
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
}

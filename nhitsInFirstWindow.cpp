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
#include <TAttFill.h>


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
void all(string path,string name);
void FillHist(TFile* file,TH2D * hist);




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

void retriggers(string filename ){
	TH1D * scaleHist= new TH1D("scaleHist","",100,-1000,1000);
	double output [4];
	vector<double> scaleList;
	vector<double> energyList;
	double scale;
	cout<<"filename : "<<filename<<endl;
	RAT::DU::DSReader dsReader( filename );
	Int_t howmany=0;	
	Int_t howmany1=0;	
	for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ ){
		const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );

		size_t count=rDS.GetEVCount();
		if(count>1) howmany++;
		if(count<2) howmany1++;
		
		
	}

cout<<"The number of second triggers : "<<howmany<<endl;
cout<<"Total number of singel triggers : "<<howmany1<<endl;
cout<<"Ratio of second to first : "<<(double) howmany/howmany1<<endl;
}



double* FindHits(string filename ){
	TH1D * scaleHist= new TH1D("scaleHist","",100,-1,1);
	double output [4];
	vector<double> scaleList;
	vector<double> energyList;
	double scale;
	cout<<"filename : "<<filename<<endl;
	RAT::DU::DSReader dsReader( filename );
	
	for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ ){
		int firstWindowHits;	
		int unbuiltHits;	
		int allHits=0;	
		const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
		const RAT::DS::MC& rmc = rDS.GetMC();
		double energy= rmc.GetScintQuenchedEnergyDeposit();

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
		//scale=firstWindowHits;
		//cout<<"ratio between 1st window/all hits : "<< scale<<endl;
		scaleHist->Fill(scale);
		scaleList.push_back(scale);
		energyList.push_back(energy);
		
	}

TCanvas *c1=new TCanvas();
scaleHist->Draw();
c1->Print((filename+".png").c_str());

output[0]=mean(energyList);
output[1]=sqrt(variance(energyList,mean(energyList)));
output[2]=mean(scaleList);
output[3]=sqrt(variance(scaleList,mean(scaleList)));
cout<<"size of scaleList "<< scaleList.size()<<endl;
return output;
}


void window(string path,string name){
	
	TH1D * scaleHist= new TH1D("scaleHist","",100,0,1);
	vector<string> poFileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po212/root/","SolarPo212");
	//vector<string> poFileList= glob(path,name);
	int n=poFileList.size();
	cout<<n<<endl;
	double x[n],y[n],ex[n],ey[n];
	
	for( int i=0; i<poFileList.size(); i++ ){
		cout<<"File List "<<poFileList[i]<<endl;
		double *out=FindHits( poFileList[i]);
		x[i]=out[0];
		ex[i]=out[1];
		y[i]=out[2];
		ey[i]=out[3];
		cout<<"This is the mean scale "<< y[i]<<endl;
		cout<<"This is the std of y "<< ey[i]<<endl;
	  }

	{
	TCanvas* c1= new TCanvas();
	TGraphErrors * gr =new TGraphErrors(n,x,y,ex,ey);
	gr->SetTitle("Number of nhits in the first window.");
	gr->GetXaxis()->SetTitle("Energy (MeV)");
	gr->GetXaxis()->SetRangeUser(0,1.5);
	gr->GetYaxis()->SetRangeUser(0,900);
	gr->GetYaxis()->SetTitle("Nhits in first window");
	gr->Draw("AC*");
	c1->Print("NhitVsEnergy.png");
	c1->Print("NhitVsEnergy.eps");
	}
// Plotting things
//	{
//	TCanvas* c1= new TCanvas();
//	c1->cd();
//	TAxis * xAxis=   scaleHist->GetXaxis();
//	TAxis * yAxis=    scaleHist->GetYaxis();
//
//	scaleHist->SetTitle(name.c_str());
//	xAxis->SetTitle("Fraction in 1st window");
//
//	scaleHist->Draw();
//	c1->Print((name+".png").c_str());
//	c1->Print((name+".eps").c_str());
//	
//	TFile fileout("window.root","UPDATE");
//	scaleHist->Write();
//	fileout.Close();
//	}








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

void FindDiff(string filestring,TH2D * E,TH2D * Erc, double lowr, double highr){
	TFile * file = TFile::Open(filestring.c_str());
	TTree* Tree = (TTree*) file->Get("output");
	Double_t energy,mcEdepQuenched,posr;
	Bool_t fitValid;

	Tree->SetBranchAddress("energy",&energy);
	Tree->SetBranchAddress("mcEdepQuenched",&mcEdepQuenched);
	Tree->SetBranchAddress("posr",&posr);
	Tree->SetBranchAddress("fitValid",&fitValid);
	Int_t n = (Int_t)Tree->GetEntries();

	for( Int_t i =0;i<n;i++){
	Tree->GetEntry(i);
	if(lowr<posr && posr< highr&& fitValid && energy>0.2){
	E->Fill(mcEdepQuenched,energy-mcEdepQuenched);
	Erc->Fill(mcEdepQuenched,(energy-mcEdepQuenched)/mcEdepQuenched);
	//E->Add(Erc,1);
		}
	}

}
void diff(string path,string name){
	
	for( int r =0; r<6;r++){
		TH2D * scaleHist;
		TH2D * scaleHist2;
		if(name=="SolarBi210"){
			scaleHist= new TH2D("scaleHist","",50,0.,1.5,100,-3,3);
			scaleHist2= new TH2D("scaleHist2","",50,0.,1.5,100,-3,3);
		}else if(name=="SolarPo212"){

			scaleHist= new TH2D("scaleHist","",50,0.,2,100,-2,2);
			scaleHist2= new TH2D("scaleHist2","",50,0.,2,100,-1.5,1.5);
		}else if(name=="SolarPo214"){

			scaleHist= new TH2D("scaleHist","",50,0,2,100,-2,2);
			scaleHist2= new TH2D("scaleHist2","",50,0.,2,100,-1.5,1.5);
		}else if(name=="SolarPo210"){

			scaleHist= new TH2D("scaleHist","",50,0.4,0.6,100,-2,2);
			scaleHist2= new TH2D("scaleHist2","",50,0.4,0.6,100,-1.5,1.5);
		}else{
		cout<<"Hist problem"<<endl;
		}

		vector<string> poFileList= glob(path,name);
		if (poFileList.size()>20){

		for( int i=0; i<20/*poFileList.size()*/; i++ ){
			FindDiff( poFileList[i],   scaleHist, scaleHist2,r*1000,(r+1)*1000);
		  }
		}else{

		for( int i=0; i<poFileList.size(); i++ ){
			FindDiff( poFileList[i],   scaleHist, scaleHist2,r*1000,(r+1)*1000);
		  }
		}

	// Plotting things
		{
		TCanvas* c1= new TCanvas();
		c1->cd();
		std::string low = SSTR(r*1000);
		std::string high = SSTR((r+1)*1000);
		TAxis * xAxis=   scaleHist->GetXaxis();
		TAxis * yAxis=    scaleHist->GetYaxis();
		scaleHist->SetTitle((name+" {"+low+"< r < "+high+"}").c_str());
		xAxis->SetTitle("Emc (MeV)");
		yAxis->SetTitle("E-Emc (MeV)");
		scaleHist->Draw();
		c1->Print((name+"_diff_"+low+"mm_"+high+"mm.png").c_str());
		c1->Print((name+"_diff_"+low+"mm_"+high+".eps").c_str());
		
		TFile fileout("diff.root","UPDATE");
		scaleHist->Write();
		fileout.Close();
		}
	
		{
		TCanvas* c2= new TCanvas();
		c2->cd();
		TAxis * xAxis=   scaleHist2->GetXaxis();
		TAxis * yAxis=    scaleHist2->GetYaxis();
		std::string low = SSTR(r*1000);
		std::string high = SSTR((r+1)*1000);
		scaleHist2->SetTitle((name+" {"+low+"< r < "+high+"}").c_str());
		xAxis->SetTitle("Emc (MeV)");
		yAxis->SetTitle("(E-Emc)/Emc");
		scaleHist2->Draw();
		c2->Print((name+"_diff_mc_"+low+"mm_"+high+"mm.png").c_str());
		c2->Print((name+"_diff_mc_"+low+"mm_"+high+".eps").c_str());
		
		TFile fileout("diff_mc.root","UPDATE");
		scaleHist->Write();
		fileout.Close();
		}
		
	}
}

void alldiff(){
	diff("/data/snoplus/OfficialProcessing/production_5_3_0/Bi210","SolarBi210");
	diff("/data/snoplus/OfficialProcessing/production_5_3_0/Po210","SolarPo210");
	diff("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po212","SolarPo212");
	diff("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po214","SolarPo214");
}

void alphas(){
	all("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po210/root","SolarPo210");
	all("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po212/root","SolarPo212");
	all("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po214/root","SolarPo214");
	all("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Be7full","SolarBe7");
	all("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/PEPfull_2","SolarPep");
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


void plotsForWriteUp(){
				{
  TFile* inputFile = new TFile("scales.root");
  inputFile->ls();
  // TH1D* Po210 = dynamic_cast<TH1D*>(inputFile->Get("scaleHist;1"));
  // TH1D* Po212 = dynamic_cast<TH1D*>(inputFile->Get("scaleHist;2"));
  // TH1D* Po214 = dynamic_cast<TH1D*>(inputFile->Get("scaleHist;3"));
  // TH1D* Be7 = dynamic_cast<TH1D*>(inputFile->Get("scaleHist;4"));
  // TH1D* Pep = dynamic_cast<TH1D*>(inputFile->Get("scaleHist;5"));
  TH1D* Po210 =(TH1D*)inputFile->FindObjectAny("scaleHist;1");
  TH1D* Po212 =(TH1D*)inputFile->FindObjectAny("scaleHist;2");
  TH1D* Po214 =(TH1D*)inputFile->FindObjectAny("scaleHist;3");
  TH1D* Be7 =(TH1D*)inputFile->FindObjectAny("scaleHist;4");
  TH1D* Pep =(TH1D*)inputFile->FindObjectAny("scaleHist;5");

  Po210->Scale(1/Po210->GetEntries()); 
  Po212->Scale(1/Po212->GetEntries()); 
  Po214->Scale(1/Po214->GetEntries()); 
  Be7->Scale(1/Be7->GetEntries()); 
  Pep->Scale(1/Pep->GetEntries()); 
	

	TH1D* hblank = new TH1D("hblank","",100,0.,1);
	TH1D* halphas = new TH1D("halphas","",100,0.,1);
	TH1D* hbetas = new TH1D("hbetas","",100,0.,1);
	TCanvas * c1= new TCanvas();
	c1->SetGrid(kTRUE);

	gStyle->SetOptStat(kFALSE);
	hblank->SetMaximum(0.3);
	hblank->GetXaxis()->SetRangeUser(0.6,1);
	hblank->SetTitle("Fraction of nhits in 1^{ st } window");
	hblank->GetXaxis()->SetTitle("Fraction of nhits in 1^{st} window");
	hblank->Draw();
	halphas->Add(Po210,1);
	halphas->Add(Po212,1);
	halphas->Add(Po214,1);
	halphas->SetFillColorAlpha(kBlue,0.1);
	halphas->DrawNormalized("same");
	hbetas->Add(Be7,1);
	hbetas->Add(Pep,1);
	hbetas->SetFillColorAlpha(kRed,0.1);
	hbetas->DrawNormalized("same");

	TLegend* leg= new TLegend(0.1,0.5,0.4,0.9);
	c1->SetGrid(kTRUE);
	leg->AddEntry(halphas,"Po210, Po212 & Po214","f");
	leg->AddEntry(hbetas,"Be7 & pep","f");
	leg->Draw();
	c1->Print("numberOfNhitsInFirstWindow.png");
	c1->Print("numberOfNhitsInFirstWindow.tex");
}
	

	TH2D * hEnergy= new TH2D("hEnergy","",100,0,2.5,100,-1,1);
	vector<string> poFileList= glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output/ntuple/","alpha_");

	for( int i=0; i<poFileList.size(); i++ ){
		TFile * file= TFile::Open(poFileList[i].c_str());	
		FillHist( file, hEnergy);
		file->Close();
	  }
// Plotting things
	{
	TCanvas* c1= new TCanvas();
	c1->cd();
	c1->SetRightMargin(1);
	c1->SetGrid(kTRUE);
	gStyle->SetOptStat(kFALSE);
	hEnergy->SetTitle("#frac{E_{MC}-E_{reco}}{E_{MC}} across alpha energy");	
	hEnergy->GetXaxis()->SetTitle("Alpha energy (MeV)");	
	// TAxis * xAxis=   scaleHist->GetXaxis();
	// TAxis * yAxis=    scaleHist->GetYaxis();
	//McVsZ->Draw("surf1");
	hEnergy->Draw("contz0");
	c1->Print("energyDiffOverEnergy.png");
	c1->Print("energyDiffOverEnergy.tex");
	c1->SaveAs("something.tex");
	c1->SaveAs("something.C");
	// c1->Print("test.eps");
	
	}
}

void FillHist(TFile* file,TH2D * hist){

				TTree* Tree = (TTree*) file->Get("output");
				Double_t para, mcEdepQuenched,energy,posr;
				Bool_t  Qfit;
				Int_t pdg1, pdg2, evIndex;
				Int_t parentpdg1,parentpdg2;

				Tree->SetBranchAddress("pdg1",&pdg1);
				Tree->SetBranchAddress("pdg2",&pdg2);
				Tree->SetBranchAddress("parentpdg1",&parentpdg1);
				Tree->SetBranchAddress("parentpdg2",&parentpdg2);

				Tree->SetBranchAddress("berkeleyAlphaBeta",&para);
				Tree->SetBranchAddress("mcEdepQuenched",&mcEdepQuenched);
				Tree->SetBranchAddress("energy",&energy);
				Tree->SetBranchAddress("posr",&posr);
				Tree->SetBranchAddress("fitValid",&Qfit);
				Tree->SetBranchAddress("evIndex",&evIndex);
				Int_t n = (Int_t)Tree->GetEntries();
				Int_t code;

				for( Int_t i =0;i<n;i++){
								Tree->GetEntry(i);

								if( Qfit && evIndex==0 && posr<6000 && energy>0.1){
												hist->Fill(energy,(mcEdepQuenched-energy)/mcEdepQuenched);
												// hist->Fill(energy,mcEdepQuenched);
								}
								}

				}





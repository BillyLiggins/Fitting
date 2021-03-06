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

void full(){
	/*This function looks at a list of Po files and finds the different between reconstructed and truth energy in radial shells.
	 *
	 *
	 */

	//TH2D * hist= new TH2D("scaleHist","",100,0,1,200,-1,1);
	TH2D * hist= new TH2D("","",100,0.8,1.2,200,-1,1);
	TH1D *reco= new TH1D("","",200,-1,1);
	TH1D *E= new TH1D("","",200,-0,5);
	TH1D *truth= new TH1D("","",200,-1,1);
	//TH1D * scaleHist= new TH1D("scaleHist","",100,0,1);
	//vector<string> poFileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po210/","SolarPo210");
	//vector<string> poFileList= glob("/data/snoplus/OfficialProcessing/production_5_3_0/Po210","SolarPo210");
	vector<string> poFileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po212","SolarPo");
	//vector<string> poFileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po214","SolarPo");

	for( int i=0; i<poFileList.size(); i++ ){
		TFile * file= TFile::Open(poFileList[i].c_str());	
		TTree* Tree = (TTree*) file->Get("output");
		Bool_t fit;
		Double_t Ereco, Etrue,x,y,z,r;
		Double_t xtrue, ytrue, ztrue;

		Tree->SetBranchAddress("energy",&Ereco);
		Tree->SetBranchAddress("mcEdepQuenched",&Etrue);
		Tree->SetBranchAddress("mcPosz",&ztrue);
		Tree->SetBranchAddress("mcPosy",&ytrue);
		Tree->SetBranchAddress("mcPosx",&xtrue);
		Tree->SetBranchAddress("posz",&z);
		Tree->SetBranchAddress("posy",&y);
		Tree->SetBranchAddress("posx",&x);
		Tree->SetBranchAddress("fitValid",&fit);
		Int_t n = (Int_t)Tree->GetEntries();

		for( Int_t i =0;i<n;i++){
		Tree->GetEntry(i);
		r=sqrt(x*x+y*y+z*z);
		if(fit && Ereco>0.2){
			hist->Fill(Etrue,(Etrue-Ereco)/Etrue);
			//three->Fill(Etrue,r,(Etrue-Ereco)/Etrue);
			//hist->Fill(Ereco,(Etrue-Ereco)/Etrue);
			reco->Fill(Etrue-Ereco);
			truth->Fill((Etrue-Ereco)/Etrue);
			E->Fill(Ereco);
		}
	}
	  }
// Plotting things
	{
	TCanvas* can= new TCanvas();
	TAxis * xAxis=	reco->GetXaxis();
	TAxis * yAxis=	reco->GetYaxis();
	xAxis->SetTitle("E_{true}-E_{reco} (MeV)");
	reco->SetTitle("Error in Reconstruction");
	reco->Draw();
	//truth->Draw("same");
	can->Print("plots/complete_DeltaE.png");
	can->Print("plots/complete_DeltaE.eps");

	xAxis=	truth->GetXaxis();
	yAxis=	truth->GetYaxis();
	xAxis->SetTitle("E_{true}-E_{reco})/E_{true}");
	truth->SetTitle("Fractional Error");
	truth->Draw();
	//truth->Draw("same");
	can->Print("plots/complete_FracE.png");
	can->Print("plots/complete_FracE.eps");

	TCanvas* c1= new TCanvas();
	c1->cd();
	xAxis=   hist->GetXaxis();
	yAxis=   hist->GetYaxis();
	yAxis->SetTitle("(E_{true}-E_{reco})/E_{true}");
	xAxis->SetTitle("E_{true} (MeV)");
	hist->SetTitle("(E_{true}-E_{reco})/E_{true} against E_{true} {E_{reco} > 0.2}  ");
	hist->Draw();
	string ss= "plots/complete_test.png";
	string ss2= "plots/complete_test.eps";
	c1->Print(ss.c_str());
	c1->Print(ss2.c_str());
	//three->Draw();
	
	TCanvas * complete = new TCanvas();
	complete->Divide(2,2,0.01,0.01);
	complete->cd(1);
	hist->Draw();
	complete->cd(2);
	truth->Draw();
	complete->cd(3);
	reco->Draw();
	complete->cd(4);
	xAxis=   E->GetXaxis();
	yAxis=   E->GetYaxis();
	xAxis->SetTitle("E_{reco} (MeV)");
	hist->SetTitle("E_{reco}");
	E->Draw();
	complete->Print("plots/complete_threePanel.png");
	}



}

void error(double rad){

	//TH2D * hist= new TH2D("scaleHist","",100,0,1,200,-1,1);
	TH2D * hist= new TH2D("scaleHist","",100,0.4,0.6,200,-1,1);
	TH3D * three= new TH3D("scaleHist","",100,0.4,0.6,100,0,6000,200,-1,1);
	TH1D *reco= new TH1D("reco","",200,-1,1);
	TH1D *truth= new TH1D("truth","",200,-1,1);
	//TH1D * scaleHist= new TH1D("scaleHist","",100,0,1);
	//vector<string> poFileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po210/","SolarPo210");
	vector<string> poFileList= glob("/data/snoplus/OfficialProcessing/production_5_3_0/Po210","SolarPo210");

	for( int i=0; i<poFileList.size(); i++ ){
		TFile * file= TFile::Open(poFileList[i].c_str());	
		TTree* Tree = (TTree*) file->Get("output");
		Bool_t fit;
		Double_t Ereco, Etrue,x,y,z,r;

		Tree->SetBranchAddress("energy",&Ereco);
		Tree->SetBranchAddress("mcEdepQuenched",&Etrue);
		Tree->SetBranchAddress("posz",&z);
		Tree->SetBranchAddress("posy",&y);
		Tree->SetBranchAddress("posx",&x);
		Tree->SetBranchAddress("fitValid",&fit);
		Int_t n = (Int_t)Tree->GetEntries();

		for( Int_t i =0;i<n;i++){
		Tree->GetEntry(i);
		r=sqrt(x*x+y*y+z*z);
		if(fit && Ereco>0.2 && rad<r && r<(rad+1000.)){
			hist->Fill(Etrue,(Etrue-Ereco)/Etrue);
			//three->Fill(Etrue,r,(Etrue-Ereco)/Etrue);
			//hist->Fill(Ereco,(Etrue-Ereco)/Etrue);
			reco->Fill(Etrue-Ereco);
			truth->Fill((Etrue-Ereco)/Etrue);
		}
	}
	  }
// Plotting things
	{
	TCanvas* can= new TCanvas();
	TAxis * xAxis=	reco->GetXaxis();
	TAxis * yAxis=	reco->GetYaxis();
	xAxis->SetTitle("#Delta E");
	reco->Draw();
	//truth->Draw("same");
	can->Print(("DeltaE_"+SSTR(rad)+"m_to_"+SSTR(rad+1000)+"m.png").c_str());
	can->Print(("DeltaE_"+SSTR(rad)+"m_to_"+SSTR(rad+1000)+"m.eps").c_str());

	xAxis=	truth->GetXaxis();
	yAxis=	truth->GetYaxis();
	xAxis->SetTitle("#Delta E/E_{true}");
	truth->Draw();
	//truth->Draw("same");
	can->Print(("FracE_"+SSTR(rad)+"m_to_"+SSTR(rad+1000)+"m.png").c_str());
	can->Print(("FracE_"+SSTR(rad)+"m_to_"+SSTR(rad+1000)+"m.eps").c_str());

	TCanvas* c1= new TCanvas();
	c1->cd();
	xAxis=   hist->GetXaxis();
	yAxis=   hist->GetYaxis();
	yAxis->SetTitle("(E_{true}-E_{reco})/E_{true}");
	xAxis->SetTitle("E_{true}");
	hist->SetTitle(("(E_{true}-E_{reco})/E_{true} against E_{true} {E_{reco} > 0.2 & "+SSTR(rad)+"mm < posr <"+SSTR(rad+1000)+"mm } ").c_str());
	hist->Draw();
	string ss= "test_"+SSTR(rad)+"m_to_"+SSTR(rad+1000)+"m.png";
	string ss2= "test_"+SSTR(rad)+"mm_to_"+SSTR(rad+1000)+"mm.eps";
	c1->Print(ss.c_str());
	c1->Print(ss2.c_str());
	//three->Draw();
	
TCanvas * complete = new TCanvas();
complete->Divide(2,2,0.01,0.01);
complete->cd(1);
hist->Draw();
complete->cd(2);
truth->Draw();
complete->cd(3);
reco->Draw();
complete->Print(("complete"+SSTR(rad)+"m_to_"+SSTR(rad+1000)+"m.png").c_str());
	}



}

void more(){
	
	for( double r =0;r<6000;r=r+1000){
		error(r);
		cout<<"From for loop r = "<<r<<endl;


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

double* FindHits(string filename ){
	TH1D * scaleHist= new TH1D("scaleHist","",100,-1000,1000);
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
		//scale=(double)firstWindowHits/(allHits+unbuiltHits);
		scale=firstWindowHits;
		cout<<"ratio between 1st window/all hits : "<< scale<<endl;
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
	vector<string> poFileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po210/root/","SolarPo210_r56_s0_p2.root");
	/* vector<string> poFileList= glob(path,name); */
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

void viva(){

	window("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output/root/","alpha_");
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
    Double_t mcEdepquenched;
    Double_t alphaBeta212;
    Double_t alphaBeta214;
    Double_t mcEdepQuenched;

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


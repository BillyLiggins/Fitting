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



void compareT(double num,TH1D * time){

	vector<string> biFileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Bi210/root/","SolarBi210");
	TH1D* hist= new TH1D("hist","",100,-100,300);
	hist->SetLineColor(kBlue);


	vector<string> poFileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po210/root/","SolarPo210");
	TH1D* hist2= new TH1D("hist2","",100,-100,300);
	hist2->SetLineColor(kBlack);

	RAT::DU::LightPathCalculator fLightPath = RAT::DU::Utility::Get()->GetLightPathCalculator();
	RAT::DU::DSReader dsreader2(biFileList[0]);
	double time_res;
	double t10=0;
	int entrycount2= dsreader2.GetEntryCount();
	vector<double> vars;
	for( int ent=0; ent</*entrycount2*/1000;ent++){

		const RAT::DS::Entry& ds=dsreader2.GetEntry(ent);

		for(int iEv=0; iEv<ds.GetEVCount(); iEv++){
			const RAT::DS::EV& rev = ds.GetEV(iEv);

			const RAT::DS::FitVertex& vertex= rev.GetFitResult("scintFitter").GetVertex(0);

			TVector3 fEventPosition= vertex.GetPosition();
			double fEventTime= vertex.GetTime();
			const RAT::DS::CalPMTs& calpmts = rev.GetCalPMTs();

			vector<double> timesorted;
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
				//cout<< "Time "<< time_res<< endl;
				timesorted.push_back(time_res);
			}
				sort(timesorted.begin(),timesorted.end());
				t10=timesorted[num+1];
//				cout<< "t10 "<< t10<< endl;
			time->Fill(t10);
//			time->Fit("gaus");
//			TF1 *fit = time->GetFunction("gaus");
//			double var=fit->GetParameter(2);
//			vars.push_back(var);
//			cout<<"vars "<<var<<endl;


		}
	}
//	double min=*min_element(vars.begin(),vars.end());

//	return min;
}

void lop(){
	TCanvas *c1=new TCanvas();
	c1->cd();
	for(double i=0;i<20;i++){
	TH1D * time=new TH1D("time","",100,-10,10);
	time->SetLineColor(i);
	time->SetMaximum(100);
	compareT(i,time);
	time->Draw("same");
	double min = time->GetRMS();
	cout<< "min : "<<min<<" the "<<i<<" entry."<<endl;
	}
}
void findTres(string filename,TH1D* hist2){

	RAT::DU::LightPathCalculator fLightPath = RAT::DU::Utility::Get()->GetLightPathCalculator();
	RAT::DU::DSReader dsreader2( filename);
	double time_res;
	double t10=0;
	int entrycount2= dsreader2.GetEntryCount();
//	TH1D * time=new TH1D("time","",200,-200,200);
	for( int ent=0; ent</*entrycount2*/2;ent++){

		const RAT::DS::Entry& ds=dsreader2.GetEntry(ent);

		for(int iEv=0; iEv<ds.GetEVCount(); iEv++){
			const RAT::DS::EV& rev = ds.GetEV(iEv);

			const RAT::DS::FitVertex& vertex= rev.GetFitResult("scintFitter").GetVertex(0);

			TVector3 fEventPosition= vertex.GetPosition();
			double fEventTime= vertex.GetTime();
			const RAT::DS::CalPMTs& calpmts = rev.GetCalPMTs();

			vector<double> timesorted;
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
				timesorted.push_back(time_res);
			}
				sort(timesorted.begin(),timesorted.end());
				t10=timesorted[10];

			for(int iPMT=0;iPMT<calpmts.GetCount();iPMT++){
				const RAT::DS::PMTCal& pmt = calpmts.GetPMT(iPMT);

				const TVector3 pmtPosition = RAT::DU::Utility::Get()->GetPMTInfo().GetPosition(pmt.GetID());
				const float pmtTime = pmt.GetTime();
				fLightPath.CalcByPosition(fEventPosition, pmtPosition);
				double distInInnerAV = fLightPath.GetDistInInnerAV();
				double distInAV = fLightPath.GetDistInAV();
				double distInWater = fLightPath.GetDistInWater();

				const double flightTime = RAT::DU::Utility::Get()->GetEffectiveVelocity().CalcByDistance(distInInnerAV, distInAV, distInWater);

				//time_res=pmtTime-flightTime-fEventTime-t10;
				time_res=pmtTime-flightTime-fEventTime;
				//cout<< "Time "<< time_res<< endl;
				hist2->Fill(time_res);
			}	

		}
	}
	//time->Draw();

}


//---------------------------Chi test statistic---------------------------
void findChi(string filename, TH1D* chiH, TH1D* exhist){
	double count=0;
	double sb1000=0;
	TH1D * true_energy=new TH1D("energy","",100,0,1);
	true_energy->SetLineColor(kRed);true_energy->SetFillColor(kRed);
	TH1D * raw_energy=new TH1D("energy","",100,0,1);
	raw_energy->SetLineColor(kBlack);
	TH1D * energy=new TH1D("energy","",100,0,1);
	energy->SetLineColor(kGreen);energy->SetFillColor(kGreen);

	RAT::DU::LightPathCalculator fLightPath = RAT::DU::Utility::Get()->GetLightPathCalculator();
	RAT::DU::DSReader dsreader2( filename);
	double time_res;
	double t10=0;
	TH1D* hist2= new TH1D("hist","",100,-100,300);
	hist2->SetLineColor(kRed);
	int entrycount2= dsreader2.GetEntryCount();
	double window=0;

	for( int ent=0; ent</*entrycount2*/1000;ent++){


		const RAT::DS::Entry& ds=dsreader2.GetEntry(ent);
		const RAT::DS::MC& mcev = ds.GetMC();
		double trueEnergy= mcev.GetScintQuenchedEnergyDeposit();

		for(int iEv=0; iEv<ds.GetEVCount(); iEv++){
			if(ds.GetEVCount()>1) window++;
			cout<<"# of EVs in each entry = "<<ds.GetEVCount()<<endl;
			const RAT::DS::EV& rev = ds.GetEV(iEv);
			true_energy->Fill(trueEnergy);
			const RAT::DS::FitVertex& vertex= rev.GetFitResult("scintFitter").GetVertex(0);
			double E=vertex.GetEnergy();
			TVector3 fEventPosition= vertex.GetPosition();
			double fEventTime= vertex.GetTime();
			const RAT::DS::CalPMTs& calpmts = rev.GetCalPMTs();

			vector<double> timesorted;
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
				timesorted.push_back(time_res);
			}
				sort(timesorted.begin(),timesorted.end());
				t10=timesorted[9];
				//cout<< "t10 "<< t10<< endl;

			for(int iPMT=0;iPMT<calpmts.GetCount();iPMT++){
				const RAT::DS::PMTCal& pmt = calpmts.GetPMT(iPMT);

				const TVector3 pmtPosition = RAT::DU::Utility::Get()->GetPMTInfo().GetPosition(pmt.GetID());
				const float pmtTime = pmt.GetTime();
				fLightPath.CalcByPosition(fEventPosition, pmtPosition);
				double distInInnerAV = fLightPath.GetDistInInnerAV();
				double distInAV = fLightPath.GetDistInAV();
				double distInWater = fLightPath.GetDistInWater();

				const double flightTime = RAT::DU::Utility::Get()->GetEffectiveVelocity().CalcByDistance(distInInnerAV, distInAV, distInWater);

				time_res=pmtTime-flightTime-fEventTime-t10;
				//cout<< "Time "<< time_res<< endl;
				hist2->Fill(time_res);
			}	

			double N= hist2->GetEntries();
			double chi=0;
			double ob, ex;
			for (int i=0; i<100;i++){
				ob= hist2->GetBinContent(i+1);	
				ex=N*( exhist->GetBinContent(i+1));	
				if (ex==0) break;
				chi=chi+(ob-ex)*(ob-ex)/ex;
			//cout<<" ob = "<<ob<<endl;
			//cout<<" ex = "<<ex<<endl;
			}
			sb1000++;
			raw_energy->Fill(E);
			if(chi>1500){
				energy->Fill(E);
				count++;
			}
			//cout<<" N = "<<N<<endl;
			//cout<<"chi = "<<chi<<endl;
			chiH->Fill(chi);
		}
	}
cout<<"should be 1000 = "<<sb1000<<endl; 
cout<<"Count is equal to "<<count<<endl; 
cout<<"Fraction of appected events = "<<count/sb1000<<endl; 

cout<<"Window = "<<window<<endl;
true_energy->Draw();
raw_energy->Draw("same");
energy->Draw("same");
}


void go(){

	int i=0;

	vector<string> biFileList= glob("/data/liggins/SummerInternship/BiPo/billSims","1MeV_electrons_labppo_test.root");
	TH1D* hist= new TH1D("hist","",100,-100,300);
	hist->SetLineColor(kRed);hist->SetLineWidth(3);

	findTres(biFileList[0],hist);

	//vector<string> poFileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Po210/root/","SolarPo210");
	//vector<string> poFileList= glob("/data/liggins/","Solar");
	vector<string> poFileList= glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Bipo212full","Solar");


	TH1D* hist2= new TH1D("hist2","",100,-100,300);
	hist2->SetLineColor(kBlue);hist2->SetLineWidth(3);

	findTres(poFileList[3],hist2);

	gStyle->SetOptStat(0);
	TLegend* t1 = new TLegend(0.6,0.7,0.89,0.88);
	TCanvas * c1= new TCanvas();
	TAxis* xaxis=hist->GetXaxis();
	xaxis->SetTitle("T_{res} (ns)");
	hist->SetTitle("Time Residual Comparison Between a 1Mev electrons and a BiPo212 Event");
	hist->DrawNormalized();	
	hist2->DrawNormalized("same");	
	t1->AddEntry(hist, "1 MeV e^{-}");
	t1->AddEntry(hist2, "BiPo212");
	t1->Draw();
	c1->Print("TimeResidual_BiPo.png");
}

//--------------------------------------------------------------

//void findTres(string filename,TH2D* hist2){
//	//you have started to try to plot nhits in the 1st window against t_res
//
//	RAT::DU::LightPathCalculator fLightPath = RAT::DU::Utility::Get()->GetLightPathCalculator();
//	RAT::DU::DSReader dsreader2( filename);
//	double time_res;
//	int entrycount2= dsreader2.GetEntryCount();
//	for( int ent=0; ent</*entrycount*/500;ent++){
//
//		const RAT::DS::Entry& ds=dsreader2.GetEntry(ent);
//
//		for(int iEv=0; iEv<ds.GetEVCount(); iEv++){
//			const RAT::DS::EV& rev = ds.GetEV(iEv);
//
//			const RAT::DS::FitVertex& vertex= rev.GetFitResult("scintFitter").GetVertex(0);
//
//			TVector3 fEventPosition= vertex.GetPosition();
//			double fEventTime= vertex.GetTime();
//			const RAT::DS::CalPMTs& calpmts = rev.GetCalPMTs();
//
//			for(int iPMT=0;iPMT<calpmts.GetCount();iPMT++){
//				const RAT::DS::PMTCal& pmt = calpmts.GetPMT(iPMT);
//
//				const TVector3 pmtPosition = RAT::DU::Utility::Get()->GetPMTInfo().GetPosition(pmt.GetID());
//				const float pmtTime = pmt.GetTime();
//				fLightPath.CalcByPosition(fEventPosition, pmtPosition);
//				double distInInnerAV = fLightPath.GetDistInInnerAV();
//				double distInAV = fLightPath.GetDistInAV();
//				double distInWater = fLightPath.GetDistInWater();
//
//				const double flightTime = RAT::DU::Utility::Get()->GetEffectiveVelocity().CalcByDistance(distInInnerAV, distInAV, distInWater);
//
//				time_res=pmtTime-flightTime-fEventTime;
//				cout<< "Time "<< time_res<< endl;
//				hist2->Fill(time_res);
//			}	
//
//		}
//                const RAT::DS::MCHits& unhits=rDS.GetMC().GetUnbuiltMCHits();
//                unbuiltHits=unhits.GetCount();  
//                scale=(double)firstWindowHits/(allHits+unbuiltHits);
//                cout<<"ratio between 1st window/all hits : "<< scale<<endl;
//		
//		for(int iEv=0; iEv<ds.GetEVCount(); iEv++){
//			const RAT::DS::EV& rev = ds.GetEV(iEv);
//
//			const RAT::DS::FitVertex& vertex= rev.GetFitResult("scintFitter").GetVertex(0);
//
//			TVector3 fEventPosition= vertex.GetPosition();
//			double fEventTime= vertex.GetTime();
//			const RAT::DS::CalPMTs& calpmts = rev.GetCalPMTs();
//
//			vector<double> timesorted;
//			for(int iPMT=0;iPMT<calpmts.GetCount();iPMT++){
//				const RAT::DS::PMTCal& pmt = calpmts.GetPMT(iPMT);
//
//				const TVector3 pmtPosition = RAT::DU::Utility::Get()->GetPMTInfo().GetPosition(pmt.GetID());
//				const float pmtTime = pmt.GetTime();
//				fLightPath.CalcByPosition(fEventPosition, pmtPosition);
//				double distInInnerAV = fLightPath.GetDistInInnerAV();
//				double distInAV = fLightPath.GetDistInAV();
//				double distInWater = fLightPath.GetDistInWater();
//
//				const double flightTime = RAT::DU::Utility::Get()->GetEffectiveVelocity().CalcByDistance(distInInnerAV, distInAV, distInWater);
//
//				time_res=pmtTime-flightTime-fEventTime;
//				//cout<< "Time "<< time_res<< endl;
//				timesorted.push_back(time_res);
//			}
//				sort(timesorted.begin(),timesorted.end());
////				for( int i =0; i< timesorted.size();i++){
////					cout<<"The is the "<<i<<"th after sorting : "<<timesorted[i]<<endl;
////				}
//				t10=timesorted[10];
//				cout<< "t10 "<< t10<< endl;
////				time->Fill(t10);
//
//			for(int iPMT=0;iPMT<calpmts.GetCount();iPMT++){
//				const RAT::DS::PMTCal& pmt = calpmts.GetPMT(iPMT);
//
//				const TVector3 pmtPosition = RAT::DU::Utility::Get()->GetPMTInfo().GetPosition(pmt.GetID());
//				const float pmtTime = pmt.GetTime();
//				fLightPath.CalcByPosition(fEventPosition, pmtPosition);
//				double distInInnerAV = fLightPath.GetDistInInnerAV();
//				double distInAV = fLightPath.GetDistInAV();
//				double distInWater = fLightPath.GetDistInWater();
//
//				const double flightTime = RAT::DU::Utility::Get()->GetEffectiveVelocity().CalcByDistance(distInInnerAV, distInAV, distInWater);
//
//				time_res=pmtTime-flightTime-fEventTime-t10;
//				//cout<< "Time "<< time_res<< endl;
//				hist2->Fill(time_res);
//			}	
//
//		}
//	}
//	//time->Draw();
//
//}


//void findTres(string filename,TH2D* hist2){
//	//you have started to try to plot nhits in the 1st window against t_res
//
//	RAT::DU::LightPathCalculator fLightPath = RAT::DU::Utility::Get()->GetLightPathCalculator();
//	RAT::DU::DSReader dsreader2( filename);
//	double time_res;
//	int entrycount2= dsreader2.GetEntryCount();
//	for( int ent=0; ent</*entrycount*/500;ent++){
//
//		const RAT::DS::Entry& ds=dsreader2.GetEntry(ent);
//
//		for(int iEv=0; iEv<ds.GetEVCount(); iEv++){
//			const RAT::DS::EV& rev = ds.GetEV(iEv);
//
//			const RAT::DS::FitVertex& vertex= rev.GetFitResult("scintFitter").GetVertex(0);
//
//			TVector3 fEventPosition= vertex.GetPosition();
//			double fEventTime= vertex.GetTime();
//			const RAT::DS::CalPMTs& calpmts = rev.GetCalPMTs();
//
//			for(int iPMT=0;iPMT<calpmts.GetCount();iPMT++){
//				const RAT::DS::PMTCal& pmt = calpmts.GetPMT(iPMT);
//
//				const TVector3 pmtPosition = RAT::DU::Utility::Get()->GetPMTInfo().GetPosition(pmt.GetID());
//				const float pmtTime = pmt.GetTime();
//				fLightPath.CalcByPosition(fEventPosition, pmtPosition);
//				double distInInnerAV = fLightPath.GetDistInInnerAV();
//				double distInAV = fLightPath.GetDistInAV();
//				double distInWater = fLightPath.GetDistInWater();
//
//				const double flightTime = RAT::DU::Utility::Get()->GetEffectiveVelocity().CalcByDistance(distInInnerAV, distInAV, distInWater);
//
//				time_res=pmtTime-flightTime-fEventTime;
//				cout<< "Time "<< time_res<< endl;
//				hist2->Fill(time_res);
//			}	
//
//		}
//			int firstWindowHits;    
//			int unbuiltHits;        
//			int allHits=0;  
//			const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
//
////			size_t count=rDS.GetMCEVCount();
////			for( int iMCEV=0; iMCEV<count; iMCEV++){
////				const RAT::DS::MCHits& hits=rDS.GetMCEV(iMCEV).GetMCHits();
////				cout <<"count : "<<iMCEV<< " hits : "<< hits.GetCount() <<endl; 
////				if (iMCEV==0){
////					firstWindowHits=hits.GetCount();
////				}
////			allHits+=hits.GetCount();
////			}
//			
//			const RAT::DS::MCHits& unhits=rDS.GetMC().GetUnbuiltMCHits();
//			unbuiltHits=unhits.GetCount();  
//			scale=(double)firstWindowHits/(allHits+unbuiltHits);
//			cout<<"ratio between 1st window/all hits : "<< scale<<endl;
//			//scaleList.push_back(scale);
//			scaleHist->Fill(scale);
//
//
//
//
//
//	}
//
//
//        vector<double> scaleList;
//        double scale;
//        cout<<"filename : "<<filename<<endl;
//        RAT::DU::DSReader dsReader( filename );
//        // Loop through entrys in rootfile
//        for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ ){
//                int firstWindowHits;    
//                int unbuiltHits;        
//                int allHits=0;  
//                const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
//
//                size_t count=rDS.GetMCEVCount();
//                for( int iMCEV=0; iMCEV<count; iMCEV++){
//                        const RAT::DS::MCHits& hits=rDS.GetMCEV(iMCEV).GetMCHits();
//                        cout <<"count : "<<iMCEV<< " hits : "<< hits.GetCount() <<endl; 
//                        if (iMCEV==0){
//                                firstWindowHits=hits.GetCount();
//                        }
//                allHits+=hits.GetCount();
//                }
//                
//                const RAT::DS::MCHits& unhits=rDS.GetMC().GetUnbuiltMCHits();
//                unbuiltHits=unhits.GetCount();  
//                scale=(double)firstWindowHits/(allHits+unbuiltHits);
//                cout<<"ratio between 1st window/all hits : "<< scale<<endl;
//                //scaleList.push_back(scale);
//                scaleHist->Fill(scale);
//        }
//
//
//}
//

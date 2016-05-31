#ifdef __BILL_UTIL__
#define __BILL_UTIL__
#define SSTR( x ) static_cast< std::ostringstream & >( \ ( std::ostringstream() << std::dec << x ) ).str()

#include <vector> 
#include <string> 
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

/* namespace util { */

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

/* }//End of namespace */
#endif

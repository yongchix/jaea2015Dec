/** \file Trace.cpp
 * \brief Implement how to do our usual tricks with traces
 */
#include <algorithm>
#include <iostream>
#include <cmath>
#include <numeric>

#include "Trace.hpp"

using namespace std;
using namespace dammIds::trace;

namespace dammIds {
    namespace trace {
    }
}

const Trace emptyTrace; ///< an empty trace for const references to point to
Plots Trace::histo(OFFSET, RANGE, "traces");

double Trace::GetTraceMax(unsigned int lo,unsigned int hi){
	double max=0,maxch=0;
	double baseline = DoBaseline(0,20);
	double max2=0;
	for(unsigned int i=0;i<hi;i++){
		if(max<at(i)){
			max=at(i);
			maxch=i;
		}
	}
  
	max2=max-baseline;
	return max2;
}


void Trace::TrapezoidalFilter(Trace &filter, const TFP &parms,
							  unsigned int lo, unsigned int hi) const {
    lo = max(lo, (unsigned int)parms.GetSize());

    filter.assign(lo, 0);

    for (unsigned int i = lo; i < hi; i++) {
        int leftSum = accumulate(begin() + i - parms.GetSize(),
                                 begin() + i - parms.GetRiseSamples()
                                 - parms.GetGapSamples(), 0);
        int rightSum = accumulate(begin() + i - parms.GetRiseSamples(),
                                  begin() + i, 0);
        filter.push_back(rightSum - leftSum);
    }
}

void Trace::ShowTraceVal(unsigned int lo, unsigned int hi){
	cout << "===================" << endl;
	for(unsigned int i=lo;i<hi;i++){
		cout << "[TraceVal in Trace.cpp]" << i << " " << at(i) << endl;
	}
	cout << "===================" << endl;
}

bool Trace::ScanPileup(unsigned int lo,unsigned int hi){
  
	double baseline = DoBaseline(2,20);
	double max1=0;
	int ch_max1=0;

	double max2=0;
	int ch_max2=0;

	int space = 15;
	double reduction=0.4;
	for(unsigned int i=lo;i<hi;i++){
		if(max1 < at(i)){
			max1 = at(i);
			ch_max1 = i;
		}
	}

	for(unsigned int i=ch_max1+space;i<hi;i++){
		if(max2 < at(i)){
			max2 = at(i);
			ch_max2 = i;
		}
	}
  
	//cout << max1*reduction << " " << max2 << endl; 
	if((max1-baseline)*reduction<(max2-baseline)){
		//  cout << "[pile up]" << ch_max1 << "  " << max1-baseline << " " << ch_max2 << " " << max2-baseline << endl; 
		return true ;
	}else{return false;
	}
  
 
}

double Trace::GetDiffMean(unsigned int lo,unsigned int hi){
  
	double t=lo;
	double integral=0;
	double mean=0;
	double sum_mean=0;				
	double ch_max=0;
	double max=0;
	double diff_mean;
	for(unsigned int i=lo;i<hi;i++){
    
		//    cout << "GetPileup " <<lo << " " << at(i) << endl;
    
		if(max<at(i)){
			max = at(i);
			ch_max = t;
		}
    
		integral  += at(i);
		sum_mean  += t*at(i); 
		t = t+1;
		//    cout << t << " " << at(i) << " " << integral << " " << sum_mean << " " << ch_max  << endl;
	}
  
	mean      = sum_mean/integral;
	diff_mean = mean-ch_max; 
	// cout << " mean " << mean << " max " << ch_max << " " << diff_mean << endl;
	return diff_mean;
}

double Trace::DoBaseline(unsigned int lo, unsigned int numBins) {
  
	if (size() < lo + numBins) {
		cerr << "Bad range in baseline calculation." << endl;
		return NAN;
	}
  
    unsigned int hi = lo + numBins;

    if (baselineLow == lo && baselineHigh == hi)
        return GetValue("baseline");

    double sum = accumulate(begin() + lo, begin() + hi, 0.0);
    double mean = sum / numBins;
    double sq_sum = inner_product(begin() + lo, begin() + hi,
                                  begin() + lo, 0.0);
    double std_dev = sqrt(sq_sum / numBins - mean * mean);

    SetValue("baseline", mean);
    SetValue("sigmaBaseline", std_dev);

    baselineLow  = lo;
    baselineHigh = hi;

    return(mean);
}

double Trace::DoDiscrimination(unsigned int lo, unsigned int numBins) {
    unsigned int high = lo+numBins;

    if(size() < high)
        return pixie::U_DELIMITER;

    int discrim = 0, max = GetValue("maxpos");
    double baseline = GetValue("baseline");

    if(size() < high)
		return pixie::U_DELIMITER;

    for(unsigned int i = max+lo; i <= max+high; i++)
		discrim += at(i)-baseline;

    InsertValue("discrim", discrim);

    return(discrim);
}

double Trace::DeduceRegression(unsigned int lo,unsigned int hi,unsigned int mode){
  
	unsigned int max_ch=0;			       
	double max=0;
	double ave[300];
	double waveform[300];
	double t[300];
	double regression[300];
	double reg;
	double reglow=-4000;
	double reghigh=-10;
	int    regcount=0;
	double regsum=0;
	int    range=10;
	double av=0,mean=0;
	double avsum=0,meansum=0;
	double av2=0,mean2=0;
	double avsum2=0,meansum2=0;
  
	double baseline = DoBaseline(2,20);
  
	for(unsigned int i=lo;i<hi;i++){
		waveform[i]=at(i)-baseline;
	}
  
	for(unsigned int i=lo+2;i<hi-2;i++){
		ave[i] = log((waveform[i-2]+waveform[i-1]+waveform[i]+waveform[i+1]+waveform[i+2])/5);
		if(max<ave[i]){
			max_ch=i;
			max=ave[i];
			t[i]=i;
		}
	}
  
  
	if(mode==0){ // fast decay component 
		for(unsigned int i = max_ch+2;i < max_ch+range;i++){  
			regression[i]=(ave[i]-ave[i-1])/(t[i]-t[i-1]);
			if(regression[i]>reglow && regression[i]<reghigh){
				regcount++;
				regsum += abs(regression[i]);
			}
			avsum+=ave[i];
			meansum+=i;
		}  
		av=avsum/(range-2);
		mean=meansum/(range-2);
	}
  
	double crosssum=0;
	double diffsum=0;
	double beta_chi=0;
	double tt;
  
	if(mode==0){
		for(unsigned int i = max_ch+2;i < max_ch+3*range;i++){  
			tt=i;
			crosssum+=(tt-mean)*(ave[i]-av);
			diffsum+=(tt-mean)*(tt-mean);
    
		}
		beta_chi = crosssum/diffsum;
		reg=beta_chi;
	}
  
	if(mode==1){ //  slow decay component 
		for(unsigned int i = max_ch+range;i < max_ch+2*range;i++){  
			regression[i]=(ave[i]-ave[i-1])/(t[i]-t[i-1]);
			if(regression[i]>reglow && regression[i]<reghigh){
				regcount++;
				regsum += abs(regression[i]);
			}
      
			avsum2+=ave[i];
			meansum2+=i;
		}
    
		av2=avsum2/(hi-2-range);
		mean2=meansum2/(hi-2-range);
	}
  

	double crosssum2=0;
	double diffsum2=0;
	double beta_chi2=0;
	double tt2;
  
	if(mode==1){
		for(unsigned int i = max_ch+range;i < max_ch+hi-2;i++){  
			tt2=i;
			crosssum2+=(tt2-mean2)*(ave[i]-av2);
			diffsum2+=(tt2-mean2)*(tt2-mean2);
    
		}
		beta_chi2 = crosssum2/diffsum2;
		reg = beta_chi2;
	}

  
	//  cout << "beta_chi " << beta_chi << " " << beta_chi2 << endl;



	//  cout << reg << endl;
    
	return reg;
}

double Trace::DoQDCSimple(unsigned int lo,unsigned int hi){
	double baseline = DoBaseline(0, 20); 
	
	//	double baseline = GetValue("baseline");
	double qdc=0;
	for(unsigned int i=lo;i<=hi;i++){
		qdc += at(i)-baseline;
		if(lo+hi >= size()) break;
	}
  
	return (qdc);
}
double Trace::DoQDC(unsigned int lo, unsigned int numBins) {
    unsigned int high = lo+numBins;

    if(size() < high)
		return pixie::U_DELIMITER;

    double baseline = GetValue("baseline");
    double qdc = 0, fullQdc = 0;
    
    for(unsigned int i = lo; i <= high; i++) {
        qdc += at(i)-baseline;
		waveform.push_back(at(i)-baseline);
    }
    
    for(unsigned int i = 0; i < size(); i++) {
		fullQdc += at(i)-baseline;
	}// QDC through out the length of a trace

    InsertValue("fullQdc", fullQdc);
    InsertValue("tqdc", qdc);
    return(qdc);
}

//-----------//
double Trace::DoQDCTail(unsigned int lo, unsigned int numBins) {
	/* This is not appropriate for 
	 * pile-up traces
	 */
	unsigned int high = lo + numBins; 
	if(size() < high)
		return pixie::U_DELIMITER; 
	
	double baseline = GetValue("baseline"); 
	double max = 0; 
	int ch_max = 0;
	
	for(unsigned int i = lo; i < high; i++) {
		if(max < at(i)) {
			max = at(i); 
			ch_max = i; 
		}
	} // found the maximum point
	
	double qdc = 0; 
	for(unsigned int i = ch_max; i <= high; i++) {
		qdc += at(i) - baseline; 
	}

	return qdc; 	
}

//---------------//
double Trace::DoPSD(unsigned int lo, unsigned int numBins) {
	unsigned int high = lo + numBins;
	if(size() < high) 
		return pixie::U_DELIMITER; 

	double baseline = DoBaseline(0, 20); 
	double pulseQDC = DoQDCSimple(lo, numBins); 
	double tailQDC = DoQDCTail(lo, numBins); 

	double psd = tailQDC/pulseQDC; 

	return psd; 
}
//-----------------//
void Trace::SumTrace(Trace& tr) {
	if(tr.size() == 0) {
		for(int i = 0; i < size(); i++) {
			tr.push_back(at(i)); 
		}
	} else if(tr.size() == size()) {
		for(int i = 0; i < size(); i++) {
			tr.at(i) += at(i); 
		}
	} else if(tr.size() > 0 && tr.size() != size()) {
		exit(EXIT_FAILURE); 
	}
}

unsigned int Trace::FindMaxInfo(unsigned int lo, unsigned int hi, unsigned int numBins) {
  
	unsigned int high = Globals::get()->traceDelay() /
        (Globals::get()->adcClockInSeconds()*1e9);
    unsigned int low = high - (Globals::get()->trapezoidalWalk() /
							   (Globals::get()->adcClockInSeconds()*1e9)) - 3;

    if(size() < high)
		return pixie::U_DELIMITER;

    Trace::const_iterator itTrace = max_element(begin()+low, end()-(size()-high));
    int maxPos = int(itTrace-begin());
    
    if(maxPos + hi > size())
		return pixie::U_DELIMITER;
    
    if(*itTrace >= 4095) {
		InsertValue("saturation", 1);
		return(-1);
    }
    
    DoBaseline(0, maxPos-lo);
    InsertValue("maxpos", maxPos);
    InsertValue("maxval", *itTrace-GetValue("baseline"));
    
    return (itTrace-begin());
}

void Trace::Plot(int id) {
    for (size_type i=0; i < size(); i++) {
        histo.Plot(id, i, 1, at(i));
    }
}

void Trace::Plot(int id, int row) {
    for (size_type i=0; i < size(); i++) {
        histo.Plot(id, i, row, at(i));
    }
}

void Trace::ScalePlot(int id, double scale) {
    for (size_type i=0; i < size(); i++) {
        histo.Plot(id, i, 1, abs(at(i)) / scale);
    }
}

void Trace::ScalePlot(int id, int row, double scale) {
    for (size_type i=0; i < size(); i++) {
        histo.Plot(id, i, row, abs(at(i)) / scale);
    }
}

void Trace::OffsetPlot(int id, double offset) {
    for (size_type i=0; i < size(); i++) {
        histo.Plot(id, i, 1, max(0., at(i) - offset));
    }
}

void Trace::OffsetPlot(int id, int row, double offset) {
    for (size_type i=0; i < size(); i++) {
        histo.Plot(id, i, row, max(0., at(i) - offset));
    }
}

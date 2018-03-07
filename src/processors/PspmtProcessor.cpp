#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <signal.h>
#include <limits.h>
#include <cmath>

#include "PspmtProcessor.hpp"
#include "DammPlotIds.hpp"
#include "Globals.hpp"
#include "Messenger.hpp"
#include "JAEACorrelator.hpp"
#include "DetectorDriver.hpp"
#include "NaIProcessor.hpp"
using namespace std;
using namespace dammIds::pspmt;

/* by Yongchi Xiao, 01/17/2018
 * Use the external variable pixelCalib as a reference for 
 * calibrating (gainmatching) energy in each pixel
 */

namespace dammIds{
    namespace pspmt{   // the offset is 1900   
		const int D_RAW1=0;
		const int D_RAW2=1;
		const int D_RAW3=2;
		const int D_RAW4=3;
		const int D_RAWD=4;
		const int D_SUM=5;
		const int DD_POS1_RAW=6;
		const int DD_POS2_RAW=7;
		const int DD_POS1=8; // 1908
		const int DD_POS2=9;
      
		const int D_ENERGY_TRACE1=10;
		const int D_ENERGY_TRACE2=11;
		const int D_ENERGY_TRACE3=12;
		const int D_ENERGY_TRACE4=13;
		const int D_ENERGY_TRACED=14;
		const int D_ENERGY_TRACESUM=15;
		const int DD_POS1_RAW_TRACE=16;
		const int DD_POS2_RAW_TRACE=17;
		const int DD_POS1_TRACE=18;
		const int DD_POS2_TRACE=19;
      
		const int D_QDC_TRACE1 = 20;
		const int D_QDC_TRACE2 = 21;
		const int D_QDC_TRACE3 = 22;
		const int D_QDC_TRACE4 = 23;
		const int D_QDC_TRACED = 24;
		const int D_ENERGY_QDCSUM = 25;
		const int DD_POS1_RAW_QDC = 26;
		const int DD_POS1_QDC     = 27;
		const int DD_DIRCAL1      = 28;

		const int DD_P1D_CHANNEL=30;
		const int DD_P1D_TRACE=31;
		const int DD_P1D_QDC=32;
		const int DD_P1D_QDCSUM=33;
      
		const int DD_P1D_IMPLANT_CHE=35;
		const int DD_P1D_DECAY_CHE=36; // 1936
		const int DD_P1D_IMPLANT_QDC=37;
		const int DD_P1D_DECAY_QDC=38;

		const int DD_MAP_IMPLANT=39;
		const int DD_MAP_DECAY=40;
        
		const int DD_MAP_IMPLANT_CHE=41;
		const int DD_MAP_DECAY_CHE=42;

		const int DD_TRACEMAX_P1D=43;
		const int DD_TRACEMAXCAL_P1D=44;
		const int DD_TRMAX_QDC=45;
      
		// decay time granx plot, 50 to 60 occupied
      
		const int DD_ENERGY_DECAY_TIME_GRANX=50;
		const int DD_ENERGY_DECAY_TIME_GRANX_TRACEE=60;

		const int DD_P1D_DECAY_TIME0=70;
		const int DD_P1D_DECAY_TIME1=71;
		const int DD_P1D_DECAY_TIME2=72;
		const int DD_P1D_DECAY_TIME3=73;
		const int DD_P1D_DECAY_TIME4=74;
		const int DD_P1D_DECAY_TIME5=75;
		const int DD_TRACE_E_QDC=76;
		const int DD_SINGLE_TRACE=77;
		const int DD_DOUBLE_TRACE=78;

		const int DD_TRACE_POS=79;
		const int DD_TRACE_DYNODE=80;
      
		const int DD_TRACE_IMPLANT_ALL=81;
		const int DD_TRACE_IMPLANT_DYNODE=82;
		const int DD_TRACE_IMPLANT_POSITION=83;
		const int DD_TRACE_IMPLANT_PILEUP=84;
     
		const int DD_TRACE_DECAY_ALL=85;
		const int DD_TRACE_DECAY_DYNODE=86;
		const int DD_TRACE_DECAY_POSITION=87;
		const int DD_TRACE_DECAY_PILEUP=88;
            
		const int DD_MWPC=89;
		const int DD_MWPC_PSPMT=90;
		const int DD_MWPC_NOPSPMT=91;
		const int DD_TRACE_SLOW=95;

		const int DD_CHE_REG=96;
		const int DD_QDC_REG=97;
		const int DD_QDC_REG_IMPLANT=99;
		const int DD_QDC_REG_DECAY=100;
		const int DD_QDC_REG2_IMPLANT=101;
		const int DD_QDC_REG2_DECAY=102;
 
      

		const int DD_REG12=103;
		const int DD_REG12_IMPLANT=104;
		const int DD_REG12_DECAY=105;
      
		const int DD_QDC_REG1=106;
		const int DD_QDC_REG2=107;

		const int DD_TRACE_HIGHENE=108;
		const int DD_TRACE_HIGHENE2=109;
		const int DD_TRACE_HIGHENE3=110;
		const int DD_ENE_QDC_DECAY=112;
		const int DD_ENE_QDC_IMPLANT=113;
		const int DD_TRACE_HIGHENE4=114;

		const int DD_P1D_QDCCAL=115;
		const int DD_QDCREG_SPECIFIC=116;

    }
}

void PspmtProcessor::PspmtData::Clear(void) {    
}

PspmtEvent::PspmtEvent(){
	ene_pspmt_  = -1.0;
	ene_trace_  = -1.0;
	time_pspmt_ = -1.0;
	time_mwpc_  = -1.0;
	is_beam_    = false;
	is_decay_   = false;
	is_veto_    = false;
}

PspmtEvent::PspmtEvent(double ene_pspmt,double ene_trace,double time_pspmt,double time_mwpc,bool is_beam,bool is_decay,bool is_veto){
	set_ene_pspmt(ene_pspmt);
	set_ene_trace(ene_trace);
	set_time_pspmt(time_pspmt);
	set_time_mwpc(time_mwpc);
	set_is_beam(is_beam);
	set_is_decay(is_decay);
	set_is_veto(is_veto);
}

PspmtProcessor::PspmtProcessor(void) : 
	EventProcessor(OFFSET, RANGE, "pspmt")
{
	associatedTypes.insert("pspmt");
	int size_x=24;
	int size_y=24;
	pixels_ = new deque<PspmtEvent>*[size_x];
	for(int i=0;i<size_x;i++){
		pixels_[i] = new deque<PspmtEvent>[size_y];
	}
}

PspmtProcessor::~PspmtProcessor(){
  
	int size_y=24;
	for(int i=0;i<size_y;i++){
		delete [] pixels_[i];
	}
	delete[] pixels_;
}

void PspmtProcessor::DeclarePlots(void) {
   
	const int posBins      = 32; 
	const int traceBins    = 256;
	const int timeBins     = 256;
	const int traceBins2   = 512;
	const int mapBins      = 1024;
	const int p1dBins      = 1024;
	const int Bins         = 4096;
	const int energyBins   = 2048;
	const int regBins      = 2048;
	const int corrBins     = 2048;
  
	// Raw 0-29
    DeclareHistogram1D(D_RAW1, energyBins, "Pspmt1 Raw");
    DeclareHistogram1D(D_RAW2, energyBins, "Pspmt2 Raw");
    DeclareHistogram1D(D_RAW3, energyBins, "Pspmt3 Raw");
    DeclareHistogram1D(D_RAW4, energyBins, "Pspmt4 Raw");
    DeclareHistogram1D(D_RAWD, energyBins, "Pspmt Dynode");
	/* commented out by YX */
	/*
    DeclareHistogram1D(D_SUM,  energyBins, "Pspmt Sum");
    DeclareHistogram2D(DD_POS1_RAW, mapBins, mapBins, "Pspmt Pos1 Raw"); // 1906
    DeclareHistogram2D(DD_POS2_RAW, mapBins, mapBins, "Pspmt Pos2 Raw"); // 1907
    DeclareHistogram2D(DD_POS1, posBins, posBins, "Pspmt Pos1"); // 1908
    DeclareHistogram2D(DD_POS2, posBins, posBins, "Pspmt Pos2"); // 1909
	*/
	DeclareHistogram2D(6, 1024, 1024, "2-D map of all signals, reproduced - I"); // 1906
	DeclareHistogram2D(7, 1024, 1024, "2-D map of all signals, reproduced - II"); // 1907
	DeclareHistogram2D(8, 32, 32, "Pixelized  2-D map, reproduced - I"); // 1908
	DeclareHistogram2D(9, 32, 32, "Pixelized 2-D map, reproduced - II"); // 1909
	
    DeclareHistogram1D(D_ENERGY_TRACE1, energyBins, "Energy1 from trace");
    DeclareHistogram1D(D_ENERGY_TRACE2, energyBins, "Energy2 from trace");
    DeclareHistogram1D(D_ENERGY_TRACE3, energyBins, "Energy3 from trace");
    DeclareHistogram1D(D_ENERGY_TRACE4, energyBins, "Energy4 from trace");
    DeclareHistogram1D(D_ENERGY_TRACED, energyBins, "EnergyD from trace");
    DeclareHistogram1D(D_ENERGY_TRACESUM,  energyBins, "Pspmt Sum");
	// pspmt pos from trace filtered energy, 1916-1919
	/*
    DeclareHistogram2D(DD_POS1_RAW_TRACE, mapBins, mapBins, "Pspmt pos Raw by Trace1");
    DeclareHistogram2D(DD_POS2_RAW_TRACE, mapBins, mapBins, "Pspmt pos Raw by Trace2");
    DeclareHistogram2D(DD_POS1_TRACE, posBins, posBins, "Pspmt pos by Trace1");
    DeclareHistogram2D(DD_POS2_TRACE, posBins, posBins, "Pspmt pos by Trace2");
	*/
	// by YX
	DeclareHistogram2D(16, 1024, 1024, "MAP by QDC"); // 1916
	DeclareHistogram2D(17, 1024, 1024, "MAP by QDC2"); // 1917
	DeclareHistogram2D(18, 32, 32, "MAP by QDC, truncated"); // 1918
	DeclareHistogram2D(19, 32, 32, "MAP by QDC2, truncated"); // 1919

    DeclareHistogram1D(D_QDC_TRACE1, energyBins, "Energy1 from QDC scaled by 10");
    DeclareHistogram1D(D_QDC_TRACE2, energyBins, "Energy2 from QDC scaled by 10");
    DeclareHistogram1D(D_QDC_TRACE3, energyBins, "Energy3 from QDC scaled by 10");
    DeclareHistogram1D(D_QDC_TRACE4, energyBins, "Energy4 from QDC scaled by 10");
    DeclareHistogram1D(D_QDC_TRACED, energyBins, "EnergyD from QDC scaled by 10");
    DeclareHistogram2D(DD_POS1_RAW_QDC, mapBins, mapBins, "Pspmt pos Raw by QDC"); // 1926
    DeclareHistogram2D(DD_POS1_QDC, posBins, posBins, "Pspmt pos by QDC"); // 1927
    DeclareHistogram2D(DD_DIRCAL1, mapBins, mapBins, "Map direction calibrated");
   
    
    // Energy resolutions 30-
    DeclareHistogram2D(DD_P1D_CHANNEL, energyBins,Bins, "Ch vs Dynode EChannel");
    DeclareHistogram2D(DD_P1D_TRACE, energyBins,Bins, "Ch vs Dynode ETrace"); // 1931
    DeclareHistogram2D(DD_P1D_QDC, energyBins,1024, "Ch vs Dynode QDC scaled by 10"); // 1932
    DeclareHistogram2D(DD_P1D_QDCSUM, energyBins,1024, "P1D vs QDCSum 4 keV/ch, Decays"); // 1933
	DeclareHistogram2D(34, energyBins, 1024, "P1D vs. QDCSum 4 keV/ch, Implants"); // 1934
    // Energy res by correlation
    DeclareHistogram2D(DD_P1D_IMPLANT_CHE, energyBins,p1dBins, "[Implant] ch vs E(ch)"); // 1935
    DeclareHistogram2D(DD_P1D_DECAY_CHE, energyBins,p1dBins,   "[Decay] ch vs E(ch)"); // 1936
    DeclareHistogram2D(DD_P1D_IMPLANT_QDC, energyBins,p1dBins, "[Implant] ch vs E(QDC)"); // 1937
    DeclareHistogram2D(DD_P1D_DECAY_QDC, energyBins,p1dBins, "[Decay] ch vs E(QDC)"); // 1938

	// by Yongchi
	DeclareHistogram2D(40, 32, 32, "Implant Map no veto"); // 1940
	DeclareHistogram1D(41, 1024, "10*log(Dt/1.e-9), Central Implants"); // 1941
	DeclareHistogram1D(42, 1024, "10*log(Dt/1.e-9), Central Decays"); // 1942	
	DeclareHistogram1D(43, 1024, "10*log(Dt/1.e-9), Border Implants"); // 1943
	DeclareHistogram1D(44, 1024, "10*log(Dt/1.e-9), Border Decays"); // 1944	
	    
    /*
    DeclareHistogram2D(DD_MAP_IMPLANT, mapBins, mapBins, "2D MAP Implant direction calib."); // 1939
    DeclareHistogram2D(DD_MAP_DECAY, mapBins, mapBins, "2D MAP Decay direction calib."); // 1940
    DeclareHistogram2D(DD_MAP_IMPLANT_CHE, mapBins, mapBins, "ChE MAP Implant direction calib."); // 1941
    DeclareHistogram2D(DD_MAP_DECAY_CHE, mapBins, mapBins, "ChE MAP Decay direction calib."); // 1942
	*/

	// by Yongchi Xiao, 01/18/2018
	DeclareHistogram2D(46, 2048, 256, "Energy (4 keV/ch) vs. Time (1 us/ch)"); // 1946 
	DeclareHistogram2D(47, 2048, 256, "Energy (4 keV/ch) vs. Time (10 us/ch)"); // 1947
	DeclareHistogram2D(48, 2048, 256, "Energy (4 keV/ch) vs. Time (100 us/ch)"); // 1948
	DeclareHistogram2D(49, 2048, 256, "Energy (4 keV/ch) vs. Time (1 ms/ch)"); // 1949
	DeclareHistogram2D(50, 2048, 256, "Energy (4 keV/ch) vs. Time (10 ms/ch)"); // 1950
	DeclareHistogram2D(51, 2048, 1024, "P1D vs. Correlated decay E. 4 keV/ch"); // 1906 -> 1951
	//
	DeclareHistogram1D(52, 128, "Pulse Shape Discrimination, Decays only"); // 1952
	DeclareHistogram1D(53, 128, "Pulse Shape Discrimination, Recoils only"); // 1953
	DeclareHistogram2D(54, 2048, 128, "Energy (4 keV/ch) vs. PSD"); // 1954

	// correlation matrix
	DeclareHistogram2D(61, 2048, 2048, "Decay Correlation Matrix, 4 keV/ch"); // 1961
    
	// commented out by YX
	/* 1950-1957
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 0, corrBins, corrBins,
					   "2nd Ty,Ex (10ns/ch)(xkeV)"); // 1950
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 1, corrBins, corrBins, 
					   "2nd Ty,Ex (100ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 2, corrBins, corrBins,
					   "2nd Ty,Ex (1us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 3, corrBins, corrBins,
					   "2nd Ty,Ex (10us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 4, corrBins, corrBins, 
					   "2nd Ty,Ex (100us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 5, corrBins, corrBins,
					   "2nd Ty,Ex (1ms/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 6, corrBins, corrBins,
					   "2nd Ty,Ex (10ms/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 7, corrBins, corrBins,
					   "2nd Ty,Ex (100ms/ch)(xkeV)"); // 1957
 
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_TRACEE + 4, corrBins, corrBins,		       "TraceMax Ty,Ex (100us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_TRACEE + 5, corrBins, corrBins, 		       "TraceMax Ty,Ex (1ms/ch)(xkeV)");
    */

  
    // Time gated Decay P1D, 1970-1975
	/*
    DeclareHistogram2D(DD_P1D_DECAY_TIME0, energyBins,p1dBins, "Decay time gated < 1  ms ");
    DeclareHistogram2D(DD_P1D_DECAY_TIME1, energyBins,p1dBins, "Decay time gated < 10 ms");
    DeclareHistogram2D(DD_P1D_DECAY_TIME2, energyBins,p1dBins, "Decay time gated < 100 ms");
    DeclareHistogram2D(DD_P1D_DECAY_TIME3, energyBins,p1dBins, "Decay time gated < 1 s");
    DeclareHistogram2D(DD_P1D_DECAY_TIME4, energyBins,p1dBins, "Decay time gated < 10 s");
    DeclareHistogram2D(DD_P1D_DECAY_TIME5, energyBins,p1dBins, "Decay time gated 5");
	*/


    // Trace Max calibration 
	/*
    DeclareHistogram2D(DD_TRACEMAX_P1D, energyBins, p1dBins,"Decay TraceMax vs P1D"); // 1943
    DeclareHistogram2D(DD_TRACEMAXCAL_P1D, energyBins, p1dBins,"Decay TraceMaxCal vs P1D"); // 1944
    DeclareHistogram2D(DD_TRMAX_QDC, energyBins, energyBins,"Decay TraceMaxCal vs QDCcal"); // 1945
	*/
    // Trace
   
    DeclareHistogram2D(DD_SINGLE_TRACE, traceBins, traceBins2,"Single traces");
    DeclareHistogram2D(DD_DOUBLE_TRACE, traceBins, traceBins2,"Pileup traces"); // 1978
    DeclareHistogram2D(DD_TRACE_POS, traceBins, traceBins2,"Position traces");
    DeclareHistogram2D(DD_TRACE_DYNODE, traceBins, traceBins2,"Dynode traces");

    DeclareHistogram2D(DD_TRACE_IMPLANT_ALL, traceBins, traceBins2,"Trace Implant All");
    DeclareHistogram2D(DD_TRACE_IMPLANT_DYNODE, traceBins, traceBins2,"Trace Implant dynode");
    DeclareHistogram2D(DD_TRACE_IMPLANT_POSITION, traceBins, traceBins2,"Trace Implant position");
    DeclareHistogram2D(DD_TRACE_IMPLANT_PILEUP, traceBins, traceBins2,"Trace Implant pileup");

    DeclareHistogram2D(DD_TRACE_DECAY_ALL, traceBins, traceBins2,"Trace Decay All");
    DeclareHistogram2D(DD_TRACE_DECAY_DYNODE, traceBins, traceBins2,"Trace Decay dynode");
    DeclareHistogram2D(DD_TRACE_DECAY_POSITION, traceBins, traceBins2,"Trace Decay position");
    DeclareHistogram2D(DD_TRACE_DECAY_PILEUP, traceBins, traceBins2,"Trace Decay pileup");

    
    DeclareHistogram2D(DD_MWPC,Bins,Bins,"MWPC in PSPMT Processor");
    DeclareHistogram2D(DD_MWPC_PSPMT,Bins,Bins,"MWPC_PSPMT in PSPMT Processor");
    DeclareHistogram2D(DD_MWPC_NOPSPMT,Bins,Bins,"MWPC_PSPMT in NoPSPMT Processor");
    DeclareHistogram2D(DD_TRACE_E_QDC, energyBins, energyBins,"TraceE vs QDC");

	/* 1996-2007
    DeclareHistogram2D(DD_CHE_REG, regBins, regBins,"ChE vs Regression");
    DeclareHistogram2D(DD_QDC_REG, regBins, regBins,"QDC vs Regression");
    
    DeclareHistogram2D(DD_QDC_REG_IMPLANT, regBins, regBins,"QDC vs Regression Implant");   
    DeclareHistogram2D(DD_QDC_REG_DECAY, regBins, regBins,"QDC vs Regression Decay");   
    DeclareHistogram2D(DD_QDC_REG2_IMPLANT, regBins, regBins,"QDC vs Regression2 Implant");   
    DeclareHistogram2D(DD_QDC_REG2_DECAY, regBins, regBins,"QDC vs Regression Decay2");   
    DeclareHistogram2D(DD_REG12, regBins, regBins,"tau1 vs tau2");   
    DeclareHistogram2D(DD_REG12_IMPLANT, regBins, regBins,"Implant tau1 vs tau2");   
    DeclareHistogram2D(DD_REG12_DECAY, regBins, regBins,"Decay tau1 vs tau2");   
    DeclareHistogram2D(DD_QDC_REG1, regBins, regBins,"QDC vs REG1");   
    DeclareHistogram2D(DD_QDC_REG2, regBins, regBins,"QDC vs REG2");   
	*/
    
    DeclareHistogram2D(DD_TRACE_HIGHENE , traceBins, traceBins2,"Trace Decay High energy"); 
    DeclareHistogram2D(DD_TRACE_HIGHENE2, traceBins, traceBins2,"Trace Decay Alpha"); // 2009
    DeclareHistogram2D(DD_TRACE_HIGHENE3, traceBins, traceBins2,"Trace Decay Beta"); // 2010
    DeclareHistogram2D(DD_TRACE_HIGHENE4, traceBins, traceBins2,"Trace Decay HighQDC/Ech ratio gated"); // 2014

    DeclareHistogram2D(DD_ENE_QDC_DECAY, energyBins, energyBins,"[Decay] Ech vs QDC");
    DeclareHistogram2D(DD_ENE_QDC_IMPLANT, energyBins, energyBins,"[Implantation] Ech vs QDC");
    
    
    DeclareHistogram2D(DD_P1D_QDCCAL, energyBins,Bins, "QDC cal");

    DeclareHistogram2D(DD_QDCREG_SPECIFIC, regBins, regBins,"REG1 vs REG2 (Specific region)");   

}


bool PspmtProcessor::PreProcess(RawEvent &event){
    if (!EventProcessor::PreProcess(event))
		return false;
    
    static const vector<ChanEvent*> &pspmtEvents = sumMap["pspmt"]->GetList();
    
    data_.Clear();
    
    double q1  =0,q2  =0,q3  =0,q4  =0,qd  =0;
    double qdc1=0,qdc2=0,qdc3=0,qdc4=0,qdcd=0;
    double tre1=0,tre2=0,tre3=0,tre4=0,tred=0;
    
    double qright=0,qleft=0,qtop=0,qbottom=0,qsum=0;
    double xright=0,xleft=0,ytop=0,ybottom=0;
    
    double qtre_r=0,qtre_l=0,qtre_t=0,qtre_b=0,qtre_s=0;
    double xtre_r=0,xtre_l=0,ytre_t=0,ytre_b=0;
    
    double qqdc_r=0,qqdc_l=0,qqdc_t=0,qqdc_b=0,qqdc_s=0;
    double xqdc_r=0,xqdc_l=0,yqdc_t=0,yqdc_b=0;
    
    double pxright=0,pxleft=0,pytop=0,pybottom=0;
    double pxtre_r=0,pxtre_l=0,pytre_t=0,pytre_b=0;
	// add by YX
	double pxqdc_r=0,pxqdc_l=0,pyqdc_t=0,pyqdc_b=0;

    double diff;
 
    // tentatively local params //
    double threshold=250;
    double slope=0.0606;
    double intercept=10.13;
    static int traceNum=-1;
    
    for (vector<ChanEvent*>::const_iterator it = pspmtEvents.begin();
         it != pspmtEvents.end(); it++) {
        
        ChanEvent *chan   = *it;
        string subtype    = chan->GetChanID().GetSubtype();
        int    ch         = chan->GetChanID().GetLocation();
        double calEnergy  = chan->GetCalEnergy();
		Trace trace       = chan->GetTrace();
        double trace_energy;
        double trace_time;
        double baseline;
        double qdc;	
	
		if(trace.HasValue("filterEnergy")){	  
			trace_time    = trace.GetValue("filterTime");
			trace_energy  = trace.GetValue("filterEnergy");
			baseline      = trace.DoBaseline(2,20);
			qdc           = trace.DoQDCSimple(20,140);
			if(ch==0){
                qdc1 = qdc/10;
                tre1 = trace_energy;
                plot(D_QDC_TRACE1,qdc1);
                plot(D_ENERGY_TRACE1,tre1);
            }else if(ch==1){
                qdc2 = qdc/10;
                tre2 = trace_energy; 
                plot(D_QDC_TRACE2,qdc2);
                plot(D_ENERGY_TRACE2,tre2);
            }else if(ch==2){
                qdc3 = qdc/10;
                tre3 = trace_energy; 
                plot(D_QDC_TRACE3,qdc3);
                plot(D_ENERGY_TRACE3,tre3);
            }else if(ch==3){
                qdc4 = qdc/10;
                tre4 = trace_energy; 	  
                plot(D_QDC_TRACE4,qdc4);
                plot(D_ENERGY_TRACE4,tre4);
            }else if(ch==4){
                qdcd = qdc;
                tred = trace_energy; 
				plot(D_QDC_TRACED,qdcd/10);
                plot(D_ENERGY_TRACED,tred);
				plot(DD_TRACE_E_QDC,tred,qdcd/10);
			}
        }
		
		
		// below is dealing with on-board energy
		/*
        if(ch==0){
            q1 = calEnergy;
            plot(D_RAW1,q1);
        }else if(ch==1){
            q2 = calEnergy;
            plot(D_RAW2,q2);
        }else if(ch==2){
            q3 = calEnergy;
            plot(D_RAW3,q3);
        }else if(ch==3){
            q4 = calEnergy;
            plot(D_RAW4,q4);
        }else if(ch==4){
            qd = calEnergy;
            plot(D_RAWD,qd);
        }
		*/
        /*
        if(q1>0 && q2>0 && q3>0 && q4>0){
		// q2     q1
		// 
		// 
		// q3     q4
		//
		// energy on the border
		    qtop    = (q1+q2)/2;
			qleft   = (q2+q3)/2;
			qbottom = (q3+q4)/2;
			qright  = (q4+q1)/2;
			// deduce the position          
            qsum    = (q1+q2+q3+q4)/2;
            xright  = (qright/qsum)*512+100;
            xleft   = (qleft/qsum)*512+100;
            ytop    = (qtop/qsum)*512+100;
            ybottom = (qbottom/qsum)*512+100;
            plot(D_SUM,qsum);
        }

        if(q1>threshold && q2>threshold && q3>threshold && q4>threshold ){
            pxleft   = trunc(slope*xleft-intercept);
            pxright  = trunc(slope*xright-intercept);
            pytop    = trunc(slope*ytop-intercept);
            pybottom = trunc(slope*ybottom-intercept);
            
            plot(DD_POS1_RAW,xright,ytop); // 1906
            plot(DD_POS2_RAW,xleft,ybottom); // 1907
            plot(DD_POS1,pxright,pytop); // 1908
            plot(DD_POS2,pxleft,pybottom); // 1909
	    
		}
		*/
		

		// below are energy from in-code filters        
		/*
        if(tre1>0 && tre2>0 && tre3>0 && tre4>0 ){
            qtre_t=(tre1+tre2)/2;
            qtre_l=(tre2+tre3)/2;
            qtre_b=(tre3+tre4)/2;
            qtre_r=(tre4+tre1)/2;
            qtre_s=(tre1+tre2+tre3+tre4)/2;
            
            xtre_r=(qtre_r/qtre_s)*512+100;
            xtre_l=(qtre_l/qtre_s)*512+100;
            ytre_t=(qtre_t/qtre_s)*512+100;
            ytre_b=(qtre_b/qtre_s)*512+100;
            
            pxtre_r = trunc(slope*xtre_r-intercept);
            pxtre_l = trunc(slope*xtre_l-intercept);
            pytre_t = trunc(slope*ytre_t-intercept);
            pytre_b = trunc(slope*ytre_b-intercept);
                        
            plot(D_ENERGY_TRACESUM,qtre_s); // 1915
	    
            if(tre1>threshold && tre2>threshold && tre3>threshold && tre4>threshold ){
				plot(DD_POS1_RAW_TRACE,xtre_r,ytre_t); // 1916
				plot(DD_POS2_RAW_TRACE,xtre_l,ytre_b);
				plot(DD_POS1_TRACE,pxtre_r,pytre_t);
				plot(DD_POS2_TRACE,pxtre_l,pytre_b);
            }    
        }
		*/
        
		// below are processing of QDC
		if(qdc1>0 && qdc2>0 && qdc3>0 && qdc4>0 ){
            qqdc_t = (qdc1+qdc2)/2;
            qqdc_l = (qdc2+qdc3)/2;
            qqdc_b = (qdc3+qdc4)/2;
            qqdc_r = (qdc4+qdc1)/2;
            qqdc_s = (qqdc_t+qqdc_l+qqdc_b+qqdc_r)/2;
			
            xqdc_r = (qqdc_r/qqdc_s)*512+100;
            xqdc_l = (qqdc_l/qqdc_s)*512+100;
            yqdc_t = (qqdc_t/qqdc_s)*512+100;
            yqdc_b = (qqdc_b/qqdc_s)*512+100;

			pxqdc_r = trunc(slope*xqdc_r-intercept);
            pxqdc_l = trunc(slope*xqdc_l-intercept);
            pyqdc_t = trunc(slope*yqdc_t-intercept);
            pyqdc_b = trunc(slope*yqdc_b-intercept);
            
			// by YX
			plot(16, xqdc_r, yqdc_t); // 1916
			plot(17, xqdc_l, yqdc_b); // 1917
			plot(18, pxqdc_r, pyqdc_t); // 1918
			plot(19, pxqdc_l, pyqdc_b); // 1919
        }
        
    } // end of one channel event
    
    EndProcess();
	return(true);
}

// for correlations
static PixelEvent implantRec[600] = {};
static PixelEvent decayRec[2][600] = {};
//
static double implantRefTime[600] = {};
static double decayRefTime[600] = {};
// for position calibrations
const double parX[4] = {2.3524e-7, -0.00018115, 0.10339, -7.1697};
const double parY[4] = {1.5449e-7, -0.00011528, 0.084563, -5.1789}; 
//
fstream outfile; 

bool PspmtProcessor::Process(RawEvent &event){
	if (!EventProcessor::Process(event))
		return false;
  
	vector<ChanEvent*> pspmtEvents = 
		event.GetSummary("pspmt:pspmt", true)->GetList();
	vector<ChanEvent*> naiEvents = 
		event.GetSummary("nai:nai", true)->GetList();
	vector<ChanEvent*> mwpcEvents = 
		event.GetSummary("mcp:1time", true)->GetList();
	vector<ChanEvent*> mwpcEEvents = 
		event.GetSummary("mcp:1position1", true)->GetList();
  
	bool has_pspmt   = false;
	bool has_mwpc    = false;
	bool has_nai     = false;
	bool has_implant = false;
	bool has_decay   = false;
	bool has_veto    = false;
	bool has_pileup  = false;
	bool has_bigqdc  = false;
	bool has_all_anode_qdc = false; 
	
	int  mult_pspmt = event.GetSummary("pspmt",true)->GetMult();
	int  mult_mwpc  = event.GetSummary("mcp",true)->GetMult();
	int  mult_nai   = event.GetSummary("nai",true)->GetMult();
	// by YX
	bool canProcess = false;
	
	if(mult_pspmt==5){
		has_pspmt=true;
	}
	if(mult_pspmt<5){
		has_pspmt=false;
	}
	if(mult_mwpc>0) {
		has_mwpc=true; // flag for mwpc
	}
	if(mult_nai>0){
		has_nai=true;
	}
	if(has_mwpc && has_pspmt){
		has_implant = true; // flag for implant
	}
	if(!has_mwpc && has_pspmt){
		has_decay   = true; // flag for decay
	}
	if(has_nai){
		has_veto    = true;  
	}
  
	double q1=0,q2=0,q3=0,q4=0,qd=0;
	double qtop=0,qleft=0,qbottom=0,qright=0;
	double qsum=0;
	double xright=0,xleft=0,ytop=0,ybottom=0;
	int px_r=-1,px_l=-1,py_t=-1,py_b=-1;
	double tre1=0,tre2=0,tre3=0,tre4=0,tred=0;
	double qdc1=0,qdc2=0,qdc3=0,qdc4=0,qdcd=0,qdcs=0;
	double qdcd_cal=0;
  
	double qdc_top=0,qdc_left=0,qdc_bottom=0,qdc_right=0;
	double yqdc_top=0,xqdc_left=0,yqdc_bottom=0,xqdc_right=0;
	double pyqdc_top=0,pxqdc_left=0,pyqdc_bottom=0,pxqdc_right=0;
	//
	double xPos, yPos, xPrePos, yPrePos; 
  
	double xcal=0,ycal=0;
	double xcal2=0,ycal2=0;
	static int traceNum=0;
	static int traceNumPosition=0;
	static int traceNumPositionDecay=0;
	static int traceNumPositionImplant=0;
	static int traceNumDynode=0;
	static int traceNumDynodeDecay=0;
	static int traceNumDynodeImplant=0;
	static int traceNumSecond=0;
	static int traceNumSecondImplant=0;
	static int traceNumSecondDecay=0;
	static int traceNumImplant=0;
	static int traceNumDecay=0;
	static int traceNumHigh=0;
	static int traceNumHigh2=0;
	static int traceNumAlpha=0;
	static int traceNumBeta=0;
	static int traceNumbigqdc=0;

	// by Yongchi Xiao, 01/18/2018
	double qdcCalib = -1; 
  
	int p1d;
	double threshold=0;
	// For position calib at PSPMT
	
	//Original SG setting
	double slope=0.0606;
	double intercept=10.13;

	// by Yongchi Xiao, for PSD purpose, 03/06/2018
	double psd = 0; 
	/*
	double slope = 5.952381e-02; 
	double intercept = 0; 
	*/
	// 
	/*
	double xslope=0.0623;
	double xoffset=-10.791;
	double yslope=0.0601;
	double yoffset=-9.780;
	*/
	/*
	double slopeX = 0.0613976395; 
	double slopeY = 0.05899345; 
	double offsetX = -4.3216135359; 
	double offsetY = -3.5612959407; 
	*/
	double offset_mirror=341;
	double offset_mirror2=365;
	double QDCHIGHGATE=2000;
	double QDCLOWGATE=700;
  
	double QDCHIGHALPHA=330;
	double QDCLOWALPHA=270;
	double QDCHIGHBETA=180;
	double QDCLOWBETA=100;
  
	double x1=0,x2=251;
	double y1=250,y2=510;
	double a_reg,b_reg;
 
	//a_reg = (y2-y1)/(x2-x1);
	//b_reg = y1;
	a_reg = 1;
	b_reg = 0;
  
  
	double static pspmttime;
	double static mwpctime;
	double mwpcene;
	int pspmtch=-1;
	int numpulse=0;
	int mwpcch;
	double regression=0;
	double regression2=0;
	int mwpcCh;
	double  mwpcE;
	Trace traceD;
	Trace trace1,trace2,trace3,trace4;
	double max_dynode;
  
	// in processor //  
	for (vector<ChanEvent*>::iterator itm = mwpcEvents.begin();
		 itm != mwpcEvents.end();
		 ++itm) {
		mwpctime = (*itm)->GetTime();
		mwpcene  = (*itm)->GetCalEnergy();
		mwpcch   = (*itm)->GetChanID().GetLocation();
	}
	for (vector<ChanEvent*>::iterator itm = mwpcEEvents.begin();
		 itm != mwpcEEvents.end();
		 ++itm) {
		mwpcE    = (*itm)->GetCalEnergy();
		mwpcCh   = (*itm)->GetChanID().GetLocation();
	}
 
	plot(DD_MWPC,mwpcene/10,mwpcE);
	if(has_pspmt){
		plot(DD_MWPC_PSPMT,mwpcene/10,mwpcE);
	}else if(!has_pspmt){
		plot(DD_MWPC_NOPSPMT,mwpcene/10,mwpcE);
	}
   
	for (vector<ChanEvent*>::const_iterator it = pspmtEvents.begin();
		 it != pspmtEvents.end(); it++) {        
		ChanEvent *chan   = *it;
		int ch            = chan->GetChanID().GetLocation();
		pspmtch           = ch;
		pspmttime         = chan->GetTime();
		double energy     = chan->GetCalEnergy();
		Trace  trace      = chan->GetTrace();        
    
		if(ch==0){
			q1 = energy;
			trace1 = chan->GetTrace();
		}else if(ch==1){
			q2 = energy;
			trace2 = chan->GetTrace();
		}else if(ch==2){
			q3 = energy;
			trace3 = chan->GetTrace();
		}else if(ch==3){
			q4 = energy;
			trace4 = chan->GetTrace();
		}else if(ch==4){
			qd = energy;
			traceD = chan->GetTrace();
		}
    
		if(trace.HasValue("filterEnergy")){          
			double trace_energy;
			double trace_time;
			double baseline;
			double qdc;
			double ch_max;
			trace_energy  = trace.GetValue("filterEnergy");
			trace_time    = trace.GetValue("filterTime");
			baseline      = trace.DoBaseline(2,20);
			qdc           = trace.DoQDCSimple(20,140);
			numpulse      = trace.GetValue("numPulses");
			
			
			if(ch==0){
				qdc1 = qdc;
				tre1 = trace_energy;
			}else if(ch==1){
				qdc2 = qdc;
				tre2 = trace_energy; 
			}else if(ch==2){
				qdc3 = qdc;
				tre3 = trace_energy; 
			}else if(ch==3){
				qdc4 = qdc;
				tre4 = trace_energy; 	  
			}else if(ch==4){
				qdcd = qdc;
				tred = trace_energy; 
				has_pileup=trace.ScanPileup(1,140);
				max_dynode=trace.GetTraceMax(1,140);
				regression    = abs(3000*trace.DeduceRegression(20,140,0));
				regression2   = abs(3000*trace.DeduceRegression(20,140,1));
				if(!has_pileup) {
					psd = trace.DoPSD(1, 90); 
					//					cout << psd << endl;
				}
			}
		}
    
    
	} // end of channel event

	////////////////////////////
	// based on qdc
	if(qdc1 > 0 && qdc2 > 0 && qdc3 > 0 && qdc4 > 0) {
		has_all_anode_qdc = true; 

	}
	
	if(has_all_anode_qdc) {
		qdc_top    = (qdc1+qdc2)/2;
		qdc_left   = (qdc2+qdc3)/2;
		qdc_bottom = (qdc3+qdc4)/2;
		qdc_right  = (qdc4+qdc1)/2;
		qdcs=(qdc1+qdc2+qdc3+qdc4)/2;

		// almost raw information
		xqdc_right  = (qdc_right/qdcs)*512 + 100;
		xqdc_left   = (qdc_left/qdcs)*512 + 100;
		yqdc_top    = (qdc_top/qdcs)*512 + 100;
		yqdc_bottom = (qdc_bottom/qdcs)*512 + 100;
		// by YX
		xPos = xqdc_right - 100; 
		yPos = yqdc_top - 100; 
		xPrePos = xPos; 
		yPrePos = yPos;
		xPos = parX[0]*pow(xPos,3) + parX[1]*pow(xPos, 2) + parX[2]*xPos + parX[3];
		yPos = parY[0]*pow(yPos,3) + parY[1]*pow(yPos, 2) + parY[2]*yPos + parY[3];
		// by YX
		pxqdc_right = trunc(xPos + 0.5);
		pyqdc_top = trunc(yPos + 0.5);

		pxqdc_left = trunc(slope*xqdc_left-intercept);
		pyqdc_bottom = trunc(slope*yqdc_bottom-intercept);

		xcal  = yqdc_top;
		ycal  = -1*(xqdc_right-offset_mirror)+offset_mirror;
		xcal2 = -1*(yqdc_bottom-offset_mirror2)+offset_mirror2;
		ycal2 = xqdc_left;
	
		p1d     = pxqdc_right + 24*pyqdc_top;
		// by YX
		if(p1d > 0 && p1d <= 576 & has_pspmt) {
			qdcCalib=(qdcs/40.*pixelCalib[p1d]); 
			canProcess = true; 
		}
		plot(DD_P1D_QDCCAL,qdcd_cal,p1d);  // QDC vs P1D
  
		/* by Yongchi Xiao
		 * make correlations beyond this point 
		 */
		/* by Yongchi Xiao, 01/26/2018
		 * blank out border and corner pixels, check correlations in the center pixels
		 * Simply cut off the first and last third pixels (1d)
		 */
		if(canProcess) { 
			// plot reproduced "raw" position, in analogy with 1916/1917
			if(qdcCalib > 0) {
				plot(6, xqdc_right  , yqdc_top  ); // 1906
				plot(7, xqdc_left  , yqdc_bottom  ); // 1907
				plot(8, pxqdc_right, pyqdc_top); // 1908
				plot(9, pxqdc_left, pyqdc_bottom); // 1909
			}		
			double timeDiffImplant = 0;
			double timeDiffDecay = 0;  
			if(has_implant && qdcCalib > 0) {
				plot(40, pxqdc_right, pyqdc_top); // 1940
				plot(34, qdcCalib, p1d); // 1934
				plot(53, psd*100); // 1953

				implantRec[p1d].AssignValues(qdcCalib, pspmttime, xPrePos, yPrePos); 
				// remove decay events
				for(int i = 0; i < 2; i++) {
					decayRec[i][p1d].Clear();
				}
			} else if(has_decay && qdcCalib > 0) {			
				plot(33, qdcCalib, p1d); // 1933, decay only
				if(!has_pileup && !has_veto) {
					plot(52, psd*100); // 1952
				}
				if(implantRec[p1d].Is_Filled() 
				   && abs(xPrePos - implantRec[p1d].GetRawPos().first) < 5.
				   && abs(yPrePos - implantRec[p1d].GetRawPos().second) < 6. ) { // the preceding ion has been found
					if(!decayRec[0][p1d].Is_Filled()) {
						decayRec[0][p1d].AssignValues(qdcCalib, pspmttime, xPrePos, yPrePos); 
						timeDiffImplant = (pspmttime - implantRec[p1d].GetTime())*Globals::get()->clockInSeconds();
						// plot 1946-1949
						if(!has_veto) {
							plot(46, qdcCalib, timeDiffImplant*1e6); 
							plot(47, qdcCalib, timeDiffImplant*1e5);
							plot(48, qdcCalib, timeDiffImplant*1e4);
							plot(49, qdcCalib, timeDiffImplant*1e3); 
							plot(50, qdcCalib, timeDiffImplant*1e2);
							plot(51, qdcCalib, p1d); // 1906 -> 1951
							// Energy vs. PSD
							plot(54, qdcCalib, psd*100); // 1954
							// output decay events falling within desired energy range
							/*
							if(abs(qdcCalib - energyCentroid) <= energyFWHM) {
								outfile.open((runName + ".scanout").c_str(), std::iostream::out | std::iostream::app); 
								outfile << pxqdc_right << "  " << pyqdc_top << "  " << p1d << "  " 
										<< qdcCalib << "  " << timeDiffImplant << endl;
								outfile.close(); 
							}
							*/
						}
					} else if(!decayRec[1][p1d].Is_Filled()) {
						decayRec[1][p1d].AssignValues(qdcCalib, pspmttime, xPrePos, yPrePos); 
						/*
						decayRec[1][p1d].time = pspmttime;
						decayRec[1][p1d].GetEnergy() = qdcCalib;
						*/
						timeDiffImplant = (pspmttime - implantRec[p1d].GetTime())*Globals::get()->clockInSeconds();
						/* Find correlated decays: 
						 * T1/2 for 109Xe: 13 ms
						 * T1/2 for 105Te: 0.62 us
						 */
						if((decayRec[1][p1d].GetTime() - decayRec[0][p1d].GetTime())*Globals::get()->clockInSeconds() < 0.5e-3
						   && (decayRec[1][p1d].GetTime() - decayRec[0][p1d].GetTime()) > 0) {
							plot(61, decayRec[0][p1d].GetEnergy(), decayRec[1][p1d].GetEnergy()); // 1961
						}
					} 
				}
			}
		} // end of gating on pixels

		/* Gate on central pixels only
		 * to see the count rate of 
		 * decays and implants
		 */
		if(canProcess) {
			if(abs(pxqdc_right - 11.5) < 4 && abs(pyqdc_top - 11.5) < 4) {// central pixels
				if(has_implant) {
					if(implantRefTime[p1d] > 0) {
						plot(41, 10*log((pspmttime - implantRefTime[p1d])*Globals::get()->clockInSeconds()/1.e-9)); // 1941
					} 
					implantRefTime[p1d] = pspmttime; 
				} else if(has_decay && !has_veto) {
					if(decayRefTime[p1d] > 0) {
						plot(42, 10*log((pspmttime - decayRefTime[p1d])*Globals::get()->clockInSeconds()/1.e-9)); // 1942
					}
					decayRefTime[p1d] = pspmttime; 
				}
			} else { // border pixels
				if(has_implant) {
					if(implantRefTime[p1d] > 0) {
						plot(43, 10*log((pspmttime- implantRefTime[p1d])*Globals::get()->clockInSeconds()/1.e-9)); // 1943
					} 
					implantRefTime[p1d] = pspmttime; 
				} else if(has_decay && !has_veto) {
					if(decayRefTime[p1d] > 0) {
						plot(44, 10*log((pspmttime- decayRefTime[p1d])*Globals::get()->clockInSeconds()/1.e-9)); // 1944
					}
					decayRefTime[p1d] = pspmttime; 
				}
			}
		}
		
	}	
	/********************************* Below are for traces **************************************/    
  
	for(vector<int>::iterator ittr = traceD.begin();ittr != traceD.end();ittr++){    
		plot(DD_SINGLE_TRACE,ittr-traceD.begin(),traceNum,*ittr);
		plot(DD_TRACE_DYNODE,ittr-traceD.begin(),traceNumDynode,*ittr);
		if(has_pileup){
			plot(DD_DOUBLE_TRACE,ittr-traceD.begin(),traceNumSecond,*ittr); // 1978, dynode double trace
		}
     
		if(has_implant){
			plot(DD_TRACE_IMPLANT_ALL,ittr-traceD.begin(),traceNumImplant,*ittr);
			plot(DD_TRACE_IMPLANT_DYNODE,ittr-traceD.begin(),traceNumDynodeImplant,*ittr);
			if(has_pileup){
				plot(DD_TRACE_IMPLANT_PILEUP,ittr-traceD.begin(),traceNumSecondImplant,*ittr);
			}
		}
		if(has_decay && !has_veto){
			plot(DD_TRACE_DECAY_ALL,ittr-traceD.begin(),traceNumDecay,*ittr);
			plot(DD_TRACE_DECAY_DYNODE,ittr-traceD.begin(),traceNumDynodeDecay,*ittr);
			if(qdcd_cal<QDCHIGHGATE && qdcd_cal>QDCLOWGATE){
				plot(DD_TRACE_HIGHENE,ittr-traceD.begin(),traceNumHigh,*ittr);
				plot(DD_QDCREG_SPECIFIC,qdcd_cal,regression);
				if(has_bigqdc){
					plot(DD_TRACE_HIGHENE4,ittr-traceD.begin(),traceNumbigqdc,*ittr);
				}
			}
       
			if(qdcd_cal<QDCHIGHALPHA && qdcd_cal>QDCLOWALPHA){
				plot(DD_TRACE_HIGHENE2,ittr-traceD.begin(),traceNumAlpha,*ittr);
			}
			if(qdcd_cal<QDCHIGHBETA && qdcd_cal>QDCLOWBETA){
				plot(DD_TRACE_HIGHENE3,ittr-traceD.begin(),traceNumBeta,*ittr);
			}
			if(has_pileup){
				plot(DD_TRACE_DECAY_PILEUP,ittr-traceD.begin(),traceNumSecondDecay,*ittr);
			}
		}
	} 
     
   
	bool has_dynode=false;
	if(traceD.HasValue("filterEnergy")){
		has_dynode=true;
	}else{
		has_dynode=false;
	}
	bool has_trace1=false;
	bool has_trace2=false;
	bool has_trace3=false;
	bool has_trace4=false;
     
	if(trace1.HasValue("filterEnergy")){
		has_trace1=true;
	}
	if(trace2.HasValue("filterEnergy")){
		has_trace2=true;
	}
	if(trace3.HasValue("filterEnergy")){
		has_trace3=true;
	}
	if(trace4.HasValue("filterEnergy")){
		has_trace4=true;
	}
        
	//////////////////////////
	// Counter for traces   // 
	//////////////////////////
     
	traceNumPosition++;
	if(has_dynode){
		traceNum++;   	  
		traceNumDynode++;
	}
	if(has_pileup){
		traceNumSecond++;
	}
	if(has_implant){
		traceNumPositionImplant++;
		traceNumImplant++;
		if(has_dynode){
			traceNumDynodeImplant++;
			if(has_pileup){
				traceNumSecondImplant++;
			}
		}
	}
	if(has_decay && !has_veto){
		traceNumPositionDecay++;
		traceNumDecay++;
		if(has_dynode){
	  
			traceNumDynodeDecay++;
	  
			if(qdcd_cal<QDCHIGHGATE && qdcd_cal>QDCLOWGATE){
				traceNumHigh++;
			}
			if(qdcd_cal<QDCHIGHBETA && qdcd_cal>QDCLOWBETA){
				traceNumBeta++;
			}
			if(qdcd_cal<QDCHIGHALPHA && qdcd_cal>QDCLOWALPHA){
				traceNumAlpha++;
			}
			if(has_pileup){
				traceNumSecondDecay++;
			}
	 
		}
       
	}

	/********************************** End of processing traces ********************************************/    
     
      
	/* Commented out by YX 
	 */
	/*
	// Correlation stuffs
	if(q1>threshold && q2>threshold && q3>threshold && q4>threshold ){
		PspmtEvent ev(qdcd_cal,trmax,pspmttime,mwpctime,has_implant,has_decay,has_veto);
		//       cout << qdcd_cal << " " << trmax << endl;
		AddEvent(ev,px_r,py_t);
	}
	*/
   
   
	EndProcess();
	return(true);
}

bool PspmtProcessor::AddEvent(PspmtEvent& event,int x,int y){
  
	int size_x=24;
	int size_y=24;
  
	if(x<0 || x>=size_x){
		return false;
	}
	if(y<0 || y>=size_y){
		return false;
	}
  
	bool has_beam  = event.get_is_beam();
	bool has_decay = event.get_is_decay();
 
	//In case of getting heavyIons, flush it and refill it
	if (has_beam){
		FlushChain(x, y);
	}
 
	pixels_[x][y].push_back(event);
 
	// In case of getting decay events, flush it
	if (has_decay){
		FlushChain(x, y);
	}
  
  
	return true;
}
bool PspmtProcessor::FlushChain(int x,int y){
  
	unsigned chain_size = pixels_[x][y].size();
	if(chain_size<1){
		pixels_[x][y].clear();
		return false;
	}
  
	int p1d=x+24*y; 
	PspmtEvent ev_first = pixels_[x][y].front();
	bool first_is_beam   = ev_first.get_is_beam();
	// bool first_is_decay  = ev_first.get_is_decay();
  
	if(!first_is_beam){
		pixels_[x][y].clear();
		return false;
	}
  
	int alphas = 0;
	double alphaE[6]={0};
	double alphaTRE[6]={0};
	double alphaTime[6]={0};
	double dt[6];
	double BeamE=0,BeamTime=0;
	int pix_counter=0;  
  

	for (deque<PspmtEvent>::iterator it = pixels_[x][y].begin();
		 it != pixels_[x][y].end();
		 ++it) { 	
    
		pix_counter++;
		bool b_ion    = (*it).get_is_beam();
		bool b_decay  = (*it).get_is_decay();
		bool b_veto   = (*it).get_is_veto();
		double energy = (*it).get_ene_pspmt();
		double energy_tr = (*it).get_ene_trace();
		double time   = (*it).get_time_pspmt()*(Globals::get()->clockInSeconds()); 
		double timecheck = (*it).get_time_pspmt();
   
    
		if(b_ion) { 
			BeamE=energy;
			BeamTime=time;
		}    	
		if(b_decay & !b_veto) {
			alphas += 1;
			if(alphas <= 6) {
				alphaE[alphas]=energy;
				alphaTRE[alphas]=energy_tr;
				alphaTime[alphas] = time; 
			}
      
			//      cout << alphaE[alphas] << " " << alphaTRE[alphas] << endl;
		}
  
    
		if(BeamTime>0 && alphaTime[1]>0){
			const unsigned int NumGranularities = 8;
			for (unsigned int i = 0; i < NumGranularities; i++) {
				const double timeResolution[NumGranularities] = 
				/// 10ns, 100ns,  1us,  10us, 100us, 1ms, 10ms, 100ms
					{10e-9, 100e-9, 1e-6, 10e-6, 100e-6, 1e-3, 10e-3, 100e-3};
				dt[i]=100000000*(alphaTime[1]-BeamTime)*(Globals::get()->clockInSeconds()/timeResolution[i]);
	
	
				//	cout << "alphaTRE " <<alphaTRE[1] << endl;
	
				plot(DD_ENERGY_DECAY_TIME_GRANX+i,alphaE[1],dt[i]);
	

				if(i==4){
					plot(DD_ENERGY_DECAY_TIME_GRANX_TRACEE+4,alphaTRE[1],dt[i]);
				}
				if(i==5){
					plot(DD_ENERGY_DECAY_TIME_GRANX_TRACEE+5,alphaTRE[1],dt[i]);
				}
	
				if(i==4){ // 100us/ch,  1ms/10ch, 10ms/100ch
					if(dt[i]<10){ // 1ms
						plot(DD_P1D_DECAY_TIME0,alphaE[1],p1d);
					}
					if(dt[i]<100){ // 10ms
						plot(DD_P1D_DECAY_TIME1,alphaE[1],p1d);
					}
					if(dt[i]<1000){ // 100ms
						plot(DD_P1D_DECAY_TIME2,alphaE[1],p1d);
					}
					if(dt[i]<10000){ // 1s
						plot(DD_P1D_DECAY_TIME3,alphaE[1],p1d);
					}
					if(dt[i]<100000){ // 10s
						plot(DD_P1D_DECAY_TIME4,alphaE[1],p1d);
					}
	  
				}
			}
		}
	}
  
	pixels_[x][y].clear();
  
	return true;
}

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
#include "DoubleTraceAnalyzer.hpp"
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
	/*
    DeclareHistogram1D(D_RAW1, energyBins, "Pspmt1 Raw");
    DeclareHistogram1D(D_RAW2, energyBins, "Pspmt2 Raw");
    DeclareHistogram1D(D_RAW3, energyBins, "Pspmt3 Raw");
    DeclareHistogram1D(D_RAW4, energyBins, "Pspmt4 Raw");
    DeclareHistogram1D(D_RAWD, energyBins, "Pspmt Dynode");
	*/
	/*
    DeclareHistogram1D(D_SUM,  energyBins, "Pspmt Sum");
    DeclareHistogram2D(DD_POS1_RAW, mapBins, mapBins, "Pspmt Pos1 Raw"); // 1906
    DeclareHistogram2D(DD_POS2_RAW, mapBins, mapBins, "Pspmt Pos2 Raw"); // 1907
    DeclareHistogram2D(DD_POS1, posBins, posBins, "Pspmt Pos1"); // 1908
    DeclareHistogram2D(DD_POS2, posBins, posBins, "Pspmt Pos2"); // 1909
	*/
	/*
	DeclareHistogram2D(6, 1024, 1024, "2-D map of all signals, reproduced - I"); // 1906
	DeclareHistogram2D(7, 1024, 1024, "2-D map of all signals, reproduced - II"); // 1907
	DeclareHistogram2D(8, 32, 32, "Pixelized  2-D map, reproduced - I"); // 1908
	DeclareHistogram2D(9, 32, 32, "Pixelized 2-D map, reproduced - II"); // 1909
	*/
	// spectra 1900s
	DeclareHistogram1D(1, 16, "Accumulation of ions vs. decays"); // 1901
	DeclareHistogram2D(2, 32, 32, "Pixelated Map of Ions"); // 1902
	DeclareHistogram2D(3, 2048, 2048, "Raw Map of Ions"); // 1903
	DeclareHistogram2D(4, 32, 32, "Pixelated Map of Decays"); // 1904
	DeclareHistogram2D(5, 2048, 2048, "Raw Map of Decays"); // 1905
	

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
   
    /*
    // Energy resolutions 30-
    DeclareHistogram2D(DD_P1D_CHANNEL, energyBins,Bins, "Ch vs Dynode EChannel"); // 1930
    DeclareHistogram2D(DD_P1D_TRACE, energyBins,Bins, "Ch vs Dynode ETrace"); // 1931
    DeclareHistogram2D(DD_P1D_QDC, energyBins,1024, "Ch vs Dynode QDC scaled by 10"); // 1932
    DeclareHistogram2D(DD_P1D_QDCSUM, energyBins,1024, "P1D vs QDCSum 4 keV/ch, Decays"); // 1933
	DeclareHistogram2D(34, energyBins, 1024, "P1D vs. QDCSum 4 keV/ch, Implants"); // 1934
    // Energy res by correlation
    DeclareHistogram2D(DD_P1D_IMPLANT_CHE, energyBins,p1dBins, "[Implant] ch vs E(ch)"); // 1935
    DeclareHistogram2D(DD_P1D_DECAY_CHE, energyBins,p1dBins,   "[Decay] ch vs E(ch)"); // 1936
    DeclareHistogram2D(DD_P1D_IMPLANT_QDC, energyBins,p1dBins, "[Implant] ch vs E(QDC)"); // 1937
    DeclareHistogram2D(DD_P1D_DECAY_QDC, energyBins,p1dBins, "[Decay] ch vs E(QDC)"); // 1938
	*/
	
	/* Spectrum 1930-1939
	 * For PSD and other properties of signals
	 */ 
	DeclareHistogram1D(31, 128, "PSD, Ions"); // 1931
	DeclareHistogram1D(32, 128, "PSD, Decays"); // 1932
	

	// by Yongchi
	    
    /*
    DeclareHistogram2D(DD_MAP_IMPLANT, mapBins, mapBins, "2D MAP Implant direction calib."); // 1939
    DeclareHistogram2D(DD_MAP_DECAY, mapBins, mapBins, "2D MAP Decay direction calib."); // 1940
    DeclareHistogram2D(DD_MAP_IMPLANT_CHE, mapBins, mapBins, "ChE MAP Implant direction calib."); // 1941
    DeclareHistogram2D(DD_MAP_DECAY_CHE, mapBins, mapBins, "ChE MAP Decay direction calib."); // 1942
	*/

	// by Yongchi Xiao, 01/18/2018
	/* for time-energy correlation
	 * damm ID = 1940-1949
	 */ 
	// for ions
	DeclareHistogram2D(40, 2048, 256, "Ion E. vs. Time (1 us/ch), Ions"); // 1940 
	DeclareHistogram2D(41, 2048, 256, "Ion E. vs. Time (10 us/ch), Ions"); // 1941
	DeclareHistogram2D(42, 2048, 256, "Ion E. vs. Time (100 us/ch), Ions"); // 1942
	DeclareHistogram2D(43, 2048, 256, "Ion E. vs. Time (1 ms/ch), Ions"); // 1943
	DeclareHistogram2D(44, 2048, 256, "Ion E. vs. Time (10 ms/ch, Ions)"); // 1944
	// for decays
	DeclareHistogram2D(45, 2048, 256, "Decay E. vs. Time (1 us/ch)"); // 1945 
	DeclareHistogram2D(46, 2048, 256, "Decay E. vs. Time (10 us/ch)"); // 1946
	DeclareHistogram2D(47, 2048, 256, "Decay E. vs. Time (100 us/ch)"); // 1947
	DeclareHistogram2D(48, 2048, 256, "Decay E. vs. Time (1 ms/ch)"); // 1948
	DeclareHistogram2D(49, 2048, 256, "Decay E. vs. Time (10 ms/ch)"); // 1949
	// 1950s, for the same purpose
	DeclareHistogram2D(50, 2048, 4096, "Decays(x) vs. Ions(y), < 256 us"); // 1950
	DeclareHistogram2D(51, 2048, 4096, "Decays(x) vs. Ions(y), < 256 10us"); // 1951
	DeclareHistogram2D(52, 2048, 4096, "Decays(x) vs. Ions(y), < 256 100us"); // 1952
	DeclareHistogram2D(53, 2048, 4096, "Decays(x) vs. Ions(y), < 256 1ms"); // 1953
	DeclareHistogram2D(54, 2048, 4096, "Decays(x) vs. Ions(y), < 256 10ms"); // 1954
	
	/* Decays only 
	 * 1960s
	 */ 
	DeclareHistogram2D(63, 4096, 1024, "Summed QDC vs. P1D, Decays"); // 1963
	DeclareHistogram2D(64, 2048, 2048, "Pile-Up signals, E1 vs. E2, same pixel"); // 1964
	DeclareHistogram2D(65, 4096, 1024, "Summed QDC vs. P1D, Ions"); // 1965
	/* trace analysis, from 1970 to 1979
	 */
	DeclareHistogram2D(70, 512, 1024, "Summed Trace, pileup All"); // 1970
	DeclareHistogram2D(71, 512, 1024, "Summed Trace, pileup Top"); // 1971
	DeclareHistogram2D(72, 512, 1024, "Summed Trace, pileup Left"); // 1972
	DeclareHistogram2D(73, 512, 1024, "Summed Trace, pileup Bottom"); // 1973
	DeclareHistogram2D(74, 512, 1024, "Summed Trace, pileup Right"); // 1974	
	//
	DeclareHistogram2D(75, 512, 1024, "Dynode Trace, pileup"); // 1975
	//
	DeclareHistogram2D(76, 512, 1024, "Anode Trace 0, pileup"); // 1976
	DeclareHistogram2D(77, 512, 1024, "Anode Trace 1, pileup"); // 1977
	DeclareHistogram2D(78, 512, 1024, "Anode Trace 2, pileup"); // 1978
	DeclareHistogram2D(79, 512, 1024, "Anode Trace 3, pileup"); // 1979


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
	/*   
    DeclareHistogram2D(DD_SINGLE_TRACE, traceBins, traceBins2,"Single traces");
    DeclareHistogram2D(DD_DOUBLE_TRACE, 512, traceBins2,"Pileup traces"); // 1978
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
	*/
    
    DeclareHistogram2D(DD_MWPC,Bins,Bins,"MWPC in PSPMT Processor");
    DeclareHistogram2D(DD_MWPC_PSPMT,Bins,Bins,"MWPC_PSPMT in PSPMT Processor");
    DeclareHistogram2D(DD_MWPC_NOPSPMT,Bins,Bins,"MWPC_PSPMT in NoPSPMT Processor");
	//    DeclareHistogram2D(DD_TRACE_E_QDC, energyBins, energyBins,"TraceE vs QDC");

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
	int  mult_pspmt = event.GetSummary("pspmt",true)->GetMult();
	int  mult_mwpc  = event.GetSummary("mcp",true)->GetMult();

	bool has_mwpc = false;
	bool has_pspmt = false;
	bool is_ion = false;
	bool is_decay = false;
	bool incomplete_ion = false; 
	bool incomplete_decay = false; 

	if(mult_mwpc > 0) has_mwpc = true; 
	if(mult_pspmt == 5) has_pspmt = true; 
	
	if(has_pspmt) {
		if(has_mwpc) is_ion = true; 
		else is_decay = true; 
	} else {
		if(has_mwpc) incomplete_ion = true; 
		else incomplete_decay = true; 
	}

	// see accumulation of different signals
	// id=1901
	if(is_ion) plot(1, 2); 
	if(is_decay) plot(1, 4); 
	if(incomplete_ion) plot(1, 6); 
	if(incomplete_decay) plot(1, 8); 


    data_.Clear();
    
    EndProcess();
	return(true);
}

// for correlations
static PixelEvent implantRecorder[600] = {};
static PixelEvent decayRecorder[2][600] = {};
// for position calibrations
const double parX[4] = {1.67e-8, -2.49e-5, 4.11e-2, -4.91}; 
const double parY[4] = {3.16e-8, -4.75e-5, 5.30e-2, -7.20}; 
//
const double parXIon[4] = {1.27e-8, -4.25e-5, 6.94e-2, -3.08e1}; 
const double parYIon[4] = {1.82e-8, -6.15e-5, 9.2e-2, -4.03e1}; 
// linear energy calibration
const double parE[2] = {0.858, 141.675}; 

static unsigned int traceNum = 0; 
fstream outfile; 

bool PspmtProcessor::Process(RawEvent &event){
	if (!EventProcessor::Process(event))
		return false;

	DetectorDriver *driver = DetectorDriver::get();   

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
	
	int  mult_pspmt = event.GetSummary("pspmt",true)->GetMult();
	int  mult_mwpc  = event.GetSummary("mcp",true)->GetMult();
	int  mult_nai   = event.GetSummary("nai",true)->GetMult();
	// by YX

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
  
	static int summedTraceNum = 0; 
	static int anodeTraceNum[4] = {0, 0, 0, 0}; 
	static int dynodeTraceNum = 0; 
	static int originalTraceNum = 0; 

	// by Yongchi Xiao, 01/18/2018
	double qdcCalib = -1; 
  	int p1d;
	double threshold=0;
	
	// by Yongchi Xiao, for PSD purpose, 03/06/2018
	double psd = 0; 

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

	/* by Yongchi Xiao
	 * If mulplicity of PSPMT is equal to 5, which
	 * means all anodes and dynodes fire,
	 * then it is time to sum over all anode traces
	 * to get a complete incoming pulse shape
	 */ 
	// by Yongchi Xiao, 03/07/2018
	Trace traceAnode[4], traceDynode; 
	double onboardEnergy[5] = {0, 0, 0, 0, 0}; 
	if(has_pspmt) { /* all anodes and dynodes fires, exactly equal to 5
					 * BUT this does not mean all five channels has 
					 * proper filtered energy or well-recorded trace
					 */ 
		for (vector<ChanEvent*>::const_iterator it = pspmtEvents.begin();
			 it != pspmtEvents.end(); it++) {        
			ChanEvent *chan = *it; 
			int ch = chan->GetChanID().GetLocation();
			Trace tr = chan->GetTrace(); 
			pspmttime = chan->GetTime(); 
			if(ch < 4) {
				onboardEnergy[ch] = chan->GetCalEnergy(); 
				traceAnode[ch] = tr; 
			}
			if(ch == 4) {
				onboardEnergy[ch] = chan->GetCalEnergy(); 
				traceDynode = tr; 
			}
		}
	}

	// verify all on-board energy
	double all_onboard_energy = 1; 
	for(int i = 0; i < 4; i++) {
		all_onboard_energy *= onboardEnergy[i]; 
	}

	// verify all trace length
	double all_trace_length = 1; 
	for(int i = 0; i < 4; i++) {
		all_trace_length *= traceAnode[i].size(); 
	}
	all_trace_length *= traceDynode.size(); 

	/******************************
	 *         FOR TRACES         *
	 *----------------------------*
	 * Traces are generally for   *
	 * decays only, but in early  *
	 * measurement, we took traces*
	 * for beams as well.         *
	 ******************************/
	bool useTraces = false; 
	PixelEvent pEventTr; 
	vector<PixelEvent> vecPileupEvents; 
	Trace traceSum, traceTop, traceBottom, traceLeft, traceRight; 		
	if(all_trace_length > 0 && has_pspmt) {
		useTraces = true;
		double sumAnode = 0; 
		for(int i = 0; i < 4; i++) { // loop over anode traces and energy
			traceAnode[i].SumTrace(traceSum); 
		}
		if(traceSum.size() == 0) exit(EXIT_FAILURE); 
		if(traceDynode.GetValue("numPulses") > 1) has_pileup = true;  
		// Prepare to deduce position information
		// top
		traceAnode[0].SumTrace(traceTop); 
		traceAnode[1].SumTrace(traceTop); 
		// left
		traceAnode[1].SumTrace(traceLeft); 
		traceAnode[2].SumTrace(traceLeft); 
		// bottom
		traceAnode[2].SumTrace(traceBottom); 
		traceAnode[3].SumTrace(traceBottom); 
		// right
		traceAnode[3].SumTrace(traceRight); 
		traceAnode[0].SumTrace(traceRight); 
		double qdcTop, qdcLeft, qdcBottom, qdcRight; 
		// by calling DoQDCSimple(lo, binNum), the baseline is re-evaluated
		qdcTop = traceTop.DoQDCSimple(1, 120); 
		qdcLeft = traceLeft.DoQDCSimple(1, 120); 
		qdcBottom = traceBottom.DoQDCSimple(1, 120); 
		qdcRight = traceRight.DoQDCSimple(1, 120); 
		double qdcSum = traceSum.DoQDCSimple(1, 120);
		// deduce position info. 
		double rx, ry; 		// position in ratio
		rx = qdcTop/qdcSum; 
		ry = qdcRight/qdcSum; 
		rx *= 1000; 
		ry *= 1000; 
		double posX, posY; // position calibration
		posX = parX[0]*pow(rx,3) + parX[1]*pow(rx, 2) + parX[2]*rx + parX[3] + 0.5; 
		posY = parY[0]*pow(ry,3) + parY[1]*pow(ry, 2) + parY[2]*ry + parY[3] + 0.5; 
		int px, py; // pixlated
		px = trunc(posX); 
		py = trunc(posY);
		int p1d = px + 24*py; // unique 1-d position
		//
		bool badCorner = false; 
		if(px >= 18 && py >= 15) badCorner = true; 
		if(p1d >= 0 && p1d < 576 && !badCorner) {
			// energy and calibration
			qdcSum /= 40.; 
			qdcSum *= pixelCalib[p1d]; 
			//			qdcSum = parE[0]*qdcSum + parE[1]; // 4 keV/ch, don't remove
			qdcSum *= 0.961; // second gain-match, don't remove		
			// create an pixel event
			pEventTr.AssignValues(qdcSum, pspmttime, rx, ry, px, py); // first event
			pEventTr.AssignTraces(traceSum, traceDynode); 
			pEventTr.AssignSignalType(has_mwpc); 
			pEventTr.AssignPileup(has_pileup); 
		} // end:reasonable position
		//
		if(has_pileup) {
			bool similarEnergy = false; 
			double pspmttime2 = traceDynode.GetValue("filterTime2"); 
			double pulse2end = pspmttime2 + 90; 
			if( (traceDynode.GetValue("filterTime2") - traceDynode.GetValue("filterTime")) > 100
				&& !has_mwpc) {
				// check if two signals come with similar energy
				double qdcDynode = traceDynode.DoQDCSimple(traceDynode.GetValue("filterTime"), 
														   traceDynode.GetValue("filterTime")+40); 
				double qdc2Dynode = traceDynode.DoQDCSimple(pspmttime2, pspmttime2+40); 
				if(abs(qdcDynode-qdc2Dynode)/((qdcDynode+qdc2Dynode)/2.)*100 < 40.) {
					similarEnergy = true; // might be successive alpha-particle
				}
				if(similarEnergy) {
				//				if(true) {
					// deduce position info. of 2nd signal
					double qdcTop2, qdcLeft2, qdcBottom2, qdcRight2, qdcSum2; 
					qdcTop2 = traceTop.DoQDCSimple(pspmttime2, pulse2end); 
					qdcLeft2 = traceLeft.DoQDCSimple(pspmttime2, pulse2end); 
					qdcBottom2 = traceBottom.DoQDCSimple(pspmttime2, pulse2end); 
					qdcRight2 = traceRight.DoQDCSimple(pspmttime2, pulse2end); 
					qdcSum2 = traceSum.DoQDCSimple(pspmttime2, pulse2end); // get total energy for 2nd signal
					double rx2, ry2;  
					rx2 = qdcTop2/qdcSum2; 
					ry2 = qdcRight2/qdcSum2; 
					rx2 *= 1000; 
					ry2 *= 1000; 
					double posX2, posY2; 
					posX2 = parX[0]*pow(rx2,3) + parX[1]*pow(rx2, 2) + parX[2]*rx2 + parX[3] + 0.5; 
					posY2 = parY[0]*pow(ry2,3) + parY[1]*pow(ry2, 2) + parY[2]*ry2 + parY[3] + 0.5; 
					int px2, py2; 
					px2 = trunc(posX2); 
					py2 = trunc(posY2);
					int p1d2; 					
					p1d2 = px2 + 24*py2; 
					bool samePixel = false; 
					//
					bool badCorner2 = false; 
					if(px2 >= 18 && py2 >= 15) badCorner2 = true; 
					if(p1d2 > 0 && p1d2 <=576 && !badCorner2) { // position is reasonable
						if(p1d2 == p1d) samePixel = true; 
						qdcSum2 /= 40.;
						qdcSum2 *= pixelCalib[p1d]; 
						//						qdcSum2 = parE[0]*qdcSum2 + parE[1]; // 4 keV/ch
						qdcSum2 *= 0.961; // second gain-match
						pspmttime2 += pspmttime;
						PixelEvent pEventTr2;
						pEventTr2.AssignValues(qdcSum2, pspmttime2, rx2, ry2, px2, py2); 
						if(samePixel) {
							vecPileupEvents.push_back(pEventTr); 
							vecPileupEvents.push_back(pEventTr2); 
						}
					} // end:position is reasonable
				} // end:similarEnergy/true				
			} // end:enough interval between pulses
		} // end:has_pileup
	}

	/******************************
	 *    FOR ONBOARD ENERGY      *
	 *----------------------------*
	 * Generally, onboard energy  *
	 * is for ions only since we  *
	 * did not record traces for  *
	 * ions in real measurement.  *
	 ******************************/
	bool useOnboardEnergy = false; 
	PixelEvent pEventE; 
	if(all_onboard_energy > 0 && has_pspmt) { 
		useOnboardEnergy = true; 
		// deal with on-board energy
		double eSum, eTop, eLeft, eBottom, eRight; 
		eSum = onboardEnergy[0] + onboardEnergy[1] + onboardEnergy[2] + onboardEnergy[3];
		eTop = onboardEnergy[0] + onboardEnergy[1]; 
		eLeft = onboardEnergy[1] + onboardEnergy[2]; 
		eBottom = onboardEnergy[2] + onboardEnergy[3]; 
		eRight = onboardEnergy[3] + onboardEnergy[0]; 
		// intial process of onboard energy
		double orx = eTop/eSum; 
		double ory = eRight/eSum; 
		orx *= 1250; 
		ory *= 1250; 
		orx += 500; 
		ory += 500; 
		double posXIon, posYIon; 
		posXIon = parXIon[0]*pow(orx, 3) + parXIon[1]*pow(orx, 2) 
			+ parXIon[2]*pow(orx, 1) + parXIon[3] + 0.5; 
		posYIon = parYIon[0]*pow(ory, 3) + parYIon[1]*pow(ory, 2) 
			+ parYIon[2]*pow(ory, 1) + parYIon[3] + 0.5; 
		int pxIon, pyIon; 
		pxIon = trunc(posXIon); 
		pyIon = trunc(posYIon); 
		int p1dIon = pxIon + 24*pyIon; 
		//
		bool badCornerIon = false; 
		if(pxIon >= 18 && pyIon >= 15) badCornerIon = true; 
		if(p1dIon >= 0 && p1dIon < 576 && !badCornerIon) {
			// create a pixel event
			pEventE.AssignValues(eSum, pspmttime, orx, ory, pxIon, pyIon); 
			pEventE.AssignSignalType(has_mwpc); 
			//			pEventE.AssignPileup(has_pileup); 
		}
	}

	/************************
	 *      IONS ONLY       *
	 ************************/
	if(has_implant && !has_veto) {
		PixelEvent pe; 
		//		if(useTraces) {pe = pEventTr;}
		if(useOnboardEnergy) {pe = pEventE; }
		if(pe.Is_Filled() && pe.Is_Ion()) {
			double orx, ory; 
			orx = pe.GetRawX(); 
			ory = pe.GetRawY(); 
			plot(3, orx, ory); // 1903
			int pxIon, pyIon; 
			pxIon = pe.GetX(); 
			pyIon = pe.GetY(); 
			plot(2, pxIon, pyIon); // 1902
			double energyIon; 
			energyIon = pe.GetEnergy(); 
			int p1dIon = pxIon + pyIon*24; 
			plot(65, energyIon, p1dIon); // 1965
			// Ion-decay correlations
			implantRecorder[p1dIon] = pe; 
			decayRecorder[0][p1dIon].Clear(); // clear recorded decays
			decayRecorder[1][p1dIon].Clear(); 
			
			bool plotTraces = false; // artificial switch
			if(useTraces && plotTraces && !has_pileup) {
				// plot all associated traces
				for(vector<int>::iterator ittr = traceSum.begin(); ittr != traceSum.end(); ittr++) {
					plot(70, ittr-traceSum.begin(), traceNum, *ittr); // 1970
				}
				// summed anode traces for positioning
				for(vector<int>::iterator it = traceTop.begin();
					it != traceTop.end(); it++) {
					plot(71, it-traceTop.begin(), traceNum, (*it)); // 1971, top
				}
				for(vector<int>::iterator it = traceLeft.begin();
					it != traceLeft.end(); it++) {
					plot(72, it-traceLeft.begin(), traceNum, (*it)); // 1972, left
				}
				for(vector<int>::iterator it = traceBottom.begin();
					it != traceBottom.end(); it++) {
					plot(73, it-traceBottom.begin(), traceNum, (*it)); // 1973, bottom
				}
				for(vector<int>::iterator it = traceRight.begin();
					it != traceRight.end(); it++) {
					plot(74, it-traceRight.begin(), traceNum, (*it)); // 1974, right
				}
				// plot dynode trace
				for(vector<int>::iterator ittr = traceDynode.begin(); ittr != traceDynode.end(); ittr++) {
					plot(75, ittr-traceDynode.begin(), traceNum, *ittr); // 1975
				}
				// Original Traces
				for(vector<int>::iterator ittr = traceAnode[0].begin(); ittr != traceAnode[0].end(); ittr++) {
					plot(76, ittr-traceAnode[0].begin(), traceNum, *ittr); // 1976
				}
				for(vector<int>::iterator ittr = traceAnode[1].begin(); ittr != traceAnode[1].end(); ittr++) {
					plot(77, ittr-traceAnode[1].begin(), traceNum, *ittr); // 1977
				}
				for(vector<int>::iterator ittr = traceAnode[2].begin(); ittr != traceAnode[2].end(); ittr++) {
					plot(78, ittr-traceAnode[2].begin(), traceNum, *ittr); // 1978
				}
				for(vector<int>::iterator ittr = traceAnode[3].begin(); ittr != traceAnode[3].end(); ittr++) {
					plot(79, ittr-traceAnode[3].begin(), traceNum, *ittr); // 1979
				}
				traceNum++; 
			} // end:plot_traces

		}
	}
	/*************************
	 *     DECAYS ONLY       *
	 *************************/
	if(has_decay && !has_veto && !has_pileup) {
		PixelEvent pe; 
		if(useTraces) {pe = pEventTr; }
		if(pe.Is_Filled() && !pe.Is_Ion()) {
			double rx, ry; 
			rx = pe.GetRawX(); 
			ry = pe.GetRawY(); 
			int px, py; 
			px = pe.GetX(); 
			py = pe.GetY(); 
			plot(5, rx, ry); // 1905
			plot(4, px, py); // 1904
			double qdcSum = pe.GetEnergy();
			int p1d = pe.GetX() + pe.GetY()*24; 
			plot(63, qdcSum, p1d); // 1963
			// Ion-decay correlations
			if(!decayRecorder[0][p1d].Is_Filled()) {
				decayRecorder[0][p1d] = pe;
				double timeDecay = decayRecorder[0][p1d].GetTime(); 
				if(implantRecorder[p1d].Is_Filled()) {
					double Dt = (timeDecay - implantRecorder[p1d].GetTime())*Globals::get()->clockInSeconds(); 
					Dt *= 2; 
					plot(45, qdcSum, Dt*1.e6); // 1945
					plot(46, qdcSum, Dt*1.e5); // 1946
					plot(47, qdcSum, Dt*1.e4); // 1947
					plot(48, qdcSum, Dt*1.e3); // 1948
					plot(49, qdcSum, Dt*1.e2); // 1949
					// plot correlated ions
					if(true) {
						//					if(implantRecorder[p1d].Has_Trace()) {
						//					if(abs(qdcSum - 410.) < 55.16) {// 109I-protons
						//						if(abs(qdcSum - 711) < 73.59) {// 109Te-alphas
						plot(40, implantRecorder[p1d].GetEnergy(), Dt*1.e6); // 1940
						plot(41, implantRecorder[p1d].GetEnergy(), Dt*1.e5); // 1941
						plot(42, implantRecorder[p1d].GetEnergy(), Dt*1.e4); // 1942
						plot(43, implantRecorder[p1d].GetEnergy(), Dt*1.e3); // 1943
						plot(44, implantRecorder[p1d].GetEnergy(), Dt*1.e2); // 1944
						// ion E. vs. decay E. 
						if(Dt < 256e-6) plot(50, qdcSum, implantRecorder[p1d].GetEnergy()); // 1950
						if(Dt < 256e-5) plot(51, qdcSum, implantRecorder[p1d].GetEnergy()); // 1951
						if(Dt < 256e-4) plot(52, qdcSum, implantRecorder[p1d].GetEnergy()); // 1952
						if(Dt < 256e-3) plot(53, qdcSum, implantRecorder[p1d].GetEnergy()); // 1953
						if(Dt < 256e-2) plot(54, qdcSum, implantRecorder[p1d].GetEnergy()); // 1954
					}// previous ion has traces
				}
				bool plotTraces = false; 
				if(has_pileup && vecPileupEvents.size() == 2) {
					// see condition for filling vecPileupEvents above
					plotTraces = true; 
				}
				plotTraces = false; // artificial switch
				if(plotTraces) {
					// plot all associated traces
					for(vector<int>::iterator ittr = traceSum.begin(); ittr != traceSum.end(); ittr++) {
						plot(70, ittr-traceSum.begin(), traceNum, *ittr); // 1970
					}
					// summed anode traces for positioning
					for(vector<int>::iterator it = traceTop.begin();
						it != traceTop.end(); it++) {
						plot(71, it-traceTop.begin(), traceNum, (*it)); // 1971, top
					}
					for(vector<int>::iterator it = traceLeft.begin();
						it != traceLeft.end(); it++) {
						plot(72, it-traceLeft.begin(), traceNum, (*it)); // 1972, left
					}
					for(vector<int>::iterator it = traceBottom.begin();
						it != traceBottom.end(); it++) {
						plot(73, it-traceBottom.begin(), traceNum, (*it)); // 1973, bottom
					}
					for(vector<int>::iterator it = traceRight.begin();
						it != traceRight.end(); it++) {
						plot(74, it-traceRight.begin(), traceNum, (*it)); // 1974, right
					}
					// plot dynode trace
					for(vector<int>::iterator ittr = traceDynode.begin(); ittr != traceDynode.end(); ittr++) {
						plot(75, ittr-traceDynode.begin(), traceNum, *ittr); // 1975
					}
					// Original Traces
					for(vector<int>::iterator ittr = traceAnode[0].begin(); ittr != traceAnode[0].end(); ittr++) {
						plot(76, ittr-traceAnode[0].begin(), traceNum, *ittr); // 1976
					}
					for(vector<int>::iterator ittr = traceAnode[1].begin(); ittr != traceAnode[1].end(); ittr++) {
						plot(77, ittr-traceAnode[1].begin(), traceNum, *ittr); // 1977
					}
					for(vector<int>::iterator ittr = traceAnode[2].begin(); ittr != traceAnode[2].end(); ittr++) {
						plot(78, ittr-traceAnode[2].begin(), traceNum, *ittr); // 1978
					}
					for(vector<int>::iterator ittr = traceAnode[3].begin(); ittr != traceAnode[3].end(); ittr++) {
						plot(79, ittr-traceAnode[3].begin(), traceNum, *ittr); // 1979
					}
					traceNum++; 
					double qdcSum2 = vecPileupEvents.at(1).GetEnergy(); 
					plot(64, qdcSum, qdcSum2); // 1964
				} // end:plot_traces

			} // end:didn't fill the first gen. of decay
		} // pe is got
	} // end:!has_mwpc and !has_veto
	
	vecPileupEvents.clear(); 
		
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

/** \file NaIProcessor.cpp
 * NaI class takes the signals from the mcp detector
 *   and calculates a 2D position based on the readout
 * SNL 2-2-08, Modified DTM 9-09
 */
#include "DammPlotIds.hpp"
#include "NaIProcessor.hpp"
#include "RawEvent.hpp"

using std::string;
using std::vector;
using namespace dammIds::nai;

/* by Yongchi Xiao
 * a notice for PspmtProcessor on 511 keV gamma-ray
 */
bool has511keV, has511keVBarrel, has511keVPlug;
const double fwhm[9] = {0, 0, 0, 0, 132.21, 187.28, 142.00, 146.5, 52.55}; 

namespace dammIds {
  namespace nai {	
    const int DD_POSENE = 0;
    const int D_ENE     = 1;
  }
}
void NaIProcessor::NaIData::Clear(void)
{
}

NaIProcessor::NaIProcessor(void) : EventProcessor(OFFSET, RANGE, "nai")
{
  associatedTypes.insert("nai");
}

void NaIProcessor::DeclarePlots(void)
{
  
  const int posBins      = S4; //    16
  const int energyBins   = 1024; // 16384
  //  DeclareHistogram2D(DD_POSENE, posBins, energyBins, "NaI Position and Energy"); // 2030
  DeclareHistogram1D(DD_POSENE+1, energyBins, "NaI 1");
  DeclareHistogram1D(DD_POSENE+2, energyBins, "NaI 2");
  DeclareHistogram1D(DD_POSENE+3, energyBins, "NaI 3");
  DeclareHistogram1D(DD_POSENE+4, energyBins, "NaI 4");
  DeclareHistogram1D(DD_POSENE+5, energyBins, "NaI 5");
  DeclareHistogram1D(DD_POSENE+6, energyBins, "NaI 6");
  DeclareHistogram1D(DD_POSENE+7, energyBins, "NaI 7");
  DeclareHistogram1D(DD_POSENE+8, energyBins, "NaI 8");
  /* by YX
   * generate a summed signal for plug-in part of NaI
   */
  DeclareHistogram1D(DD_POSENE+9, energyBins, "Plug-in signal"); // 2039
}

bool NaIProcessor::PreProcess(RawEvent &event) {
	if (!EventProcessor::Process(event)) 
		return false;

	has511keV = false;
	has511keVBarrel = false;
	has511keVPlug = false;
	
	static const vector<ChanEvent*> &naiEvents = sumMap["nai"]->GetList();
	data.Clear();
	
	double sumPlugEnergy = 0; 
	for (vector<ChanEvent*>::const_iterator it = naiEvents.begin();
		 it != naiEvents.end(); it++) {
		ChanEvent *chan = *it;
      
		string subtype   = chan->GetChanID().GetSubtype();
		int number       = chan->GetChanID().GetLocation();
		double calEnergy = chan->GetCalEnergy();
		double naiTime   = chan->GetTime();
		using std::cout;
		using std::endl; 
		
		// check 511 keV gamma
		if(abs(calEnergy - 511) < fwhm[number]/2.) {
			has511keV = true; 
			has511keVBarrel = true; 
		}

		if(number < 4 ) {
			sumPlugEnergy += calEnergy; 
		} 
	}

	if(sumPlugEnergy > 0) {
		sumPlugEnergy *= 0.2687; // calibrated to 1 keV/ch
		plot(DD_POSENE+9, sumPlugEnergy); // 2039
		// check 511 keV gamma
		if(abs(sumPlugEnergy - 511) < fwhm[8]/2.) {
			has511keV = true;
			has511keVPlug = true; 
		}
	}


	EndProcess();
	return(true);
	
}


bool NaIProcessor::Process(RawEvent &event)
{
  if (!EventProcessor::Process(event))
    return false;

  static const vector<ChanEvent*> &naiEvents = sumMap["nai"]->GetList();

  data.Clear();
  for (vector<ChanEvent*>::const_iterator it = naiEvents.begin();
       it != naiEvents.end(); it++) {
      ChanEvent *chan = *it;

      //      std::cout << "========NaIProcessor " << std::endl;
      
      string subtype   = chan->GetChanID().GetSubtype();
      int number       = chan->GetChanID().GetLocation();
      double calEnergy = chan->GetCalEnergy();
      double naiTime   = chan->GetTime();
      using std::cout;
      using std:: endl; 
      
      if(number<13){
	plot(DD_POSENE,number,calEnergy);
      }
      if(number==0){
	plot(DD_POSENE+1,calEnergy);
      }
      if(number==1){
	plot(DD_POSENE+2,calEnergy);
      }
      if(number==2){
	plot(DD_POSENE+3,calEnergy);
      }
      if(number==3){
	plot(DD_POSENE+4,calEnergy);
      }
      if(number==4){
	plot(DD_POSENE+5,calEnergy);
      }
      if(number==5){
	plot(DD_POSENE+6,calEnergy);
      }
      if(number==6){
	plot(DD_POSENE+7,calEnergy);
      }
      if(number==7){
	plot(DD_POSENE+8,calEnergy);
      }
      
      


  }
      
      
      EndProcess();
      return true;
}

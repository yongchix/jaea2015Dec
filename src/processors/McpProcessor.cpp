/** \file McpProcessor.cpp
 * Mcp class takes the signals from the mcp detector
 *   and calculates a 2D position based on the readout
 * SNL 2-2-08, Modified DTM 9-09
 */
#ifdef useroot
#include <TTree.h>
#endif

#include "DammPlotIds.hpp"
#include "McpProcessor.hpp"
#include "RawEvent.hpp"

using std::string;
using std::vector;
using namespace dammIds::mcp;
using namespace std;

namespace dammIds {
  namespace mcp {	
    const int D_CATHODE   = 1;
    const int D_TAC       = 2;
    const int DD_POSXY    = 3;
    const int D_LEFT = 4;
    const int D_RIGHT = 5;
    const int DD_XE_YE = 6;
  }
}


void McpProcessor::McpData::Clear(void)
{
  for (size_t i=0; i < nPos; i++)
    raw[i] = 0;
  xpos = ypos = 0.;

  mult = 0;
}

McpProcessor::McpProcessor(void) : EventProcessor(OFFSET, RANGE, "mcp")
{
  associatedTypes.insert("mcp");
}

void McpProcessor::DeclarePlots(void)
{
  
  const int posBins   = SE; //
  const int posBins2  = SE;
  const int posBins2D = S9; // 512
  const int posBins3  = 2048;
  const int posBins4  = 512;
  DeclareHistogram1D(D_CATHODE, posBins, "Cathode");  // 2001
  DeclareHistogram1D(D_TAC, posBins2, "TAC");// 2002
  DeclareHistogram2D(DD_POSXY, posBins4, posBins3, "MCP Left 2D"); //2003


}


bool McpProcessor::Process(RawEvent &event)
{
  
  if (!EventProcessor::Process(event))
    return false;

  double qTotal, qRight, qTop, Cathode;
  double tac;
  static const vector<ChanEvent*> &mcpEvents = sumMap["mcp"]->GetList();
  
  // static const vector<ChanEvent*> &pspmtEvents = sumMap["pspmt"]->GetList();
  
  //int mult_pspmt = event.GetSummary("pspmt",true)->GetMult();
  //cout << "pspmt "<< mult_pspmt << endl;
  data.Clear();
  
  for (vector<ChanEvent*>::const_iterator it = mcpEvents.begin();
       it != mcpEvents.end(); it++) {
      ChanEvent *chan = *it;

      string subtype   = chan->GetChanID().GetSubtype();
      double calEnergy = chan->GetCalEnergy();
      double mcpTime   = chan->GetTime();
      

   
      if(subtype == "1time") {
	tac=calEnergy;
	
	//	cout << "tac " << tac << endl;  
      } else if (subtype == "1position1"){
	Cathode = calEnergy;
	data.raw[0] = calEnergy;
	data.time[0] = mcpTime;
	data.mult++;
      } else if (subtype == "1position2"){
	data.raw[1] = calEnergy;
	data.mult++;
	data.time[1] = mcpTime;
      } else if (subtype == "1position3"){
	data.raw[2] = calEnergy;
	data.mult++;
	data.time[2] = mcpTime;
      } else if (subtype == "1position4"){
	data.raw[3] = calEnergy;
	data.time[3] = mcpTime;
	data.mult++;
      }
  }
    
  // calculation of position from charge collected on the four corners
  // magic numbers here
  data.raw[0] *= 1.3;
    
  qTotal = data.raw[0] + data.raw[1] + data.raw[2] + data.raw[3];
  qRight = data.raw[3] + data.raw[0];
  // qLeft   = data.raw[2] + data.raw[1];  
  qTop   = data.raw[0] + data.raw[1];
  // qBottom = data.raw[2] + data.raw[3];
    
  data.xpos = (qRight / qTotal) * 512. - 75; //horizontal MCP pos
  data.ypos = (qTop   / qTotal) * 512. - 75; //vertical MCP pos
  // qLeft, qBottom not used
  
  if (data.mult >= 2) {
    using namespace dammIds::mcp;

double tDiff=  data.time[1] - data.time[0]+20;
//cout << tDiff << endl;
    plot(D_CATHODE, Cathode);    
    plot(D_TAC, tac);
    data.raw[2]=data.raw[2]/5.;
    plot(DD_POSXY,tac/10,Cathode);
    
  }



  
  
  EndProcess();
  return (data.mult == 3);
}

#ifdef useroot
bool McpProcessor::AddBranch(TTree *tree)
{
  if (tree) {
    TBranch *mcpBranch = 
      tree->Branch(name.c_str(), &data, "raw[4]/D:xpos:ypos:mult/I");

    return (mcpBranch != NULL);
  } 

  return false;
}

void McpProcessor::FillBranch(void)
{
  if (!HasEvent())
    data.Clear();  
}
#endif //useroot

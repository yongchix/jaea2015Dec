/** \file SHECorrelator.hpp
 *
 */

#ifndef __PSPMTCORRELATOR_HPP_
#define __PSPMTCORRELATOR_HPP_

#include <vector>
#include <deque> 
#include <sstream>

#include "Plots.hpp"
#include "Globals.hpp"
#include "ChanEvent.hpp"
#include "Messenger.hpp"
#include "WalkCorrector.hpp"
#include "Calibrator.hpp"
#include "DammPlotIds.hpp"

class PspmtEvent {
public:
  PspmtEvent();
  ~PspmtEvent(){}
  double get_energy();
  double get_time();
  void set_energy(double energy);
  void set_time(double time);
private:
  int x;
  int y;
  double energy_;
  double time_;
};



#endif

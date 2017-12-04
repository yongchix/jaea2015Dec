/*! \file SheCorrelator.cpp
 *
 * The SheCorrelator is designed to recontruct dssd events in a SHE experiment
 * and to correlate chains of alphas in dssd pixels
 */

#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>

#include "DammPlotIds.hpp"
#include "PspmtProcessor.hpp"
#include "PspmtCorrelator.hpp"
#include "DetectorDriver.hpp"
#include "Exceptions.hpp"
#include "Notebook.hpp"

using namespace std;


PspmtEvent::PspmtEvent(){
  energy_ = -1;
  time_   = -1;
}

double PspmtEvent::get_energy(){
  return energy_;
}
double PspmtEvent::get_time(){
  return time_;
}
void PspmtEvent::set_energy(double energy){
  energy_ = energy;
}
void PspmtEvent::set_time(double time){
  time_ = time;
}
  

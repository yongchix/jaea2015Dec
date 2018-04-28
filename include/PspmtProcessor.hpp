/** \file PspmtProcessor.hpp
 *  \brief A processor to handle pixelated PMTs
 *  \author Shintaro Go
 *  \date November 16, 2015
 */

#ifndef __PSPMTPROCESSOR_HPP__
#define __PSPMTPROCESSOR_HPP__

#include "RawEvent.hpp"
#include "EventProcessor.hpp"
#include "Globals.hpp"
#include "Trace.hpp"
#include <vector>

using namespace std; 
  
class PspmtEvent{
public:
  PspmtEvent();
  PspmtEvent(double ene_pspmt,double ene_trace,double time_pspmt, double time_mwpc,bool is_beam, bool is_decay, bool is_veto);
  ~PspmtEvent(){}
  void set_ene_pspmt(double ene_pspmt){
    ene_pspmt_ = ene_pspmt;
  }
  void set_ene_trace(double ene_trace){
    ene_trace_ = ene_trace;
  }
  void set_time_pspmt(double time_pspmt){
    time_pspmt_ = time_pspmt;
  }
  void set_time_mwpc(double time_mwpc){
    time_mwpc_ = time_mwpc;
  }
  void set_is_beam(bool is_beam){
    is_beam_   = is_beam;
  }
  void set_is_decay(bool is_decay){
    is_decay_  = is_decay;
  }
  void set_is_veto(bool is_veto){
    is_veto_   = is_veto;
  }
  double get_ene_pspmt(){
    return ene_pspmt_;
  }
  double get_ene_trace(){
    return ene_trace_;
  }
  double get_time_pspmt(){
    return time_pspmt_;
  }
  double get_time_mwpc(){
    return time_mwpc_;
  }
  bool get_is_beam() const{
    return is_beam_;
  }
  bool get_is_decay() const{
    return is_decay_;
  }
  bool get_is_veto() const{
    return is_veto_;
  }
  
private:
  double ene_pspmt_;
  double ene_trace_;
  double time_pspmt_;
  double time_mwpc_;
  bool is_beam_;
  bool is_decay_;
  bool is_veto_;
};

class PspmtProcessor : public EventProcessor {
public:
  PspmtProcessor();
  ~PspmtProcessor();
  bool AddEvent(PspmtEvent& event,int x,int y);
  bool FlushChain(int x,int y);
  virtual void DeclarePlots(void);
  virtual bool PreProcess(RawEvent &event);
  virtual bool Process(RawEvent &event);
private:
  std::deque<PspmtEvent>** pixels_;
 
  struct PspmtData {
    void Clear(void);
  } data_;
};

// by YX
class PixelEvent {
private: 
	double energy;
	double time;
	double rx, ry; 
	int x, y;  
	//
	Trace sumAnodeTrace, dynodeTrace; 
	bool is_ion; 
	bool is_pileup; 
public:
	PixelEvent() {
		energy = -1;
		time = -1;
		rx = -1; ry = -1; 
		x = -1; y = -1; 
		is_ion = false; 
		is_pileup = false; 
	}
	~PixelEvent() {
		energy = -1;
		time = -1;
		rx = -1; ry = -1; 
		x = -1; y = -1; 
		is_ion = false; 
		is_pileup = false; 

	}
	void Clear() {
		energy = -1;
		time = -1;
		rx = -1; ry = -1; 
		x = -1; y = -1; 
		is_ion = false; 
		is_pileup = false; 
	}
	//
	bool Is_Filled() {return (time > 0);}
	bool Has_Trace() {return (dynodeTrace.size() > 0); }
	bool Is_Pileup() {return is_pileup; }
	bool Is_Ion() {return is_ion; }
	//
	void AssignValues(double e, double t, double rx_, double ry_, int x_, int y_) {
		energy = e; 
		time = t; 
		rx = rx_; ry = ry_; 
		x = x_; y = y_; 
	}
	void AssignTraces(Trace &trAnode, Trace &trDynode) {
		sumAnodeTrace = trAnode; 
		dynodeTrace = trDynode; 
	}
	void AssignSignalType(bool &type_) {is_ion = type_; }
	void AssignPileup(bool &pileup_) {is_pileup = pileup_; }
	//
	Trace GetSummedTrace() {return sumAnodeTrace; }
	Trace GetDynodeTrace() {return dynodeTrace; }
	double GetTime() {return time; }
	double GetEnergy() {return energy; }
	double GetRawX() {return rx; }
	double GetRawY() {return ry; }
	int GetX() {return x;}
	int GetY() {return y;}
};

#endif // __PSPMTPROCESSOR_HPP__

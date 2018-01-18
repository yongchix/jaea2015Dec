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

struct PixelEvent {
	double energy;
	double time;
	PixelEvent() {
		energy = -1;
		time = -1;
	}
	~PixelEvent() {
		energy = -1;
		time = -1;
	}
	void Clear() {
		energy = -1;
		time = -1;
	}
	bool Is_Filled() {
		return (time > 0); 
	}
};

#endif // __PSPMTPROCESSOR_HPP__

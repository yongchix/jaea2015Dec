/** \file NaIProcessor.hpp
 * 
 * 
 */

#ifndef __NAIPROCESSOR_HPP_
#define __NAIPROCESSOR_HPP_

#include "EventProcessor.hpp"

extern bool has511keV, has511keVBarrel, has511keVPlug; 

class NaIProcessor : public EventProcessor
{
public:
	NaIProcessor();
	virtual void DeclarePlots(void);
	virtual bool PreProcess(RawEvent &rEvent); 
	virtual bool Process(RawEvent &rEvent);
private:
	struct NaIData {
		void Clear(void);
	} data;
};

#endif // __NAIPROCESSOR_HPP_

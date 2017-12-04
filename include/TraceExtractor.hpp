/** \file TraceExtractor.hpp
 *  \brief Header file for the TraceExtractor class
 *
 *  \author David Miller
 *  \date January 2011
 */
#ifndef __TRACEEXTRACTOR_HPP_
#define __TRACEEXTRACTOR_HPP_
#include <string>

#include "TraceAnalyzer.hpp"

#include "DammPlotIds.hpp"

class Trace;

//! \brief A class to extract traces from events
class TraceExtractor : public TraceAnalyzer {
public:
    /** Default Constructor */
    TraceExtractor(){};

    /** Constructor taking the type and subtype to plot
    * \param [in] aType : a type to plot the traces for
    * \param [in] aSubtype : a subtype to plot the traces for */
    TraceExtractor(const std::string &aType, const std::string &aSubtype);

    /** Default Destructor */
    ~TraceExtractor() {};

    /** Declare the plots for the analyzer */
    virtual void DeclarePlots(void);

    /** The main analysis driver
    * \param [in] trace : the trace to analyze
    * \param [in] aType : the type being analyze
    * \param [in] aSubtype : the subtype begin analyzed */
    virtual void Analyze(Trace &trace, const std::string &aType,
                         const std::string &aSubtype);
protected:
    static const int traceBins; //!< The number of bins for the trace length
    static const int numTraces; //!< The number of traces to analyze

    std::string type; //!< the detector type
    std::string subtype; //!< The detector subtype
};
#endif // __TRACEEXTRACTOR_HPP_

/* File TriggerFrequency.hpp
 *
 * Created       : Mon Sep  4 16:43:39 CDT 2006
 * Author        : ddboline
 *
 * Purpose       : 
 *
 * Last modified : 
 * Comments      : 
 */



#ifndef TriggerFrequency_HPP_
#define TriggerFrequency_HPP_


#include "cafe/Processor.hpp"
#include "cafe/Event.hpp"

using namespace cafe;

    class TriggerFrequency : public cafe::Processor
{
	/**
    This gives the number of times each trigger fired, nothing else really
         */
    public:
        TriggerFrequency(const char *name);
        ~TriggerFrequency();

        bool processEvent(cafe::Event& event);
        void finish();

        ClassDef(TriggerFrequency, 0);
    private:
        std::map< std::string , int > _trigger_frequencies;
        std::map< std::pair< std::string , std::string > , int > _trigger_correlations;
};

#endif


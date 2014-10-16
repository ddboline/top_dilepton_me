/* File TriggerFrequency.cpp
 *
 * Created       : Mon Sep  4 16:43:39 CDT 2006
 * Author        : ddboline
 *
 * Purpose       : Template to build a new processor in package "top_cafe"
 *
 * Last modified :
 * Comments      :
*/

#include <stdexcept>

#include "cafe/Processor.hpp"
#include "cafe/Collection.hpp"
#include "cafe/Config.hpp"

#include "top_dilepton_me/TriggerFrequency.hpp"

using namespace std;

TriggerFrequency::TriggerFrequency(const char *name) 
    : cafe::Processor(name)
{
    Config config(name);
    vector<string> temp = config.getVString("Triggers", " ,");
    for( vector<string>::iterator st_it = temp.begin() ; st_it != temp.end() ; ++st_it )
    {
        _trigger_frequencies[*st_it] = 0;
    }

    out() << "Trigger[" << name << "] Triggers = ";
    for( std::map< std::string , int >::iterator it = _trigger_frequencies.begin() ; it != _trigger_frequencies.end() ; ++it )
    {
        out() << it->first << " " ;
        for( std::map< std::string , int >::iterator it_2 = _trigger_frequencies.begin() ; it_2 != _trigger_frequencies.end() ; ++it_2 )
        {
            if( it->first == it_2->first ) 
                continue;
            pair< string , string > temp( it->first , it_2->first );
            _trigger_correlations[temp] = 0;
        }
    }
    out() << endl;
}

    TriggerFrequency::~TriggerFrequency()
{
}
 
    bool TriggerFrequency::processEvent(cafe::Event& event)
{
    using namespace std;

    Collection<TMBTrigger> triggers = event.getTriggers();
    vector<string> trigs;

    //
	// This is proportional to the number of triggers in the event
	// and about log(numUserTriggers).
    //
    for(Collection<TMBTrigger>::const_iterator it = triggers.begin(); it != triggers.end(); ++it) 
    {
        if( _trigger_frequencies.find( (*it).getTrgName() ) != _trigger_frequencies.end() )
        {
            _trigger_frequencies[(*it).getTrgName()] += 1;
            trigs.push_back( (*it).getTrgName() );
        }
    }

    for( vector<string>::iterator it = trigs.begin() ; it != trigs.end() ; ++it )
    {
        for( vector<string>::iterator it_2 = trigs.begin() ; it_2 != trigs.end() ; ++it_2 )
        {
            pair< string , string > temp( *it , *it_2 );
            _trigger_correlations[temp] += 1;
        }
    }
    return true;
}

    void TriggerFrequency::finish() 
{ 
    printf("============================================\n");
    printf("  TriggerFrequency[ %s ] SUMMARY\n\n",name().c_str());

    out() << " there are " << _trigger_frequencies.size() << " seperate triggers in this list " << endl;

    for( std::map< std::string , int >::iterator it = _trigger_frequencies.begin() ; it != _trigger_frequencies.end() ; ++it )
    {
        if( it->second > 0 )
            out() << it->first << " " << it->second << endl;
    }

    for( map< pair< string, string> , int >::iterator it = _trigger_correlations.begin() ; it != _trigger_correlations.end() ; ++it )
    {
        if( it->second > 0 )
            out() << it->first.first << " " << it->first.second << " " << it->second << endl;
    }

    printf("\n============================================\n");
}

ClassImp(TriggerFrequency); ;

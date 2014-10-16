#ifndef top_cafe_Probability_DILEP_HPP_
#define top_cafe_Probability_DILEP_HPP_

//Inherit probProcessor virtual abstracts
#include "caf_trigger/probProcessor.hpp"
#include "cafe/Stat.hpp"
#include "cafe/StatSample.hpp"

struct TrigStruct {
    std::string version ;
    std::string trigger_name  ; // L1 or L2 or L3 name 
    std::string object ; // muon, electron, l1atxx, ...
    std::string level ;
    const cafe::StatSample* sample ;
    std::string name ; // full name (including term name, version, level, ...) 

// important to have L1, L2, L3 in name,  see objectProbability.cpp
    TrigStruct(std::string v, std::string t, std::string o, std::string l, std::string type="", const cafe::StatSample* s=0)
    : version(v), trigger_name(t), object(o), level(l), sample(s), name("Dilepton Version: " + v + " Trigger: " + t + " Type " + type + " Object: " + o + " Level: L" + l) 
    {
//                 if (sample) name = sample->name() + " " + name ;
    }
    bool operator< (const TrigStruct& trig) const {return name < trig.name ;}
} ;

struct Eff {
    double eff ;
    double err ;
    double n ;
    double x2 ;
    Eff() : eff(0) , err(0), n(0), x2(0) {}
    Eff& operator += (double value) {
        n++ ;
        eff += (value - eff)/n ;
        x2 += value*value ;
        if (n>1) {
            err = (x2 - n*eff*eff)/(n-1) ;
            err = err > 0 ? sqrt(err) : 0 ;
        }
    //    std::cout << "A " << value << " " << n << " " << x2 << " " << eff << " " << err 
    //	      << " " << x2 - n*eff*eff<< std::endl ;
        return *this ;
    }
    double uncert() const {return n>0 ? err/sqrt(n) : 0 ;}
} ;

class Probability_DILEP : public probProcessor {
    public:
        std::map<TrigStruct, Eff> _triggers ; // average trigger efficiencies for each trigger term
        std::map<std::string,Eff> _versions ; // average efficiency for each trigger term for different versions
        std::map<std::string,Eff> _globtrigs; // average combine efficiencies L1xL2xL3 for the different global trigger versions
//         cafe::Event* _event ;
        std::string _outname ;
        std::set<std::string> _single_triggers;
        std::map<std::string,double> _single_trigger_prescale;
        std::set<std::string> _overlap;
        std::set<std::string> _overlap2;

        int nevents;

        cafe::Stat* STAT ;

    public:
        Probability_DILEP( const char * name );
        ~Probability_DILEP();

        bool processEvent( cafe::Event & event );

          //Define your effInfo objects here
        void defineEffInfo(std::map< std::string, eff_utils::EffInfo > &effInfo);

        ////////////////////////////////////////////////////////////
        //These methods define the probability calculations
        double calcProb(std::string version);

        ////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////

        double _Sigma;

        void finish() ;

    public:
        ClassDef(Probability_DILEP,0);
};

inline
std::string tex_convert(const std::string& init) {

  using namespace std ;
  string name ;
  for (unsigned int i=0; i<init.size(); i++) {
    string ch ;
    if (init.substr(i,1) == "_") ch = "\\_" ;
    else if (init.substr(i,2) == ">=") {
      ch = "$\\geq$" ;
      i++ ;
    }
    else if (init.substr(i,2) == "<=") {
      ch = "$\\leq$" ;
      i++ ;
    }
    else if (init.substr(i,1) == "<") ch = "$<$" ;
    else if (init.substr(i,1) == ">") ch = "$>$" ;
    else ch = init.substr(i,1) ;
    name = name + ch ;
  }
  return name ;
}

#endif // top_cafe_Probability_DILEP_HPP_


#include "top_dilepton_me/Probability_DILEP.hpp"
#include "cafe/Event.hpp"
#include <stdexcept>
#include <fstream>
#include <sstream>

using namespace std;

struct TableLine {
    std::string version ;
    std::string trigger_name  ;
    std::vector<std::string> eff ;
    std::string lumi ;

  TableLine(std::string v, std::string t, std::string l) :
          version(v), trigger_name(t), lumi(l) { }

  std::ostream& out(ostream& os) {
      os << version << " & "
              << tex_convert(trigger_name) << " & " ;
      for (std::vector<std::string>::const_iterator it = eff.begin();
           it != eff.end() ; it++) os << *it << " & " ;
      os << lumi << "\\\\ \\hline" ;
      return os ;
  }
} ;

struct TableLine1 {
    std::string version ;
    std::string trigger_name  ;
    std::vector<std::string> tname ;
    std::vector<std::string> eff ;

  TableLine1(std::string v, std::string t) :
          version(v), trigger_name(t) { }

  std::ostream& out(ostream& os) {
      os << version << " & "
              << tex_convert(trigger_name) << " & " ;
      os << tex_convert(*(tname.begin())) << " & " << *(eff.begin()) << " \\\\ " ;
      for (unsigned int i = 1; i < tname.size() ; i++) 
          os << " & & " << tex_convert(*(tname.begin()+i)) << " & " << *(eff.begin()+i) << " \\\\ " ;
      os << " \\hline" ;

      return os ;
  }
} ;

Probability_DILEP::~Probability_DILEP()
{
}

Probability_DILEP::Probability_DILEP( const char * name )
    : probProcessor(name)
//  , _event( 0 )
{
    nevents = 0;
    STAT = 0;
}

bool Probability_DILEP::processEvent( cafe::Event & event )
{
//     _event = &event ;
    nevents++;

    cafe::StatPointer stat ;
    event.get("StatPointer", stat) ;
    STAT = stat.pointer();
    probProcessor::_saveweights = false;

    return probProcessor::processEvent(event) ;
}

void Probability_DILEP::defineEffInfo( std::map< std::string, eff_utils::EffInfo > & effInfo )
{
    std::vector<const cafe::StatSample*> samples ;
    if (STAT) samples = STAT->get_samples()  ;
    else samples.push_back(new StatSample(name())) ;

    _channel = "TriggerProbability";

  //We grab terms from the CAFe configuration file here
    cafe::Config config(name());


    _outname = config.get("OutputFileName", "table_trigger");
    TString _systematic = config.get( "Systematic" , "nominal" );
    _Sigma = 0.0;
    if( _systematic == "plus" ) _Sigma = 1.0;
    if( _systematic == "minus" ) _Sigma = -1.0;


    vector<string> global_versions = config.getVString("GlobalTriggerVersions", " ,");
    if ( global_versions.size() == 0 )
        throw runtime_error("Probability_EMU: At least one trigger version must be specified!") ;

    for ( vector<string>::const_iterator it = global_versions.begin(); it != global_versions.end() ; it++)
    {
        string version = *it ;
        for (vector<const StatSample*>::const_iterator st = samples.begin(); st != samples.end(); st++)
            _versions[(*st)->name()+ ":" + version] = Eff() ;
        cafe::Config version_config(version);

        vector<string> overlap = version_config.getVString("Overlap", " ");
        vector<string> overlap2 = version_config.getVString( "Overlap2" , " ");
        vector<string> single_triggers = version_config.getVString("Triggers", " ");
        vector<string> triggers;

        for ( vector<string>::const_iterator tt = overlap.begin(); tt != overlap.end() ; tt++)
        {
            _overlap.insert(*tt);
            triggers.push_back(*tt);
        }

        for ( vector<string>::const_iterator tt = overlap2.begin(); tt != overlap2.end() ; tt++)
        {
            _overlap2.insert(*tt);
            triggers.push_back(*tt);
        }

        for ( vector<string>::const_iterator tt = single_triggers.begin(); tt != single_triggers.end() ; tt++) 
        {
            _single_triggers.insert(*tt);
            triggers.push_back(*tt);
            string trigger_name = *tt ;
            for (vector<const StatSample*>::const_iterator st = samples.begin(); st != samples.end(); st++)
            {
                _globtrigs[(*st)->name()+ ":" + version + ":" + trigger_name] = Eff() ;
            }
        }

        for ( vector<string>::const_iterator tt = triggers.begin(); tt != triggers.end() ; tt++)
        {
            string trigger_name = *tt ;
            cafe::Config trigname_config(trigger_name);

            double prescale_fact = trigname_config.get( "Prescale" , 1.0 );
            _single_trigger_prescale[trigger_name] = prescale_fact;

            for (int l = 1; l <= 3 ; l++) 
            {
                string level = "1" ;
                if (l==2) level = "2" ;
                else if (l==3) level = "3" ;

                vector<string> objects = trigname_config.getVString("ObjectsL"+level, " ,");
                for ( vector<string>::const_iterator ot = objects.begin(); ot != objects.end() ; ot++)
                {
                    string object = *ot ;
                    cafe::Config obj_config(object);

                    string effName = obj_config.get("EffName" , "eff");
                    string effType = obj_config.get("EffType" , "Binned");
                    vector<string> effVarNames = obj_config.getVString("EffVarNames" , " ,");
                    string objRelativeTo = obj_config.get("ObjRelativeTo" , "");
                    string objType = obj_config.get("ObjType" , "");
                    string objQuality  = obj_config.get("ObjQuality", "");
                    string trigVersionLow  = obj_config.get("TriggerVersionLow", "");
                    string trigVersionHigh  = obj_config.get("TriggerVersionHigh", "");
                    int runRangeLow  = obj_config.get("RunRangeLow", -1);
                    int runRangeHigh  = obj_config.get("RunRangeHigh", -1);
                    string trigVersion  = obj_config.get("TriggerVersion", "");

                    if (objType.size() ==0) continue ;

                    for (vector<const StatSample*>::const_iterator st = samples.begin(); st != samples.end(); st++)
                    {
                        TrigStruct trig(version, trigger_name, object, level, objType , *st) ;
//                         if( _triggers.find( trig ) != _triggers.end() )
//                             continue;
                        _triggers[trig] = Eff() ;
                        string effInfoName = trig.name ;
                        effInfo[effInfoName.c_str()].EffName(effName);
                        effInfo[effInfoName.c_str()].EffVarNames(effVarNames);
                        effInfo[effInfoName.c_str()].EffType(effType);
                        effInfo[effInfoName.c_str()].ObjType(objType);
                        if (objRelativeTo != "")
                            effInfo[effInfoName.c_str()].ObjRelativeTo(objRelativeTo);
                        effInfo[effInfoName.c_str()].ObjQuality(objQuality);
                        if (trigVersionLow != "")
                            effInfo[effInfoName.c_str()].TriggerVersionLow(trigVersionLow);
                        if (trigVersionHigh != "")
                            effInfo[effInfoName.c_str()].TriggerVersionHigh(trigVersionHigh);
                        if (trigVersion != "")
                            effInfo[effInfoName.c_str()].TriggerVersion(trigVersion);
                        if (runRangeLow >= 0 )
                            effInfo[effInfoName.c_str()].RunRangeLow(runRangeLow);
                        if (runRangeHigh >= 0)
                            effInfo[effInfoName.c_str()].RunRangeHigh(runRangeHigh);
                    }
                }
            }
        }
    }
}

double Probability_DILEP::calcProb( std::string version )
{
    // combine L1xL2xL3 efficiency for this event
    map<string,float> triggs ;

    // efficiency for each trigger term in this event
    map<TrigStruct,float> trig_term ;

    string sname = "";

//     cout << " size " << _triggers.size() << endl;

    for(map<TrigStruct, Eff>::iterator tt = _triggers.begin() ; tt != _triggers.end(); tt++) 
    {
        const TrigStruct* trig = &(tt->first) ;
//         cout << " trigger " << trig->version << " " << version << " " << ( trig->version == version ) << endl;
//         cout << " sample " << trig->sample->tagged(_event) << endl;
        if (trig->version != version) continue ;
//         if (!trig->sample->tagged(_event) ) continue ;

        sname = trig->sample->name() ;

        const objectProbabilities *prob = this->probObject(trig->name) ;

        const EffInfo* info = prob->Info() ;

        float eff = 1.0 ;
//         cout << " obj type " << info->ObjType() << endl;
        if (info->ObjType().find("muon") != string::npos || info->ObjType().find("Muon") != string::npos || info->ObjType().find("MUON") != string::npos)
        {
            if(info->ObjQuality().find("l3trk") != string::npos || info->ObjQuality().find("l1trk") != string::npos || info->ObjQuality().find("ttk") != string::npos || info->ObjQuality().find("TRK") != string::npos || info->ObjQuality().find("pt4") != string::npos || info->ObjQuality().find("TLM") != string::npos )
            {
                for(unsigned int iobj = 0; iobj < MU.size(); ++iobj)
                {
                    const TMBTrack* track = MU[iobj].GetChargedTrack() ;
                    if (!track)
                    {
//                         cout << "Probability_EMU: Muon track pointer 0!" << endl ;
                        continue ;
                    }
                    if( track->getChi2Ndf() < 4.0 )
                    {
                        double stat = _Sigma ;
                        double temp_eff = prob->getProbability(*track, stat);
//                     cout << " temp_eff muon track "<< info->ObjQuality() << " " << track->z() << " " << track->det_etaCFT() << " "  << temp_eff << endl;
                        eff *= (1.0-temp_eff) ;
                    }
                }
                for(unsigned int iobj = 0; iobj < EM.size(); ++iobj)
                {
//                     const TMBTrack* track = EM[iobj].getPtrChp() ;
                    const TMBTrack* track = EM[iobj].GetChargedTrack() ;
                    if (!track) 
                    {
//                         cout << "Probability_EMU: Electron track pointer 0!" << endl ;
                        continue ;
                    }
                    if( track->getChi2Ndf() < 4.0 )
                    {
                        double stat = _Sigma ;
                        double temp_eff = prob->getProbability(*track, stat);
//                     cout << " temp_eff electron track " << info->ObjQuality() << " " << track->z() << " " << track->det_etaCFT() << " "  << temp_eff << endl;
                        eff *= (1.0-temp_eff) ;
                    }
                }
                for(unsigned int iobj = 0 ; iobj < TRK.size() ; iobj++ )
                {
                    if( TRK[iobj].getChi2Ndf() < 4.0 )
                    {
                        double stat = _Sigma ;
                        double temp_eff = prob->getProbability(TRK[iobj], stat);
//                         cout << " temp_eff electron track " << info->ObjQuality() << " " << TRK[iobj].z() << " " << TRK[iobj].det_etaCFT() << " "  << temp_eff << endl;
                        eff *= (1.0-temp_eff) ;
                    }
                }
                eff = 1.0 - eff ;
            }
            else
            {
                for(unsigned int iobj = 0; iobj < MU.size(); ++iobj)
                {
                    double stat = _Sigma;
                    double temp_eff = prob->getProbability( MU[iobj] , stat );
//                     cout << " temp_eff muon " << info->ObjQuality() << " " << MU[iobj].Eta() << " " << MU[iobj].Phi() << " " << temp_eff << endl;
                    eff *= ( 1.0 - temp_eff );
                }
                eff = 1.0 - eff;
            }
        }
        else if (info->ObjType().find("Electron") != string::npos || info->ObjType().find("electron") != string::npos || info->ObjType().find("ELECTRON") != string::npos || info->ObjType().find("EM") != string::npos || info->ObjType().find("em") != string::npos) 
        {
            for(unsigned int iobj = 0; iobj < EM.size(); ++iobj)
            {
                double stat = _Sigma;
                double temp_eff = prob->getProbability( EM[iobj] , stat );
//                 cout << " temp_eff electron " << info->ObjQuality() << " " << temp_eff << endl;
                eff *= ( 1.0 - temp_eff );
            }
            eff = 1.0 - eff;
        }
        else if( info->ObjType().find("jet") != string::npos || info->ObjType().find("Jet") != string::npos || info->ObjType().find("JET") != string::npos )
        {
            for(unsigned int iobj = 0; iobj < JET.size(); ++iobj)
            {
                double stat = _Sigma;
                double temp_eff = prob->getProbability( JET[iobj] , stat );
//                 cout << " temp_eff jet " << info->ObjQuality() << " " << temp_eff << endl;
                eff *= ( 1.0 - temp_eff );
            }
            eff = 1.0 - eff;
        }
//         else throw runtime_error("Probability_EMU: Unknown object") ;

    // including sample name (ttbar, ztautau, ..), trigger version and 
    // full trigger name
        string tname = sname + ":" + version + ":" + trig->trigger_name ;
//         cout << " ProbEMU: add to triggs map: " << tname << endl;
        map<string, float>::iterator it = triggs.find(tname) ;
        if (it == triggs.end()) triggs[tname] = eff ;
        else it->second *= eff ; // L1xL2xL3 efficiency

        tt->second += eff ; // average efficiency for each trigger term
        trig_term[tt->first] = eff ; // current efficiency for each trigger term
    }

//     cout << " triggs " << triggs.size() << endl;

  // ORing
    float probability = 0.0, probplus = 0.0 , probneg = 0.0;
    float eff_MUEM1 = 0.0, eff_MUEM2 = 0.0;

  // several sample in the same event
    unsigned int pos = triggs.begin()->first.find(":") ;
    string sample = triggs.begin()->first.substr(0,pos);

  // Sum the efficiency for single triggers, substract overlap efficiency
    for (map<string,float>::const_iterator it = triggs.begin() ; it != triggs.end() ; it++)
    {
        if (sample != it->first.substr(0,it->first.find(":"))) continue ;
        unsigned int pos2 = it->first.rfind(":");
        string name = it->first.substr(pos2+1,string::npos);

        set<string>::const_iterator itplus = _single_triggers.find(name);
        set<string>::const_iterator itneg = _overlap.find(name);
        set<string>::const_iterator itplus2 = _overlap2.find(name);
        double prescale = 1.0;
        if( itplus!=_single_triggers.end() )
            prescale = _single_trigger_prescale[*itplus];
        if( itneg!=_overlap.end() )
            prescale = _single_trigger_prescale[*itneg];
        if( itplus2!=_overlap2.end() )
            prescale = _single_trigger_prescale[*itplus2];
//         if( prescale != 1.0 )
//             cout << " prescale_fact " << prescale << endl;
        if (itplus!=_single_triggers.end()) 
        {
            probplus += it->second * prescale;
//             if (name=="MUEM1_LEL12_TRK5") eff_MUEM1 = it->second;
//             if (name=="MUEM2_LEL12_TRK5") eff_MUEM2 = it->second; 
//      if (name=="MUEM1_SH12_TRK5") eff_MUEM1 = it->second;
//      if (name=="MUEM2_SH12_TRK5") eff_MUEM2 = it->second;
        }
        else if (itneg != _overlap.end())
        {
            probneg += it->second * prescale;
        }
        else if (itplus2 != _overlap2.end() )
        {
            probplus += it->second * prescale;
        }
        else
        {
            cout << " ERROR in Probability_EMU Oring: could not find efficiency" << " for trigger: " << name << endl;
        }
    }

    probability = probplus - probneg;

    for (map<string,float>::const_iterator it = triggs.begin() ; it != triggs.end() ; it++)
    {

    // compare sample name ("ttll", "ztautau", etc)
        if (sample != it->first.substr(0,it->first.find(":"))) continue ;

        unsigned int pos2 = it->first.rfind(":");
        string name = it->first.substr(pos2+1,string::npos);
        set<string>::const_iterator its = _single_triggers.find(name);
        if (its!=_single_triggers.end())
        {
            map<string, Eff>::iterator jt = _globtrigs.find(it->first);
            if (jt == _globtrigs.end())
            {
                cout << "Internal error 1 in Probability_EMU: " << name << endl ;
            }
            else
            {
                jt->second += probability ;
            }
        } // its!=_single_triggers.end
    } // for triggs

    if ( debug() > 1)
        cout << "EMU probability for version " << sname << ":" << version << " = " << probability << endl;
    if (probability < 0 )
        probability = 0 ;
    else if ( probability > 1.0 )
        probability = 1.0 ;

    _versions[sname+":"+version] += probability ;

    return probability;
}


void Probability_DILEP::finish()
{
    probProcessor::finish() ;

    if( nevents == 0 )
        return;

    if (_outname.empty()) return ;

    ofstream tex((_outname + ".tex").c_str());
    tex.setf(ios::fixed);
    int pr = 2 ;
    tex.precision(pr) ;

    vector<TableLine> lines ;
    vector<TableLine1> lines1 ;
    float ilumi = 0 ;

    for(map<string, Eff>::const_iterator tt = _versions.begin() ; tt != _versions.end(); tt++)
    {
        ostringstream eff ;
        eff.setf(ios::fixed);
        eff.precision(pr) ;
        eff << "$ " <<  100.0*tt->second.eff << " \\pm " 
                << 100.0*tt->second.uncert()  << " $" ;
        string name = tt->first ;
        int pos = name.find(":") ;
        string version = name.substr(pos+1, name.size()-pos-1) ;

//             if(name.find("jet") != string::npos) continue ;

//         if(name.find("ttll") == string::npos &&
//            name.find("ww") == string::npos &&
//            name.find("ztautau_60_130") == string::npos 
//           ) continue ;

//         if(name.find("ttll") != string::npos) {
// @@       name.find("ttll_") == string::npos) {
        std::map<std::string, float>::const_iterator lt = _mapVersionLumi.find(version);

        if (lt == _mapVersionLumi.end()) 
        {
            cout << "Probability_EMU: Version [" << version << "] does not found in the Lumi map" << endl ;
            continue ;
        }
        ostringstream lm ;
        ilumi += lt->second ;
        if (version == "v14") 
        {
            float lumi = 81.65 ;
            lm << "\\begin{minipage}{3.2cm} " << lt->second - lumi 
                    << "(refixed skim) + " << lumi << "(p17.09 skim) \\end{minipage}" ;
        }
        else
            lm << "$ " <<  lt->second << " $" ;
        if (version == "v8") 
        {
            lines.insert(lines.begin(), TableLine(version,"",lm.str())) ;
            lines.begin()->eff.push_back(eff.str()) ;
        } else {
            lines.push_back(TableLine(version,"",lm.str())) ;
            lines.back().eff.push_back(eff.str()) ;
        }
//         }
//         else
//         {
//             for (vector<TableLine>::iterator lt = lines.begin() ; lt != lines.end() ; lt++) 
//             {
//                 if (lt->version == version) 
//                 {
//                     lt->eff.push_back(eff.str()) ;
//                     break ;
//                 }
//             }
//         }
    }

    for(map<string, Eff>::const_iterator tt = _globtrigs.begin() ; tt != _globtrigs.end(); tt++) 
    {
        string name = tt->first ;
        int pos = name.find(":") ;
        int pos1 = name.find(":",pos+1) ;
        string version = name.substr(pos+1, pos1-pos-1) ;
//         if(name.find("ttll") != string::npos && name.find("ttll_") == string::npos) break ;
        for (vector<TableLine>::iterator lt = lines.begin() ; lt != lines.end() ; lt++) 
        {
            if (lt->version != version) continue ;
            string trgname = name.substr(pos1+1, name.size()-pos1-1) ;
            if (trgname.find("_v11") != string::npos)
                trgname.replace(trgname.find("_v11"),4,"") ;
            else if (trgname.find("_v8") != string::npos)
                trgname.replace(trgname.find("_v8"),3,"") ;
            lt->trigger_name = trgname ;
            if (lt->version == "v8") lt->version = "v8 -- v11" ;
            else if (lt->version == "v13") lt->version = "v13 -- v13.3" ;
            else if (lt->version == "v13.3") lt->version = "v13.3 -- v14" ;
        }
    }

    for(map<TrigStruct, Eff>::const_iterator tt = _triggers.begin() ; tt != _triggers.end(); tt++) 
    {
        string name = tt->first.sample->name() ;
//         if(name.find("ttll") == string::npos) continue ;
//         if(name.find("jet") != string::npos) continue ;

        string trg = tt->first.trigger_name ;

        ostringstream eff ;
        eff.setf(ios::fixed);
        eff.precision(pr) ;
        eff << "$ " <<  100.0*tt->second.eff << " \\pm " 
                << 100.0*tt->second.uncert()  << " $" ;

        string obj = tt->first.object ;
        if (obj.find("_v8") != string::npos)
            obj.replace(obj.find("_v8"),3,"") ;
        else if (obj.find("_v1") != string::npos)
            obj.replace(obj.find("_v1"),4,"") ;
    
        vector<TableLine1>::iterator lt = lines1.begin();
        for (; lt != lines1.end(); lt++)
            if ( lt->version == tt->first.version && lt->trigger_name == trg ) break ;
    
        if (lt != lines1.end()) 
        {
            lt->tname.push_back(obj) ;
            lt->eff.push_back(eff.str()) ;
        }
        else 
        {
            if (tt->first.version == "v8") 
            {
                lines1.insert(lines1.begin(), TableLine1(tt->first.version,trg)) ;
                lines1.begin()->tname.push_back(obj) ;
                lines1.begin()->eff.push_back(eff.str()) ;
            }
            else 
            {
                lines1.push_back(TableLine1(tt->first.version,trg)) ;
                lines1.back().tname.push_back(obj) ;
                lines1.back().eff.push_back(eff.str()) ;
            }
        }
    }

    for (vector<TableLine1>::iterator lt = lines1.begin() ; lt != lines1.end() ; lt++) 
    {
        if (lt->version == "v8") lt->version = "v8 -- v11" ;
        else if (lt->version == "v13") lt->version = "v13 -- v13.3" ;
        else if (lt->version == "v13_3") lt->version = "v13.3 -- v14" ;
        if (lt->trigger_name.find("_v11") != string::npos) 
            lt->trigger_name.replace(lt->trigger_name.find("_v11"),4,"") ;
        else if (lt->trigger_name.find("_v8") != string::npos) 
            lt->trigger_name.replace(lt->trigger_name.find("_v8"),3,"") ;
    }


    cout << "========================================= " << endl ;
    for(map<TrigStruct, Eff>::const_iterator tt = _triggers.begin() ; tt != _triggers.end(); tt++)
        cout << tt->first.name << ": " 
                << 100.0*tt->second.eff << " +-" << 100.0*tt->second.uncert()  << endl ;

    cout << "========================================= " << endl ;
    for(map<string, Eff>::const_iterator tt = _globtrigs.begin() ; tt != _globtrigs.end(); tt++)
        cout << tt->first << ": " 
                << 100.0*tt->second.eff << " +-" << 100.0*tt->second.uncert()  << endl ;
    cout << "========================================= " << endl ;

    for(map<string, Eff>::const_iterator tt = _versions.begin() ; tt != _versions.end(); tt++)
        cout << tt->first << ": " 
                << 100.0*tt->second.eff << " +-" << 100.0*tt->second.uncert()  << endl ;
    cout << "========================================= " << endl ;

    tex << "" << endl;
    tex << "% Trigger Efficiencies / Lumi table %" << endl;
    tex << endl;
    tex << "\\begin{table}[p]" << endl;
    tex << "\\begin{center}" << endl;
    tex << "\\begin{tabular}{l|l|r|r|r|r}" << endl;
    tex << "  \\hline  \\hline" << endl;
    tex << "Trigger Version & Triggers Names & "
            << "\\multicolumn{3}{c|}{Efficiencies [\\%]} & Luminosity, pb$^{-1}$ \\\\ \\cline{3-5}" 
            << endl;      
    tex << "&& $t\\bar{t}$ MC &  $WW$ MC & $Z\\to \\tau\\tau$ MC & \\\\ \\hline" 
            << endl;      

    for (vector<TableLine>::iterator lt = lines.begin() ;
         lt != lines.end() ; lt++) {
             lt->out(tex) ;
             tex << endl ;
         }

         tex << " \\hline Total & & " ; 
         if (STAT) {
             vector<const StatSample*> samples = STAT->get_samples() ;
             for (vector<const StatSample*>::const_iterator st = samples.begin(); st != samples.end(); st++) 
             {
                 string name = (*st)->name() ;
//                                        if(name.find("jet") != string::npos) continue ;
//                                        if(name.find("ttll") == string::npos &&
//                                           name.find("ww") == string::npos &&
//                                           name.find("ztautau_60_130") == string::npos 
//                                          ) continue ;
                 if( (*st)->size() == 0 )
                     continue;
                 const StatWeight*  trig = (*st)->eventWeight("TriggerProbability") ;

                 if (!trig) {
                     cout << "Probability_EMU:  No trigger weight found in Stat!" << endl ;
                     continue ;
                 }
                 if (trig->nevents()>0)
                     tex << " $ " << 100.0*trig->weight_average() << " \\pm " << 100.0*trig->err() << " $ " ;
                 tex << " & " ;
             }
         }

         tex << ilumi << " \\\\ \\hline " << endl ;
         tex << " \\end{tabular}" << endl;
         tex << " \\caption{\\label{tab:emu_trigg} The triggers used in the analysis with corresponding luminosities and efficiencies.}" << endl;
         tex << "\\end{center}" << endl;
         tex << "\\end{table}" << endl; 
        
         tex << endl;
         tex << "% L1,L2,L3 Trigger Efficiencies table %" << endl;
         tex << endl;
         tex << "\\begin{table}[p]" << endl;
         tex << "\\begin{center}" << endl;
         tex << "\\begin{tabular}{|l|l|l|r|r|r|}" << endl;
         tex << "\\hline" << endl;
         tex << "Trigger Version & Global Triggers Names & L1, L2, L3 Trigger Term& "
                 << "Efficiencies [\\%] \\\\ \\hline" 
                 << endl;      
        
         for (vector<TableLine1>::iterator lt = lines1.begin() ;
              lt != lines1.end() ; lt++) {
                  lt->out(tex) ;
                  tex << endl ;
              }
        
              tex << " \\end{tabular}" << endl;
              tex << " \\caption{\\label{tab:emu_trigg_terms} Individual triggers efficiency}" << endl;
              tex << "\\end{center}" << endl;
              tex << "\\end{table}" << endl; 
        
              tex.close() ;
}

ClassImp(Probability_DILEP);



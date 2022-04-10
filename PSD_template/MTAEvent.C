#include "MTAEvent.h"
#include "TWfd4Fold.h"
#include "TTPCEventProperties.h"
#include "TTPCUpstreamEventProperties.h"

#include "TTPCIsland.h"

#include "TTPCMiniPulse.h"
#include "TTPCMiniTPulse.h"
#include "TTPCRawPulse.h"
#include "TTPCGaussPulse.h"
#include "TTPCTOTPulse.h"
#include "TTPCQuadPulse.h"
#include "TGondolaCluster.h"
#include "TTPCIslandPulsesWrapper.h"

#include "Parameters.h"

extern TPC_PARAM tpc_parameters;

using std::vector;

TClonesArray* MTAEvent::fgWfd4Fold = 0;

MTAEvent::MTAEvent()
  : fTPCFusionIsland(NULL)
{
  printf("Creating MTAEvent object.\n");
  fMusunMTAEvent = new TMusunMTAEvent();

  if (!fgWfd4Fold) {
    fgWfd4Fold = new TClonesArray("TWfd4Fold", 1000,kTRUE);
    fgWfd4Fold->BypassStreamer(kFALSE);
  }
  fWfd4Fold = fgWfd4Fold;
  fNWfd4Fold = 0;
  fIsTrigger=0;
  fSkimEvent=0;
  fSucceededInPrereadingMuEntrances = false;

  fTPCEventProperties = new TTPCEventProperties();
  fTPCRawPulseEventProperties = new TTPCEventProperties();
  fTPCTOTPulseEventProperties = new TTPCEventProperties();
  fTPCTOTRoadTrackEventProperties = new TTPCEventProperties();
  fTPCTOTRoadTrackStopThresholdEventProperties = new TTPCEventProperties();
  fTPCUpstreamEventProperties = new TTPCUpstreamEventProperties();
}

MTAEvent::~MTAEvent()
{
  delete fMusunMTAEvent;
  if(fTPCEventProperties) { 
    delete fTPCEventProperties;
    fTPCEventProperties = NULL;
  }
  if(fTPCRawPulseEventProperties) { 
    delete fTPCRawPulseEventProperties;
    fTPCRawPulseEventProperties = NULL;
  }
  if(fTPCTOTPulseEventProperties) { 
    delete fTPCTOTPulseEventProperties;
    fTPCTOTPulseEventProperties = NULL;
  }
  if(fTPCTOTRoadTrackEventProperties) { 
    delete fTPCTOTRoadTrackEventProperties;
    fTPCTOTRoadTrackEventProperties = NULL;
  }
  if(fTPCTOTRoadTrackStopThresholdEventProperties) { 
    delete fTPCTOTRoadTrackStopThresholdEventProperties;
    fTPCTOTRoadTrackStopThresholdEventProperties = NULL;
  }
  if(fTPCUpstreamEventProperties) { 
    delete fTPCUpstreamEventProperties;
    fTPCUpstreamEventProperties = NULL;
  }
}

void MTAEvent::SetMusunEvent(TMusunEvent *event)
{
  fMusunEvent = event;
}

void MTAEvent::SetMusunMTAEvent(TMusunMTAEvent *musunmtaevent)
{
  fMusunMTAEvent = musunmtaevent;
}

void MTAEvent::Clear(Option_t *option)
{
  fNWfd4Fold = 0;
  if(fTPCEventProperties) { 
    delete fTPCEventProperties;
    fTPCEventProperties = NULL;
  }
  if(fTPCRawPulseEventProperties) { 
    delete fTPCRawPulseEventProperties;
    fTPCRawPulseEventProperties = NULL;
  }
  if(fTPCTOTPulseEventProperties) { 
    delete fTPCTOTPulseEventProperties;
    fTPCTOTPulseEventProperties = NULL;
  }
  if(fTPCTOTRoadTrackEventProperties) { 
    delete fTPCTOTRoadTrackEventProperties;
    fTPCTOTRoadTrackEventProperties = NULL;
  }
  if(fTPCTOTRoadTrackStopThresholdEventProperties) { 
    delete fTPCTOTRoadTrackStopThresholdEventProperties;
    fTPCTOTRoadTrackStopThresholdEventProperties = NULL;
  }
  if(fTPCUpstreamEventProperties) { 
    delete fTPCUpstreamEventProperties;
    fTPCUpstreamEventProperties = NULL;
  }

  if(fTPCFusionIsland) { 
    delete fTPCFusionIsland;
    fTPCFusionIsland=NULL;
  }
}

void MTAEvent::Reset(Option_t *option)
{
  delete fgWfd4Fold; fgWfd4Fold = 0;
}

void MTAEvent::SetTPCEventProperties(TTPCEventProperties* ep)
{
  if(fTPCEventProperties) { 
    delete fTPCEventProperties;
  }
  fTPCEventProperties = ep;
}

void MTAEvent::SetTPCRawPulseEventProperties(TTPCEventProperties* ep)
{
  if(fTPCRawPulseEventProperties) { 
    delete fTPCRawPulseEventProperties;
  }
  fTPCRawPulseEventProperties = ep;
}

void MTAEvent::SetTPCTOTPulseEventProperties(TTPCEventProperties* ep)
{
  if(fTPCTOTPulseEventProperties) { 
    delete fTPCTOTPulseEventProperties;
  }
  fTPCTOTPulseEventProperties = ep;
}

void MTAEvent::SetTPCTOTRoadTrackEventProperties(TTPCEventProperties* ep)
{
  if(fTPCTOTRoadTrackEventProperties) { 
    delete fTPCTOTRoadTrackEventProperties;
  }
  fTPCTOTRoadTrackEventProperties = ep;
}

void MTAEvent::SetTPCTOTRoadTrackStopThresholdEventProperties(TTPCEventProperties* ep)
{
  if(fTPCTOTRoadTrackStopThresholdEventProperties) { 
    delete fTPCTOTRoadTrackStopThresholdEventProperties;
  }
  fTPCTOTRoadTrackStopThresholdEventProperties = ep;
}

void MTAEvent::SetTPCUpstreamEventProperties(TTPCUpstreamEventProperties* ep)
{
  if(fTPCUpstreamEventProperties) { 
    delete fTPCUpstreamEventProperties;
  }
  fTPCUpstreamEventProperties = ep;
}

void MTAEvent::SetPointers(TTree* tree)
{
  SetRunInfoEvent(0);  // initialization to zero is important
  TBranch *br_runinfo = tree->GetBranch("RunInfoEvent");
  if (!br_runinfo) printf("did not find the RunInfoEvent branch\n");
  br_runinfo->SetAddress(&fRunInfoEvent);
}

void MTAEvent::SetRunInfoEvent(TRunInfoEvent *event)
{
  fRunInfoEvent = event; 
}

TWfd4Fold* MTAEvent::AddWfd4Fold()
{
   TClonesArray &Wfd4Folds = *fWfd4Fold;
   TWfd4Fold *Wfd4Fold = new(Wfd4Folds[fNWfd4Fold++]) TWfd4Fold();
   return Wfd4Fold;
}

vector<TTPCPulseCluster*>* MTAEvent::GetTPCPulseClusters()
{
  return &fTPCPulseClusters;
}

vector<TTPCPulseCluster*>* MTAEvent::GetTPCRawPulseClusters()
{
  return &fTPCRawPulseClusters;
}

vector<TTPCPulseCluster*>* MTAEvent::GetTPCTOTPulseClusters()
{
  return &fTPCTOTPulseClusters;
}

vector<const TMusunNeutronPulse*>* MTAEvent::GetAllNeutronPulses()
{
  return &fAllNeutron;
}

vector<const TMusunNeutronPulse*>* MTAEvent::GetNeutronPulses()
{
  return &fNeutron;
}

vector<const TMusunNeutronPulse*>* MTAEvent::GetNeutronMustopPulses()
{
  return &fNeutronMustop;
}

//**Get electron stuff
vector<TGondolaCluster*>* MTAEvent::GetGondolaClusters()
{
  return &fGondolaClusters;
}

vector<TePCCluster*>* MTAEvent::GetePCClusters()
{
  return &fePCClusters;
}

vector<TGlobalElectronTrack*>* MTAEvent::GetGlobalElectronTracks()
{
  return &fGlobalElectronTracks;
}

vector<TGlobalElectronTrack_SingleTube*>* MTAEvent::GetGlobalElectronTracks_SingleTube()
{
  return &fGlobalElectronTracks_SingleTube;
}

//subselections
vector<TGlobalElectronTrack*> MTAEvent::GetTDCCathodeORTracks()
{
  vector<TGlobalElectronTrack*> tracks;
  vector<TGlobalElectronTrack*>* alltracks = GetGlobalElectronTracks();
  if(alltracks->size() < 1) return tracks;
  for(int i = 0; i<alltracks->size(); i++)
  {
    if(alltracks->at(i)->CathodeOR_TDC()) tracks.push_back(alltracks->at(i));
  }
  return tracks;
}

vector<TGlobalElectronTrack*> MTAEvent::GetWFDCathodeORTracks()
{
  vector<TGlobalElectronTrack*> tracks;
  vector<TGlobalElectronTrack*>* alltracks = GetGlobalElectronTracks();
  if(alltracks->size() < 1) return tracks;
  for(int i = 0; i<alltracks->size(); i++)
  {
    if(alltracks->at(i)->CathodeOR_WFD() ) tracks.push_back(alltracks->at(i));
  }
  return tracks;
}

vector<TGlobalElectronTrack*> MTAEvent::GetTDCCathodeANDTracks()
{
  vector<TGlobalElectronTrack*> tracks;
  vector<TGlobalElectronTrack*>* alltracks = GetGlobalElectronTracks();
  if(alltracks->size() < 1) return tracks;
  for(int i = 0; i<alltracks->size(); i++)
  {
    if(alltracks->at(i)->CathodeAND_TDC() ) tracks.push_back(alltracks->at(i));
  }
  return tracks;
}

vector<TGlobalElectronTrack*> MTAEvent::GetWFDCathodeANDTracks()
{
  vector<TGlobalElectronTrack*> tracks;
  vector<TGlobalElectronTrack*>* alltracks = GetGlobalElectronTracks();
  if(alltracks->size() < 1) return tracks;
  for(int i = 0; i<alltracks->size(); i++)
  {
    if(alltracks->at(i)->CathodeAND_WFD() ) tracks.push_back(alltracks->at(i));
  }
  return tracks;
}

vector<TGlobalElectronTrack*> MTAEvent::GetTDCAnodeTracks()
{
  vector<TGlobalElectronTrack*> tracks;
  vector<TGlobalElectronTrack*>* alltracks = GetGlobalElectronTracks();
  if(alltracks->size() < 1) return tracks;
  for(int i = 0; i<alltracks->size(); i++)
  {
    if(alltracks->at(i)->Anode_TDC() ) tracks.push_back(alltracks->at(i));
  }
  return tracks;
}

vector<TGlobalElectronTrack*> MTAEvent::GetWFDAnodeTracks()
{
  vector<TGlobalElectronTrack*> tracks;
  vector<TGlobalElectronTrack*>* alltracks = GetGlobalElectronTracks();
  if(alltracks->size() < 1) return tracks;
  for(int i = 0; i<alltracks->size(); i++)
  {
    if(alltracks->at(i)->Anode_WFD() ) tracks.push_back(alltracks->at(i));
  }
  return tracks;
}

vector<TGlobalElectronTrack*> MTAEvent::GetGlobalElectronTracks(TGondolaCluster* cluster)
{
  vector<TGlobalElectronTrack*> tracks;
  vector<TGlobalElectronTrack*>* alltracks = GetGlobalElectronTracks();
  if(alltracks->size() < 1) return tracks;
  for(int i = 0; i<alltracks->size(); i++)
  {
    TGondolaCluster* clusteroftrack = alltracks->at(i)->GetGondolaCluster();
    if(clusteroftrack == cluster) tracks.push_back(alltracks->at(i));
  }
  return tracks;  
}

vector<TGondolaCluster*> MTAEvent::GetG4foldsTDC()
{
  vector<TGondolaCluster*> vout; 
  vector<TGondolaCluster*>* allclusters = &fGondolaClusters; 
  if(allclusters->size() < 1) return vout;
  for(int i = 0; i < allclusters->size(); i++)
  {
    if( allclusters->at(i)->IsTDCCluster() &&  allclusters->at(i)->GetGondolaClusterSize() == 4) { vout.push_back(allclusters->at(i));}
  }
  return vout;
}

vector<TGondolaCluster*> MTAEvent::GetG4foldsWFD()
{
  vector<TGondolaCluster*> clusters;
  vector<TGondolaCluster*>* allclusters = GetGondolaClusters();
  if(allclusters->size() < 1) return clusters;
  for(int i = 0; i<allclusters->size(); i++)
  {
    if(allclusters->at(i)->IsWFDCluster() &&  allclusters->at(i)->GetGondolaClusterSize() == 4) clusters.push_back(allclusters->at(i));
  }
  return clusters;
}



vector<const TMusunNeutronPulse*>* MTAEvent::GetGammaPulses()
{
  return &fGamma;
}
vector<const TMusunNeutronPulse*>* MTAEvent::GetFusionNeutronPulses()
{
  return &fFusionNeutron;
}

vector<const TMusunNeutronPulse*>* MTAEvent::GetDelayedElectronNeutronPulses()
{
  return &fDelayedElectronNeutron;
}

vector<const TMusunNeutronPulse*>* MTAEvent::GetCaptureNeutronPulses()
{
  return &fCaptureNeutron;
}


/** A more general function to flag events with muon tracks. 
 *  At the time of writing this comment, this function is equivalent to the
 *  HasContPadTrack function, but this might have chnaged with the addition
 *  of new track definitions.
 */
bool MTAEvent::HasMuonTrack()
{
  return false;
}

vector<TTPCGenericPulse*>* MTAEvent::GetTPCMiniPulses()
{
  return &GetMusunEvent()->GetTPCMiniPulses_HG();
}

vector<TTPCGenericPulse*>* MTAEvent::GetTPCMiniTPulses()
{
  return &GetMusunEvent()->GetTPCMiniTPulses_HG();
}

vector<TTPCGenericPulse*>* MTAEvent::GetTPCGaussPulses()
{
  return &GetMusunEvent()->GetTPCGaussPulses_HG();
}

vector<TTPCGenericPulse*>* MTAEvent::GetTPCRawPulses()
{
  return &GetMusunEvent()->GetTPCRawPulses_HG();
}

vector<TTPCGenericPulse*>* MTAEvent::GetTPCTOTPulses()
{
  return &GetMusunEvent()->GetTPCTOTPulses_HG();
}

vector<TTPCGenericPulse*>* MTAEvent::GetTPCQuadPulses()
{
  return &GetMusunEvent()->GetTPCQuadPulses_HG();
}


/** Fixes the mistake in the wiring map introduced in the 
  * run4 analysis. This bug will also be fixed upstream,
  * so this function will be obsolete on the next production
  * pass.
  * 
  * This function swaps pads 7 and 8 (one-indexed) on the TPC
  * pulse and island objects, reversing the wiring mistake noted
  * in elog:
  *
  * https://muon.npl.washington.edu/elog/musun/analysis-run4/464
  *
  * This function first checks the wiring map in the TPC_PARAM object
  * loaded from the OdbTree in the tree#####.root file. If the wiring
  * for channel 7 (zero-indexed) matches the fixed version (WFD channel 3),
  * this function returns without doing anything. For reference:
  * 
  * In the INCORRECT wiring map, we had:
  *   - TPC channel 6 (zero-indexed), which is "/Analyzer/TPC/PAD_01_07/"
  *     in the ODB, was mapped to WFD crate 2, slot 16, channel 2.
  *   - TPC channel 7 (zero-indexed), which is "/Analyzer/TPC/PAD_01_08/"
  *     in the ODB, was mapped to WFD crate 2, slot 16, channel 3.
  *
  * In the CORRECT wiring map, we have:
  *   - TPC channel 6 (zero-indexed), which is "/Analyzer/TPC/PAD_01_07/"
  *     in the ODB, is mapped to WFD crate 2, slot 16, channel 3.
  *   - TPC channel 7 (zero-indexed), which is "/Analyzer/TPC/PAD_01_08/"
  *     in the ODB, is mapped to WFD crate 2, slot 16, channel 2.
  *
  */
void MTAEvent::FixTPCWiringMistakeRun4()
{
  // Check if it needs to be fixed
  TWFDHardwareMap* wfd = tpc_parameters.GetWFDHardwareMap(7,2);
  if(wfd->fChannel==2) return;

  vector< vector<TTPCGenericPulse*>*> vec_all_pulses;
  vec_all_pulses.push_back(&GetMusunEvent()->GetTPCAllPulses_HG());
  vec_all_pulses.push_back(&GetMusunEvent()->GetTPCAllPulses_LG());

  for(size_t i=0; i<vec_all_pulses.size(); i++) {
    for(size_t j=0; j<vec_all_pulses[i]->size(); j++) {
      TTPCGenericPulse* pulse = vec_all_pulses.at(i)->at(j);
      if(pulse->GetPadZeroIndexed() == 6) {
        pulse->SetPadZeroIndexed(7);
      } else if (pulse->GetPadZeroIndexed() == 7) {
        pulse->SetPadZeroIndexed(6);
      }
    }
  }

  vector< vector<TTPCIsland*>* > vec_all_islands;
  vec_all_islands.push_back(&GetMusunEvent()->GetTPCIslands_LG());
  vec_all_islands.push_back(&GetMusunEvent()->GetTPCIslands_HG());

  for(size_t i=0; i<vec_all_islands.size(); i++) {
    for(size_t j=0; j<vec_all_islands[i]->size(); j++) {
      TTPCIsland* island = vec_all_islands.at(i)->at(j);
      if(island->GetPadZeroIndexed() == 6) {
        island->SetPadZeroIndexed(7);
      } else if(island->GetPadZeroIndexed() == 7) {
        island->SetPadZeroIndexed(6);
      }
    }
  }

}

/** Returns TPC pulses of the given class. 
  *
  * Multiplexes into the GetTPC*Pulses() methods. If pulse_class is NULL, returns NULL.
  * If pulse_class isn't a known pulse class, prints an error and returns NULL.
  *
  * Note that if a new pulse class is added, this method must be updated.
  */
vector<TTPCGenericPulse*>* MTAEvent::GetTPCPulses(TClass* pulse_class)
{
       if(pulse_class==NULL) { return NULL; }
  else if(pulse_class==TTPCMiniPulse::Class()) { return GetTPCMiniPulses(); }
  else if(pulse_class==TTPCMiniTPulse::Class()) { return GetTPCMiniTPulses(); }
  else if(pulse_class==TTPCRawPulse::Class()) { return GetTPCRawPulses(); }
  else if(pulse_class==TTPCTOTPulse::Class()) { return GetTPCTOTPulses(); }
  else if(pulse_class==TTPCGaussPulse::Class()) { return GetTPCGaussPulses(); }
  else if(pulse_class==TTPCQuadPulse::Class()) { return GetTPCQuadPulses(); }
  else {
    printf("Error: MTAEvent::GetTPCPulses() : Invalide pulse class specified : %s. "
           "Please update this method or check the pulse class requested.\n",
           pulse_class->GetName());
    return NULL;
  }

}

vector<TTPCIsland*>* MTAEvent::GetTPCIslands()
{
  return &GetMusunEvent()->GetTPCIslands_HG();
}

void MTAEvent::SetTPCFusionIsland(TTPCIsland* island)
{
  if(fTPCFusionIsland) { 
    delete fTPCFusionIsland;
    fTPCFusionIsland=NULL;
  }
  fTPCFusionIsland = island;
}

TTPCIslandPulsesWrapper* MTAEvent::GetIslandPulseWrapper(int blockTime, int pad)
{
  vector<TTPCIslandPulsesWrapper*>& wrappers = GetMusunEvent()->GetTPCIslandPulsesWrappers_HG();
  if(wrappers.size()==0) { return NULL;  }

  int j = -1; 
  for(int i=0; i < wrappers.size(); i++)
  {
    if(wrappers.at(i)->GetBlockTime() == blockTime && (wrappers.at(i)->GetPadZeroIndexed()+1) == pad ) j = i;     
  }
  if(j > -1) return wrappers.at(j); else return NULL;
}


TTPCIslandPulsesWrapper* MTAEvent::GetIslandPulseWrapper(double time, int pad)
{
  vector<TTPCIslandPulsesWrapper*>& wrappers = GetMusunEvent()->GetTPCIslandPulsesWrappers_HG();
  if(wrappers.size()==0) { return NULL;  }

  int j = -1; 
  for(int i=0; i < wrappers.size(); i++)
  {
    if(wrappers.at(i)->GetEndTime() > time && wrappers.at(i)->GetBlockTime_ns() < time && (wrappers.at(i)->GetPadZeroIndexed()+1) == pad ) j = i;     
  }
  if(j > -1) return wrappers.at(j); else return NULL;
}



TTPCIslandPulsesWrapper* MTAEvent::GetIslandPulseWrapper(TTPCGenericPulse* pulse)
{
  int bTime = pulse->GetBlockTime();
  int pad = pulse->GetPadOneIndexed();
  return  GetIslandPulseWrapper(bTime,pad);  
}

bool MTAEvent::IsCleanBlock()
{
  bool clean = false;
  TMusunEvent* musunEvent = GetMusunEvent();

  if( !musunEvent->HasMuSCMatchError() && !musunEvent->HasCaenError() && !musunEvent->HasBlockError() )
  {
    if( !musunEvent->GetSparkInBlock() )
    {
      clean = true;
    }
  }  
  return clean;
}

bool MTAEvent::IsCleanEvent()
{
  bool clean = false;
  TMusunEvent* musunEvent = GetMusunEvent();
  if( IsCleanBlock() && !musunEvent->TDCForcedTrigger() && !musunEvent->TDCForcedTriggerPulser() )
  {
    clean = true;
  }  
}

/** Gets muon tracks based on the specified algorithm.
  * Your choices are the following
  *
  *  - "BasicMini"
  *  - "BasicTOT"
  *  - "RoadTOT"
  *  - "RoadStopThreshTOT"
  *
  * Since this uses string comparison to figure out which vector to return,
  * you shouldn't use this for performance-intensive lookups. It's better
  * to get the TTPCEventProperties that you want directly from the appropriate
  * MTAEvent method. This method is supposed to be used for the event display.
  */
vector<TTPCGenericMuonTrack*> MTAEvent::GetMuonTracks(const char* track_name)
{
  TString name(track_name);
  name.ToLower();

  TTPCEventProperties* ep = NULL;
       if(name == "basicmini") ep = GetTPCEventProperties();
  else if(name == "basictot") ep = GetTPCTOTPulseEventProperties();
  else if(name == "roadtot") ep = GetTPCTOTRoadTrackEventProperties();
  else if(name == "roadstopthreshtot") ep = GetTPCTOTRoadTrackStopThresholdEventProperties();
  else ep = NULL;

  vector<TTPCGenericMuonTrack*> ret;
  if(ep == NULL) { ret = vector<TTPCGenericMuonTrack*>(); }
  else {
    ret = ep->GetAllMuonTracks();
  }

  return ret;
}

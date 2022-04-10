#ifndef ROOT_MTAEvent
#define ROOT_MTAEvent

/** A wrapper around TMusunEvent, providing a global storage for
 * MTA similar to gData in MU. Also provides access to the members
 * of TMusunEvent (either through direct methods or indirect:
 * GetMusunEvent()->GetTPCPulses()).
 */

#include "TMucapEvent.h"
#include "TMusunEvent.h"
#include "TKickerEvent.h"
#include "TMusunMTAEvent.h"
#include "TRunInfoEvent.h"
#include "TGondolaCluster.h"
#include "TGlobalElectronTrack.h" 
#include "TGlobalElectronTrack_SingleTube.h" 
#include "TePCCluster.h"

#include <TTree.h>
#include <TObject.h>
#include <vector>

class TWfd4Fold;
class TTPCPulseCluster;
class TMusunNeutronPulse;
class TTPCEventProperties;
class TTPCUpstreamEventProperties;
class TGondolaCluster;
class TePCCluster;
class TGlobalElectronTrack;
class TGlobalElectronTrack_SingleTube;
class TTPCGenericMuonTrack;

class MTAEvent
{

 private:
  TMusunEvent *fMusunEvent;
  TMusunMTAEvent *fMusunMTAEvent;
  TRunInfoEvent *fRunInfoEvent;  
  Int_t entryMusunEvent;
  Int_t          fNWfd4Fold; //Number of WFD 4-folds
  TClonesArray  *fWfd4Fold; //array of Wfd 4-folds
  static TClonesArray *fgWfd4Fold;
  bool fIsTrigger;
  bool fSkimEvent;
  bool fSucceededInPrereadingMuEntrances;
   
  std::vector<TTPCPulseCluster*> fTPCPulseClusters;
  std::vector<TTPCPulseCluster*> fTPCRawPulseClusters;
  std::vector<TTPCPulseCluster*> fTPCTOTPulseClusters;

  std::vector<const TMusunNeutronPulse*> fAllNeutron;
  std::vector<const TMusunNeutronPulse*> fFusionNeutron;
  std::vector<const TMusunNeutronPulse*> fDelayedElectronNeutron;
  std::vector<const TMusunNeutronPulse*> fNeutron;
  std::vector<const TMusunNeutronPulse*> fNeutronMustop;
  std::vector<const TMusunNeutronPulse*> fGamma;
  std::vector<const TMusunNeutronPulse*> fCaptureNeutron;
  TTPCIsland* fTPCFusionIsland;

  TTPCEventProperties* fTPCEventProperties;
  TTPCEventProperties* fTPCRawPulseEventProperties;
  TTPCEventProperties* fTPCTOTPulseEventProperties;
  TTPCEventProperties* fTPCTOTRoadTrackEventProperties;
  TTPCEventProperties* fTPCTOTRoadTrackStopThresholdEventProperties;
  TTPCUpstreamEventProperties* fTPCUpstreamEventProperties;

  std::vector<TGondolaCluster*>       fGondolaClusters;
  std::vector<TePCCluster*>           fePCClusters;
  std::vector<TGlobalElectronTrack*>  fGlobalElectronTracks;
  std::vector<TGlobalElectronTrack_SingleTube*>  fGlobalElectronTracks_SingleTube;

 public:
  MTAEvent();
  virtual ~MTAEvent();

  void SetPointers(TTree* tree);
  void SetMusunEvent(TMusunEvent *event);
  TMusunEvent *GetMusunEvent() { return fMusunEvent; }

  void SetMusunMTAEvent(TMusunMTAEvent *mmtaevent);
  TMusunMTAEvent *GetMusunMTAEvent() { return fMusunMTAEvent; }

  void SetRunInfoEvent(TRunInfoEvent *event);
  TRunInfoEvent *GetRunInfoEvent() { return fRunInfoEvent; }
  int GetPSIRunNr() { return fRunInfoEvent->GetPsiRunNum();}

  void Clear(Option_t *option="");
  static void Reset(Option_t *option= "");

  void SetEntry(Int_t value) { entryMusunEvent = value; }
  Int_t GetEntry() { return entryMusunEvent; }

  Int_t         GetNWfd4Fold() const { return fNWfd4Fold; }
  TClonesArray *GetWfd4Folds() const { return fWfd4Fold; }

  TMuEntrance* GetMuEntrance() { return GetMusunEvent()->GetMuEntrance(); }

  //TPC Pulse Stuff****************
  std::vector<TTPCPulseCluster*>* GetTPCPulseClusters();
  std::vector<TTPCPulseCluster*>* GetTPCRawPulseClusters();
  std::vector<TTPCPulseCluster*>* GetTPCTOTPulseClusters();

  std::vector<TTPCGenericPulse*>* GetTPCMiniPulses();
  std::vector<TTPCGenericPulse*>* GetTPCMiniTPulses();
  std::vector<TTPCGenericPulse*>* GetTPCGaussPulses();
  std::vector<TTPCGenericPulse*>* GetTPCRawPulses();
  std::vector<TTPCGenericPulse*>* GetTPCTOTPulses();
  std::vector<TTPCGenericPulse*>* GetTPCPulses(TClass* pulse_class=NULL);
  std::vector<TTPCGenericPulse*>* GetTPCQuadPulses();
  std::vector<TTPCGenericPulse*>* GetTPCQuadPulses(int pad, double tLow, double tHigh);


  std::vector<TTPCIsland*>* GetTPCIslands();

  TTPCIslandPulsesWrapper* GetIslandPulseWrapper(int blockTime, int pad); // 'pad' one indexed
  TTPCIslandPulsesWrapper* GetIslandPulseWrapper(double time, int pad);  //return wrapper if it covers 'time'
  TTPCIslandPulsesWrapper* GetIslandPulseWrapper(TTPCGenericPulse* pulse);

  //Get electrons ***************
  std::vector<TGondolaCluster*>*        GetGondolaClusters();
  std::vector<TePCCluster*>*            GetePCClusters();
  std::vector<TGlobalElectronTrack*>*   GetGlobalElectronTracks();
  std::vector<TGlobalElectronTrack_SingleTube*>*   GetGlobalElectronTracks_SingleTube();
  //sub selections
  std::vector<TGlobalElectronTrack*>   GetTDCCathodeORTracks();
  std::vector<TGlobalElectronTrack*>   GetWFDCathodeORTracks();
  std::vector<TGlobalElectronTrack*>   GetTDCCathodeANDTracks();
  std::vector<TGlobalElectronTrack*>   GetWFDCathodeANDTracks();
  std::vector<TGlobalElectronTrack*>   GetTDCAnodeTracks(); //don`t look at the Cathodes, but they might be there
  std::vector<TGlobalElectronTrack*>   GetWFDAnodeTracks();
  std::vector<TGlobalElectronTrack*>   GetGlobalElectronTracks(TGondolaCluster* cluster); //returns global electron tracks returning this cluster
  std::vector<TGondolaCluster*> GetG4foldsTDC();
  std::vector<TGondolaCluster*> GetG4foldsWFD();



  //Neutron Stuff*************
  std::vector<const TMusunNeutronPulse *>* GetAllNeutronPulses();
  std::vector<const TMusunNeutronPulse *>* GetNeutronPulses();
  std::vector<const TMusunNeutronPulse *>* GetNeutronMustopPulses();
  std::vector<const TMusunNeutronPulse *>* GetGammaPulses();
  std::vector<const TMusunNeutronPulse *>* GetFusionNeutronPulses();
  std::vector<const TMusunNeutronPulse *>* GetDelayedElectronNeutronPulses();
  std::vector<const TMusunNeutronPulse *>* GetCaptureNeutronPulses();

  // The event properties objects are a summary of all TPC information
  TTPCEventProperties*
    GetTPCEventProperties() {return fTPCEventProperties;}

  TTPCEventProperties*
    GetTPCRawPulseEventProperties() {return fTPCRawPulseEventProperties;}

  TTPCEventProperties*
    GetTPCTOTPulseEventProperties() {return fTPCTOTPulseEventProperties;}

  TTPCEventProperties*
    GetTPCTOTRoadTrackEventProperties() {return fTPCTOTRoadTrackEventProperties;}

  TTPCEventProperties*
    GetTPCTOTRoadTrackStopThresholdEventProperties() {return fTPCTOTRoadTrackStopThresholdEventProperties;}

  TTPCUpstreamEventProperties*
    GetTPCUpstreamEventProperties() {return fTPCUpstreamEventProperties;}

  std::vector<TTPCGenericMuonTrack*> GetMuonTracks(const char* track_name);

  void SetTPCEventProperties(TTPCEventProperties* p);
  void SetTPCRawPulseEventProperties(TTPCEventProperties* p);
  void SetTPCTOTPulseEventProperties(TTPCEventProperties* p);
  void SetTPCTOTRoadTrackEventProperties(TTPCEventProperties* p);
  void SetTPCTOTRoadTrackStopThresholdEventProperties(TTPCEventProperties* p);
  void SetTPCUpstreamEventProperties(TTPCUpstreamEventProperties* p);

  bool HasMuonTrack();

  TWfd4Fold    *AddWfd4Fold();
  void          SetNWfd4Fold(Int_t n) {fNWfd4Fold = n; }

  void SetTrigger() {fIsTrigger =1;}
  void ResetTrigger() {fIsTrigger =0;}
  bool IsTrigger() {return fIsTrigger;}

  void FixTPCWiringMistakeRun4();

  //skimmer options
  void SetSkimEventCaptureRecoilAnalysis(){ SetSkimEvent(); fMusunEvent->SetskimCaptureRecoilAnalysis(true);}
  void SetSkimEventMusunNeutronAnalysis(){ SetSkimEvent(); fMusunEvent->SetskimMusunNeutronAnalysis(true);}
  void SetSkimEventTPCClusterAnalysis(){ SetSkimEvent(); fMusunEvent->SetskimTPCClusterAnalysis(true);}
  void SetSkimEventTPCFusionAnalysis(){ SetSkimEvent(); fMusunEvent->SetskimTPCFusionAnalysis(true);}
  void SetSkimEventXRayMuStopAnalysis(){ SetSkimEvent(); fMusunEvent->SetskimXRayMuStopAnalysis(true);}

  void SetSkimEvent() {fSkimEvent = 1;}
  void ResetSkimEvent() {fSkimEvent = 0;}
  bool GetSkimEvent() {return fSkimEvent;}

  void SetTPCFusionIsland(TTPCIsland* island);
  TTPCIsland* GetTPCFusionIsland(){return fTPCFusionIsland;}


  bool IsCleanBlock();
  bool IsCleanEvent();

  void SetSucceededInPrereadingMuEntrances(bool b) {
    fSucceededInPrereadingMuEntrances = b; }
  bool SucceededInPrereadingMuEntrances() {
    return fSucceededInPrereadingMuEntrances; }



};

#endif

#ifndef MusunNeutronHists_h
#define MusunNeutronHists_h 1

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>
#include "MTAEvent.h"
#include "FillHistBase.h"
#include "TMusunMTAEvent.h"

//#include "MusunNeutronBlockAnalysis.h"
/** MuSun MTA module histogramming many neutron properties.
 *
 * Does everything with the neutron detectors for the MTA level analysis.
 */
class MusunNeutronHists : public FillHistBase
{
 public:
  MusunNeutronHists();
  virtual ~MusunNeutronHists();
  int ProcessEntry(MTAEvent *mtaEvent);

  TH3 *hCapture;
  //TH3 *hPSD;
  TH3 *hNeutronCS;
  TH3 *hNeutronMustop;
  TH3 *hGamma;
  TH3 *hDelay;
  TH3 *hFusion;
  
 private:
  void MakeHist();
  void FillHist();


};

#endif

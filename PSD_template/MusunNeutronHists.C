
/*
  MODULE DESCRIPTION:
 */

#include <vector>
#include <fstream>
#include <string>
//#include <stdio.h>

#include "MusunNeutronHists.h"
#include "MusunNeutronAnalysis.h"
#include "TMusunWfdPulse.h"
#include "TMusunNeutronPulse.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TTPCEventProperties.h"
#include "TTPCPulseCluster.h"

using std::vector;
//#include "MusunNeutronBlockAnalysis.h"

using std::vector;

//extern MUONEVENT_PARAM muonevent_parameters;
extern MTAEvent* g_mtaEvent;

MusunNeutronHists::~MusunNeutronHists()
{
/*
  delete hCapture;
  delete hNeutronCS;
  delete hNeutronMustop;
  delete hGamma;
  delete hDelay;
  delete hFusion;
  // delete hPSD;
 */
}

MusunNeutronHists::MusunNeutronHists(){
  //printf("Contructor NEUTRON HISTS -------------------\n");
  MakeHist();

}

int MusunNeutronHists::ProcessEntry(MTAEvent *mtaEvent){
  //printf("Process entry NEUTRON HISTS -------------------\n");
  FillHist();
  return 0;
}

void MusunNeutronHists::MakeHist(){

  printf("MAKING NEUTRON HISTS -------------------\n");


  hNeutronCS = new TH3F("hNeutronCS",
		   "Time and Energy spectra with neutron cut",
		   1600, -40000.0, 40000., 260, -1000., 25000.,8,-0.5,7.5);
  hFusion= new TH3D("hFusion",
                    "Time between mu entrance and Neutron pulses and energy for fusions",
                    1600, -40000.0, 40000.,500, 0., 25000.,8, -0.5, 7.5);

  printf("NEUTRON HISTS MADE -------------------\n");

}


void MusunNeutronHists::FillHist(){
  
  TMuEntranceBase* muent = g_mtaEvent->GetMuEntrance();
  Double_t muEntTime = muent->GetTime() ;
  //printf("FILLING NEUTRON HISTS -------------------\n");
  if(!muent->HasBestEntrance() || g_mtaEvent->GetMusunEvent()->HasBlockError()) return;  //Pile up protection

  TMusunEvent *event = g_mtaEvent->GetMusunEvent();
  TClonesArray &nC = *(event->GetMusunNeutronPulse());

  for(int y=0; y<event->GetNmusunneutronpulse(); y++){
    TMusunNeutronPulse *neu = (TMusunNeutronPulse*)(nC[y]);
    double neutronpeaktime = neu->GetPeakTime();
    double neutronchannel = neu->GetChannel();
    double brentArea = neu->GetBrentArea();
    if(neu->IsNeutronChiSquare()) hNeutronCS->Fill(neutronpeaktime-muEntTime,brentArea,neutronchannel);
    if(neu->IsFusionNeutronChiSquare()) hFusion->Fill(neutronpeaktime-muEntTime,brentArea,neutronchannel);
  }
}

#ifndef __MTAGLOBALS__
#define __MTAGLOBALS__
#include <TFile.h>


struct MODSWITCHES {
  int TPCEnergyTokeV;
  int MuEntrancePileupCorrection;
  int MusunNeutronAnalysis;
  int MusunNeutronHists;
  int MusunNeutronSTAnalysis;
  int MusunNeutronSTHists;
  int TPCClusterAnalysis;
  int TPCClusterHists;
  int TPCMiniStopHists;
  int TPCRawPulseClusterAnalysis;
  int TPCRawPulseClusterHists;
  int TPCRawMiniComparison;
  int TPCTOTPulseClusterAnalysis;
  int TPCTOTPulseClusterHists;
  int TPCTOTStopHists;
  int TPCTOTRoadStopHists;
  int TPCTOTRoadStopThreshStopHists;
  int TPCTOTMiniComparison;
  int TPCTOTMCComparison;
  int TPCTOTRoadTrackAnalysis;
  int TPCFusionAnalysis;
  int TPCPTAnalysis;
  int XRayMuStopAnalysis;
  int CaptureRecoilAnalysis;
  int TPCmuPC_Correlation;
  int TPCUpstreamAnalysis;
  int TPCUpstreamHists;
  int LifetimeVsImpactPar;
  int LifetimeVsXYZ;
  int LifetimeVsStopEnergy;
  int LifetimeVsEntranceDetectors;
  int LifetimeVsMCTruth;
  int LifetimeHists;
  int MuEPairTreeBuilder;
  int ElectronTomographyAnalysis;
  int MCTruthCheck;
  int TPCUpstreamCalibration;
  int McResponseCalibration;
  int MichelElectronStudy;
  int GlobalElectronAnalysis;
  int GondolaClusterHistograms;
  int GlobalElectronGondOnlyLifetime;
  int GlobalElectronTracksLifetimeSingleTube;
  int GlobalElectronTracksLifetime;
  int DelayedEventTreeBuilder;
  int TPCUpstreamEnergies;
  int ElectronBackgroundAnalysis;
  int TPCRoadAndClusterStopComparison;
  int GondolaWFDPedestalAnalysis;
  int TPCTrackPatternAnalysis;
  int EntranceEfficiencyAnalysis;
};

struct MUE_MODSWITCHES {
};

#endif

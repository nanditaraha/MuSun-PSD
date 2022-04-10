 #include "MTAGlobals.h"
TFile *fptr;

// main modules
MODSWITCHES g_modSwitches = {
  1, // int TPCEnergyTokeV;
  0, // int MuEntrancePileupCorrection;
  0, // int MusunNeutronAnalysis;
  0, // int MusunNeutronHists;
  1, // int MusunNeutronSTAnalysis
  0, // int MusunNeutronSTHists
  1, // int TPCClusterAnalysis
  1, // int TPCClusterHists
  1, // int TPCMiniStopHists
  0, // int TPCRawPulseClusterAnalysis
  0, // int TPCRawPulseClusterHists
  0, // int TPCRawMiniComparison
  1, // int TPCTOTPulseClusterAnalysis
  1, // int TPCTOTPulseClusterHists
  1, // int TPCTOTStopHists
  1, // int TPCTOTRoadStopHists
  1, // int TPCTOTRoadStopThreshStopHists
  0, // int TPCTOTMiniComparison
  0, // int TPCTOTMCComparison
  1, // int TPCTOTRoadTrackAnalysis
  0, // int TPCFusionAnalysis
  0, // int TPCPTAnalysis
  0, // int XRayMuStopAnalysis
  0, // int CaptureRecoilAnalysis
  0, // int TPCmuPC_Correlation
  0, // int TPCUpstreamAnalysis
  0, // int TPCUpstreamHists
  1, // int LifetimeVsImpactPar
  1, // int LifetimeVsXYZ
  1, // int LifetimeVsStopEnergy
  1, // int LifetimeVsEntranceDetectors
  0, // int LifetimeVsMCTruth
  1, // int LifetimeHists
  1, // int MuEPairTreeBuilder
  0, // int ElectronTomographyAnalysis
  0, // int MCTruthCheck
  0, // int TPCUpstreamCalibration
  0, // int McResponseCalibration
  0, // int MichelElectronStudy
  1, // int GlobalElectronAnalysis
  0, // int GondolaClusterHistograms
  0, // int GlobalElectronGondOnlyLifetime
  1, // int GlobalElectronTracksLifetimeSingleTube
  1, // int GlobalElectronTracksLifetime
  1, // int DelayedEventTreeBuilder
  0, // int TPCUpstreamEnergies 
  1, // int ElectronBackgroundAnalysis
  0, // int TPCRoadAndClusterStopComparison
  0, // int GondolaWFDPedestalAnalysis
  1, // int TPCTrackPatternAnalysis
  1, // int EntranceEfficiencyAnalysis

};

MUE_MODSWITCHES g_mue_modSwitches = {
};

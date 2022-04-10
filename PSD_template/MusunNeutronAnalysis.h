#ifndef MusunNeutronAnalysis_h
#define MusunNeutronAnalysis_h 1

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>
#include "MTAEvent.h"
#include "FillHistBase.h"
#include "TMusunMTAEvent.h"
#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include "TPulse.h"
//#include "MusunNeutronBlockAnalysis.h"
/** MuSun MTA module histogramming many neutron properties.
  * 
  * Does everything with the neutron detectors for the MTA level analysis.
  */
class MusunNeutronAnalysis : public FillHistBase
{
 public:
  MusunNeutronAnalysis(int run=37027);
  //MusunNeutronAnalysis();
  virtual ~MusunNeutronAnalysis();
  int ProcessEntry(MTAEvent *mtaEvent);
  
 private:
  bool isNeutronST;
  bool isNoise;
  bool isNeutron1;
  bool mustopintpc;
  bool isDelayedElectron;
  bool isGammaST;
  bool isNeutronCS;
  bool isGammaCS;

  bool isNeutronCS_1;
  bool isGammaCS_1;

  bool notGamma;
  bool notNeutron;
  bool isNeutronBr;
  bool isGammaBr;

  int i;
  int padHits[48];
  TTree *fmtaNeutronEventTree;
  TBranch *fMusunMTANeutronEventBranch;
  TMusunMTAEvent *fMusunMTANeutronEvent;
   

  TFile *fin,*fin1;
  TH1 *hFitAvg;

  
  TH1 *hTrueAllG[8];
  TH1 *hTrueAllN[8];
  
  TH1 *hFitAllG[8];
  TH1 *hFitAllN[8];
  

  TGraph *grTime;
  TH1 *hBinFill;
  TH1 *hIslandNo;

  TH1 *hTrueTimeAllG[8];
  TH1 *hTrueTimeAllN[8];
  
  TH3 *hTime;
  TH3 *hArea;
  TH3 *hArea_Br;
  TH3 *hArea_Br1;
  TH3 *hPed;
  TGraph *brentGraph;
  TF1 *f;
  TF1 *f_ped;
  int block_channel[1500][8];
  int block, channel, peak;
  double peakTime, halfTime, trueTime;
  double tTimeAllN[8], tTimeAllG[8];
  TCanvas *c1;
  int count;
  double func_min;
  int sample;
  double avg_ped;
  Int_t pSize; //pulse size
  Double_t fArea_low, fTime_low, fPed_low, rArea_low, rTime_low, rPed_low;
  double chiN, chiG, chiN_low, chiG_low, chiN_dof, chiG_dof;


  TH3 *hXRay1;
  TH3 *hXRay2;
  TH3 *hXRay3;
  TH3 *hXRay;

  TH3 *hBinCorrelationN[8];
  TH3 *hBinCorrelationNLow;
  TH3 *hBinCorrelationNHigh;

  TH3 *hBinCorrelationG[8];
  TH3 *hBinCorrelationN5;
  TH3 *hBinCorrelationG5;
  TH3 *hBinCorrelationN10;
  TH3 *hBinCorrelationG10;

  TH3 *hBinCorrelationNTail;
  TH3 *hBinCorrelationNTailLow;
  TH3 *hBinCorrelationNTailHigh;

  TH3 *hResidueBinAllN;
  TH3 *hResidueAdcAllN;
  TH3 *hResidueBinN;
  TH3 *hResidueBinLowN;
  TH3 *hResidueBinHighN;
  TH3 *hResidueBinG;
  TH3 *hResidueAdcG;

  TH3 *hNeutronElectronCS;
  TH3 *hFusionElectronCS;


  TH1 *hNeutronPulseTemplate[8]; 
  TH1 *hGammaPulseTemplate[8]; 

  TH1 *hNeutronPulseTemplateAll[8];
  TH1 *hGammaPulseTemplateAll[8];

  TH1 *hNeutronPulseTemplateLow[8];
  TH1 *hGammaPulseTemplateLow[8];
  TH1 *hNeutronPulseTemplateHigh[8];
  TH1 *hGammaPulseTemplateHigh[8];

  TH1 *hTrueTime[8];
  TH1 *hSingleNeutronPulseShape[8];
  TH1 *hSingleGammaPulseShape[8];

  TH1 *hNeutronSum[8];
  TH1 *hNeutronAvg[8];
  TH1 *hNeutronSigma[8];
  TH1 *hGammaSum[8];
  TH1 *hGammaAvg[8];
  TH1 *hGammaSigma[8];
  TH1 *hNeutronReadSum[8];
  TH1 *hGammaReadSum[8];
  TH1 *hNeutronReadSigma[8];
  TH1 *hGammaReadSigma[8];
  TH1 *hNeutronReadN;
  TH1 *hGammaReadN;
  TH1 *hNeutronN;
  TH1 *hGammaN;

  TH1 *hPseudoTimeBinGammaAll[8];
  TH1 *hPseudoTimeBinNeutronAll[8];
  TH1 *hPseudoTimeBinGammaLow[8];
  TH1 *hPseudoTimeBinNeutronLow[8];
  TH1 *hPseudoTimeBinGamma[8];
  TH1 *hPseudoTimeBinNeutron[8];
  TH1 *hPseudoTimeBinGammaHigh[8];
  TH1 *hPseudoTimeBinNeutronHigh[8];

  TH1 *hNeutronNeutronChiSq;
  TH1 *hNeutronGammaChiSq;
  TH1 *hGammaNeutronChiSq;
  TH1 *hGammaGammaChiSq;

  TH1 *hMuentTime;
  TH1 *hTrueTimeBinNeutron;
  TH1 *hTrueTimeBinGamma;
  TH1 *hTrueBin;

  TH2 *hPeakTimeNSamples;

  TH2 *hTpcWTPulseFitAmplitudeVsPad;
  TH2 *hOneTpcWTPulseFitAmplitudeVsPad;
  
  TH1 *hTpcVetoAmplitude;
  TH1 *hTpcEntAmplitude;
  TH1 *hTpcStopAmplitude;
  
  TH1 *hTpcVetoFusionAmplitude;
  TH1 *hTpcEntFusionAmplitude;
  TH1 *hTpcStopFusionAmplitude;

  TH1 *hTpcVetoHits;
  TH1 *hTpcEntHits;
  TH1 *hTpcStopHits;
  TH1 *hTpcAllHits;
  TH2 *hTpcVetoAmplitudeVsPad1;
  TH2 *hTpcEntAmplitudeVsPad1;
  TH2 *hTpcStopAmplitudeVsPad1;


  TH1 *hTpcVetoHitsVsPad1;
  TH1 *hTpcEntHitsVsPad1;
  TH1 *hTpcStopHitsVsPad1;

  TH1 *hTpcIntegral;
  TH1 *hTimeDifferenceNoMuonStopTpc;
  TH1 *hDriftTime;
  TH1 *hTpcAmp_mta;
  TH1 *hTpcIntegralStop;
  TH1 *hTpcWfd_NPulses;
  //TH1 *hTpcWfd_MustopAmp;
  //TH1 *hTpcAmpNoStop;
  //TH1 *hMuStopDist; 
  //TH1 *hFastMuStop;
  TH1 *hMyMuStop; 

  TH1 *hPulseEnergy;
  TH1 *hPulseAmp;
  TH2 *hPadN; 
  TH2 *hPadMultiplicity;
  TH2 *hPadMultiplicityFusion;
  TH2 *hFitIntegral_1_vs_2;

  TH1 *hTimeDiff_nDet_eSC_all;
  TH1 *hTimeDiff_nDet_eSC_Neutron;
  TH1 *hTimeDiff_muSC_eSC;

  TH3 *hTimeDifferenceNeutronStop;
 
  TH3 *hTpcAmplitude;
  TH3 *hTpcAllAmplitude;
  TH3 *hTpcNeutronAmplitude;
  TH3 *hTpcNeutronStopAmplitude;
  TH3 *hTpcFusionAmplitude;
  TH3 *hTpcFusionStopAmplitude;
  TH3 *hNdetAllPulses;
  TH2 *hTpcFusionAmplitudeVsPad;
  TH2 *hTpcFusionNostopAmplitudeVsPad;
  TH3 *hTpcCaptureAmplitude;
  TH3 *hTpcCaptureStopAmplitude;
  TH3 *hTpcCaptureNostopAmplitude;

  TH3 *hTpcHits;
  //TH3 *hTpcAllHits;
  TH3 *hTpcNeutronHits;
  TH3 *hTpcNeutronStopHits;
  TH3 *hTpcFusionHits;
  TH3 *hTpcCaptureHits;
  

  //PSD histograms
  TH3 *hPSD; 
  TH3 *hPSD_chiSqDiff;
  TH3 *hPSD_chiSqDiff_dof;
  TH3 *hPSD_chiSqDiff_Low;
  TH3 *hPSD_chiSqSum;
  TH3 *hPSD_chiSqRatio;
  TH3 *hPSD_chiSq;
  TH3 *hPSD_chiN;// Neutron chiSq - no cuts
  TH3 *hPSD_chiG;// Gamma chiSq - no cuts
  TH3 *hPSD_chiN_dof;
  TH3 *hPSD_chiG_dof;
  TH3 *hPSD_chiSqGT;
  TH3 *hPSD_ratio; 
  TH3 *hPSD_ratio_neutron;
  TH3 *hPSD_ratio_neutron_g;
  TH3 *hPSD_Neutron;// Neutron chiSq - neutron cut
  TH3 *hPSD_Gamma;
  TH3 *hPSD_Neutron_dof;
  TH3 *hPSD_Gamma_dof;
  

  TH3 *hCaptureCS;
  TH3 *hNeutron;
  TH3 *hGamma;
  TH3 *hNeutronMustop;
 
  TH3 *hCapture_g;
  TH3 *hNeutron_g;
  TH3 *hGamma_g;
  TH3 *hNeutronMustop_g;
  TH3 *hDelay_g;
  TH3 *hFusion_g;

  TH3 *hGammaMustop;
  TH3 *hDelayCS;
  TH3 *hDelayFineCS;

  TH3 *hDelay_e_mu;
  TH3 *hDelay_e_mu_cut;
  TH3 *hDelay_e_mu_cut1;
  TH3 *hDelay_e_mu_cut2;
  TH3 *hFusion;
  
  TH2 *hTotalvsSlow30Chan7;
  TH2 *hTotalvsSlow30HalfChan7;
  TH2 *hTotalvsSlow30Chan7_new;
  TH2 *hTimeDifferenceNeutron_Cut;
  void NewNeutron(TMusunNeutronPulse *raw, int Block_no);
  void SumNeutron(TMusunNeutronPulse *raw);
  void SigmaNeutron(TMusunNeutronPulse *raw);
  void readFadcFile(int run);

  void makeSumHist();
  void makeSigmaHist();

  void makePseudoHist();
  void makeTempHist(int run);
  void makeFitHist();
  void readTempFile(int run);
  void readAvgFile();
  void readSigmaFile();
  void makeTpcHist();
  void makeImpHist();
  void makePSD_hist();

};

#endif

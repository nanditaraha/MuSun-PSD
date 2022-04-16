//pc_timeB, tpc_time;
double tpc_fit;
double tpc_integral, tpc_integral1,tpc_integral2,tpc_amp,tpc_amp1,tpc_amp2,tpc_amp3 ; // store tpc wfd pulses integral & amplitude
double fusionAmp = 0.;
double fusionStopAmp = 0.;
int pulse_energy, pulse_amp;
int padN=0, padV=0, padE=0, padS=0;
int myStop=0, fastStop=0, muStop=0;

Int_t neutronchannel;

int neutronPulses[8]={0};
for (int i=0; i<48; i++) padHits[i]=0;

int padHitsFusion[48]={0};

int vetoHit[48]={1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1};
int entHit[48]= {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
int stopHit[48]={0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0};


int vetoHits=0;
int entHits=0;
int stopHits=0;
double tpc_amp_veto = 0;
double tpc_amp_ent = 0;
double tpc_amp_stop = 0;

int veto=0;
int ent=0;
int stop=0;


const int NDETECTORS  = 8;
float gain[NDETECTORS]={0.0};
float gain_before_38156[NDETECTORS]={7.2,5.9,6.6,7.31,9.1,9.12,4.69,7.59};
float gain_after_38156[NDETECTORS ]={6.12,6.21,6.21,6.18,6.61,7.47,4.32,4.23};
float fusEnergy[NDETECTORS]={0};
if(run<38156)  for (int i = 0; i<8; i++) gain[i] = gain_before_38156[i];    //Before sheilding - initial
 else for (int i = 0; i<8; i++) gain[i] = gain_after_38156[i];//After sheilding

//float fusEnergy[8]={5096.,4944.,5008.,4768.,6488.,4032.,3552.,3528.};//gain*800
for (int i=0; i<8; i++)   fusEnergy[i]=gain[i]*800;
//printf("Run %d, gain %f fus ch 4 %f\n",run, gain[3],fusEnergy[3]);
float highEnergy[8]={25000.,25000.,25000.,12200.,25000.,25000.,13500.,15000.};



//THESE ARE THE NEW RUN 4 VALUES  //
/*
  float slowTot_Pulses[8][4]={{0.29,0.25,0.25,0.25},{0.23,0.24,0.25,0.27},{0.28,0.28,0.3,0.3},{0.27,0.28,0.3,0.3},
        {0.21,0.21,0.22,0.24},{0.21,0.22,0.24,0.24},{0.28,0.28,0.28,0.28},{0.25,0.24,0.23,0.23}}; // Array to find after pulses

  float slowTotNeutron[8][4]={{0.136, 0.124, 0.096, 0.069},{0.12, 0.112, 0.1,0.071},{0.136, 0.122,0.096,0.078},{0.159, 0.117, 0.104,0.115},
        {0.091, 0.083, 0.074,0.060},{0.14, 0.146, 0.104,0.09},{0.161,0.143,0.112,0.124},{0.144, 0.131, 0.116,0.102}};// Array to find neutronshXRay->GetZaxis()->SetRange(8,8)

  double tot_min[8][5]={{470.,870.,1200.,3200.,25000.},{420.,700.,1060.,1600.,20000.},{500.,650.,1690.,1500.,25000.},{550.,1000.,3000.,11100.,15500.},
  {900.,1200.,1600.,2500.,25000.},{1400.,1620.,2300.,3250.,25000.},{520.,900.,2000.,13500.,17300.},{510.,810.,1200.,1600.,15000.}};
*/
//TMuEntrance *muent = onEvent.at(muN).GetMuEntrance();
/*
  float slowTot_Pulses[8][2]={{0.262,0.228},{0.222,0.185},{0.281,0.222},{0.262,0.23},
                              {0.222,0.222},{0.262,0.217},{0.262,0.222},{0.262,0.222}}; // Array to find after puls\
es                                                                                                                  

  float slowTotNeutron[8][2]={{0.103, 0.077},{0.079,0.05},{0.122,0.085},{0.122,0.117},
                              {0.079,0.058},{0.127,0.095},{0.122,0.118},{0.109,0.105}};

 double tot_min[8][3]={{1300.,3600.,25000.},{1600.,3700,25000.},{1300.,3700,25000.},{1500.,3700.,16000.},
        {1600.,3700.,25000.},{1750.,3700.,25000.},{1400.,3700.,15000.},{1350.,3700.,15500. }};


*/


float slowTot_Pulses[8][3]={{0.222,0.18,0.15},{0.18,0.16,0.15},{0.215,0.2,0.2},{0.2,0.18,0.165},
			    {0.177,0.16,0.15},{0.213,0.2,0.19},{0.212,0.19,0.18},{0.2,0.18,0.16}}; // Array to find after pulses

float slowTotNeutron[8][3]={{0.103, 0.077,0.077},{0.09,0.05,0.05},{0.122,0.085,0.087},{0.122,0.117,0.117},
			    {0.079,0.058,0.059},{0.127,0.127,0.126},{0.122,0.118,0.118},{0.14,0.105,0.105}};

double tot_min[8][4]={{1300.,3000.,6000.,25000.},{1600.,3700,6000.,25000.},{1300.,3700,6000.,25000.},{1500.,2800.,6000.,15000.},
		      {1600.,3000.,6000.,25000.},{1750.,2800.,6000.,25000.},{1400.,2800.,6000.,15000.},{1350.,2800.,6000.,15000. }};

// TMusunMTAEvent *fMusunMTANeutronEvent = mtaEvent->GetMusunMTAEvent();
TRunInfoEvent *fRunInfoEvent = mtaEvent->GetRunInfoEvent();
//printf("\n\n");


fMusunMTANeutronEventBranch->SetAddress(&fMusunMTANeutronEvent);
fMusunMTANeutronEvent->Clear();
TMusunEvent *event = mtaEvent->GetMusunEvent();
TTPCEventProperties *tpcEventProperties = mtaEvent->GetTPCEventProperties();

Int_t evnum = event->GetEventNumber();
TMuEntrance *muent = event->GetMuEntrance();
Double_t muEntTime = muent->GetT() ;
//printf("------- MusunNeutronAnalysis Process Entry *** Muent time -------------\n");
//hMuentTime->Fill(muEntTime);
//CODE TO GENERATE AUTO CORRELATION FOR MUONS
Int_t n_muent = muent->GetNmusc();

//printf("# of Muons entered %d\n",n_muent);
/*
  for (int i = 0; i<n_muent; i++)
    {
       TMuEntrance *muent1 = event->GetMuEntrance();
       printf("Muon time %d\n",muEntTime);
       for (int i = 0; i<n_muent; i++)
        {
	   Double_t muEntTime1 = muent1->GetT();
	      printf("Correlated Muon time %d\n",muEntTime1);
	         hMuonAutocorrelation->Fill(muEntTime-muEntTime1);
		      }
		      }*/
//TClonesArray &wC = *(event->GetMusunWfdPulse());
TClonesArray &nC = *(event->GetMusunNeutronPulse());
TClonesArray &eC = *(event->GetGondCaen());

mustopintpc = false;

//  if (muent->HasMuStopInTPC())fastStop = 1;
mustopintpc = tpcEventProperties->HasMuonStop() ;
//hMuStop->Fill(mustopintpc);

//CODE FOR TPC PULSES(just to check Drift Time tpc pulses created wrt to mu ent)
/*
  vector<TTPCGenericPulse*>& tpcVector = 
    mtaEvent->GetMusunEvent()->GetTPCMiniPulses();
  
  Int_t nwfdpulses = tpcVector.size();
  `// hTpcWfd_NPulses->Fill(nwfdpulses);
  
  
  for(int y=0; y<nwfdpulses; y++){    
    TTPCGenericPulse* tpc = tpcVector[y];
    padN = tpc->GetPadOneIndexed();
    padHits[padN-1] = padHits[padN-1]+1;
    tpc_amp = tpc->GetAmplitude();
    tpc_time = tpc->GetCenterTime();
    hDriftTime->Fill(tpc_time - muEntTime);
    hTpcWTPulseFitAmplitudeVsPad->Fill(tpc->GetPadOneIndexed(),tpc_amp); 
    
  }
  //CODE TO STUDY AND DERIVE MY OWN MUSTOP DEFINITION
  for (int i=1; i<49; i++){
         
    for(int y=0; y<nwfdpulses; y++){ 
      TTPCGenericPulse* tpc = tpcVector[y];
      if (tpc->GetPadOneIndexed()==i)
      {
        tpc_amp_veto = tpc_amp_veto + vetoHit[i-1]*tpc->GetAmplitude();
	  tpc_amp_ent = tpc_amp_ent + entHit[i-1]*tpc->GetAmplitude();
	    tpc_amp_stop = tpc_amp_stop + stopHit[i-1]*tpc->GetAmplitude();
	    }
      
    }
    
    veto = veto + vetoHit[i-1]*padHits[i-1];
    ent = ent + entHit[i-1]*padHits[i-1];  
    stop = stop + stopHit[i-1]*padHits[i-1];
    if (padHits[i-1]==1){
      for(int y=0; y<nwfdpulses; y++){ 
      
      TTPCGenericPulse* tpc = tpcVector[y];
      if (tpc->GetPadOneIndexed() == i) 
        {
	    tpc_amp = tpc->GetAmplitude();
	        hOneTpcWTPulseFitAmplitudeVsPad->Fill(tpc->GetPadOneIndexed(),tpc_amp);
		  }
      }//loop for # of tpc pulses
      
    }//loop of hits = 1
  }// loop over all pads
    
  hTpcVetoHits->Fill(veto);
  hTpcEntHits->Fill(ent);
  hTpcStopHits->Fill(stop);

  hTpcAllHits->Fill(nwfdpulses);

  hTpcHits->Fill(veto,ent,stop );
    
  hTpcVetoAmplitude->Fill(tpc_amp_veto);
  hTpcEntAmplitude->Fill(tpc_amp_ent);
  hTpcStopAmplitude->Fill(tpc_amp_stop);
  
  hTpcAmplitude->Fill(tpc_amp_veto,tpc_amp_ent,tpc_amp_stop );

  for (int i=1; i<49; i++)
    hPadMultiplicity->Fill(padHits[i-1],i); 

  if(tpc_amp_veto<12 && tpc_amp_ent<25 && tpc_amp_stop>60 && tpc_time - muEntTime>4000. && tpc_time - muEntTime<14000. ) 
    {
      myStop=1;
      mustopintpc = true;
    }
  hMyMuStop->Fill(myStop);

  muStop = myStop + 2*fastStop;

  for ( int gc = 0 ; gc < event->GetNgondCaen(); gc++ ) {
    TGondOnly *gon = (TGondOnly*)(eC[gc]);
    eSC_time = gon->GetT();
    hTimeDiff_muSC_eSC->Fill(muEntTime-eSC_time);
  }

  // hMuStopDist->Fill(muStop);
  


  
    for(int y=0; y<event->GetNwfdpulse(); y++){     
      TMusunWfdPulse *tpc = (TMusunWfdPulse*)(wC[y]);
      //printf("No. of tpc wfd pulses %d\n",event->GetNwfdpulse());
      tpc_timeC = tpc.GetPulseCenterTime();
      tpc_timeB = tpc.GetBlockTime();
      tpc_fit = tpc.GetFitTime();
      //printf("\t Tpc Block time %f\t Pulse Centre time %f\t Fit time %f\t MUSC time %f\t Tpc_Time & MUSC time diff %f\n\n",tpc_timeB, tpc_timeC, tpc_fit, muEntTime,tpc_timeB - muEntTime);
             
      tpc_integral = tpc.GetFitIntegral();
      // printf("Pulse energy %d\t",tpc_integral);
      tpc_amp = tpc.GetAmplitude();
      padN = tpc.GetPadOneIndexed();
      if (y<2)
      {
      if (y==0) 
        {
	    tpc_integral1 = tpc.GetFitIntegral();
	        // printf("Pulse energy for pulse1 %d\t",tpc_integral1);
		  }
		  if (y==1) 
		    {
		        tpc_integral2 = tpc.GetFitIntegral();
			    //printf("Pulse energy for pulse2 %d\n",tpc_integral2);
			      }
			      hFitIntegral_1_vs_2->Fill(tpc_integral1,tpc_integral2);
			       }
			              pulse_energy = tpc.GetPulseEnergy();
				             pulse_amp = tpc.GetPulseAmplitude();
					            hTimeDifferenceMuontpc.Fill(tpc_timeB - muEntTime);
						           hTpcIntegral->Fill(tpc_integral);
							          hPulseEnergy->Fill(pulse_energy);
								         hPulseAmp->Fill(pulse_amp);
									        hTpcAmp_mta->Fill(tpc_amp);
										       hPadN->Fill(padN);     
										              if(mustopintpc){
											       hTimeDifferenceMuonStoptpc.Fill(tpc_timeB - muEntTime);
											        hTpcIntegralStop->Fill(tpc_integral);
												 hTpcWfd_MustopAmp->Fill(tpc_amp);
												        }
													       // printf(" :::::::::Mustops in tpc at mta = %d\n", muent->HasMuStopInTPC() ); 
													              
													              if(!mustopintpc){
														       hTimeDifferenceNoMuonStopTpc->Fill(tpc_timeB - muEntTime);
														        hTpcAmpNoStop->Fill(tpc_amp);
															 }
															        
    }
      
  }//Tpc Wfd loop ends ------------CODE FOR TPC ENDS HERE
*/

//Mustop condition from tpc
//if(tpcEventProperties->HasMuonStop()) mustopintpc=true;

/***************NEUTRON LOOP STARTS*********************/

//  printf("#of electrons in one muent %d\n",event->GetNgondCaen());
// printf("Moun Time at MTA level %f\n",muEntTime );
// printf("# of neu %d\n",event->GetNmusunneutronpulse());
//vector<const TMusunNeutronPulse*>& neuVector = muonEvent.at(muN).GetNeutronPulses();
// vector<TMusunNeutronPulse*> neuVector =  mtaEvent->GetMusunEvent()->GetNeutronPulseVector();
// vector<const TMusunNeutronPulse*> neuVector = *mtaEvent->GetNeutronPulses();
int a;
int triggerSample1=0,triggerSample2=0 ;
//double sum[8][100]={0};
//  for (int i=0; i<100; i++)
//printf("Chan %d bin %d Avg %f and sigma %f\n",neutronchannel,i,avg[neutronchannel][i],sigma_bin[neutronchannel][i]);
//This is the main neutron loop that does all things
neutron_multi = event->GetNmusunneutronpulse();
for(int y=0; y<event->GetNmusunneutronpulse(); y++){
  //printf("Pulse main ////////////////\n");
  TMusunNeutronPulse *neu = (TMusunNeutronPulse*)(nC[y]);
  //      const TMusunNeutronPulse *neu = neuVector[i];
  neutronchannel = neu->GetChannel();
  const vector <UShort_t> &samples = neu->GetSamples();
  underflow_pulse = false;
  for (int i=0; i<samples.size(); i++){
    if (samples[i]==8191) {
      underflow_pulse = true;
      break;
    }
  }

  //printf("Dynamic Range %f\n***********\n",TMusunNeutronPulse::GetMaxEnergyForPSD(1));
  //printf("Neutron Channels are %d\n",neutronchannel);
  neuN++;
  /******** FINDING THRESHOLDS
	          Int_t thresSample = -999;
		        const vector <UShort_t> &samples = neu->GetSamples();
			      for (int k=0; k<samples.size(); k++){
			            //Threshold of 1420 for run 38854 was found from its odb file//
				          if(samples[k]>1420) {
					        thresSample = k;
						      break;
						            }
							          }
								        
  */
  //Without applying gains!!!

  double neutronTotalintegral = neu->GetTotalIntegral();
  double neutronSlowintegral30 = neu->GetSlowIntegral30();

  // With gain calibrations!!
  /* double neutronTotalintegral = neu->GetTotalIntegral()/gain[neutronchannel]; 
       double neutronSlowintegral30 = neu->GetSlowIntegral30()/gain[neutronchannel];
  */
  double slowTot = neutronSlowintegral30/neutronTotalintegral;
  double neutronpeaktime1,neutronpeaktime2;

  double neutronpeaktime = neu->GetPeakTime();
  int neutronblocktime = neu->GetBlockTime();
  double neutronpeaktimeOnly; //= neu->GetPeakTimeOnly();

  //printf("Block no %d\t channel %d\t peak slip %d MTA LEVEL\n",event->GetBlockNumber(),neutronchannel,block_channel[event->GetBlockNumber()-1][neutronchannel]);
  // printf("Time diff b/w moun & neutrons %f***** MTA LEVEl\n",neutronpeaktime-muEntTime);
  //neutronPulses[neutronchannel]++;

  /****** THE CORRECTION TO FADC CLOCK SLIP IS APPLIED HERE*****/
  //if(block_channel[event->GetBlockNumber()][neutronchannel]==1)
  neutronpeaktime1=neutronpeaktime-block_channel[event->GetBlockNumber()][neutronchannel]*5.88;
  //if(block_channel[event->GetBlockNumber()][neutronchannel]==2)
  // neutronpeaktime1=neutronpeaktime;
  //if(block_channel[event->GetBlockNumber()][neutronchannel]==0) continue;

  /*** CORRECTION FOR OVERLOADS***/
  if(neutronTotalintegral>9000.)
    neutronpeaktime=neutronpeaktime-(-8.4+0.000914*neutronTotalintegral-1.184e-8*neutronTotalintegral*neutronTotalintegral);
  //The constant term -8.4 is the extrapolation to zero subtracted from the constant p0.

  /********************* CODE FOR NEUTRON/GAMMA CUTS**********/

  isNeutronST = false; //ST stands for slow & total integral method - did not remove as it is required to make templates
  isGammaST = false;
  isFusionCS = false;
  isFusionBr = false;

  if(neu->GetPeakTime()*170./1000.-neu->GetBlockTime()<30 && samples.size()>45) //takes care of all noisy pulses - peaks > 30 or so are ignored
    {
      for (int i=0; i<4; i++)
	{
	  if (slowTot<slowTot_Pulses[neutronchannel][i]
	                      && slowTot>slowTotNeutron[neutronchannel][i] && neutronTotalintegral>tot_min[neutronchannel][i]
	      && neutronTotalintegral<tot_min[neutronchannel][i+1])
	    isNeutronST=true;
	  //printf("****************Neutron %d or Gamma %d************\n", isNeutronST, isGammaST);

	}
      for (int i=0; i<4; i++)
	{
	  if (slowTot<slowTot_Pulses[neutronchannel][i] && !isGammaST
	                      && slowTot<slowTotNeutron[neutronchannel][i]-0.06 && neutronTotalintegral>tot_min[neutronchannel][i]
	      && neutronTotalintegral<tot_min[neutronchannel][i+1])
	    isGammaST=true;
	}
    }



  /****This is to call to the function that corrects/defines the peak, half and true times - function transported from mu level*****/

  iHi=55; //High Fit range initialized here so that every neutron pulse has a default max. fit range of 30
  pSize = samples.size()-1;

  //NOTE: While making configuration text files for FADC offsets I comment NewNeutron call to spped up -
  //MUST be uncommented after the configuration files made
  if(neu->GetPeakTime()*170./1000.-neu->GetBlockTime()<pSize &&
     neu->GetPeakTime()*170./1000.-neu->GetBlockTime()<iHi && muent->HasBestEntrance()) //Prevents noise
    {
      //SumNeutron(neu);
      NewNeutron(neu, event->GetBlockNumber());
    }

  //Electron loop starts - to find delayed electons

  for ( int gc = 0 ; gc < event->GetNgondCaen(); gc++ ) {
    //printf("Electron No. %d\n",event->GetNgondCaen());
    TGondOnly *gon = (TGondOnly*)(eC[gc]);
    //tdiff_eSC_nue = eSC_time - neutronpeaktime;
    double eSC_time = gon->GetT();
    if (eSC_time > neutronpeaktime + 200  && eSC_time < neutronpeaktime + 5000 && !isDelayedElectron)
      isDelayedElectron = true;
    //printf("Delayed e found ***************\n");
    if (eSC_time < neutronpeaktime + 200  && eSC_time > neutronpeaktime + 5000 && isDelayedElectron)
      isDelayedElectron = false;

  }//electron loop ends

  if (isDelayedElectron && mustopintpc && isNeutronCS && muent->HasBestEntrance())
    isFusionCS = true;

  if (isDelayedElectron && mustopintpc && isNeutronBr && muent->HasBestEntrance())
    isFusionBr = true;

  neu->fIsFusionNeutronChiSquare = isFusionCS;
  neu->fIsFusionNeutronTemplateSlowTotal = isFusionBr;

  if (event->GetNgondCaen()>0 && isNeutronCS && mustopintpc && neu->GetTotalIntegral()<highEnergy[neutronchannel] ){
    //hNeutronElectronCS->Fill(neutronpeaktime-muent->GetT(),neu->GetTotalIntegral(),neutronchannel);
    if (neu->GetTotalIntegral()<fusEnergy[neutronchannel])
      hFusionElectronCS->Fill(neutronpeaktime-muent->GetT(),neu->GetTotalIntegral(),neutronchannel);
    //printf("N time %f, energy %f, channel %d\n", neutronpeaktime,neu->GetTotalIntegral(),neutronchannel);
  }


 }//End of Neutron main loop
fmtaNeutronEventTree->Delete("all");
return 0;
} //End of process entry

/*-- All FUNCTIONS DEFINED HERE -------------------------------------------------*/

void MusunNeutronAnalysis::readFadcFile(int run){
  int data_block;
  char name[1000];
  sprintf(name,"/scratch/01490/nanditar/offset_baxkup/ds_mixed/offset_%d.txt",run);
  ifstream infile;
  infile.open(name);
  for(int i=0;i<1500;i++){
    for(int j=0;j<8;j++) {
      block_channel[i][j]=0;
    }
  }

  while (!infile.eof())
    {
      infile>>data_block;
      block=data_block;

      infile>>data_block;
      peak=data_block;

      infile>>data_block;
      channel=data_block;

      block_channel[block-1][channel]=peak;
    }

  infile.close();
}

static double pedestal(int nsamples, const vector <UShort_t> &samples)
{
  double d = 0;
  int ns = 0;

  for(int i = 0; i < 6; i++) {
    if(i < nsamples) {
      d += samples[i];
      ns++;
    }
  }

  return d / ns;
}

double pedestalDouble(int nsamples, double *samples)
{

  int max = 0;
  int maxpoint = 0;
  if(samples[0]>max)
    {
      max=samples[0];
      maxpoint=0;
    }
  if(samples[1]>max)
    {
      max=samples[1];
      maxpoint=1;
    }
  if(samples[2]>max)
    {
      max=samples[2];
      maxpoint=2;
    }
  if(maxpoint==0)
    return samples[1]>samples[2]?samples[1]:samples[2];
  if(maxpoint==1)
    return samples[0]>samples[2]?samples[0]:samples[2];
  if(maxpoint==2)
    return samples[0]>samples[1]?samples[0]:samples[1];

}

static bool isOverflow(int nsamples,const vector <UShort_t> &samples)
{
  for(int i = 0; i < nsamples; i++) {
    if(samples[i] & 0x1000) return true;
  }

  return false;
}

static int peakBinShort(int nsamples, const vector <UShort_t> &samples)
{
  int peakBin = 0;
  short peak = -1;
  for(int i = 0; i < nsamples; i++) {
    if(samples[i] > peak) {
      peak = samples[i];
      peakBin = i;
    }
  }
  return peakBin;
}

static int peakBinDouble(int nsamples, double *samples)
{
  int peakBin = 0;
  double peak = -1;
  for(int i = 0; i < nsamples; i++) {
    if(samples[i] > peak) {
      peak = samples[i];
      peakBin = i;
    }
  }
  return peakBin;
}

// This is my attempt to take the bin corresponding to 50% of the peak sample instead of the bin for peak sample

static int halfBinDouble(int nsamples, double *samples, double ped, double peak)
{
  /*To find the halfway bin, we need the half the value of peak sample wrt ped i.e. 0.5(peak-ped)
    But ped has to be added to this half value as integralFromPeakTime function subtracts the ped ultimately.
    We end up with 
     0.5(peak-ped)+ ped = 0.5(peak+ped)
  */

  double halfway = 0.5*(peak+ped);
  double half = -1;
  int halfBin = 0;
  for(int i = 0; i < nsamples; i++)
    {
      if(samples[i] > halfway) {
	halfway=samples[i];
	halfBin = i;
	break;
      }
    }
  // printf("Block no. = % d channel %d block_channel = %d\n",Block_no,channel, block_channel[Block_no-1][channel] );

  return halfBin;
}


static double integralFromPeakTime(int nsamples, double *samples, double binWidth,
				   int pb, double ped, double t1, double t2)
{
  int bin1 = pb + (int) (t1/binWidth);
  int bin2 = pb + (int) (t2/binWidth);

  if(bin1 < 0) bin1 = 0;
  if(bin2 >= nsamples) bin2 = nsamples - 1;

  double integral = 0;
  for(int i = bin1; i <= bin2; i++) {
    integral += samples[i];
  }
  return 1.0*(integral - ped*(bin2-bin1+1));
}

static double integralFromTrueTime(int nsamples, double *samples, double binWidth,
				   double pb, double ped, double t1, double t2)
{
  int bin1 = (int)(pb + t1/binWidth);
  int bin2 = (int)(pb + t2/binWidth);

  if(bin1 < 0) {
    //printf("Off left edge\n");
    bin1 = 0;
  }
  if(bin2 >= nsamples) {
    //printf("Off right edge\n");
    bin2 = nsamples - 1;
  }

  double integral = 0;
  for(int i = bin1; i <= bin2; i++) {
    integral += samples[i];
  }
  return 1.0*(integral - ped*(bin2-bin1+1));
}


void MusunNeutronAnalysis::makeTempHist(int run){
  for(int i=0;i<8;i++){
    char name[80], title[120];

    // for 2nd run of MTA to make template histograms
    sprintf(name, "hNeutronPulseTemplateAll%d",i+1);
    sprintf(title, "Neutron Pulse Template channel %d",i+1);
    hNeutronPulseTemplateAll[i] = new TH1D(name,title,800, 0,80);

    sprintf(name, "hGammaPulseTemplateAll%d",i+1);
    sprintf(title, "Gamma Pulse Template channel %d",i+1);
    hGammaPulseTemplateAll[i] = new TH1D(name,title, 800, 0,80);

    sprintf(name, "hTrueAllG%d",i+1);
    sprintf(title, "Gamma true time channel %d",i+1);
    hTrueAllG[i] = new TH1D(name,title,10, 0,1);

    sprintf(name, "hTrueAllN%d",i+1);
    sprintf(title, "Neutron true time channel %d",i+1);
    hTrueAllN[i] = new TH1D(name,title,10, 0,1);
  }

  // 2nd run of mta to save true times and make templates
  char name[30];
  //sprintf(name,"hist_try%d.root",run);
  //fin = TFile::Open(name);
  fin = TFile::Open("hist_try41061.root");
  printf("NOt called:: yet Run number is %d\n",run);

  for(int i=0;i<8;i++){
    char name1[81];
    sprintf(name1,"hTrueTimeGammaAll%d",i+1);
    hTrueAllG[i] = (TH1*)fin->Get(name1);
    sprintf(name1,"hTrueTimeNeutronAll%d",i+1);
    hTrueAllN[i] = (TH1*)fin->Get(name1);

  }
}



void MusunNeutronAnalysis::NewNeutron(TMusunNeutronPulse *neu, int Block_no){

  isNeutronCS = false; //CS stands for chi squared for pulse temp fit
  notGamma = false; // this has an additional cut based on gamma chi2 distribution
  isGammaCS = false;

  isNeutronCS_1 = false; //stricter definitions based on isNeutronCS and notGamma
  isGammaCS_1 = false;

  isNeutronBr = false; //Br stands for Brent area from fit
  isGammaBr = false;

  isNoise = false;
  nOver=0;

  countN++;
  int t0 = neu->GetBlockTime();
  const vector <UShort_t> &samples = neu->GetSamples();

  nsamples = samples.size();

  int channel = neu->GetChannel();

  double binWidth = 1000.0/170; //clock spped is 170 MHz

  double ped = neu->GetPed();
  //printf("Pedestal %f\n",ped);
  //double ped = pedestal(nsamples, samples);
  //printf("Pedestal %f\n",ped);

  //int pb = peakBinShort(nsamples, samples);
  int peakSample = neu->GetPeakSample();

  /************ To deal with pulses after resampling and splining ************/
  int resampleFactor = 20;
  int numResamples = nsamples*resampleFactor;
  double resamples[numResamples];
  resamplePulse(samples, nsamples, resamples, numResamples);

  double ped2 = pedestalDouble(numResamples, resamples);

  int pb2 = peakBinDouble(numResamples, resamples);

  int peakSample2 = resamples[pb2];

  int hb = halfBinDouble(numResamples, resamples, ped2, peakSample2 );

  if(block_channel[Block_no][channel]==1)
    peakTime = ((double)t0 + ((double)pb2/resampleFactor))*binWidth-binWidth;
  else  peakTime = ((double)t0 + ((double)pb2/resampleFactor))*binWidth;

  if(block_channel[Block_no][channel]==1)
    halfTime = ((double)t0 + ((double)hb)/resampleFactor)*binWidth-binWidth;
  else halfTime = ((double)t0 + ((double)hb)/resampleFactor)*binWidth;

  /************ To deal with pulses without splining and resampling for pseudo times************/
  int resampleFactor1 = 1;
  int numResamples1 = nsamples*resampleFactor1;
  double resamples1[numResamples1];
  resamplePulse(samples, nsamples, resamples1, numResamples1);
  //printf("Resample factor = 1 ....\n");

  int a;
  int pb1 = peakBinDouble(numResamples1, resamples1);
  if(block_channel[Block_no][channel]==1)
    pb1=pb1+1;
  else pb1 = pb1;

  double ptb;
  // if(neu->GetTotalIntegral()>2000 && neu->GetTotalIntegral()<7000){
  //if(isGammaST)
  ptb = pseudoBinDouble(numResamples1, resamples1);
  //}
  //h1->GetXaxis()->SetRangeUser(ptb-0.25,ptb);
  //if(neu->GetTotalIntegral()>2000 && neu->GetTotalIntegral()<7000)
  // tTime1= h1->GetBinContent(h1->FindBin(ptb,0,0));//This is the time within the peak bin for a correction
  double pseudoTime = ((double)pb1/resampleFactor1 + (double)ptb/resampleFactor1)*binWidth;
  double peakTime1 =  ((double)pb1/resampleFactor1)*binWidth;
  /*
  if(block_channel[Block_no][channel]==1)
    trueTime = ((double)t0+tTime1[channel])*binWidth + peakTime1 - binWidth; //true time as it takes the fluctuations within the bin (resoles a clock tick) + the peak bin
  else trueTime = ((double)t0+tTime1[channel])*binWidth + peakTime1;
  */
  double tTime;


  /////////////// 2nd MTA run ----- never delete if you want to make templates very imp to make them//////////////////////////

  //This part reads pseudo times to give true times and creates pseudo time histos - Comment for 3rd run and use for
  //1st and 2nd run as commented
  /*
  if(neu->GetTotalIntegral()>2000 && neu->GetTotalIntegral()<7000 && samples.size()>=60){
    if(isNeutronST && samples.size()>=60){
      tTimeAllN[channel] = hTrueAllN[channel]->GetBinContent(hTrueAllN[channel]->FindBin(ptb)); //2nd MTA
      //hPseudoTimeBinNeutronAll[channel]->Fill(ptb); //for 1st run only
    }
  
    if(isGammaST && samples.size()>=60){
      tTimeAllG[channel] = hTrueAllG[channel]->GetBinContent(hTrueAllG[channel]->FindBin(ptb));//2nd MTA
      //hPseudoTimeBinGammaAll[channel]->Fill(ptb);//for 1st run only
    }
    
  }
  
  
  ////// This part creates templates - 2nd run.. Comment for 3rd run of mta//////  
  
  int sum_samples_NeutronAll[8]={0};
  int sum_samples_GammaAll[8]={0};

  //const vector <UShort_t> &samples = neu->GetSamples();
  //sum_samples_Neutron value ranges from 14000 to 21000
  if(neu->GetTotalIntegral()>2000 && neu->GetTotalIntegral()<7000 && isNeutronST && samples.size()>=60){
    for(int i=0;i<samples.size()-1;i++)
      sum_samples_NeutronAll[channel]+=samples[i]-neu->GetPed();
  }
  if(neu->GetTotalIntegral()>2000 && neu->GetTotalIntegral()<7000 && isGammaST && samples.size()>=60){
    for(int i=0;i<samples.size()-1;i++)
      sum_samples_GammaAll[channel]+=samples[i]-neu->GetPed();
}
     
  if(neu->GetTotalIntegral()>2000 && neu->GetTotalIntegral()<7000 && isNeutronST && samples.size()>=60){
    for(int i=0;i<samples.size()-1;i++)
      hNeutronPulseTemplateAll[channel]->Fill(i+(14-pb1-tTimeAllN[channel]),(double)(samples[i]-neu->GetPed())/sum_samples_NeutronAll[channel]);

  }
  
  if(neu->GetTotalIntegral()>2000 && neu->GetTotalIntegral()<7000 && isGammaST && samples.size()>=60){
    for(int i=0;i<samples.size()-1;i++)
      hGammaPulseTemplateAll[channel]->Fill(i+(14-pb1-tTimeAllG[channel]),(double)(samples[i]-neu->GetPed())/sum_samples_GammaAll[channel]);
  }

  */

  //End of template creation 2nd run- Rest is only for 3rd run to fit templates

  //FOR CHECKING SINGLE PULSES and FITTING WITH TEMPLATES //
  // MAKE WITH EVENTDISPLAY after make clean////

  //------------------Neutrons fitted with Neutron templates---------------//
  //----------------------------------------------------------------------//

  char name[82];
  pSize = samples.size()-1;
  Ped=neu->GetPed();
  //To get rid of noises
  if(fTime>pSize && fTime>iHi) return;
  hSingleNeutronPulseShape[neu->GetChannel()]->Add(hSingleNeutronPulseShape[neu->GetChannel()],-1);

  for (int i=0; i<samples.size(); i++){
    hSingleNeutronPulseShape[neu->GetChannel()]->Fill(i,samples[i]);
  }

  if(iHi>pSize) iHi=pSize;

  hPulse=hSingleNeutronPulseShape[neu->GetChannel()];

  for (int i=iLo; i<iHi; i++)
    if(hPulse->GetBinContent(i)>4095) nOver++;

  double resN=0, resN1=0, resN2=0, resG=0, resG1=0, resG2=0;

  hFitAllN[neu->GetChannel()]->Scale(1/(hFitAllN[neu->GetChannel()]->Integral(10*iLo,10*iHi)));
  h3=hFitAllN[neu->GetChannel()];


  f_avg= new TF1("f_avg",fitTemp,iLo,iHi,3);
  //f_avg->SetLineColor(kYellow);
  Area = hSingleNeutronPulseShape[channel]->Integral()-pSize*neu->GetPed();
  if(Area<0) {
    Area = hSingleNeutronPulseShape[channel]->Integral(iLo,iHi)-(iHi-iLo)*neu->GetPed();
    hFitAllN[neu->GetChannel()]->GetXaxis()->SetRange(iLo*10,iHi*10);
    hFitAllN[neu->GetChannel()]->Scale(1/(hFitAllN[neu->GetChannel()]->Integral(iLo*10,iHi*10)));
    h3=hFitAllN[neu->GetChannel()];

  }

  f_avg->SetParNames("Area", "Peak time","Pedestal");
  brent_count=0;


  min((double)pb1);//Original function with weights=1
  chiN=BrentMin();
  neu->SetChiN(chiN);
  //printf("Neutron chiSq set %f \t chin %f Neutron %d gamma %d\n",neu->GetChiN(),chiN, isNeutronST, isGammaST);

  chiN_dof=chiN/(iHi-iLo-nOver-1);


  fAreaN = AreaBr;
  f_avg->SetParameters(AreaBr,fTime,Ped);
  //IMP - Use this section for finding chi2, but for residue plots neutrons should be fitted to neutron templates only & so is true for gammas
  //This fit is done again to fit neutrons with neutron templetes below

  // if(isNeutronST && samples.size()>55 && neu->GetOverflow()==0){

  //hNeutronElectronCS->Fill(neu->GetPeakTime(),neu->GetTotalIntegral(),channel);
  if(neu->GetTotalIntegral()>2000 && neu->GetTotalIntegral()<7000 && isNeutronST && samples.size()>=60){
    for (int i=1; i<samples.size()-1; i++){
      resN = hPulse->GetBinContent(i)- f_avg->Eval(hPulse->GetBinCenter(i)); //residue ith bin
      //if((i+15-pb1)<samples.size()){
      resN1 = hPulse->GetBinContent(i)- f_avg->Eval(hPulse->GetBinCenter(i)); //residue ith bin wrt to pb1
      resN2 = hPulse->GetBinContent(i+1)- f_avg->Eval(hPulse->GetBinCenter(i+1));//residue i+1th bin
      //}

      double sigmaN, sigmaN1;

      sigmaN = sqrt(hNeutronReadSigma[channel]->GetBinContent(i+(15-pb1)));
      sigmaN1 = sqrt(hNeutronReadSigma[channel]->GetBinContent(i+(16-pb1)));

      //hNeutronSigma[channel]->Fill(i+(14-pb1),resN1*resN1/hNeutronReadN->GetBinContent(channel+1));//2nd run after making templates

      //if(i<82 && chiN<3700){
      if(hNeutronReadN->GetBinContent(channel+1)>0 && TMath::Abs(resN1/sigmaN)<=5.00){//last run
	hBinCorrelationN[channel]->Fill(resN1/sigmaN,resN2/sigmaN1,i+14-pb1); //1st make hNeutronSigma -last run
	hResidueBinN->Fill(i+14-pb1,resN/sigmaN,channel);//i+14-pb1 aligns all peaks on bin 14 - last run
      }
      //printf("Bin %d, peak %d, res1 %f, bin content %f, channel %d, func %f, sigma %f res/sigma %f\n", i,pb1,resN1, hPulse->GetBinContent(i), channel, f_avg->Eval(i) , sigmaN, resN1/sigmaN);
      //}

    }
  }




  //------------------Pulses fitted with Gamma template------------------//
  //----------------------------------------------------------------------//

  hSingleGammaPulseShape[channel]->Add(hSingleGammaPulseShape[channel],-1);
  for (int i=0; i<samples.size(); i++){
    hSingleGammaPulseShape[channel]->Fill(i,samples[i]);

  }

  //hPulse=hSingleGammaPulseShape[neu->GetChannel()];
  //hFitG[channel] = (TH1*)fin1->Get("hGammaPulseTemplate[channel]");
  hFitAllG[neu->GetChannel()]->Scale(1/(hFitAllG[neu->GetChannel()]->Integral(10*iLo,10*iHi)));
  h3=hFitAllG[neu->GetChannel()];

  Area = hSingleGammaPulseShape[channel]->Integral()-pSize*neu->GetPed();

  if(Area<0){
    Area = hSingleGammaPulseShape[channel]->Integral(iLo,iHi)-(iHi-iLo)*neu->GetPed();
    hFitAllG[neu->GetChannel()]->GetXaxis()->SetRange(iLo*10,iHi*10);
    hFitAllG[neu->GetChannel()]->Scale(1/(hFitAllG[neu->GetChannel()]->Integral(iLo*10,iHi*10)));
    h3=hFitAllG[neu->GetChannel()];
  }

  //printf("pb1 %d\n",pb1);
  min((double)pb1);
  //printf("Later: pb1 %d\n",pb1);
  brent_count=0;

  chiG=BrentMin();
  //printf("Gamma chiSq %f\n",chiG);
  chiG_dof=chiG/(iHi-iLo-nOver-1);

  fAreaG = AreaBr;
  f_avg->SetParameters(AreaBr,fTime,Ped);
  f_avg->SetLineColor(kBlue);

  //if(isGammaST && samples.size()>55 && neu->GetOverflow()==0){
  if(neu->GetTotalIntegral()>2000 && neu->GetTotalIntegral()<7000 && isGammaST && samples.size()>=60){
    //printf("Pulse area %f, Brent Area %f, chiN %f, chiG %f Overflow %d\n",neu->GetTotalIntegral(), brentArea, chiN, chiG, neu->GetOverflow());

    hCaptureCS->Fill(neu->GetPeakTime(),neu->GetTotalIntegral(),channel);
    for (int i=1; i<samples.size()-1; i++){
      resG = hPulse->GetBinContent(i)- f_avg->Eval(hPulse->GetBinCenter(i));
      //if((i+15-pb1)<samples.size()){
      resG1 = hPulse->GetBinContent(i)- f_avg->Eval(hPulse->GetBinCenter(i));
      resG2 = hPulse->GetBinContent(i+1) - f_avg->Eval(hPulse->GetBinCenter(i+1));
      //}
      //printf("Peak %d sample %f function %f #G %f\n",  pb1,hPulse->GetBinContent(i),f_avg->Eval(hPulse->GetBinCenter(i)),hGammaReadN->GetBinContent(channel+1));

      double sigmaG = sqrt(hGammaReadSigma[channel]->GetBinContent(i+(15-pb1)));
      double sigmaG1 = sqrt(hGammaReadSigma[channel]->GetBinContent(i+(16-pb1)));
      //hGammaSigma[channel]->Fill(i+(14-pb1),resG1*resG1/hGammaReadN->GetBinContent(channel+1));//2nd run

      //if(i<82 && chiG<3700){
      if(hGammaReadN->GetBinContent(channel+1)>0 && TMath::Abs(resG1/sigmaG)<=5.0)//3rd run or last run
	{
	  hBinCorrelationG[channel]->Fill(resG1/sigmaG,resG2/sigmaG1,i+14-pb1); //1st make hGammaSigma - last run
	  hResidueBinG->Fill(i+14-pb1,resG/sigmaG,channel);//1st make hGammaSigma - last run
	}
      //}
      //printf("Bin %d\t Residue %f\n", i, hPulse->GetBinContent(i)-f_avg->Eval(i));

    }

    //printf("Energy %f\t chi2 %f\n", neu->GetTotalIntegral(), chiG );

#if 0  //For Event display to see pulses set to 1...

    hSingleGammaPulseShape[channel]->Draw();
    f_avg->Draw("same");
    c1->Update();
    a=getchar();
    printf("Next Pulse\n");
    #endif

  }


  //printf("Pulse area %f, Brent Area %f, chiN %f, chiG %f Overflow %d\n",neu->GetTotalIntegral(), brentArea, chiN, chiG, neu->GetOverflow());

  //double brentArea;
  if(chiN<chiG) brentArea=fAreaN;
  if(chiN>chiG) brentArea=fAreaG;
  //if(neu->GetTotalIntegral()>11000 && neu->GetTotalIntegral()<12500) printf("Brent area %f\t total area %f\t Slow area %f\n",brentArea, neu->GetTotalIntegral(),neu->GetSlowIntegral30());
  if(fTime<pSize && fTime<iHi){

    for(int i=0;i<2;i++){
      if((chiN-chiN_Cuts[channel][i])/(brentArea-area_nCuts[channel][i])<slope_nCuts[channel][i]
	 && brentArea>area_nCuts[channel][i] && brentArea<area_nCuts[channel][i+1] && neu->GetSlowIntegral30()/brentArea > slowTotalBr[i])
	isNeutronCS=true;

      if((chiN-chiN_Cuts[channel][i])/(brentArea-area_nCuts[channel][i])>slope_nCuts[channel][i]
	 && brentArea>area_nCuts[channel][i] && brentArea<area_nCuts[channel][i+1])
	isGammaCS=true;
    }


    // for(int i=0;i<1;i++){
    if ((chiG-chiG_Cuts[channel])/(brentArea-area_gCuts[channel][0])>slope_gCuts[channel]
	&& brentArea>area_gCuts[channel][0] && brentArea<area_gCuts[channel][1])
      notGamma=true;
    if(brentArea>area_gCuts[channel][1] && isNeutronCS) notGamma=true;

    if (brentArea>area_gCuts[channel][0] && brentArea<area_gCuts[channel][1]
	&& (chiG-chiG_Cuts[channel])/(brentArea-area_gCuts[channel][0])<slope_gCuts[channel])
      notNeutron=true;
    if(brentArea>area_gCuts[channel][1] && !isNeutronCS) notNeutron=true;
    //}

    for(int i=0;i<2;i++){
      if((neu->GetSlowIntegral30()-slowBr_Cuts[channel][i])/(brentArea-area_BrCuts[channel][i])>slope_BrCuts[channel][i]
	 && brentArea>area_BrCuts[channel][i] && brentArea<area_BrCuts[channel][i+1])
	isNeutronBr=true;

      if((neu->GetSlowIntegral30()-slowBr_Cuts[channel][i])/(brentArea-area_BrCuts[channel][i])<slope_BrCuts[channel][i]
	 && brentArea>area_BrCuts[channel][i] && brentArea<area_BrCuts[channel][i+1])
	isGammaBr=true;
    }
    if(notGamma && isNeutronCS)
      isNeutronCS_1 = true;

    if(notNeutron && isGammaCS)
      isGammaCS_1 = true;
  }
  neu->SetBrentArea(brentArea);
  neu->fIsNeutronChiSquare=isNeutronCS;
  neu->fIsNeutronTemplateSlowTotal=isNeutronBr;
  hPSD_chiSqDiff->Fill(chiN-chiG,brentArea,channel); //Only for my thesis - remove later
  hPSD_chiSq->Fill(neu->GetSlowIntegral30()/brentArea,brentArea,channel); //Only for my thesis - remove later
  hPSD->Fill(neu->GetSlowIntegral30()/neu->GetTotalIntegral(),neu->GetTotalIntegral(),channel); //Only for my thesis - remove later

#if 0  //For Event display to see pulses set to 1...
  hSingleGammaPulseShape[channel]->Draw();
  //hSingleNeutronPulseShape[channel]->Fit("f_avg","RWW");
  f_avg->Draw("same");
  c1->Update();
  a=getchar(); //wait for a keyboard prompt
  printf("Pulse area %f, Brent Area %f, chiN %f, chiG %f\n",neu->GetTotalIntegral(), brentArea, chiN, chiG);
  printf("Next Pulse\n");
  #endif

  //printf("Pulse area %f, brentAr %f, slow area %f, chiN %f, chiG %f, channel %d\n", neu->GetTotalIntegral(),brentArea,
  //neu->GetSlowIntegral30(), chiN, chiG, channel);
  //hPSD_chiN->Fill(chiN,brentArea,channel);

#if 0 //For PSD/chi2 plots using templates
      //hPSD_chiSqDiff->Fill(chiN-chiG,brentArea,channel);
      //hPSD_chiSqDiff_dof->Fill(chiN_dof-chiG_dof,brentArea,channel);
  hPSD_chiN->Fill(chiN,brentArea,channel);
  hPSD_chiG->Fill(chiG,brentArea,channel);
  hPSD_chiN_chiG->Fill(chiN,chiG,channel);
  //hPSD_chiN_dof->Fill(chiN_dof,brentArea,channel);
  //hPSD_chiG_dof->Fill(chiG_dof,brentArea,channel);
  hPSD_ratio->Fill(brentArea,neu->GetSlowIntegral30()/brentArea,channel);
  //printf("Slow total %d \t chiN %d \n",isNeutronST,isNeutronCS);
  if(isNeutronST && !isNeutronCS){
    //printf("Slow total %d \t chiN %d \t brentarea %f \t ratio %f \t",isNeutronST,isNeutronCS,brentArea,neu->GetSlowIntegral30()/brentArea);
    hPSD_notChiN_ST->Fill(neu->GetSlowIntegral30()/brentArea,brentArea,channel);
  }
  if(!isNeutronST && isNeutronCS)
    hPSD_ChiN_notST->Fill(neu->GetSlowIntegral30()/brentArea,brentArea,channel);

  //if(isNeutronST) hPSD_Neutron->Fill(chiN,brentArea,channel);

  //if(isGammaST) hPSD_Gamma->Fill(chiG,brentArea,channel);
  if(isNeutronBr) {
    hPSD_Brent_neutron->Fill(neu->GetSlowIntegral30()/brentArea,brentArea,channel);
    hPSD_chiN_neutron->Fill(chiN,brentArea,channel);
  }


  if(isNeutronCS) {
    hPSD_ratio_neutron->Fill(brentArea,neu->GetSlowIntegral30()/brentArea,channel);
    hPSD_chiN_neutronCS->Fill(chiN,brentArea,channel);
  }

  if(isNeutronCS_1) {
    hPSD_ratio_neutron_g->Fill(brentArea,neu->GetSlowIntegral30()/brentArea,channel);
    hPSD_chiN_neutron_g->Fill(chiN,brentArea,channel);
  }

  //hPSD_chiG_neutron_g->Fill(chiG,brentArea,channel);
  //hPSD_chiSq->Fill(chiN,chiG,channel);
  //hPSD_ratio->Fill(brentArea,neu->GetSlowIntegral30()/brentArea,channel);
  //hPSD->Fill(brentArea,neu->GetSlowIntegral30(),channel);
  //hArea_Br->Fill(neu->GetTotalIntegral(),brentArea,neu->GetChannel());

  #endif

  //printf("Brent area %f isNeutronChi2 %d and isNeutronbr %d\n", neu->GetBrentArea(),neu->IsNeutronChiSquare(),neu->IsNeutronTemplateSlowTotal());

#if 0 //Redefinded total/slow etc. for various PSD methods @ MTA level


  double integralCutoff, totalIntegralHalf, slowIntegralHalf20, slowIntegralHalf25, slowIntegralHalf30, slowIntegralHalf35, slowIntegralHalf40,
    totalIntegral, slowIntegral20, slowIntegral25, slowIntegral30, slowIntegral35, slowIntegral40, totalIntegralTrue, slowIntegralTrue20, slowIntegralTrue25, slowIntegralTrue30, slowIntegralTrue35, slowIntegralTrue40;

  integralCutoff = 20*binWidth;


  //Hre is the original code

      totalIntegral =
	integralFromPeakTime(numResamples, resamples, binWidth/resampleFactor, pb2, ped, -3*binWidth, integralCutoff)/resampleFactor;
          slowIntegral20 =
	    integralFromPeakTime(numResamples, resamples, binWidth/resampleFactor, pb2, ped, 3*binWidth, integralCutoff)/resampleFactor;
	      slowIntegral25 =
		integralFromPeakTime(numResamples, resamples, binWidth/resampleFactor, pb2, ped, 4*binWidth, integralCutoff)/resampleFactor;
	          slowIntegral30 =
		    integralFromPeakTime(numResamples, resamples, binWidth/resampleFactor, pb2, ped, 5*binWidth, integralCutoff)/resampleFactor;
		      slowIntegral35 =
			integralFromPeakTime(numResamples, resamples, binWidth/resampleFactor, pb2, ped, 6*binWidth, integralCutoff)/resampleFactor;
		          slowIntegral40 =
			    integralFromPeakTime(numResamples, resamples, binWidth/resampleFactor, pb2, ped, 7*binWidth, integralCutoff)/resampleFactor;


			      totalIntegralHalf =
				integralFromPeakTime(numResamples, resamples, binWidth/resampleFactor, hb, ped, -(3-1.33)*binWidth, integralCutoff+1.33)/resampleFactor;
			          slowIntegralHalf20 =
				    integralFromPeakTime(numResamples, resamples, binWidth/resampleFactor, hb, ped, (3+1.33)*binWidth, integralCutoff+1.33)/resampleFactor;
				      slowIntegralHalf25 =
					integralFromPeakTime(numResamples, resamples, binWidth/resampleFactor, hb, ped, (4+1.33)*binWidth, integralCutoff+1.33)/resampleFactor;
				          slowIntegralHalf30 =
					    integralFromPeakTime(numResamples, resamples, binWidth/resampleFactor, hb, ped, (5+1.33)*binWidth, integralCutoff+1.33)/resampleFactor;
					      slowIntegralHalf35 =
						integralFromPeakTime(numResamples, resamples, binWidth/resampleFactor, hb, ped, (6+1.33)*binWidth, integralCutoff+1.33)/resampleFactor;
					          slowIntegralHalf40 =
						    integralFromPeakTime(numResamples, resamples, binWidth/resampleFactor, hb, ped, (7+1.33)*binWidth, integralCutoff+1.33)/resampleFactor;

						  //To be used for true times only
						      totalIntegralTrue =
							integralFromTrueTime(numResamples1, resamples1, binWidth/resampleFactor1, (double)pb1+tTime, ped, -3*binWidth, integralCutoff)/resampleFactor1;
						          slowIntegralTrue20 =
							    integralFromTrueTime(numResamples1, resamples1, binWidth/resampleFactor1,(double) pb1+tTime, ped, 3*binWidth, integralCutoff)/resampleFactor1;
							      slowIntegralTrue25 =
								integralFromTrueTime(numResamples1, resamples1, binWidth/resampleFactor1, (double)pb1+tTime, ped, 4*binWidth, integralCutoff)/resampleFactor1;
							          slowIntegralTrue30 =
								    integralFromTrueTime(numResamples1, resamples1, binWidth/resampleFactor1, (double)pb1+tTime, ped, 5*binWidth, integralCutoff)/resampleFactor1;
								      slowIntegralTrue35 =
									integralFromTrueTime(numResamples1, resamples1, binWidth/resampleFactor1, (double)pb1+tTime, ped, 6*binWidth, integralCutoff)/resampleFactor1;
								          slowIntegralTrue40 =
									    integralFromTrueTime(numResamples1, resamples1, binWidth/resampleFactor1, (double)pb1+tTime, ped, 7*binWidth, integralCutoff)/resampleFactor1;

									  #endif
}



void MusunNeutronAnalysis::SumNeutron(TMusunNeutronPulse *neu){

  const vector <UShort_t> &samples = neu->GetSamples();

  int channel = neu->GetChannel();
  //printf("Nuetron # %d Gmma # %d multi %d channel %d energy %f samples_size %d gamma %d\n", pulse_N[channel], pulse_G[channel],neutron_multi, channel,neu->GetTotalIntegral(),samples.size(), isGammaST );

  if(neu->GetTotalIntegral()>2000 && neu->GetTotalIntegral()<7000 && isNeutronST && samples.size()>=60){
    //if(isNeutronST && samples.size()>=55 && neu->GetOverflow()==0){
    hNeutronN->Fill((double)channel);
    //printf("Peak %d sumN %d sumG %d\n",  pb1,sum_samples_NeutronAll[channel], sum_samples_GammaAll[channel]);
  }

  if(neu->GetTotalIntegral()>2000 && neu->GetTotalIntegral()<7000 && isGammaST && samples.size()>=60){
    //if(isGammaST && samples.size()>=55 && neu->GetOverflow()==0){
    //printf("Gamma: Peak %d sumN %d sumG %d\n",  pb1,sum_samples_NeutronAll[channel], sum_samples_GammaAll[channel]);
    hGammaN->Fill((double)channel);

  }
}


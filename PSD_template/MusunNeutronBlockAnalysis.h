#ifndef ROOT_MusunNeutronBlockAnalysis_h
#define ROOT_MusunNeutronBlockAnalysis_h 1
#include <vector>
#include <fstream>
#include <TH2.h>
#include "FillHistBase.h"
using std::vector;
//#include "Parameters.h"

class MusunNeutronBlockAnalysis : public FillHistBase
{
 private:
  std::ofstream fout;
 public:
  int left,right;
  int blockNo; int peak;
  MusunNeutronBlockAnalysis();
  void blockAnalysis(int run, int Block, int Channel);  
  virtual ~MusunNeutronBlockAnalysis();
  };

#endif

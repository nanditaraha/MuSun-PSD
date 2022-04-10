#include "MusunNeutronBlockAnalysis.h"
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>

using namespace std;

MusunNeutronBlockAnalysis::~MusunNeutronBlockAnalysis(){

}


MusunNeutronBlockAnalysis::MusunNeutronBlockAnalysis(){
  blockNo=-1;
  left = 0; right = 0;
  peak=-1;
  //printf("MusunNeutronBlockAnalysis .... called****\n");

}

void MusunNeutronBlockAnalysis::blockAnalysis(int run, int Block, int Channel){
  if(blockNo<Block){
    char name[82];
    sprintf(name,"offset_%d.txt",run);
    fout.open(name,ios::app);
    //printf("opend file %d\n",run);
    blockNo = Block;
    if (right>left) peak=1;
    if (left>right) peak = 0;
    if(left==right) peak=-1;
    fout<<blockNo<<"\n"<<peak<<"\n"<<Channel<<"\n";
    
    left=0;
    right=0;
  }
  fout.close();
  if(blockNo>Block)
    printf("WARNING!!! Blocks are not sequential \n");
}


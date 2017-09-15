
#include "TRuntimeObjects.h"

#include <iostream>
#include <map>
#include <cstdio>

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TRandom.h>
#include <TObject.h>
#include <TLine.h>

#include "TGretina.h"
#include "TS800.h"
#include "TBank29.h"
#include "TS800.h"
#include "TFastScint.h"
#include "GCutG.h"

#include "TChannel.h"
#include "GValue.h"


// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.

std::vector<GCutG*> incoming_cuts = {0};
std::vector<GCutG*> outgoing_cuts = {0};
std::vector<GCutG*> time_energy_cuts = {0};

int gates_loaded=0;


std::map<int,int> detMap = {{26,1}, {30,2}, {34,3}, {38,4}, {25,5}, {29,6}, {33,7}, {37,8}, {27,9}, {31,10}, {35,11}, {39,12},
			    {24,13}, {28,14}, {32,15}, {36,16}, {63,17}, {71,18}, {79,19}, {59,20}, {67,21}, {58,22}, {66,23},			             {60,24}, {68,25}, {76,26}, {62,27}, {70,28}, {78,29}, {56,30}, {64,31}, {57,32}, {65,33}, {61,34},
			    {69,35}, {77,36}};

//std::map<int,int> detMap = {{26,0}, {30,1}, {34,2}, {38,3}, {25,4}, {29,5}, {33,6}, {37,7}, {27,8}, {31,9}, {35,10}, {39,11},
//			    {24,12}, {28,13}, {32,14}, {36,15}, {59, 16}, {63, 17}, {67, 18}, {71, 19}, {79, 20}, {58, 21},
//			    {60, 22}, {66, 23}, {68, 24}, {76, 25}, {56, 26}, {62, 27}, {64, 28}, {70, 29}, {78, 30}, {57, 31},
//			    {61, 32}, {65, 33}, {69, 34}, {77, 35}, {57, 36}, {61, 37}, {65, 38},  {69, 39}};

std::map<int,int> quadMap = {{5,1}, {6,2}, {7,3}, {8,4}, {14,5}, {16,6}, {18,7}, {13,8}, {15,9}};

std::map<int,int> crysThetaMap = {{26,1}, {30,1}, {34,1}, {38,1}, {25,2}, {29,2}, {33,2}, {37,2}, {27,3}, {31,3}, {35,3}, {39,3},
				  {24,4}, {28,4}, {32,4}, {36,4}, {63,5}, {71,5}, {79,5}, {59,6}, {67,6}, {58,7}, {66,7}, {60,8},
				  {68,8}, {76,8}, {62,9}, {70,9}, {78,9},
				  {56,10}, {64,10}, {57,11}, {65,11}, {61,12}, {69,12}, {77,12}};
/*
detMap[26]=1;
detMap[30]=2;
detMap[34]=3;
detMap[38]=4;

detMap[25]=5;
detMap[29]=6;
detMap[33]=7;
detMap[37]=8;

detMap[27]=9;
detMap[31]=10;
detMap[35]=11;
detMap[39]=12;

detMap[24]=13;
detMap[28]=14;
detMap[32]=15;
detMap[36]=16;

detMap[63]=17;
detMap[71]=18;
detMap[79]=19;

detMap[59]=20;
detMap[67]=21;

detMap[58]=22;
detMap[66]=23;

detMap[60]=24;
detMap[68]=25;
detMap[76]=26;

detMap[62]=27;
detMap[70]=28;
detMap[78]=29;

detMap[56]=30;
detMap[64]=31;

detMap[57]=32;
detMap[65]=33;

detMap[61]=34;
detMap[69]=35;
detMap[77]=36;
*/

bool HandleTiming(TRuntimeObjects &obj, TCutG *incoming, TCutG* outgoing) {
  TS800 *s800  = obj.GetDetector<TS800>();
  //TBank29 *bank29  = obj.GetDetector<TBank29>();
  
  if(!s800)
    {return false;}

  int E1UpSize = s800->GetMTof().E1UpSize();
  int XfpSize = s800->GetMTof().XfpSize();
  int ObjSize = s800->GetMTof().ObjSize();
  std::string dirname = "Timing";


  ///////////////////////////////////////////////////////
  //incoming gates on Xfp - E1 vs. OBJ - E1 Uncorrected//
  ///////////////////////////////////////////////////////
  if(incoming) {
    {dirname=Form("Timing_%s",incoming->GetName());}
    bool passed = false;  
    if(incoming->IsInside(s800->GetMTof().GetCorrelatedObjE1(),s800->GetMTof().GetCorrelatedXfpE1())) {
          passed=true;
        }
    if(!passed)
      return false;
  }

  obj.FillHistogram(dirname,"E1UpSize",10,0,10,E1UpSize);

  obj.FillHistogram(dirname,"XfpSize",10,0,10,XfpSize);

  obj.FillHistogram(dirname,"ObjSize",10,0,10,ObjSize);

  obj.FillHistogram(dirname,"Register",10,0,10,s800->GetReg());
  

  obj.FillHistogram(dirname,"E1",4000,0,64000,s800->GetMTof().GetCorrelatedE1Up());

  obj.FillHistogram(dirname,"Xfp",4000,0,64000,s800->GetMTof().GetCorrelatedXfp());

  obj.FillHistogram(dirname,"Obj",4000,0,64000,s800->GetMTof().GetCorrelatedObj());

  obj.FillHistogram(dirname,"Xfp-E1",4000,-10000,10000,s800->GetMTof().GetCorrelatedXfpE1());

  obj.FillHistogram(dirname,"Obj-E1",4000,-10000,10000,s800->GetMTof().GetCorrelatedObjE1());

  obj.FillHistogram(dirname,"Obj-Xfp",2000,-6000,-4000
		                     ,s800->GetMTof().GetCorrelatedObj() - s800->GetMTof().GetCorrelatedXfp());

  obj.FillHistogram(dirname,"Xfp-Obj",2000,4000,6000
		                     ,s800->GetMTof().GetCorrelatedXfp() - s800->GetMTof().GetCorrelatedObj());  

  if(s800->GetReg() == 2) {
   obj.FillHistogram(dirname,"E1_Coin_raw",4000,0,64000,s800->GetMTof().fE1Up);

   obj.FillHistogram(dirname,"Xfp_Coin_raw",4000,0,64000,s800->GetMTof().fXfp);

   obj.FillHistogram(dirname,"Obj_Coin_raw",4000,0,64000,s800->GetMTof().fObj);
  }

  if(s800->GetReg() == 1) {
   obj.FillHistogram(dirname,"E1_Sing_raw",4000,0,64000,s800->GetMTof().fE1Up);

   obj.FillHistogram(dirname,"Xfp_Sing_raw",4000,0,64000,s800->GetMTof().fXfp);

   obj.FillHistogram(dirname,"Obj_Sing_raw",4000,0,64000,s800->GetMTof().fObj);
  }
     
  obj.FillHistogram(dirname,"E1_raw",4000,0,64000,s800->GetMTof().fE1Up);

  obj.FillHistogram(dirname,"Xfp_raw",4000,0,64000,s800->GetMTof().fXfp);
 
  obj.FillHistogram(dirname,"Obj_raw",4000,0,64000,s800->GetMTof().fObj);
  
  for(int i=0; i<E1UpSize; i++) {
   for(int j=0; j<XfpSize;j++) {  
     obj.FillHistogram(dirname,"Xfp-E1_raw",4000,2000,6000,s800->GetMTof().fXfp.at(j) - s800->GetMTof().fE1Up.at(i));
   }
  }
  
  for(int i=0; i<E1UpSize; i++) {
   for(int j=0; j<ObjSize;j++) { 
     obj.FillHistogram(dirname,"Obj-E1_raw",4000,-2000,0,s800->GetMTof().fObj.at(j) - s800->GetMTof().fE1Up.at(i));
   }
  }

  if(!outgoing)
    {return false;}
  if(!outgoing->IsInside(s800->GetMTofObjE1(),s800->GetIonChamber().Charge()))
    {return false;}

  obj.FillHistogram(dirname,Form("Obj-Xfp_%s",outgoing->GetName()),2000,-6000,-4000
		    ,s800->GetMTof().GetCorrelatedObj() - s800->GetMTof().GetCorrelatedXfp());

  obj.FillHistogram(dirname,Form("Xfp-Obj_%s",outgoing->GetName()),2000,4000,6000
		    ,s800->GetMTof().GetCorrelatedXfp() - s800->GetMTof().GetCorrelatedObj()); 
            
  return true;
}

bool HandleS800(TRuntimeObjects &obj) {
  TS800 *s800  = obj.GetDetector<TS800>();
  
  if(!s800)
    {return false;}
  
  return true;
}

bool s800Singles(TRuntimeObjects &obj,TCutG *incoming, TCutG* outgoing) {
  TS800 *s800  = obj.GetDetector<TS800>();  

  if(!s800)
    {return false;}
  
  std::string dirname;

  ///////////////////////////////////////////////////////
  //incoming gates on Xfp - E1 vs. OBJ - E1 Uncorrected//
  ///////////////////////////////////////////////////////
  
  if(incoming) {
    bool passed = false;  
    if(incoming->IsInside(s800->GetMTof().GetCorrelatedObjE1(),s800->GetMTof().GetCorrelatedXfpE1())) {
          passed=true;
        }
    if(!passed)
      return false;
  }

  ///////////////////////////////////////////////////////
  //incoming gates on Xfp - E1 vs. OBJ - E1 (Corrected)//
  ///////////////////////////////////////////////////////
  /*
  if(incoming) {
    bool passed = false;  
    if(incoming->IsInside(s800->GetMTofObjE1(),s800->GetMTofXfpE1())) {
          dirname=Form("S800_%s",incoming->GetName());
          passed=true;
        }
    if(!passed)
      return false;
  }*/

  //s800 Singles Only
  dirname = "s800Singles";
  if(incoming)
    {dirname=Form("s800Singles_%s",incoming->GetName());}

  if(s800->GetReg() == 1) {

   obj.FillHistogram(dirname,"Ion_v_Obj-E1",1400,-1500,-800,s800->GetMTofObjE1()
                                           ,1600,20000,35000,s800->GetIonChamber().Charge());

   obj.FillHistogram(dirname,"Ion_v_Obj-E1_Uncorrected",700,-1500,-800,s800->GetMTof().GetCorrelatedObjE1()
                                                       ,1600,20000,35000,s800->GetIonChamber().Charge());

   obj.FillHistogram(dirname,"IonChamber",1600,20000,35000,s800->GetIonChamber().Charge());

   obj.FillHistogram(dirname, "CRDC1X", 1200, -300, 300, s800->GetCrdc(0).GetDispersiveX());

   obj.FillHistogram(dirname,"CRDC1Y_v_CRDC1X",1200,-300,300,s800->GetCrdc(0).GetDispersiveX()
                                              ,800,-200,200,s800->GetCrdc(0).GetNonDispersiveY());

   obj.FillHistogram(dirname,"CRDC1X_v_Ion",1200, -300, 300, s800->GetCrdc(0).GetDispersiveX()
		                           ,1500,20000,35000,s800->GetIonChamber().Charge());

   obj.FillHistogram(dirname,"XFP_v_Xfp-E1",1400,3000,4400,s800->GetMTofXfpE1()
                                           ,1200,-300,300,s800->GetXFP());

   obj.FillHistogram(dirname,"XFP_v_Xfp-E1_Uncorrected",700,3000,4400,s800->GetMTof().GetCorrelatedXfpE1()
                                                       ,600,-300,300,s800->GetXFP());

   obj.FillHistogram(dirname,"XFP_v_Obj-E1",1200,-1400,-800,s800->GetMTofObjE1()
                                           ,1200,-300,300,s800->GetXFP());

   obj.FillHistogram(dirname,"XFP_v_Obj-E1_Uncorrected",600,-1400,-800,s800->GetMTof().GetCorrelatedObjE1()
                                                       ,600,-300,300,s800->GetXFP());

   obj.FillHistogram(dirname,"Xfp-E1_v_Obj-E1",1400,-1500,-700,s800->GetMTofObjE1()
        	                              ,2800,3000,4400,s800->GetMTofXfpE1());

   obj.FillHistogram(dirname,"Xfp-E1_v_Obj-E1_Uncorrected",800,-1500,-700,s800->GetMTof().GetCorrelatedObjE1()
		                                          ,1400,3000,4400,s800->GetMTof().GetCorrelatedXfpE1());

   if(!outgoing)
    {return false;}
  
  if(!outgoing->IsInside(s800->GetMTofObjE1(),s800->GetIonChamber().Charge()))
    {return false;}

   obj.FillHistogram(dirname,Form("Ion_v_Obj-E1_%s",outgoing->GetName()),1400,-1500,-800,s800->GetMTofObjE1()
                                                                        ,1600,20000,35000,s800->GetIonChamber().Charge());

   obj.FillHistogram(dirname,Form("Ion_v_Obj-E1_Uncorrected_%s",outgoing->GetName())
		                                               ,700,-1500,-800,s800->GetMTof().GetCorrelatedObjE1()
		                                               ,1600,20000,35000,s800->GetIonChamber().Charge());

   obj.FillHistogram(dirname,Form("IonChamber_%s",outgoing->GetName()),1600,20000,35000,s800->GetIonChamber().Charge());

   obj.FillHistogram(dirname, Form("CRDC1X_%s",outgoing->GetName()), 1200, -300, 300, s800->GetCrdc(0).GetDispersiveX());

   obj.FillHistogram(dirname,Form("CRDC1Y_v_CRDC1X_%s",outgoing->GetName()),1200,-300,300,s800->GetCrdc(0).GetDispersiveX()
                                                                           ,800,-200,200,s800->GetCrdc(0).GetNonDispersiveY());

   obj.FillHistogram(dirname,Form("CRDC1X_v_Ion_%s",outgoing->GetName())
		                                   ,1200,-300,300,1500,20000,35000,s800->GetIonChamber().Charge()
		                                   ,s800->GetCrdc(0).GetDispersiveX());

   obj.FillHistogram(dirname,Form("XFP_v_Xfp-E1_%s",outgoing->GetName()),1400,3000,4400,s800->GetMTofXfpE1()
                                                                        ,1200,-300,300,s800->GetXFP());

   obj.FillHistogram(dirname,Form("XFP_v_Xfp-E1_Uncorrected_%s",outgoing->GetName())
		                                               ,700,3000,4400,s800->GetMTof().GetCorrelatedXfpE1()
                                                               ,600,-300,300,s800->GetXFP());

   obj.FillHistogram(dirname,Form("XFP_v_Obj-E1_%s",outgoing->GetName()),1200,-1400,-800,s800->GetMTofObjE1()
                                                                        ,1200,-300,300,s800->GetXFP());

   obj.FillHistogram(dirname,Form("XFP_v_Obj-E1_Uncorrected_%s",outgoing->GetName())
		                                               ,600,-1400,-800,s800->GetMTof().GetCorrelatedObjE1()
                                                               ,600,-300,300,s800->GetXFP());

   obj.FillHistogram(dirname,Form("Xfp-E1_v_Obj-E1_%s",outgoing->GetName()),1400,-1500,-700,s800->GetMTofObjE1()
        	                                                           ,2800,3000,4400,s800->GetMTofXfpE1());

   obj.FillHistogram(dirname,Form("Xfp-E1_v_Obj-E1_Uncorrected_%s",outgoing->GetName())
		                                                  ,800,-1500,-700,s800->GetMTof().GetCorrelatedObjE1()
		                                                  ,1400,3000,4400,s800->GetMTof().GetCorrelatedXfpE1());
 
  }

  return true;
  
}

bool HandleS800_Gated(TRuntimeObjects &obj,TCutG *incoming, TCutG* outgoing) {
  TS800 *s800  = obj.GetDetector<TS800>();  

  if(!s800)
    {return false;}
  
  std::string dirname = "S800";

  ///////////////////////////////////////////////////////
  //incoming gates on Xfp - E1 vs. OBJ - E1 Uncorrected//
  ///////////////////////////////////////////////////////
  
  if(incoming) {
    bool passed = false;  
    if(incoming->IsInside(s800->GetMTof().GetCorrelatedObjE1(),s800->GetMTof().GetCorrelatedXfpE1())) {
          dirname=Form("S800_%s",incoming->GetName());
          passed=true;
        }
    if(!passed)
      return false;
  }

  ///////////////////////////////////////////////////////
  //incoming gates on Xfp - E1 vs. OBJ - E1 (Corrected)//
  ///////////////////////////////////////////////////////
  /*
  if(incoming) {
    bool passed = false;  
    if(incoming->IsInside(s800->GetMTofObjE1(),s800->GetMTofXfpE1())) {
          dirname=Form("S800_%s",incoming->GetName());
          passed=true;
        }
    if(!passed)
      return false;
  }*/
                                 
  //Incoming PID Plots
  dirname = "IncPID";
  if(incoming)
    {dirname=Form("IncPID_%s",incoming->GetName());}

  obj.FillHistogram(dirname,"Xfp-E1_v_Obj-E1",1400,-1500,-700,s800->GetMTofObjE1()
        	                             ,2800,3000,4400,s800->GetMTofXfpE1());

  obj.FillHistogram(dirname,"Xfp_v_Obj_Uncorrected",1500,11000,14000,s800->GetMTof().GetCorrelatedObj()
		                                   ,1250,16000,18500,s800->GetMTof().GetCorrelatedXfp());

  obj.FillHistogram(dirname,"Xfp-E1_v_Obj-E1_Uncorrected",800,-1500,-700,s800->GetMTof().GetCorrelatedObjE1()
		                                         ,1400,3000,4400,s800->GetMTof().GetCorrelatedXfpE1());
  /* 
  for(int j=0; j < XfpSize; j++) {
   for(int i=0; i < ObjSize; i++) {
     
     obj.FillHistogram(dirname,"Xfp_v_Obj_raw",4000,11000,15000,s800->GetMTof().fObj.at(i)
                                              ,3000,16000,19000,s800->GetMTof().fXfp.at(j));
      
   }
  }
  
  for(int k=0; k < E1UpSize; k++) {
   for(int j=0; j < XfpSize; j++) {
    for(int i=0; i < ObjSize; i++) {
      
      obj.FillHistogram(dirname,"Xfp-E1_v_Obj-E1_raw",800,-1500,-700,s800->GetMTof().fObj.at(i)-s800->GetMTof().fE1Up.at(k)
                                                     ,1600,3000,4600,s800->GetMTof().fXfp.at(j)-s800->GetMTof().fE1Up.at(k));
     
    } 
   }
  }
  */

  //Outgoing PID Plots
  dirname = "OutPID";
  if(incoming)
    {dirname=Form("OutPID_%s",incoming->GetName());}
  
  obj.FillHistogram(dirname,"Ion_v_Obj-E1",1400,-1500,-800,s800->GetMTofObjE1()
                                          ,1600,20000,35000,s800->GetIonChamber().Charge());

  obj.FillHistogram(dirname,"Ion_v_Obj-E1_Uncorrected",700,-1500,-800,s800->GetMTof().GetCorrelatedObjE1()
                                                      ,1600,20000,35000,s800->GetIonChamber().Charge());

  obj.FillHistogram(dirname,"IonChamber",1600,20000,35000,s800->GetIonChamber().Charge());

  
  //CRDC Positions
  dirname = "CRDCPos";
  if(incoming)
    {dirname=Form("CRDCPos_%s",incoming->GetName());}
  
  obj.FillHistogram(dirname, "CRDC1X", 1200, -300, 300, s800->GetCrdc(0).GetDispersiveX());

  obj.FillHistogram(dirname, "CRDC1Y", 800, -200, 200, s800->GetCrdc(0).GetNonDispersiveY());

  obj.FillHistogram(dirname, "CRDC2X", 1200, -300, 300, s800->GetCrdc(1).GetDispersiveX());

  obj.FillHistogram(dirname, "CRDC2Y", 800, -200, 200, s800->GetCrdc(1).GetNonDispersiveY());

  obj.FillHistogram(dirname, "CRDC1Y - CRDC2Y", 800, -200, 200,
		    s800->GetCrdc(0).GetNonDispersiveY() - s800->GetCrdc(1).GetNonDispersiveY()); 

  obj.FillHistogram(dirname,"CRDC1Y_v_CRDC1X",1200,-300,300,s800->GetCrdc(0).GetDispersiveX()
                                             ,800,-200,200,s800->GetCrdc(0).GetNonDispersiveY());

  //obj.FillHistogram(dirname,"CRDC1X_v_CRDC1Y_Uncal",600,0,300,s800->GetCrdc(0).GetDispersiveX()
  //                                           ,2000,0,2000,s800->GetCrdc(0).GetNonDispersiveY());
   
  obj.FillHistogram(dirname,"CRDC2Y_v_CRDC2X",1200,-300,300,s800->GetCrdc(1).GetDispersiveX()
                                             ,800,-200,200,s800->GetCrdc(1).GetNonDispersiveY());

  //obj.FillHistogram(dirname,"CRDC2X_v_CRDC2Y_Uncal",600,0,300,s800->GetCrdc(1).GetDispersiveX()
  //                                           ,2000,0,2000,s800->GetCrdc(1).GetNonDispersiveY());

  obj.FillHistogram(dirname,"CRDC1X_v_Ion",1200, -300, 300, s800->GetCrdc(0).GetDispersiveX()
		                          ,1500,20000,35000,s800->GetIonChamber().Charge());

  
  //Focal Plane, Target
  dirname = "FPandTA";
  if(incoming)
    {dirname=Form("FPandTA_%s",incoming->GetName());}
  
  obj.FillHistogram(dirname,"AFP",1000,-0.1,0.1,s800->GetAFP());

  obj.FillHistogram(dirname,"BFP",1000,-0.15,0.15,s800->GetBFP());

  obj.FillHistogram(dirname,"ATA",1000,-0.1,0.1,s800->GetAta());

  obj.FillHistogram(dirname,"YTA",1000,-50,50,s800->GetYta());

  obj.FillHistogram(dirname,"BTA",1000,-0.1,0.1,s800->GetBta());

  obj.FillHistogram(dirname,"DTA",1000,-0.5,0.5,s800->GetDta());

  obj.FillHistogram(dirname,"YFP_v_BFP",1000,-0.1,0.1,s800->GetBFP(),400,-200,200,s800->GetYFP());

  obj.FillHistogram(dirname,"AFP_v_BFP",1000,-0.1,0.1,s800->GetBFP(),1000,-0.1,0.1,s800->GetAFP());

  obj.FillHistogram(dirname,"ATA_v_BTA",1000,-0.1,0.1,s800->GetBta(),1000,-0.1,0.1,s800->GetAta());

  obj.FillHistogram(dirname,"Azita",400,-6.3,6.3,s800->Azita());

  
  //Focal Plane vs ToF
  dirname = "FPvToF";
  if(incoming)
    {dirname=Form("FPvToF_%s",incoming->GetName());}
  
  obj.FillHistogram(dirname,"XFP_v_Obj-E1",1200,-1400,-800,s800->GetMTofObjE1()
                                          ,1200,-300,300,s800->GetXFP());

  obj.FillHistogram(dirname,"XFP_v_Obj-E1_Uncorrected",600,-1400,-800,s800->GetMTof().GetCorrelatedObjE1()
                                                      ,600,-300,300,s800->GetXFP());

  obj.FillHistogram(dirname,"AFP_v_Obj-E1",800,-1300,-900,s800->GetMTofObjE1()
                                          ,1000,-0.05,0.05,s800->GetAFP());

  obj.FillHistogram(dirname,"AFP_v_Obj-E1_Uncorrected",400,-1300,-900,s800->GetMTof().GetCorrelatedObjE1()
                                                      ,1000,-0.05,0.05,s800->GetAFP());

  obj.FillHistogram(dirname,"XFP_v_Xfp-E1",1400,3000,4400,s800->GetMTofXfpE1()
                                          ,1200,-300,300,s800->GetXFP());

  obj.FillHistogram(dirname,"XFP_v_Xfp-E1_Uncorrected",700,3000,4400,s800->GetMTof().GetCorrelatedXfpE1()
                                                      ,600,-300,300,s800->GetXFP());

  obj.FillHistogram(dirname,"AFP_v_Xfp-E1",1400,3000,4400,s800->GetMTofXfpE1()
                                          ,1000,-0.05,0.05,s800->GetAFP());

  obj.FillHistogram(dirname,"AFP_v_Xfp-E1_Uncorrected",700,3000,4400,s800->GetMTof().GetCorrelatedXfpE1()
                                                      ,500,-0.05,0.05,s800->GetAFP());

  
  ///////////////////////////////////////////////////////////////////
  //OUTGOING GATES on Ion Chamber Charge vs Object - E1 (Corrected)//
  ///////////////////////////////////////////////////////////////////
 
  if(!outgoing)
    {return false;}
  if(!outgoing->IsInside(s800->GetMTofObjE1(),s800->GetIonChamber().Charge()))
    {return false;}

  //ToF vs EvtNum
  dirname = "ToFvEvtNum";
  if(incoming)
    {dirname=Form("ToFvEvtNum_%s",incoming->GetName());}

  obj.FillHistogram(dirname,Form("Obj-E1_v_EvtNum_Uncorrected_%s",outgoing->GetName()),5000,0,5000000,s800->GetEventCounter()
		                                                 ,2000,-4000,2000,s800->GetMTof().GetCorrelatedObjE1());

  obj.FillHistogram(dirname,Form("Xfp-E1_v_EvtNum_Uncorrected_%s",outgoing->GetName()),5000,0,5000000,s800->GetEventCounter()
		                                                 ,2000,1000,7000,s800->GetMTof().GetCorrelatedXfpE1());

  obj.FillHistogram(dirname,Form("Obj-E1_v_EvtNum_%s",outgoing->GetName()),5000,0,5000000,s800->GetEventCounter()
		                                                          ,2000,-4000,2000,s800->GetMTofObjE1());

  obj.FillHistogram(dirname,Form("Xfp-E1_v_EvtNum_%s",outgoing->GetName()),5000,0,5000000,s800->GetEventCounter()
		                                                          ,2000,1000,7000,s800->GetMTofXfpE1());

  obj.FillHistogram(dirname,Form("EvtNum_%s",outgoing->GetName()),5000,0,5000000,s800->GetEventCounter());
    
  //Incoming PID Plots
  dirname = "IncPID";
  if(incoming)
    {dirname=Form("IncPID_%s",incoming->GetName());}
  
  obj.FillHistogram(dirname,Form("Xfp-E1_v_Obj-E1_%s",outgoing->GetName()),1400,-1500,-700,s800->GetMTofObjE1()
        	                                                          ,2800,3000,4400,s800->GetMTofXfpE1());

  obj.FillHistogram(dirname,Form("Xfp_v_Obj_Uncorrected_%s",outgoing->GetName())
		                                           ,1500,11000,14000,s800->GetMTof().GetCorrelatedObj()
		                                           ,1250,16000,18500,s800->GetMTof().GetCorrelatedXfp());

  obj.FillHistogram(dirname,Form("Xfp-E1_v_Obj-E1_Uncorrected_%s",outgoing->GetName())
		                                                ,400,-1500,-700,s800->GetMTof().GetCorrelatedObjE1()
		                                                ,1400,3000,4400,s800->GetMTof().GetCorrelatedXfpE1());
  /* 
  for(int j=0; j < XfpSize; j++) {
   for(int i=0; i < ObjSize; i++) {
    obj.FillHistogram(dirname,Form("Xfp_v_Obj_raw_%s",outgoing->GetName()),4000,11000,15000,s800->GetMTof().fObj.at(i)
                                                                          ,3000,16000,19000,s800->GetMTof().fXfp.at(j));
   }
  }
  
  for(int k=0; k < E1UpSize; k++) { 
   for(int j=0; j < XfpSize; j++) {
    for(int i=0; i < ObjSize; i++) {
      obj.FillHistogram(dirname,Form("Xfp-E1_v_Obj-E1_raw_%s",outgoing->GetName())
			             ,800,-1500,-700,s800->GetMTof().fObj.at(i)-s800->GetMTof().fE1Up.at(k)
			             ,600,3400,4000,s800->GetMTof().fXfp.at(j)-s800->GetMTof().fE1Up.at(k));
    } 
   }
  }
  */
  
 //Outgoing PID Plots
 dirname = "OutPID";
 if(incoming)
   {dirname=Form("OutPID_%s",incoming->GetName());}
  
 obj.FillHistogram(dirname,Form("Ion_v_Obj-E1_%s",outgoing->GetName()),1400,-1500,-800,s800->GetMTofObjE1()
                                                                      ,1600,20000,35000,s800->GetIonChamber().Charge());

 obj.FillHistogram(dirname,Form("Ion_v_Obj-E1_Uncorrected_%s",outgoing->GetName())
		                                             ,700,-1500,-800,s800->GetMTof().GetCorrelatedObjE1()
                                                             ,1600,20000,35000,s800->GetIonChamber().Charge());

 obj.FillHistogram(dirname,Form("IonChamber_%s",outgoing->GetName()),1600,20000,35000,s800->GetIonChamber().Charge());

 //CRDC Positions
 dirname = "CRDCPos";
 if(incoming)
    {dirname=Form("CRDCPos_%s",incoming->GetName());}
  
 obj.FillHistogram(dirname,Form("CRDC1X_%s",outgoing->GetName()), 1200, -300, 300, s800->GetCrdc(0).GetDispersiveX());

 obj.FillHistogram(dirname,Form("CRDC1Y_%s",outgoing->GetName()), 800, -200, 200, s800->GetCrdc(0).GetNonDispersiveY());

 obj.FillHistogram(dirname,Form("CRDC1Y_v_CRDC1X_%s",outgoing->GetName()),1200,-300,300,s800->GetCrdc(0).GetDispersiveX()
                                                                         ,800,-200,200,s800->GetCrdc(0).GetNonDispersiveY());

 obj.FillHistogram(dirname,Form("CRDC1X_v_Ion_%s",outgoing->GetName()),1500,20000,35000,s800->GetIonChamber().Charge()
		                                                      ,1200,-300,300,s800->GetCrdc(0).GetDispersiveX());

 //Focal Plane, Target
 dirname = "FPandTA";
 if(incoming)
    {dirname=Form("FPandTA_%s",incoming->GetName());}
 
 obj.FillHistogram(dirname,Form("AFP_%s",outgoing->GetName()),1000,-0.1,0.1,s800->GetAFP());

 obj.FillHistogram(dirname,Form("BFP_%s",outgoing->GetName()),1000,-0.15,0.15,s800->GetBFP());

 obj.FillHistogram(dirname,Form("ATA_%s",outgoing->GetName()),1000,-0.1,0.1,s800->GetAta());

 obj.FillHistogram(dirname,Form("YTA_%s",outgoing->GetName()),1000,-50,50,s800->GetYta());

 obj.FillHistogram(dirname,Form("BTA_%s",outgoing->GetName()),1000,-0.1,0.1,s800->GetBta());

 obj.FillHistogram(dirname,Form("DTA_%s",outgoing->GetName()),1000,-0.5,0.5,s800->GetDta());

 obj.FillHistogram(dirname,Form("YFP_v_BFP_%s",outgoing->GetName()),1000,-0.1,0.1,s800->GetBFP()
                                                                   ,400,-200,200,s800->GetYFP());

 obj.FillHistogram(dirname,Form("AFP_v_BFP_%s",outgoing->GetName()),1000,-0.1,0.1,s800->GetBFP()
		                                                   ,1000,-0.1,0.1,s800->GetAFP());

 obj.FillHistogram(dirname,Form("ATA_v_BTA_%s",outgoing->GetName()),1000,-0.1,0.1,s800->GetBta()
                                                                   ,1000,-0.1,0.1,s800->GetAta());

 obj.FillHistogram(dirname,Form("Azita_%s",outgoing->GetName()),400,-6.3,6.3,s800->Azita());


 //Focal Plane vs ToF
 dirname = "FPvToF";
 if(incoming)
    {dirname=Form("FPvToF_%s",incoming->GetName());}

 obj.FillHistogram(dirname,Form("XFP_v_Obj-E1_%s",outgoing->GetName()),1200,-1400,-800,s800->GetMTofObjE1()
                                                                      ,1200,-300,300,s800->GetXFP());

 obj.FillHistogram(dirname,Form("XFP_v_Obj-E1_Uncorrected_%s",outgoing->GetName())
		                                             ,1200,-1400,-800,s800->GetMTof().GetCorrelatedObjE1()
                                                             ,600,-300,300,s800->GetXFP());

 obj.FillHistogram(dirname,Form("AFP_v_Obj-E1_%s",outgoing->GetName()),800,-1300,-900,s800->GetMTofObjE1()
                                                                      ,1000,-0.05,0.05,s800->GetAFP());

 obj.FillHistogram(dirname,Form("AFP_v_Obj-E1_Uncorrected_%s",outgoing->GetName())
		                                             ,400,-1300,-900,s800->GetMTof().GetCorrelatedObjE1()
                                                             ,1000,-0.05,0.05,s800->GetAFP());

 obj.FillHistogram(dirname,Form("XFP_v_Xfp-E1_%s",outgoing->GetName()),1400,3000,4400,s800->GetMTofXfpE1()
                                                                      ,1200,-300,300,s800->GetXFP());

 obj.FillHistogram(dirname,Form("XFP_v_Xfp-E1_Uncorrected_%s",outgoing->GetName())
		                                             ,700,3000,4400,s800->GetMTof().GetCorrelatedXfpE1()
                                                             ,600,-300,300,s800->GetXFP());

 obj.FillHistogram(dirname,Form("AFP_v_Xfp-E1_%s",outgoing->GetName()),1400,3000,4400,s800->GetMTofXfpE1()
                                                                      ,1000,-0.05,0.05,s800->GetAFP());

 obj.FillHistogram(dirname,Form("AFP_v_Xfp-E1_Uncorrected_%s",outgoing->GetName())
		                                             ,700,3000,4400,s800->GetMTof().GetCorrelatedXfpE1()
                                                             ,500,-0.05,0.05,s800->GetAFP());
  
  return true;
    
}

bool HandleGretina(TRuntimeObjects &obj) {
  TGretina *gretina  = obj.GetDetector<TGretina>(); 
   
  if(!gretina)
    {return false;}
  std::string dirname = "Gretina_Ungated";
  for(unsigned int i=0;i<gretina->Size();i++) {
      TGretinaHit hit = gretina->GetGretinaHit(i);

      if(hit.GetCoreEnergy() > 100)
	{

         obj.FillHistogram(dirname,"GretinaCoreEnergy",2000,0,4000,hit.GetCoreEnergy());

         obj.FillHistogram(dirname,"GretinaSummarySpectrum",200,0,200,hit.GetCrystalId()
			                                   ,2000,0,4000,hit.GetCoreEnergy());
                                                 
         obj.FillHistogram(dirname,"GretinaPositionSpectrum",226,0,3.2,hit.GetTheta()
                                                            ,452,0,6.3,hit.GetPhi());
	 

	 obj.FillHistogram(dirname,"CrystalId",100,0,100,hit.GetCrystalId());

         obj.FillHistogram(dirname,"QuadNum",100,0,100,hit.GetHoleNumber());

	 obj.FillHistogram(dirname,"DetMap",100,0,100,detMap[hit.GetCrystalId()]);

	 obj.FillHistogram(dirname,"CrysThetaMap",100,0,100,crysThetaMap[hit.GetCrystalId()]);

         obj.FillHistogram(dirname,"QuadMap",100,0,100,quadMap[hit.GetHoleNumber()]);
	 
	 
	 obj.FillHistogram(dirname,"Theta_v_CrysId",100,0,100,hit.GetCrystalId()
			                           ,226,0,3.2,hit.GetTheta());

	 obj.FillHistogram(dirname,"Phi_v_CrysId",100,0,100,hit.GetCrystalId()
			                         ,452,0,6.3,hit.GetPhi());

	 obj.FillHistogram(dirname,"Theta_v_QuadNum",100,0,100,hit.GetHoleNumber()
			                            ,226,0,3.2,hit.GetTheta());

	 obj.FillHistogram(dirname,"Phi_v_QuadNum",100,0,100,hit.GetHoleNumber()
			                          ,452,0,6.3,hit.GetPhi());

	 obj.FillHistogram(dirname,"HitTheta_v_DetMap",100,0,100,detMap[hit.GetCrystalId()]
			                              ,226,0,3.2,hit.GetTheta());

	 obj.FillHistogram(dirname,"HitPhi_v_DetMap",100,0,100,detMap[hit.GetCrystalId()]
			                            ,452,0,6.3,hit.GetPhi());

	 obj.FillHistogram(dirname,"HitTheta_v_CrysThetaMap",100,0,100,crysThetaMap[hit.GetCrystalId()]
			                                    ,226,0,3.2,hit.GetTheta());

	 obj.FillHistogram(dirname,"HitPhi_v_CrysThetaMap",100,0,100,crysThetaMap[hit.GetCrystalId()]
			                                  ,452,0,6.3,hit.GetPhi());

	 obj.FillHistogram(dirname,"HitTheta_v_QuadMap",100,0,100,quadMap[hit.GetHoleNumber()]
			                               ,226,0,3.2,hit.GetTheta());

	 obj.FillHistogram(dirname,"HitPhi_v_QuadMap",100,0,100,quadMap[hit.GetHoleNumber()]
			                             ,452,0,6.3,hit.GetPhi());

	 obj.FillHistogram(dirname,"CrystalTheta_v_DetMap",100,0,100,detMap[hit.GetCrystalId()]
			                                  ,226,0,3.2,(gretina->GetCrystalPosition(hit.GetCrystalId())).Theta());

	 obj.FillHistogram(dirname,"CrystalPhi_v_DetMap",100,0,100,detMap[hit.GetCrystalId()]
			   ,452,0,6.3,(gretina->GetCrystalPosition(hit.GetCrystalId())).Phi());

	 //obj.FillHistogram(dirname,"CoreEnergy_v_DetMap",38,0,38,detMap[hit.GetCrystalId()]
	 //		                                ,2000,0,4000,hit.GetCoreEnergy());
	 
                                                             
         //plots with s800
         //TS800 *s800 = obj.GetDetector<TS800>();
         //if(!s800) {return false;}
      
        } 
    }


  return true;

}

bool Gretina_PIDOUTGatedSpectra(TRuntimeObjects &obj,TGretinaHit hit,TS800* s800,TBank29* bank29,
				std::string dirname, TCutG* incoming, TCutG* outgoing, TCutG* time_energy=NULL) {

  if(hit.GetCoreEnergy() > 100) {

    if(incoming)
      {dirname+=Form("_%s",incoming->GetName());}

    if(time_energy)
      {dirname+=Form("_%s",time_energy->GetName());}

      TVector3 track = s800->Track();
      double yta = s800->GetYta();

      obj.FillHistogram(dirname,Form("HitTheta_v_DetMap_%s",outgoing->GetName()),100,0,100,crysThetaMap[hit.GetCrystalId()]
			                                                        ,226,0,3.2,hit.GetTheta());

      obj.FillHistogram(dirname,Form("HitPhi_v_DetMap_%s",outgoing->GetName()),100,0,100,crysThetaMap[hit.GetCrystalId()]
			                                                      ,452,0,6.3,hit.GetPhi());

      obj.FillHistogram(dirname,Form("HitTheta_v_CrysThetaMap_%s",outgoing->GetName()),100,0,100,crysThetaMap[hit.GetCrystalId()]
			                                                              ,226,0,3.2,hit.GetTheta());

      obj.FillHistogram(dirname,Form("HitPhi_v_CrysThetaMap_%s",outgoing->GetName()),100,0,100,crysThetaMap[hit.GetCrystalId()]
			                                                            ,452,0,6.3,hit.GetPhi());
      
      //Gamma Ray Spectra
      obj.FillHistogram(dirname,Form("Gamma(Beta)_%s",outgoing->GetName()),1000,0,4000,hit.GetDoppler(GValue::Value("BETA")));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_%s",outgoing->GetName())
			                                   ,1000,0,4000,hit.GetDoppler(GValue::Value("BETA"),&track));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA)_%s",outgoing->GetName())
			,2000,0,4000,hit.GetDopplerYta(GValue::Value("BETA"),yta,&track));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA&DTA)_%s",outgoing->GetName())
			,2000,0,4000,hit.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")),yta,&track));
      

      //Gretina Timing
      obj.FillHistogram(dirname,Form("Gamma_v_Time_%s",outgoing->GetName())
			                              ,1200,-400,800,bank29->Timestamp() - hit.GetTime()
                                                      ,2000,0,8000,hit.GetDoppler(GValue::Value("BETA"),&track));

      obj.FillHistogram(dirname,Form("Gamma_v_TimeStamp_%s",outgoing->GetName())
		                                           ,600,-400,800,bank29->Timestamp() - hit.Timestamp()
                                                           ,2000,0,8000,hit.GetDoppler(GValue::Value("BETA"),&track));

      if(!incoming || !time_energy)
	{return false;}  

      dirname = Form("GretAngCorrs_%s_%s",incoming->GetName(),time_energy->GetName());
      
      TVector3 incBeam = TVector3(0.0,0.0,1.0);
      TVector3 reacPlane = track.Cross(incBeam);
      TVector3 detPlane = (hit.GetPosition()).Cross(incBeam);

      //simple angle between the vectors
      double angle = TMath::ACos((reacPlane.Dot(detPlane))/(reacPlane.Mag()*detPlane.Mag()));

      double reac_phi = reacPlane.Phi();
      if(reac_phi < 0)
	{reac_phi += TMath::TwoPi();}

      double det_phi = detPlane.Phi();
      if(det_phi < 0)
	{det_phi += TMath::TwoPi();}

      double phi1 = reac_phi - det_phi;
      if(phi1 < 0)
	{phi1 += TMath::TwoPi();}

      double phi2 = (2.0*TMath::Pi()-s800->Azita())-TMath::Abs(hit.GetPhi());
      if(phi2 < 0)
	{phi2 += TMath::TwoPi();}
      
      //phi1 and phi2 should be equivalent
      //phi and angle are no longer equivalent (from above); phi runs from 0 to 2pi, angle runs from 0 to pi
      //Phi, Angle, and Azita
      obj.FillHistogram(dirname,Form("Phi2-Phi1_%s",outgoing->GetName()),400,-6.3,6.3,phi2-phi1);
      
      //Phi1
      obj.FillHistogram(dirname,Form("Phi1_%s",outgoing->GetName()),400,-6.3,6.3,phi1);

      obj.FillHistogram(dirname,Form("Gamma(Beta)_v_Phi1_%s",outgoing->GetName())
			                                    ,400,-6.3,6.3,phi1
                                                            ,2000,0,4000,hit.GetDoppler(GValue::Value("BETA")));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_v_Phi1_%s",outgoing->GetName()),400,-6.3,6.3,phi1
			                                          ,2000,0,4000,hit.GetDoppler(GValue::Value("BETA"),&track));
      
      //Phi2
      obj.FillHistogram(dirname,Form("Phi2_%s",outgoing->GetName()),400,-6.3,6.3,phi2);

      
      obj.FillHistogram(dirname,Form("Gamma(Beta)_v_Phi2_%s",outgoing->GetName())
			                                    ,400,-6.3,6.3,phi2
			                                    ,2000,0,4000,hit.GetDoppler(GValue::Value("BETA")));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_v_Phi2_%s",outgoing->GetName())
			                                          ,400,-6.3,6.3,phi2
			                                          ,2000,0,4000,hit.GetDoppler(GValue::Value("BETA"),&track));

      //Angle
      obj.FillHistogram(dirname,Form("Angle_%s",outgoing->GetName()),400,-6.3,6.3,angle);

      obj.FillHistogram(dirname,Form("Gamma(Beta)_v_Angle_%s",outgoing->GetName())
			                                     ,400,-6.3,6.3,angle
                                                             ,2000,0,4000,hit.GetDoppler(GValue::Value("BETA")));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_v_Angle_%s",outgoing->GetName()),400,-6.3,6.3,angle
      			                                           ,2000,0,4000,hit.GetDoppler(GValue::Value("BETA"),&track));
      
      //Azita
      obj.FillHistogram(dirname,Form("Azita_%s",outgoing->GetName()),400,-6.3,6.3,s800->Azita());

      obj.FillHistogram(dirname,Form("Gamma(Beta)_v_Azita_%s",outgoing->GetName())
			                                     ,400,-6.3,6.3,s800->Azita()
                                                             ,2000,0,4000,hit.GetDoppler(GValue::Value("BETA")));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_v_Azita_%s",outgoing->GetName())
			                                           ,400,-6.3,6.3,s800->Azita()
			                                           ,2000,0,4000,hit.GetDoppler(GValue::Value("BETA"),&track));
      
      
      //Quad, Crystal, and Theta-Phi Summaries
      dirname = Form("GretinaSummaries_%s_%s",incoming->GetName(),time_energy->GetName());

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_v_CrystalId_%s",outgoing->GetName())
			,100,0,100,hit.GetCrystalId()
			,2000,0,4000,hit.GetDoppler(GValue::Value("BETA"),&track));
      
      obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_v_QuadNum_%s",outgoing->GetName())
			                                             ,100,0,100,hit.GetHoleNumber()
			                                             ,2000,0,4000,hit.GetDoppler(GValue::Value("BETA"),&track));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_v_Quad(T-P)_%s",outgoing->GetName())
		        ,100,0,100,quadMap[hit.GetHoleNumber()]
			,2000,0,4000,hit.GetDoppler(GValue::Value("BETA"),&track));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA)_v_Quad(T-P)_%s",outgoing->GetName())
			,100,0,100,quadMap[hit.GetHoleNumber()]
			,2000,0,4000,hit.GetDopplerYta(GValue::Value("BETA"),yta,&track));

      for(int i=0; i<30; i++) {
	
	obj.FillHistogram(dirname,Form("Gamma(Beta_%f)_v_Crys(Theta-Phi)_%s",(0.384+(i*.001)),outgoing->GetName())
			                                                    ,38,0,38,detMap[hit.GetCrystalId()]
			                                                    ,2000,0,4000,hit.GetDoppler(0.384+(i*.001)));
      }

      for(int i=0; i<30; i++) {
	
	obj.FillHistogram(dirname,Form("Gamma(Beta_%f&Track)_v_Crys(Theta-Phi)_%s",(0.384+(i*.001)),outgoing->GetName())
			  ,38,0,38,detMap[hit.GetCrystalId()],2000,0,4000,hit.GetDoppler(0.384+(i*.001),&track));
      }

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_v_Crys(Theta-Phi)_%s",outgoing->GetName())
		        ,38,0,38,detMap[hit.GetCrystalId()]
			,2000,0,4000,hit.GetDoppler(GValue::Value("BETA"),&track));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA)_v_Crys(Theta-Phi)_%s",outgoing->GetName())
			,38,0,38,detMap[hit.GetCrystalId()]
			,2000,0,4000,hit.GetDopplerYta(GValue::Value("BETA"),yta,&track));
  
      obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA&DTA)_v_Crys(Theta-Phi)_%s",outgoing->GetName())
      			,38,0,38,detMap[hit.GetCrystalId()]
			,2000,0,4000,hit.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")),yta,&track));


      //Detector (Crystal) Level
      dirname=Form("GammaDetLevel_%s_%s",incoming->GetName(),time_energy->GetName());

      obj.FillHistogram(dirname,Form("Gamma(Beta)_Crys(T-P)%i_%s",detMap[hit.GetCrystalId()],outgoing->GetName())
			,2000,0,4000,hit.GetDoppler(GValue::Value("BETA")));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_Crys(T-P)%i_%s",detMap[hit.GetCrystalId()],outgoing->GetName())
			,2000,0,4000,hit.GetDoppler(GValue::Value("BETA"),&track));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA)_Crys(T-P)%i_%s",detMap[hit.GetCrystalId()],outgoing->GetName())
			,2000,0,4000,hit.GetDopplerYta(GValue::Value("BETA"),yta,&track));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA&DTA)_Crys(T-P)%i_%s",detMap[hit.GetCrystalId()],outgoing->GetName())
			,2000,0,4000,hit.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")),yta,&track));
      
      
      //Azita Correlations
      dirname=Form("AzitaDetLevelCorrs_%s_%s",incoming->GetName(),time_energy->GetName());

      obj.FillHistogram(dirname,Form("Gamma(Beta)_v_Azita_Crys(T-P)%i_%s",detMap[hit.GetCrystalId()],outgoing->GetName())
			                                                 ,1000,0,4000,hit.GetDoppler(GValue::Value("BETA"))
			                                                 ,400,-6.3,6.3,s800->Azita());
			                                                       

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_v_Azita_Crys(T-P)%i_%s",detMap[hit.GetCrystalId()],outgoing->GetName())
			,1000,0,4000,hit.GetDoppler(GValue::Value("BETA"),&track)
			,400,-6.3,6.3,s800->Azita());
      
      //Plane Angle (phi) Correlations
      dirname=Form("PhiDetLevelCorrs_%s_%s",incoming->GetName(),time_energy->GetName());

      obj.FillHistogram(dirname,Form("Gamma(Beta)_v_Phi1_Crys(T-P)%i_%s",detMap[hit.GetCrystalId()],outgoing->GetName())
			                                                ,1000,0,4000,hit.GetDoppler(GValue::Value("BETA"))
			                                                ,400,-6.3,6.3,phi1);

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_v_Phi1_Crys(T-P)%i_%s",detMap[hit.GetCrystalId()],outgoing->GetName())
			,1000,0,4000,hit.GetDoppler(GValue::Value("BETA"),&track)
			,400,-6.3,6.3,phi1);
      
      //YTA Correlations
      dirname=Form("YTADetLevelCorrs_%s_%s",incoming->GetName(),time_energy->GetName());

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_v_YTA_Crys(T-P)%i_%s",detMap[hit.GetCrystalId()],outgoing->GetName())
			,60,-15,15,yta
			,1000,0,4000,hit.GetDoppler(GValue::Value("BETA"),&track));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA)_v_YTA_Crys(T-P)%i_%s",detMap[hit.GetCrystalId()],outgoing->GetName())
			,60,-15,15,yta
			,1000,0,4000,hit.GetDopplerYta(GValue::Value("BETA"),yta,&track));
      
      //DTA Correlations
      dirname=Form("DTADetLevelCorrs_%s_%s",incoming->GetName(),time_energy->GetName());

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA)_v_DTA_Crys(T-P)%i_%s",detMap[hit.GetCrystalId()],outgoing->GetName())
	                ,500,-0.5,0.5,s800->GetDta(),1000,0,4000,hit.GetDopplerYta(GValue::Value("BETA"),yta,&track));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA&DTA)_v_DTA_Crys(T-P)%i_%s",detMap[hit.GetCrystalId()]
			,outgoing->GetName()),500,-0.5,0.5,s800->GetDta()
			,1000,0,4000,hit.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")),yta,&track));
      

      //Quad (Hole) Level 
      dirname=Form("GammaQuadLevel_%s_%s",incoming->GetName(),time_energy->GetName());

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_Quad(T-P)%i_%s",quadMap[hit.GetHoleNumber()],outgoing->GetName())
			,2000,0,4000,hit.GetDoppler(GValue::Value("BETA"),&track));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA)_Quad(T-P)%i_%s",quadMap[hit.GetHoleNumber()],outgoing->GetName())
			,2000,0,4000,hit.GetDopplerYta(GValue::Value("BETA"),yta,&track));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA&DTA)_Quad(T-P)%i_%s",quadMap[hit.GetHoleNumber()],outgoing->GetName())
			,2000,0,4000,hit.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")),yta,&track));

      //Azita Correlations
      dirname=Form("AzitaQuadLevelCorrs_%s_%s",incoming->GetName(),time_energy->GetName());

      obj.FillHistogram(dirname,Form("Gamma(Beta)_v_Azita_Quad(T-P)%i_%s",quadMap[hit.GetHoleNumber()],outgoing->GetName())
			                                                 ,1000,0,4000,hit.GetDoppler(GValue::Value("BETA"))
			                                                 ,400,-6.3,6.3,s800->Azita());
			                                                       

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_v_Azita_Quad(T-P)%i_%s",quadMap[hit.GetHoleNumber()],outgoing->GetName())
			,1000,0,4000,hit.GetDoppler(GValue::Value("BETA"),&track)
		        ,400,-6.3,6.3,s800->Azita());
      
      //Plane Angle (phi) Correlations
      dirname=Form("PhiQuadLevelCorrs_%s_%s",incoming->GetName(),time_energy->GetName());

      obj.FillHistogram(dirname,Form("Gamma(Beta)_v_Phi1_Quad(T-P)%i_%s",quadMap[hit.GetHoleNumber()],outgoing->GetName())
		        ,1000,0,4000,hit.GetDoppler(GValue::Value("BETA"))
			,400,-6.3,6.3,phi1);

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_v_Phi1_Quad(T-P)%i_%s",quadMap[hit.GetHoleNumber()],outgoing->GetName())
			,1000,0,4000,hit.GetDoppler(GValue::Value("BETA"),&track)
			,400,-6.3,6.3,phi1);
      
      //YTA Correlations
      dirname=Form("YTAQuadLevelCorrs_%s_%s",incoming->GetName(),time_energy->GetName());
      
      obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_v_YTA_Quad(T-P)%i_%s",quadMap[hit.GetHoleNumber()],outgoing->GetName())
			,60,-15,15,yta
			,1000,0,4000,hit.GetDoppler(GValue::Value("BETA"),&track));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA)_v_YTA_Quad(T-P)%i_%s",quadMap[hit.GetHoleNumber()]
			,outgoing->GetName())
			,60,-15,15,yta
			,1000,0,4000,hit.GetDopplerYta(GValue::Value("BETA"),yta,&track));
      
      //DTA Correlations
      dirname=Form("DTAQuadLevelCorrs_%s_%s",incoming->GetName(),time_energy->GetName());

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA)_v_DTA_Quad(T-P)%i_%s",quadMap[hit.GetHoleNumber()]
			,outgoing->GetName())
	                ,500,-0.5,0.5,s800->GetDta()
			,1000,0,4000,hit.GetDopplerYta(GValue::Value("BETA"),yta,&track));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA&DTA)_v_DTA_Quad(T-P)%i_%s",quadMap[hit.GetHoleNumber()]
			,outgoing->GetName())
      	     		,500,-0.5,0.5,s800->GetDta()
			,1000,0,4000,hit.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")),yta,&track));

      //Crys Theta Level
      dirname=Form("GammaCrysThetaLevel_%s_%s",incoming->GetName(),time_energy->GetName());

      obj.FillHistogram(dirname,Form("Gamma(Beta)_Crys(Theta)%i_%s",crysThetaMap[hit.GetCrystalId()],outgoing->GetName())
			,2000,0,4000,hit.GetDoppler(GValue::Value("BETA")));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_Crys(Theta)%i_%s",crysThetaMap[hit.GetCrystalId()],outgoing->GetName())
			,2000,0,4000,hit.GetDoppler(GValue::Value("BETA"),&track));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA)_Crys(Theta)%i_%s",crysThetaMap[hit.GetCrystalId()]
			,outgoing->GetName())
			,2000,0,4000,hit.GetDopplerYta(GValue::Value("BETA"),yta,&track));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA&DTA)_Crys(Theta)%i_%s",crysThetaMap[hit.GetCrystalId()]
			,outgoing->GetName())
			,2000,0,4000,hit.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")),yta,&track));
      
      
      //Azita Correlations
      dirname=Form("AzitaCrysThetaLevelCorrs_%s_%s",incoming->GetName(),time_energy->GetName());

      obj.FillHistogram(dirname,Form("Gamma(Beta)_v_Azita_Crys(Theta)%i_%s",crysThetaMap[hit.GetCrystalId()],outgoing->GetName())
			                                                   ,1000,0,4000,hit.GetDoppler(GValue::Value("BETA"))
			                                                   ,400,-6.3,6.3,s800->Azita());
			                                                       

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_v_Azita_Crys(Theta)%i_%s",crysThetaMap[hit.GetCrystalId()]
		        ,outgoing->GetName())
			,1000,0,4000,hit.GetDoppler(GValue::Value("BETA"),&track)
			,400,-6.3,6.3,s800->Azita());
      
      //Plane Angle (phi) Correlations
      dirname=Form("PhiCrysThetaLevelCorrs_%s_%s",incoming->GetName(),time_energy->GetName());

      obj.FillHistogram(dirname,Form("Gamma(Beta)_v_Phi1_Crys(Theta)%i_%s",crysThetaMap[hit.GetCrystalId()],outgoing->GetName())
			                                                  ,1000,0,4000,hit.GetDoppler(GValue::Value("BETA"))
			                                                  ,400,-6.3,6.3,phi1);

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_v_Phi1_Crys(Theta)%i_%s",crysThetaMap[hit.GetCrystalId()]
			,outgoing->GetName())
			,1000,0,4000,hit.GetDoppler(GValue::Value("BETA"),&track)
			,400,-6.3,6.3,phi1);
      
      //YTA Correlations
      dirname=Form("YTACrysThetaLevelCorrs_%s_%s",incoming->GetName(),time_energy->GetName());

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track)_v_YTA_Crys(Theta)%i_%s",crysThetaMap[hit.GetCrystalId()]
			,outgoing->GetName())
			,60,-15,15,yta
			,1000,0,4000,hit.GetDoppler(GValue::Value("BETA"),&track));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA)_v_YTA_Crys(Theta)%i_%s",crysThetaMap[hit.GetCrystalId()]
			,outgoing->GetName())
			,60,-15,15,yta
			,1000,0,4000,hit.GetDopplerYta(GValue::Value("BETA"),yta,&track));
      
      //DTA Correlations
      dirname=Form("DTACrysThetaLevelCorrs_%s_%s",incoming->GetName(),time_energy->GetName());

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA)_v_DTA_Crys(Theta)%i_%s",crysThetaMap[hit.GetCrystalId()]
		        ,outgoing->GetName())
	                ,500,-0.5,0.5,s800->GetDta()
			,1000,0,4000,hit.GetDopplerYta(GValue::Value("BETA"),yta,&track));

      obj.FillHistogram(dirname,Form("Gamma(Beta&Track&YTA&DTA)_v_DTA_Crys(Theta)%i_%s",crysThetaMap[hit.GetCrystalId()]
			,outgoing->GetName()),500,-0.5,0.5,s800->GetDta()
			,1000,0,4000,hit.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")),yta,&track));

    }//end if(hit.CoreEnergy > 100) 

  return true;

}

bool HandleGretina_Gated(TRuntimeObjects &obj,TCutG *incoming, TCutG* outgoing, TCutG* time_energy) {
  TGretina *gretina  = obj.GetDetector<TGretina>(); 
  TS800 *s800 = obj.GetDetector<TS800>();
  TBank29 *bank29 = obj.GetDetector<TBank29>();
  
  if(!gretina || !s800 || !bank29)
    {return false;}

  std::string dirname = "Gretina";

  ///////////////////////////////////////////////////
  //incoming gates on Xfp-E1 vs. OBJ-E1 Uncorrected//
  ///////////////////////////////////////////////////
  
  if(incoming) {
    bool passed = false;  
    if(incoming->IsInside(s800->GetMTof().GetCorrelatedObjE1(),s800->GetMTof().GetCorrelatedXfpE1())) {  
          passed=true;
        }
    if(!passed)
      return false;
  }

  ///////////////////////////////////////////////////////
  //incoming gates on Xfp - E1 vs. OBJ - E1 (Corrected)//
  ///////////////////////////////////////////////////////
  /*
  if(incoming) {
    bool passed = false;  
    if(incoming->IsInside(s800->GetMTofObjE1(),s800->GetMTofXfpE1())) {
          dirname=Form("S800_%s",incoming->GetName());
          passed=true;
        }
    if(!passed)
      return false;
  }*/
    
  ///////////////////////////////////////////////////////////////////
  //OUTGOING GATES on Ion Chamber Charge vs Object - E1 (Corrected)//
  ///////////////////////////////////////////////////////////////////
  
  if(!outgoing)
    {return false;}
  
  if(!outgoing->IsInside(s800->GetMTofObjE1(),s800->GetIonChamber().Charge()))
    {return false;}

  if(!time_energy) {
   for(unsigned int i=0;i<gretina->Size();i++) {

    TGretinaHit hit = gretina->GetGretinaHit(i);
    Gretina_PIDOUTGatedSpectra(obj,hit,s800,bank29,dirname,incoming,outgoing);

   }
   /*
   for(int i=0;i<gretina->AddbackSize();i++) {

    TGretinaHit hit = gretina->GetAddbackHit(i);
    Gretina_PIDOUTGatedSpectra(obj,hit,s800,bank29,"GretinaAddback",incoming,outgoing);

   }*/
  }

  //////////////////////////////////////////////////////////////////////////
  //Time-Energy Gate on Gamma Enery vs Banke29 TimeStamp - GretinaHit Time//
  //////////////////////////////////////////////////////////////////////////
  
  if(time_energy) {
   for(unsigned int i=0;i<gretina->Size();i++) {
     
    TGretinaHit hit = gretina->GetGretinaHit(i);
    TVector3 track = s800->Track();
    if(!time_energy->IsInside(bank29->Timestamp() - hit.GetTime(),hit.GetDoppler(GValue::Value("BETA"),&track)))
      {return false;}
    
    Gretina_PIDOUTGatedSpectra(obj,hit,s800,bank29,dirname,incoming,outgoing,time_energy);  
    
   }

   /*
   for(int i=0;i<gretina->AddbackSize();i++) {
     
    TGretinaHit hit = gretina->GetAddbackHit(i);
    TVector3 track = s800->Track();
    if(!time_energy->IsInside(bank29->Timestamp() - hit.GetTime(),hit.GetDoppler(GValue::Value("BETA"),&track)))
      {return false;}

      Gretina_PIDOUTGatedSpectra(obj,hit,s800,bank29,"GretinaAddback",incoming,outgoing,time_energy);  
    
   }*/
  }
   
  return true;

}

extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
  //TGretina *gretina = obj.GetDetector<TGretina>();
  //TBank29  *bank29  = obj.GetDetector<TBank29>();
  //TS800    *s800    = obj.GetDetector<TS800>();

  TList    *list    = &(obj.GetObjects());
  int numobj = list->GetSize();  
   
  TList *gates = &(obj.GetGates());

  if(gates_loaded!=gates->GetSize()) {
    TIter iter(gates);
    while(TObject *obj = iter.Next()) {
      GCutG *gate = (GCutG*)obj;
      std::string tag = gate->GetTag();
      if(!tag.compare("incoming")) {
        incoming_cuts.push_back(gate);
        std::cout << "incoming: << " << gate->GetName() << std::endl;
      } else if(!tag.compare("outgoing")) {
        outgoing_cuts.push_back(gate);
        std::cout << "outgoing: << " << gate->GetName() << std::endl;
      } else if(!tag.compare("time_energy")) {
        time_energy_cuts.push_back(gate);
        std::cout << "time_energy: << " << gate->GetName() << std::endl;
      }
      gates_loaded++;
    }
  }

  
  for(size_t i=0;i<incoming_cuts.size();i++) {
   for(size_t j=0;j<outgoing_cuts.size();j++) {
     HandleTiming(obj,incoming_cuts.at(i),outgoing_cuts.at(j));
   }
  }
  
  //HandleS800(obj);
  /*
  for(size_t i=0;i<incoming_cuts.size();i++) {
    for(size_t j=0;j<outgoing_cuts.size();j++) {
      s800Singles(obj,incoming_cuts.at(i),outgoing_cuts.at(j));
    }
  }
  */  
  
  for(size_t i=0;i<incoming_cuts.size();i++) {
   for(size_t j=0;j<outgoing_cuts.size();j++) {
     HandleS800_Gated(obj,incoming_cuts.at(i),outgoing_cuts.at(j));
   }
  }
  
  HandleGretina(obj);
  
  for(size_t i=0;i<incoming_cuts.size();i++) {
   for(size_t j=0;j<outgoing_cuts.size();j++) {
    for(size_t k=0;k<time_energy_cuts.size();k++) {
     HandleGretina_Gated(obj,incoming_cuts.at(i),outgoing_cuts.at(j),time_energy_cuts.at(k));
    }
   }
  }
   

//      obj.FillHistogram(dirname,histname,
//          600,-600,600,bank29->Timestamp()-hit.GetTime(),
//          2000,0,4000,hit.GetCoreEnergy());


  

//  TList *gates = &(obj.GetGates());
//  if(gretina) {
//    for(unsigned int i=0;i<gretina->Size();i++) {
//      TGretinaHit hit = gretina->GetGretinaHit(i);
//      histname = "Gretina_Bank29_time";
//      obj.FillHistogram(dirname,histname,
//          600,-600,600,bank29->Timestamp()-hit.Timestamp(),
//          2000,0,4000,hit.GetCoreEnergy());
//      histname = "Gretina_t0_Bank29_time";
//      obj.FillHistogram(dirname,histname,
//          600,-600,600,bank29->Timestamp()-hit.GetTime(),
//          2000,0,4000,hit.GetCoreEnergy());
//    }
//  }

  if(numobj!=list->GetSize())
    list->Sort();
}





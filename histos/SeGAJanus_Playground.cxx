#include "TRuntimeObjects.h"

#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <limits>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TObject.h"
#include "TRandom.h"

#include "GValue.h"
#include "TJanusDDAS.h"
#include "TNSCLScalers.h"
#include "TSega.h"
#include "GCutG.h"

#include "TVector3.h"
#include "TReaction.h"
#include "TSRIM.h"
#include "TRandom.h"
#include "TNucleus.h"

static const bool BEFORE = GValue::Value("BEFORE");

static std::map<int,std::vector<std::string>> system_map = {
                                                            {0,{"87Rb","208Pb","rb87_in_pb208","pb208_in_pb208"}},
                                                            {1,{"106Cd","208Pb","cd106_in_pb208","pb208_in_pb208"}},
                                                            {2,{"106Cd","48Ti","cd106_in_ti48","ti48_in_ti48"}},
                                                            {3,{"106Cd","184W","cd106_in_ti48","w184_in_ti48"}}
                                                           };

static std::map<int,int> D1_RingMap = {
                                       {119,1},{118,2},{117,3},{116,4},{115,5},{114,6},{113,7},{112,8},{96,9},{97,10},{98,11},
				       {99,12},{100,13},{101,14},{102,15},{103,16},{104,17},{105,18},{106,19},{107,20},{108,21},
				       {109,22},{110,23},{111,24}
                                      };

static std::map<int,std::array<int,2>> nMap = {{1,{2,8}},{2,{1,3}},{3,{2,4}},{4,{3,5}},{5,{4,6}},{6,{5,7}},
					       {7,{6,8}},{8,{7,1}},{9,{10,16}},{10,{9,11}},{11,{10,12}},
					       {12,{11,13}},{13,{12,14}},{14,{13,15}},{15,{14,16}},
					       {16,{15,9}}};

static int min_level = 1;
static int max_level = 1;

static long prev_ts = 0;
static long prev_chan_ts = 0;

class BinCenters {
    public:
        BinCenters(int bins, double low, double high)
            : bins(bins), low(low), high(high) { }

        class iterator {
            public:
                iterator(BinCenters& axis, int binnum)
                    : axis(axis), binnum(binnum) { }

                double operator*() const;

                bool operator==(const iterator& other) const {
                    return
                        (&axis == &other.axis) &&
                        (binnum == other.binnum);
                }

                bool operator!=(const iterator& other) const {
                    return !(*this == other);
                }

                iterator& operator++() {
                    binnum++;
                    return *this;
                }

            private:
                BinCenters& axis;
                int binnum;
        };
        friend class BinCenters::iterator;

        iterator begin() {
            return iterator(*this, 0);
        }

        iterator end() {
            return iterator(*this, bins);
        }

    private:
        int bins;
        double low;
        double high;
};

double BinCenters::iterator::operator*() const {
    return axis.low + (binnum+0.5) * (axis.high-axis.low)/axis.bins;
}

void MakeSega(TRuntimeObjects& obj, TSega& sega) {

  std::string dirname = "SeGA";

  obj.FillHistogram(dirname,"Mult",20,0,20,sega.Size());

  double seg_sum_charge = 0;
  double seg_sum_energy = 0;

  //for(auto hit : sega) {
  for(unsigned int i=0;i<sega.Size();i++) {

    TSegaHit hit = sega.GetSegaHit(i);

    obj.FillHistogram(dirname,"GammaEnergy",4000,0,4000,hit.GetEnergy());
    
    //obj.FillHistogram("SeGA_Detectors",Form("Energy_det%02d",hit.GetDetnum()),4000,0,4000,hit.GetEnergy());
    //obj.FillHistogram("SeGA_Detectors",Form("Charge_det%02d",hit.GetDetnum()),9000,0,36000,hit.Charge());

    obj.FillHistogram(dirname,"SeGA_Energy_Summary",18,0,18,hit.GetDetnum(),4000,0,4000,hit.GetEnergy());
    
    //long cc_ts = hit.Timestamp();

    double num_segs = hit.GetNumSegments();
    seg_sum_charge = 0;
    seg_sum_energy = 0;

    //obj.FillHistogram("SeGA_Detectors",Form("SegMult_det%02d",hit.GetDetnum()),33,0,33,num_segs);

    //obj.FillHistogram("SeGA_Detectors",Form("GammaEnergy_v_SegMult_det%02d",hit.GetDetnum()),
    //		      33,0,33,num_segs,
    //		      2000,0,4000,hit.GetEnergy());

    if(num_segs == 0) {
      obj.FillHistogram(dirname,"GammaEnergy_NoSegs",4000,0,4000,hit.GetEnergy());
      obj.FillHistogram("SeGA_Detectors_NoSegs",Form("GammaEnergy_det%02d",hit.GetDetnum()),4000,0,4000,hit.GetEnergy());
    }
    else if(num_segs > 0) {
      obj.FillHistogram(dirname,"GammaEnergy_YesSegs",4000,0,4000,hit.GetEnergy());
      obj.FillHistogram("SeGA_Detectors_YesSegs",Form("GammaEnergy_det%02d",hit.GetDetnum()),4000,0,4000,hit.GetEnergy());
    }
    
    for(int j=0;j<num_segs;j++) {

      TSegaSegmentHit seg = hit.GetSegment(j);
      //long seg_ts = seg.Timestamp();
       seg_sum_charge += seg.Charge();
       seg_sum_energy += seg.GetEnergy();

       //obj.FillHistogram(dirname,"Core_Seg_tdiff",2000,-20000,20000,cc_ts-seg_ts);

       //obj.FillHistogram(dirname,"CoreEn_v_Core_Seg_tdiff",2000,-20000,20000,cc_ts-seg_ts,
       //		 2000,0,4000,hit.GetEnergy());

    }//end segment loop

    if(hit.Charge() <= 0) {

      obj.FillHistogram("SeGA_Detectors_NoCoreCharge",Form("SegMult_det%02d",hit.GetDetnum()),33,0,33,num_segs);
      obj.FillHistogram("SeGA_Detectors_NoCoreCharge",Form("SegSumCharge_det%02d",hit.GetDetnum()),33000,0,33000,seg_sum_charge);
      
      obj.FillHistogram("SeGA_Detectors_NoCoreCharge",Form("SegSumEnergy_det%02d",hit.GetDetnum()),4000,0,4000,seg_sum_energy);
      obj.FillHistogram(dirname,"SegSumEnergy_NoCoreCharge",4000,0,4000,seg_sum_energy);
      
    }
    /*
    obj.FillHistogram("SeGA_Core_Seg",Form("SegSum_v_CoreCharge_det%02d",hit.GetDetnum()),
		      2048,0,32768,hit.Charge(),2048,0,32768,seg_sum_charge);

    obj.FillHistogram("SeGA_Core_Seg",Form("SegSum_v_CoreEnergy_det%02d",hit.GetDetnum()),
		      2000,0,4000,hit.GetEnergy(),3000,0,30000,seg_sum_energy);

    
    if(sega.Size() == 1) {

      obj.FillHistogram("SeGA_Core_Seg",Form("SegSum_v_CoreCharge_Mult1_det%02d",hit.GetDetnum()),
		        2048,0,32768,hit.Charge(),2048,0,32768,seg_sum_charge);

      obj.FillHistogram("SeGA_Core_Seg",Form("SegSum_v_CoreEnergy_Mult1_det%02d",hit.GetDetnum()),
		      2000,0,4000,hit.GetEnergy(),2500,0,20000,seg_sum_energy);
    }
    */

    for(unsigned int j=0;j<sega.Size();j++) {

      if(i==j) {
	continue;
      }
      
      TSegaHit hit1 = sega.GetSegaHit(j);

      obj.FillHistogram(dirname,"GammaGammaMatrix",
		        2500,0,4000,hit.GetEnergy(),
		        2500,0,4000,hit1.GetEnergy());

      if(hit.GetDetnum()!=16 && hit1.GetDetnum()!=16) {
	obj.FillHistogram(dirname,"GammaGammaMatrix_NoDet16",
			  2500,0,4000,hit.GetEnergy(),
			  2500,0,4000,hit1.GetEnergy());
      }
      
    } //end second hit loop
   
  }//end hit loop
}

void BasicJanusSpectra(TRuntimeObjects& obj, TJanusDDASHit& hit, std::string dirname, std::string h_suf) {

  obj.FillHistogram(dirname,Form("%sCharge_Summary_Hit",h_suf.c_str()),
                    128,0,128,hit.GetFrontChannel(),
                    3000,0,30000,hit.Charge());
  obj.FillHistogram(dirname,Form("%sCharge_Summary_Hit",h_suf.c_str()),
                    128,0,128,hit.GetBackChannel(),
                    3000,0,30000,hit.GetBackHit().Charge());

  obj.FillHistogram(dirname,Form("%sRing_det%d",h_suf.c_str(),hit.GetDetnum()),26,0,26,hit.GetRing());
  obj.FillHistogram(dirname,Form("%sSector_det%d",h_suf.c_str(),hit.GetDetnum()),34,0,34,hit.GetSector());

  obj.FillHistogram(dirname,Form("%sRing_v_Sector_det%d",h_suf.c_str(),hit.GetDetnum()),
		    26,0,26,hit.GetRing(),34,0,34,hit.GetSector());

  double theta = hit.GetPosition(BEFORE).Theta()*TMath::RadToDeg();
  double phi = hit.GetPosition(BEFORE).Phi();
  double phi1 = hit.GetPosition(BEFORE).Phi();
  if(phi < 0) {
    phi+=TMath::TwoPi();
  }
  phi*=TMath::RadToDeg();
  double perp = hit.GetPosition(BEFORE).Perp();

  double shift = 0.5*TMath::TwoPi()/32.0;
  obj.FillHistogram(dirname,Form("%sPerp_v_Phi_det%d",h_suf.c_str(),hit.GetDetnum()),
		    32,-TMath::Pi()+shift,TMath::Pi()+shift,phi1,24,1.1,3.5,perp);

  obj.FillHistogram(dirname,Form("%sTheta_det%d",h_suf.c_str(),hit.GetDetnum()),360,0,180,theta);
  obj.FillHistogram(dirname,Form("%sPhi_det%d",h_suf.c_str(),hit.GetDetnum()),108,-360,720,phi);
  
  obj.FillHistogram(dirname,Form("%sTheta_v_Phi_det%d",h_suf.c_str(),hit.GetDetnum()),
		    108,-360,720,phi,
		    360,0,180,theta);

  obj.FillHistogram(dirname,Form("%sRing_v_Phi_det%d",h_suf.c_str(),hit.GetDetnum()),
		    108,-360,720,phi,
		    26,0,26,hit.GetRing());

  obj.FillHistogram(dirname,Form("%sY_v_X_det%d",h_suf.c_str(),hit.GetDetnum()),
		    200,-4,4,hit.GetPosition(BEFORE).X(),
		    200,-4,4,hit.GetPosition(BEFORE).Y());

  obj.FillHistogram(dirname,Form("%sSectorCharge_v_Ring_Hit_det%d",h_suf.c_str(),hit.GetDetnum()),
		    26,0,26,hit.GetRing(),
		    3000,0,30000,hit.Charge());

  obj.FillHistogram(dirname,Form("%sRingCharge_v_Ring_Hit_det%d",h_suf.c_str(),hit.GetDetnum()),
		    26,0,26,hit.GetRing(),
		    3000,0,30000,hit.GetBackHit().Charge());

  obj.FillHistogram(dirname,Form("%sSectorCharge_v_RingCharge_Hit_det%d",h_suf.c_str(),hit.GetDetnum()),
		    2000,0,30000,hit.GetBackHit().Charge(),
		    2000,0,30000,hit.Charge());

  obj.FillHistogram(dirname,Form("%sSectorCharge_v_Theta_Hit_det%d",h_suf.c_str(),hit.GetDetnum()),
		    500,0,180,theta,
		    3000,0,30000,hit.Charge());

  obj.FillHistogram(dirname,Form("%sRingCharge_v_Theta_Hit_det%d",h_suf.c_str(),hit.GetDetnum()),
		    500,0,180,theta,
		    3000,0,30000,hit.GetBackHit().Charge());
  
  auto tdiff = hit.GetBackHit().Timestamp() - hit.Timestamp();
   
  obj.FillHistogram(dirname,Form("%sFrontBack_Tdiff_Hit_det%d",h_suf.c_str(),hit.GetDetnum()),2750,-11000,11000,tdiff);

  auto cdiff = hit.GetBackHit().Charge() - hit.Charge();
  
  obj.FillHistogram(dirname,Form("Cdiff_Hit_det%d",hit.GetDetnum()),2000,0,20000,cdiff);

  obj.FillHistogram(dirname,Form("%sCdiff_v_Tdiff_Hit_det%d",h_suf.c_str(),hit.GetDetnum()),
		    2000,-10000,10000,tdiff,
		    2000,0,20000,cdiff);

  obj.FillHistogram(dirname,Form("%sRingCharge_v_Tdiff_Hit_det%d",h_suf.c_str(),hit.GetDetnum()),
		    2000,-10000,10000, tdiff,
		    2000,0,30000,hit.GetBackHit().Charge());

  obj.FillHistogram(dirname,Form("%sRingCharge_v_Cdiff_Hit_det%d",h_suf.c_str(),hit.GetDetnum()),
		    2000,0,20000,cdiff,
		    2000,0,30000,hit.GetBackHit().Charge());

  obj.FillHistogram(dirname,Form("%sSectorCharge_v_Tdiff_Hit_det%d",h_suf.c_str(),hit.GetDetnum()),
		    2000,-10000,10000,tdiff,
		    2000,0,30000,hit.Charge());

  obj.FillHistogram(dirname,Form("%sSectorCharge_v_Cdiff_Hit_det%d",h_suf.c_str(),hit.GetDetnum()),
		    2000,0,20000,cdiff,
		    2000,0,30000,hit.Charge());
 
}

void SegmentTaggedJanusSpectra(TRuntimeObjects& obj, TJanusDDASHit& hit, std::string d_pref,
			       std::string d_suf,std::string h_suf) {
  /*
  obj.FillHistogram(Form("%sJanus_Rings_%s",d_pref.c_str(),d_suf.c_str()),
		    Form("%sCharge_Summary_Hit_det%d_Ring%02d",h_suf.c_str(),hit.GetDetnum(),hit.GetRing()),
		    128,0,128,hit.GetFrontChannel(),
		    3000,0,30000,hit.Charge());
  
  obj.FillHistogram(Form("%sJanus_Rings_%s",d_pref.c_str(),d_suf.c_str()),
		    Form("%sCharge_Summary_Hit_det%d_Ring%02d",h_suf.c_str(),hit.GetDetnum(),hit.GetRing()),
		    128,0,128,hit.GetBackChannel(),
		    3000,0,30000,hit.GetBackHit().Charge());

  obj.FillHistogram(Form("%sJanus_Sectors_%s",d_pref.c_str(),d_suf.c_str()),
		    Form("%sSectorCharge_v_RingCharge_det%d_Sect%02d",h_suf.c_str(),hit.GetDetnum(),hit.GetSector()),
		    2000,0,30000,hit.GetBackHit().Charge(),
		    2000,0,30000,hit.Charge());

  double phi = hit.GetPosition(BEFORE).Phi();
  double perp = hit.GetPosition(BEFORE).Perp();
  
  double shift = 0.5*TMath::TwoPi()/32.0;
  obj.FillHistogram(Form("%sJanus_Rings_%s",d_pref.c_str(),d_suf.c_str()),
		    Form("%sPerp_v_Phi_det%d_Ring%02d",h_suf.c_str(),hit.GetDetnum(),hit.GetRing()),
		    32,-TMath::Pi()+shift,TMath::Pi()+shift,phi,24,1.1,3.5,perp);
  

  obj.FillHistogram(Form("%sJanus_Sectors_%s",d_pref.c_str(),d_suf.c_str()),
		    Form("%sSectorCharge_v_Ring_det%d_Sect%02d",
			 h_suf.c_str(),hit.GetDetnum(),hit.GetSector()),
		    26,0,26,hit.GetRing(),
		    2000,0,30000,hit.Charge());
  */
  
  obj.FillHistogram(Form("%sJanus_Rings_%s",d_pref.c_str(),d_suf.c_str()),
		    Form("%sSector_Hit_det%d_Ring%02d",h_suf.c_str(),hit.GetDetnum(),hit.GetRing()),
		    34,0,34,hit.GetSector());
}

void LeveledJanusSpectra(TRuntimeObjects& obj, TJanusDDAS& janus, int lvl) {

  std::string d_pref = "";
  
  if(lvl==0)
    {d_pref="Total_";}
  else if(lvl==1)
    {d_pref="Primary_";}
  else if(lvl==2)
    {d_pref="Secondary_";}
  else if(lvl==3)
    {d_pref="Tertiary_";}

  std::string dirname = d_pref + "Janus_Ungated";

  //if(janus.GetDetectorChannels(1).size()) {
  //obj.FillHistogram(dirname,"UsableChanMult_det1",100,0,100,janus.GetDetectorChannels(1).size());
  obj.FillHistogram(dirname,"All_HitMult_det1",50,0,50,janus.GetDetectorHits(1,lvl).size());   
    //}

    //if(janus.GetDetectorChannels(0).size()) {
  //obj.FillHistogram(dirname,"UsableChanMult_det0",100,0,100,janus.GetDetectorChannels(0).size());
  obj.FillHistogram(dirname,"All_HitMult_det0",50,0,50,janus.GetDetectorHits(0,lvl).size());   
    //}
  
  for(size_t i=0; i<janus.GetUnusedChannels(lvl).size(); i++) {
    auto &chan = janus.GetUnusedChannels(lvl).at(i);
    obj.FillHistogram(dirname,"Charge_Summary_UnusedChans",
		      128,0,128,chan.GetFrontChannel(),
		      3500,0,35000,chan.Charge());
    
    if(chan.GetDetnum() == 1 && chan.IsRing()) {
      obj.FillHistogram(dirname,"Charge_Summary_UnusedD1RingChans",
			26,0,26,D1_RingMap[chan.GetFrontChannel()],
			3500,0,35000,chan.Charge());
    }
    
    for(size_t j=i+1;j<janus.GetUnusedChannels(lvl).size();j++) {
      auto &chan1 = janus.GetUnusedChannels(lvl).at(j);
      
      if(chan1.GetDetnum() == chan.GetDetnum()) {
        obj.FillHistogram(dirname,Form("Tdiff_UnusedChans_det%d",chan.GetDetnum()),
	  		  2750,-11000,11000,chan.Timestamp()-chan1.Timestamp());
      }
        
    }//end second janus chan loop 
  }//end janus chan loop

  for(size_t i=0;i<janus.GetUnusedSectors(lvl).size();i++) {
    for(size_t j=0;j<janus.GetUnusedRings(lvl).size();j++) {
      auto &sect = janus.GetUnusedSectors(lvl).at(i);
      auto &ring = janus.GetUnusedRings(lvl).at(j);

      if(ring.GetDetnum() == sect.GetDetnum()) {
        obj.FillHistogram(dirname,Form("SectorCharge_v_RingCharge_UnusedChans_det%d",sect.GetDetnum()),
		          2000,0,30000,ring.Charge(),
		          2000,0,30000,sect.Charge());

	obj.FillHistogram(dirname,Form("Sector_Ring_TDiff_UnusedChans_det%d",sect.GetDetnum()),
			  2750,-11000,11000,sect.Timestamp()-ring.Timestamp());
      }
    }
  }

  
  for(size_t i=0;i<janus.GetHits(lvl).size();i++) {
    TJanusDDASHit& hit = janus.GetHits(lvl).at(i);
    BasicJanusSpectra(obj,hit,dirname,"All_");  

    if(hit.GetDetnum() == 1) {

      SegmentTaggedJanusSpectra(obj,hit,d_pref,"Ungated","All_");
      
      obj.FillHistogram(dirname,Form("All_SectorCharge_v_Ring_Hit_det1_Mult%d",(int)janus.GetDetectorHits(1,lvl).size()),
			26,0,26,hit.GetRing(),
			3000,0,30000,hit.Charge());
    }
    
    for(size_t j=i+1;j<janus.GetHits(lvl).size();j++) {
      TJanusDDASHit& hit1 = janus.GetHits(lvl).at(j);

      if(hit.GetDetnum() == hit1.GetDetnum()) {

	obj.FillHistogram(dirname,Form("All_MultiHit_TDiff_det%d",hit.GetDetnum()),
			  2750,-11000,11000,hit.Timestamp()-hit1.Timestamp());

	obj.FillHistogram(dirname,Form("All_MultiHit_TDiff_det%d_Mult%d",hit.GetDetnum(),
				       (int)janus.GetDetectorHits(hit.GetDetnum(),lvl).size()),
			  2750,-11000,11000,hit.Timestamp()-hit1.Timestamp());

      }
    } 
  }

  for(size_t i=0;i<janus.GetDoubleHits(lvl).size();i++) {
    TJanusDDASHit& hit = janus.GetDoubleHits(lvl).at(i);
    BasicJanusSpectra(obj,hit,dirname,"TwoHit_");
  }

  for(size_t i=0;i<janus.GetSingleHits(lvl).size();i++) {
    TJanusDDASHit& hit = janus.GetSingleHits(lvl).at(i);
    BasicJanusSpectra(obj,hit,dirname,"Singles_");    
  }

  for(size_t i=0;i<janus.GetAddbackHits(lvl).size();i++) {
    TJanusDDASHit& hit = janus.GetAddbackHits(lvl).at(i);
    BasicJanusSpectra(obj,hit,dirname,"Addback_");  
  }
  
}

void MakeJanus(TRuntimeObjects& obj, TJanusDDAS& janus) {
  
  std::string dirname = "Janus_Ungated";
  /*
  int size1 = janus.GetDetectorHits(1,1).size();
  int size2 = janus.GetDetectorHits(1,2).size();
  int size3 = janus.GetDetectorHits(1,3).size();
  
  int code = 0;
  if(size1 && !size2 && !size3)
    {code=1;}
  else if(!size1 && size2 && !size3)
    {code=2;}
  else if(!size1 && !size2 && size3)
    {code=3;}
  else if(size1 && size2 && !size3)
    {code=4;}
  else if(size1 && !size2 && size3)
    {code=5;}
  else if(!size1 && size2 && size3)
    {code=6;}
  else if(size1 && size2 && size3)
    {code=7;}

  obj.FillHistogram(dirname,"Hit_Level_Code_det1",10,0,10,code);

  int chan_diff = janus.GetAllChannels().size() - janus.GetChannels().size() - janus.GetBadChargeChannels().size();
  
  obj.FillHistogram(dirname,"Total_Channel_Counting",100,-50,50,chan_diff);

  int used_chans1 = janus.GetChannels().size() - janus.GetUnusedChannels1().size();
  int num_chans1 = janus.GetSingleHits1().size()*2.0 + janus.GetDoubleHits1().size()*(3.0/2.0)
                   + janus.GetAddbackHits(1).size()*3.0;
  
  obj.FillHistogram(dirname,"Primary_Hit_Chan_Counting",100,-50,50,used_chans1-num_chans1);

  int used_chans2 = janus.GetUnusedChannels1().size() - janus.GetUnusedChannels2().size();
  int num_chans2 = janus.GetSingleHits2().size()*2.0 + janus.GetDoubleHits2().size()*(3.0/2.0)
                   + janus.GetAddbackHits(2).size()*3.0;
  
  obj.FillHistogram(dirname,"Secondary_Hit_Chan_Counting",100,-50,50,used_chans2-num_chans2);

  int used_chans3 = janus.GetUnusedChannels2().size() - janus.GetUnusedChannels3().size();
  int num_chans3 = janus.GetSingleHits3().size()*2.0 + janus.GetDoubleHits3().size()*(3.0/2.0)
                   + janus.GetAddbackHits(3).size()*3.0;
  
  obj.FillHistogram(dirname,"Tertiary_Hit_Chan_Counting",100,-50,50,used_chans3-num_chans3);
  */

  if(prev_ts > 0) {
    obj.FillHistogram(dirname,"Janus_Event_TDiff_zoom",10000,0,10000000,janus.Timestamp()-prev_ts);
    obj.FillHistogram(dirname,"Janus_Event_TDiff",100000,0,100000000,janus.Timestamp()-prev_ts);
    obj.FillHistogram(dirname,"Janus_Event_TDiff_expand",100000,0,1000000000,janus.Timestamp()-prev_ts);
    obj.FillHistogram(dirname,"Janus_Event_TDiff_expand_2",100000,0,10000000000,janus.Timestamp()-prev_ts);
  }
  prev_ts = janus.Timestamp();
  
  for(size_t i=0; i<janus.GetAllChannels().size(); i++) {
    auto &chan = janus.GetJanusChannel(i);

    obj.FillHistogram(dirname,"Charge_Summary_AllChans",
		      128,0,128,chan.GetFrontChannel(),
		      3500,0,35000,chan.Charge());

    if(chan.GetDetnum() == 1 && chan.IsRing()) {
      obj.FillHistogram(dirname,"Charge_Summary_AllD1RingChans",
			26,0,26,D1_RingMap[chan.GetFrontChannel()],
			3500,0,35000,chan.Charge());
    }  
  }//end all channel loop

  for(size_t i=0; i<janus.GetChannels().size(); i++) {
    auto &chan = janus.GetChannels().at(i);

    if(prev_chan_ts > 0) {
      obj.FillHistogram(dirname,"Janus_Chan_TDiff_zoom",10000,0,10000000,chan.Timestamp()-prev_chan_ts);
      obj.FillHistogram(dirname,"Janus_Chan_TDiff",100000,0,100000000,chan.Timestamp()-prev_chan_ts);
      obj.FillHistogram(dirname,"Janus_Chan_TDiff_expand",100000,0,1000000000,chan.Timestamp()-prev_chan_ts);
      obj.FillHistogram(dirname,"Janus_Chan_TDiff_expand_2",100000,0,10000000000,chan.Timestamp()-prev_chan_ts);
    }
    prev_chan_ts = chan.Timestamp();

    obj.FillHistogram(dirname,"Charge_Summary_UsableChans",
		      128,0,128,chan.GetFrontChannel(),
		      3500,0,35000,chan.Charge());

    if(chan.GetDetnum() == 1 && chan.IsRing()) {
      obj.FillHistogram(dirname,"Charge_Summary_UsableD1RingChans",
			26,0,26,D1_RingMap[chan.GetFrontChannel()],
			3500,0,35000,chan.Charge());
    }

    for(size_t j=i+1;j<janus.GetChannels().size();j++) {
      auto &chan1 = janus.GetChannels().at(j);

      if(chan1.GetDetnum() == chan.GetDetnum()) {
        obj.FillHistogram(dirname,Form("Tdiff_UsableChans_det%d",chan.GetDetnum()),
	  		  2750,-11000,11000,chan.Timestamp()-chan1.Timestamp());
      }
      
    }//end second janus chan loop 
  }//end janus chan loop

  for(size_t i=0;i<janus.GetSectors().size();i++) {
    for(size_t j=0;j<janus.GetRings().size();j++) {
      
      auto &sect = janus.GetSectors().at(i);
      auto &ring = janus.GetRings().at(j);

      if(ring.GetDetnum() == sect.GetDetnum()) {
        obj.FillHistogram(dirname,Form("SectorCharge_v_RingCharge_UsableChans_det%d",sect.GetDetnum()),
		          2000,0,30000,ring.Charge(),
		          2000,0,30000,sect.Charge());

	auto tdiff = ring.Timestamp()-sect.Timestamp();

	obj.FillHistogram(dirname,Form("Sector_Ring_TDiff_UsableChans_det%d",sect.GetDetnum()),
			  2750,-11000,11000,tdiff);

	obj.FillHistogram(dirname,Form("SectorCharge_v_TDiff_UsableChans_det%d",sect.GetDetnum()),
		          2750,-11000,11000,tdiff,
		          2000,0,30000,sect.Charge());
      }
      
    }//end ring loop
  }//end sector loop
  
  //LeveledJanusSpectra(obj,janus,0);
  LeveledJanusSpectra(obj,janus,1);
  //LeveledJanusSpectra(obj,janus,2);
  //LeveledJanusSpectra(obj,janus,3); 
  
}

void RingRuthNorm(TRuntimeObjects& obj, TJanusDDASHit& hit, std::string dirname, std::string h_suf, int part, bool max_ring,
		  bool sol2) {

  static TSRIM b_srim(system_map[GValue::Value("REAC_KEY")].at(2).c_str());

  // Factors of 1e3 are because TNucleus and TReaction use MeV, while TSRIM uses keV.
  // (mg/cm^2) / (mg/cm^3) * (10^4 um/cm) distance energy loss calc
  static TReaction reac = TReaction(std::make_shared<TNucleus>(system_map[GValue::Value("REAC_KEY")].at(0).c_str()),
				    std::make_shared<TNucleus>(system_map[GValue::Value("REAC_KEY")].at(1).c_str()),
				    std::make_shared<TNucleus>(system_map[GValue::Value("REAC_KEY")].at(0).c_str()),
				    std::make_shared<TNucleus>(system_map[GValue::Value("REAC_KEY")].at(1).c_str()),
				    b_srim.GetAdjustedEnergy(GValue::Value("BEAM_ENERGY")*1000,
				    (GValue::Value("TARGET_THICKNESS")/GValue::Value("TARGET_DENSITY"))*10000.0/2.0)/1000.0
				    ); 

  double theta = hit.GetPosition(BEFORE).Theta();
  if(sol2)
    {theta*=-1;}

  double theta_cm = reac.ConvertThetaLabToCm(theta,part);
  double dRuth = 1.0/TMath::Power(TMath::Sin(theta_cm/2.0),4);
  //double solid_angle = hit.GetCmSolidAngle(reac,2,BEFORE,max_ring); //sr
  double solid_angle = TJanusDDAS::CmSolidAngle(hit.GetDetnum(),hit.GetRing(),reac,part,BEFORE,max_ring,sol2);

  obj.FillHistogram(dirname,Form("%sRingNorm_det%d",h_suf.c_str(),hit.GetDetnum()),
		    26,0,26,hit.GetRing(),1.0/(dRuth*solid_angle));
}

void MakeGatedJanus(TRuntimeObjects& obj, TJanusDDAS& janus, GCutG* gate) {

  if(!gate)
    {return;}
  
  std::string d_pref = "";
  std::string dirname = "";
 
  for(int lvl=min_level;lvl<max_level+1;lvl++) {

    if(lvl==0)
      {d_pref="Total_";}
    else if(lvl==1)
      {d_pref="Primary_";}
    else if(lvl==2)
      {d_pref="Secondary_";}
    else if(lvl==3)
      {d_pref="Tertiary_";}
    
    dirname = d_pref + Form("Janus_%s",gate->GetName());

    int detNum = 1;
    int part = 0;
    std::string gate_name = gate->GetName();
    if(!gate_name.compare("CdUS")) {
      detNum=0;
      part=2;
    }
    
    if(!gate_name.compare("Cd") || !gate_name.compare("CdS")) {
      part = 2;
    }
    if(!gate_name.compare("TiS") || !gate_name.compare("TiSBig") || !gate_name.compare("WS") || !gate_name.compare("Pb")
       || !gate_name.compare("PbBig")) {
      part = 3;
    }

    bool max_ring = false;
    if(!gate_name.compare("CdS")) {
      max_ring = true;
    }

    /*
    bool norm = true;
    if(!gate_name.compare("TiR") || !gate_name.compare("CdR") || !gate_name.compare("WR")) {
      norm = false;
    }
    */
    bool norm = false;

    int j_count = 0;
    for(size_t i=0;i<janus.GetHits(lvl).size();i++) {
      TJanusDDASHit& hit = janus.GetHits(lvl).at(i);
      
      if(hit.GetDetnum() != detNum || !gate->IsInside(hit.GetRing(),hit.Charge())) 
	{continue;}
      j_count++;

      bool sol2 = false;
      if(hit.GetRing() < 3 && hit.Charge() < 3000 && !gate_name.compare("CdS")) {
        sol2 = true;
      }

      if(norm) {
        RingRuthNorm(obj,hit,dirname,"All_",part,max_ring,sol2);
      }

      BasicJanusSpectra(obj,hit,dirname,"All_");
      SegmentTaggedJanusSpectra(obj,hit,d_pref,gate->GetName(),"All_");

      for(size_t j=i+1;j<janus.GetHits(lvl).size();j++) {
	TJanusDDASHit& hit1 = janus.GetHits(lvl).at(j);
        if(hit1.GetDetnum() != detNum || !gate->IsInside(hit1.GetRing(),hit1.Charge())) 
	  {continue;}

	obj.FillHistogram(dirname,Form("All_MultiHit_TDiff_det%d",hit.GetDetnum()),
			  2750,-11000,11000,hit.Timestamp() - hit1.Timestamp());

	//Wrong!
        //obj.FillHistogram(dirname,Form("All_MultiHit_TDiff_det%d_Mult%d",hit.GetDetnum(),j_count),
	//		  2750,-11000,11000,hit.Timestamp() - hit1.Timestamp());
	
      } //end second hit loop	
    } //end hit loop
    obj.FillHistogram(dirname,Form("All_GateMult_det%d",detNum),100,0,100,j_count);

    j_count = 0;
    for(size_t i=0;i<janus.GetSingleHits(lvl).size();i++) {
      TJanusDDASHit& hit = janus.GetSingleHits(lvl).at(i);
      
      if(hit.GetDetnum() != detNum || !gate->IsInside(hit.GetRing(),hit.Charge())) 
	{continue;}
      j_count++;

      bool sol2 = false;
      if(hit.GetRing() < 3 && hit.Charge() < 3000 && !gate_name.compare("CdS")) {
        sol2 = true;
      }

      if(norm) {
        RingRuthNorm(obj,hit,dirname,"Singles_",part,max_ring,sol2);
	RingRuthNorm(obj,hit,dirname,"SinglesAndTwoHit_",part,max_ring,sol2);
      }

      //BasicJanusSpectra(obj,hit,dirname,"Singles_");
      //SegmentTaggedJanusSpectra(obj,hit,d_pref,gate->GetName(),"Singles_");

      //BasicJanusSpectra(obj,hit,dirname,"SinglesAndTwoHit_");

      for(size_t j=i+1;j<janus.GetSingleHits(lvl).size();j++) {
	TJanusDDASHit& hit1 = janus.GetSingleHits(lvl).at(j);
        if(hit1.GetDetnum() != detNum || !gate->IsInside(hit1.GetRing(),hit1.Charge())) 
	  {continue;}

	obj.FillHistogram(dirname,Form("Singles_MultiHit_TDiff_det%d",hit.GetDetnum()),
			  2750,-11000,11000,hit.Timestamp()-hit1.Timestamp());

	//Still wrong
        //obj.FillHistogram(dirname,Form("Singles_MultiHit_TDiff_det%d_Mult%d",hit.GetDetnum(),j_count),
	//		  2750,-11000,11000,hit.Timestamp()-hit1.Timestamp());
	
      } //end second singles hit loop     
    } //end singles hit loop
    obj.FillHistogram(dirname,Form("Singles_GateMult_det%d",detNum),100,0,100,j_count);

    j_count=0;
    for(size_t i=0;i<janus.GetDoubleHits(lvl).size();i++) {
      TJanusDDASHit& hit = janus.GetDoubleHits(lvl).at(i);
      
      if(hit.GetDetnum() != detNum || !gate->IsInside(hit.GetRing(),hit.Charge())) 
	{continue;}
      j_count++;

      bool sol2 = false;
      if(hit.GetRing() < 3 && hit.Charge() < 3000 && !gate_name.compare("CdS")) {
        sol2 = true;
      }

      //BasicJanusSpectra(obj,hit,dirname,"TwoHit_");

      if(norm) {
        RingRuthNorm(obj,hit,dirname,"SinglesAndTwoHit_",part,max_ring,sol2);
      }
      
      //BasicJanusSpectra(obj,hit,dirname,"SinglesAndTwoHit_");

      for(size_t j=i+1;j<janus.GetDoubleHits(lvl).size();j++) {
	TJanusDDASHit& hit1 = janus.GetDoubleHits(lvl).at(j);
        if(hit1.GetDetnum() != detNum || !gate->IsInside(hit1.GetRing(),hit1.Charge())) 
	  {continue;}

	obj.FillHistogram(dirname,Form("TwoHit_MultiHit_TDiff_det%d",hit.GetDetnum()),
			  2750,-11000,11000,hit.Timestamp()-hit1.Timestamp());

	//Still wrong
        //obj.FillHistogram(dirname,Form("TwoHit_MultiHit_TDiff_det%d_Mult%d",hit.GetDetnum(),j_count),
	//		  2750,-11000,11000,hit.Timestamp()-hit1.Timestamp());
      } //end second doubles hit loop    
    } //end doubles hit loop
    obj.FillHistogram(dirname,Form("TwoHit_GateMult_det%d",detNum),100,0,100,j_count);

    j_count=0;
    for(auto& hit : janus.GetAddbackHits(lvl)) {
      if(hit.GetDetnum() != detNum || !gate->IsInside(hit.GetRing(),hit.Charge())) 
	{continue;}
      j_count++;

      bool sol2 = false;
      if(hit.GetRing() < 3 && hit.Charge() < 3000) {
        sol2 = true;
      }

      if(norm) {
        RingRuthNorm(obj,hit,dirname,"Addback_",part,max_ring,sol2);
      }
      
      SegmentTaggedJanusSpectra(obj,hit,d_pref,gate->GetName(),"Addback_");
      BasicJanusSpectra(obj,hit,dirname,"Addback_");
      
    } //end addback hit loop
    obj.FillHistogram(dirname,Form("Addback_GateMult_det%d",detNum),100,0,100,j_count);
    
  } //end level loop
}

void MakeDoubleGatedJanus(TRuntimeObjects& obj, TJanusDDAS& janus, GCutG* gate1, GCutG* gate2) {

  if(!gate1 || !gate2)
    {return;}
  
  for(int lvl=min_level;lvl<max_level+1;lvl++) {

    if(janus.GetHits(lvl).size() != 2) {
      continue;
    }

    int detNum = 1;
    TJanusDDASHit& hit = janus.GetHits(lvl).at(0);
    TJanusDDASHit& hit1 = janus.GetHits(lvl).at(1);
    
    if(hit.GetDetnum() != detNum ||
       (!gate1->IsInside(hit.GetRing(),hit.Charge()) && !gate2->IsInside(hit.GetRing(),hit.Charge()))) 
      {return;}

    if(hit1.GetDetnum() != detNum ||
       (!gate1->IsInside(hit1.GetRing(),hit1.Charge()) && !gate2->IsInside(hit1.GetRing(),hit1.Charge()))) 
      {return;}

    if(gate1->IsInside(hit.GetRing(),hit.Charge()) && gate1->IsInside(hit1.GetRing(),hit1.Charge()))
      {return;}

    if(gate2->IsInside(hit.GetRing(),hit.Charge()) && gate2->IsInside(hit1.GetRing(),hit1.Charge()))
      {return;}

    std::string d_pref = "";
    if(lvl==0)
      {d_pref="Total_";}
    else if(lvl==1)
      {d_pref="Primary_";}
    else if(lvl==2)
      {d_pref="Secondary_";}
    else if(lvl==3)
      {d_pref="Tertiary_";}

    std::string dirname = d_pref + Form("Janus_%sand%s",gate1->GetName(),gate2->GetName());
    std::string dsuf = Form("%sand%s",gate1->GetName(),gate2->GetName());
    BasicJanusSpectra(obj,hit,dirname,"All_");
    SegmentTaggedJanusSpectra(obj,hit,d_pref,dsuf,"All_");
    BasicJanusSpectra(obj,hit1,dirname,"All_");
    SegmentTaggedJanusSpectra(obj,hit1,d_pref,dsuf,"All_");
    
    bool hit_g1 = false;
    if(gate1->IsInside(hit.GetRing(),hit.Charge())) {
      hit_g1 = true;
    }

    std::string dirname1 = d_pref + Form("Janus_%sand%s_%s",gate1->GetName(),gate2->GetName(),
					 gate1->GetName());
    std::string dirname2 = d_pref + Form("Janus_%sand%s_%s",gate1->GetName(),gate2->GetName(),
					 gate2->GetName());
    std::string dsuf1 = Form("%sand%s_%s",gate1->GetName(),gate2->GetName(),gate1->GetName());
    std::string dsuf2 = Form("%sand%s_%s",gate1->GetName(),gate2->GetName(),gate2->GetName());
    
    if(hit_g1) {
      BasicJanusSpectra(obj,hit,dirname1,"All_");
      SegmentTaggedJanusSpectra(obj,hit,d_pref,dsuf1,"All_");

      BasicJanusSpectra(obj,hit1,dirname2,"All_");
      SegmentTaggedJanusSpectra(obj,hit1,d_pref,dsuf2,"All_");
    }
    else {
      BasicJanusSpectra(obj,hit,dirname2,"All_");
      SegmentTaggedJanusSpectra(obj,hit,d_pref,dsuf2,"All_");

      BasicJanusSpectra(obj,hit1,dirname1,"All_");
      SegmentTaggedJanusSpectra(obj,hit1,d_pref,dsuf1,"All_");
    }
    
  } //end level loop
  
  return;
}

void MakeSeGAJanus(TRuntimeObjects& obj, TSega& sega, TJanusDDAS& janus) {

   std::string dirname = "SeGAJanus";
   int lvl = 1;

   int d1_mult = janus.GetDetectorHits(1,lvl).size();
   int d0_mult = janus.GetDetectorHits(0,lvl).size();
   
   obj.FillHistogram(dirname,"sega_mult",20,0,20,sega.Size());

   obj.FillHistogram(dirname,"All_jMult_det1",40,0,40,d1_mult);

   obj.FillHistogram(dirname,"All_jMult_det0",40,0,40,d0_mult);

   for(auto& s_hit : sega) {

     double s_energy = s_hit.GetEnergy();

     obj.FillHistogram(dirname,"GammaEnergy",4000,0,4000,s_energy);
     obj.FillHistogram(dirname,"GammaEnergy_summary",20,0,20,s_hit.GetDetnum(),4000,0,4000,s_energy);

     if(d1_mult && !d0_mult) {
       obj.FillHistogram(dirname,"GammaEnergy_jDet01Only",4000,0,4000,s_energy);
       obj.FillHistogram(dirname,"GammaEnergy_summary_jDet01Only"
			 ,20,0,20,s_hit.GetDetnum(),4000,0,4000,s_energy);
     }
     
   }

   for(auto& j_hit : janus.GetHits(lvl)) {
     for(auto& s_hit : sega) {

       double s_energy = s_hit.GetEnergy();

       auto tdiff_sec = (s_hit.Timestamp() - j_hit.Timestamp());

       obj.FillHistogram(dirname,Form("GammaEnergy_v_tdiff_det%d",j_hit.GetDetnum()),
			 1000,-1000,3000,tdiff_sec,
			 2500,0,5000,s_energy);

       if(j_hit.GetDetnum() == 1) {
         obj.FillHistogram(dirname,Form("GammaEnergy_v_tdiff_det1_jMult%d",d1_mult),
			   1000,-1000,3000,tdiff_sec,
			   2500,0,5000,s_energy);
       }

       obj.FillHistogram(dirname,Form("ParticleCharge_v_GammaEnergy_jdet%d",j_hit.GetDetnum()),
			 2500,0,5000,s_energy,
			 3000,0,30000,j_hit.Charge());
       
       obj.FillHistogram(dirname,"P_Gam_Tdiff_Hit",2750,-11000,11000,tdiff_sec);
       
     }//end sega hit loop
   }//end janus hit loop

   
}

void MakeScatteredSeGAJanus(TRuntimeObjects& obj, TSega& sega, TJanusDDAS& janus, GCutG* time_gate, GCutG* scatter_gate) {
  
  if(!scatter_gate || !time_gate) 
    {return;}
  std::string gate_name = scatter_gate->GetName();
  
  static TSRIM b_srim(system_map[GValue::Value("REAC_KEY")].at(2).c_str());
  static TSRIM t_srim(system_map[GValue::Value("REAC_KEY")].at(3).c_str());
  
  // Factors of 1e3 are because TNucleus and TReaction use MeV, while TSRIM uses keV.
  // (mg/cm^2) / (mg/cm^3) * (10^4 um/cm) distance energy loss calc
  static TReaction reac = TReaction(std::make_shared<TNucleus>(system_map[GValue::Value("REAC_KEY")].at(0).c_str()),
				    std::make_shared<TNucleus>(system_map[GValue::Value("REAC_KEY")].at(1).c_str()),
				    std::make_shared<TNucleus>(system_map[GValue::Value("REAC_KEY")].at(0).c_str()),
				    std::make_shared<TNucleus>(system_map[GValue::Value("REAC_KEY")].at(1).c_str()),
				    b_srim.GetAdjustedEnergy(GValue::Value("BEAM_ENERGY")*1000,
				    (GValue::Value("TARGET_THICKNESS")/GValue::Value("TARGET_DENSITY"))*10000.0/2.0)/1000.0
				   );
  int detNum = -1;
  int part = -1;
  
  if(!gate_name.compare("Cd") || !gate_name.compare("CdS") || !gate_name.compare("CdS_FW")) {
    detNum = 1;
    part = 2;
  }
  
  else if(!gate_name.compare("TiS") || !gate_name.compare("TiSBig") || !gate_name.compare("WS") || !gate_name.compare("UnK")) {
    detNum = 1;
    part = 3;
  }
  
  else if(!gate_name.compare("CdUS")) {
    detNum = 0;
    part = 2;
  }

  std::string d_pref = "";
  std::string dirname = "";
  std::string jCh_dirname;

  for(int lvl=min_level;lvl<max_level+1;lvl++) {

    if(lvl==0)
      {d_pref="Total_";}
    else if(lvl==1)
      {d_pref="Primary_";}
    else if(lvl==2)
      {d_pref="Secondary_";}
    else if(lvl==3)
      {d_pref="Tertiary_";}
    
    dirname = d_pref + "SeGAJanus_Scattered" + scatter_gate->GetName();
    std::string j_dirname = dirname;
    jCh_dirname = d_pref + "SeGAJanus_jChargeGate_Scattered" + scatter_gate->GetName();
    if(time_gate) {
      dirname += Form("_%s",time_gate->GetName());
      jCh_dirname += Form("_%s",time_gate->GetName());
    }
    
    std::vector<int> times_used;
    times_used.resize(sega.Size());
    std::fill(times_used.begin(),times_used.end(),0);

    std::vector<TJanusDDASHit> good_jHits;

    std::string h_suf = "All_";
    /*
    if(time_gate) {

      std::vector<TJanusDDASHit> very_good_jHits;
      
      for(unsigned int i=0;i<janus.GetHits(lvl).size();i++) {
	auto& j_hit = janus.GetHits(lvl).at(i);
      
	if(j_hit.GetDetnum() != detNum || !scatter_gate->IsInside(j_hit.GetRing(),j_hit.Charge())) 
	  {continue;}
      
	for(unsigned int j=0;j<sega.Size();j++) {
	  auto& s_hit = sega.GetSegaHit(j);

	  if(times_used.at(j))
	    {continue;}

	  auto tdiff = (s_hit.Timestamp() - j_hit.Timestamp());
	  double s_energy = s_hit.GetEnergy();
       
	  if(!time_gate->IsInside(tdiff,s_energy) || s_energy < 600) 
	    {continue;}
	  
	  times_used.at(j)++;
	  very_good_jHits.push_back(j_hit);
	  
	  break;

	}
      }

      std::string tmp_dirname = dirname + "_sEnGate";
      for(auto& j_hit : very_good_jHits) {
	BasicJanusSpectra(obj,j_hit,tmp_dirname,h_suf);
	obj.FillHistogram(Form("%sSeGAJanus_Scattered%s_jRings_TE_sEnGate",d_pref.c_str(),
			       scatter_gate->GetName()),
			  Form("%sSector_Hit_det%d_Ring%02d",h_suf.c_str(),j_hit.GetDetnum(),
			       j_hit.GetRing()),
			  34,0,34,j_hit.GetSector());
      }
      
    }
    */
    times_used.clear();
    times_used.resize(sega.Size());
    std::fill(times_used.begin(),times_used.end(),0);

    int j_count = 0;
    for(auto& j_hit : janus.GetHits(lvl)) {
      if(j_hit.GetDetnum() != detNum || !scatter_gate->IsInside(j_hit.GetRing(),j_hit.Charge())) 
	{continue;}
      good_jHits.push_back(j_hit);
      j_count++;

      BasicJanusSpectra(obj,j_hit,j_dirname,h_suf);
      obj.FillHistogram(Form("%sSeGAJanus_Scattered%s_jRings",d_pref.c_str(),scatter_gate->GetName()),
			Form("%sSector_Hit_det%d_Ring%02d",h_suf.c_str(),j_hit.GetDetnum(),j_hit.GetRing()),
			34,0,34,j_hit.GetSector());

      //int s_count = 0;
      //int s_count_total
      for(unsigned int i=0;i<sega.Size();i++) {
	auto& s_hit = sega.GetSegaHit(i); 

	//s_count_total++;
	double s_energy = s_hit.GetEnergy();
	auto tdiff = (s_hit.Timestamp() - j_hit.Timestamp());

	obj.FillHistogram(dirname,"GammaEnergy_v_tdiff_NoTimeGate",
			  2000,-10000,10000,tdiff,
			  2500,0,5000,s_energy);

	if(time_gate) {
	  if(!time_gate->IsInside(tdiff,s_energy)) 
	    {continue;}
	}
	if(times_used.at(i))
	  {continue;}
	times_used.at(i)++;
	//s_count++;

	obj.FillHistogram(dirname,"GammaEnergy_v_tdiff",
			  1000,-1000,3000,tdiff,
			  2500,0,5000,s_energy);

	double num_segs = s_hit.GetNumSegments();

	if(num_segs == 0) {
	  obj.FillHistogram(dirname,"GammaEnergy_NoSegs",4000,0,4000,s_energy);
	}
	else if(num_segs > 0) {
	  obj.FillHistogram(dirname,"GammaEnergy_YesSegs",4000,0,4000,s_energy);
	}

	auto j_pos = j_hit.GetPosition(BEFORE);
	
	//This will fix the doppler correction of the JANUS events with switched sectors
	//Histograms with the suffix jChargeGate show this nicely
	auto j_pos_rot = j_hit.GetPosition(BEFORE);
	j_pos_rot.SetPhi(j_pos_rot.Phi() + TMath::Pi());
	
	auto s_pos = s_hit.GetPosition();
	auto angle = s_pos.Angle(j_pos)*TMath::RadToDeg();

	bool sol2 = false;
	if(!gate_name.compare("CdS") && j_hit.GetRing() < 3 && j_hit.Charge() < 3000)
	  {sol2 = true;}

	double beta;
	//double beta_noS2;
	if(part == 2) {
	  beta = TJanusDDAS::GetBeta(j_pos.Theta(),reac,part,b_srim,sol2);
          //beta_noS2 = TJanusDDAS::GetBeta(j_pos.Theta(),reac,part,b_srim,0);
	}
	else {
          beta = TJanusDDAS::GetBeta(j_pos.Theta(),reac,part,t_srim,sol2);
          //beta_noS2 = TJanusDDAS::GetBeta(j_pos.Theta(),reac,part,t_srim,0);
	}	  
	
	obj.FillHistogram(dirname,"Beta_v_Ring",26,0,26,j_hit.GetRing(),100,0,0.1,beta);
	obj.FillHistogram(dirname,"Angle_v_Ring",26,0,26,j_hit.GetRing(),90,0,180,angle);
	  
	double s_energy_dc_z = s_hit.GetDoppler(beta,j_pos);
	double s_energy_dc_z_rot = s_hit.GetDoppler(beta,j_pos_rot);
	//double s_en_dc_z_noS2 = s_hit.GetDoppler(beta_noS2,j_pos);

	
	obj.FillHistogram(dirname,"GammaEnergyDCz_v_tdiff",
			  1000,-1000,3000,tdiff,
			  2500,0,5000,s_energy_dc_z);

	obj.FillHistogram(dirname,"GammaEnergy",4000,0,4000,s_energy);
	obj.FillHistogram(dirname,"GammaEnergy_summary",20,0,20,s_hit.GetDetnum(),4000,0,4000,s_energy);
	  
	obj.FillHistogram(dirname,"GammaEnergyDCz",12000,0,4000,s_energy_dc_z);
	obj.FillHistogram(dirname,"GammaEnergyDCz_summary",20,0,20,s_hit.GetDetnum(),8000,0,4000,s_energy_dc_z);

	for(int i=1;i<101;i++) {
	  double tmp_zoff = -5.0 + 0.1*i;
	  obj.FillHistogram(dirname,"sZ_Scan",100,-5.0,5.1,tmp_zoff,
			    2000,0,4000,s_hit.GetDoppler(beta,j_pos,TVector3(0,0,tmp_zoff)));
	}

		
	obj.FillHistogram(Form("%sScattered%s_jRings_%s",d_pref.c_str(),scatter_gate->GetName(),time_gate->GetName()),
			  Form("GammaEnergy_R%02d",j_hit.GetRing()),2000,0,4000,s_energy);
	  
	obj.FillHistogram(Form("%sScattered%s_jRings_%s",d_pref.c_str(),scatter_gate->GetName(),time_gate->GetName()),
			  Form("GammaEnergyDCz_R%02d",j_hit.GetRing()),2000,0,4000,s_energy_dc_z);
	
	
	obj.FillHistogram(dirname,Form("ParticleCharge_v_GammaEnergy_jdet%d",j_hit.GetDetnum()),
			  2500,0,2500,s_energy,
			  1000,0,30000,j_hit.Charge());

 	//obj.FillHistogram(dirname,Form("ParticleCharge_v_GammaEnergyDCz_noS2_jdet%d",j_hit.GetDetnum()),
 	//		  2500,0,2500,s_en_dc_z_noS2,
	//		  3000,0,30000,j_hit.Charge());
	
 	obj.FillHistogram(dirname,Form("ParticleCharge_v_GammaEnergyDCz_jdet%d",j_hit.GetDetnum()),
			  2500,0,2500,s_energy_dc_z,
			  3000,0,30000,j_hit.Charge());

	/*
	if(!gate_name.compare("Cd")) {
	  for(double b : BinCenters(120,0.06,0.12)) {
	    obj.FillHistogram(Form("%sScattered%s_jRings_BetaScan_%s",d_pref.c_str(),
				   scatter_gate->GetName(),time_gate->GetName()),
			      Form("GammaEnergyDCz_BetaScan_R%02d",j_hit.GetRing()),
			      120,0.06,0.12,b,
			      2000,0,4000,s_hit.GetDoppler(b,j_pos,TVector3(0,0,GValue::Value("Sega_Z_offset"))));
	  }
	}
	else if(!gate_name.compare("CdUS")) {
	  for(double b : BinCenters(80,0.01,0.05)) {
	    obj.FillHistogram(Form("%sScattered%s_jRings_BetaScan_%s",d_pref.c_str(),scatter_gate->GetName(),time_gate->GetName()),
			      Form("GammaEnergyDCz_BetaScan_R%02d",j_hit.GetRing()),
			      80,0.01,0.05,b,
			      2000,0,4000,s_hit.GetDoppler(b,j_pos,TVector3(0,0,GValue::Value("Sega_Z_offset"))));
	  }
	}

	obj.FillHistogram(Form("%sScattered%s_jRings_%s",d_pref.c_str(),scatter_gate->GetName(),time_gate->GetName()),
 			  Form("ParticleCharge_v_GammaEnergy_R%02d",j_hit.GetRing()),
 			  1000,0,2000,s_energy,
       	         	  3000,0,30000,j_hit.Charge());

 	obj.FillHistogram(Form("%sScattered%s_jRings_%s",d_pref.c_str(),scatter_gate->GetName(),time_gate->GetName()),
 			  Form("ParticleCharge_v_GammaEnergyDCz_R%02d",j_hit.GetRing()),
			  1000,0,2000,s_energy_dc_z,
			  3000,0,30000,j_hit.Charge());
	*/

	obj.FillHistogram(dirname,"GammaEnergy_v_Angle",4000,0,4000,s_energy,90,0,180,angle);
	obj.FillHistogram(dirname,"GammaEnergyDCz_v_Angle",8000,0,4000,s_energy_dc_z,180,0,180,angle);

	obj.FillHistogram(dirname,"GammaEnergy_v_CosAngle",4000,0,4000,s_energy,
			  200,-1.1,1.1,TMath::Cos(angle*TMath::DegToRad()));
	obj.FillHistogram(dirname,"GammaEnergyDCz_v_CosAngle",8000,0,4000,s_energy_dc_z,
			  200,-1.1,1.1,TMath::Cos(angle*TMath::DegToRad()));

	/*
	obj.FillHistogram(Form("%sScattered%s_jRings_%s_ThetaCorr",d_pref.c_str(),scatter_gate->GetName(),time_gate->GetName()),
			  Form("GammaEnergy_v_Angle_R%02d",j_hit.GetRing()),
			  4000,0,4000,s_energy,
			  180,0,180,angle);
	
	obj.FillHistogram(Form("%sScattered%s_jRings_%s_ThetaCorr",d_pref.c_str(),scatter_gate->GetName(),time_gate->GetName()),
			  Form("GammaEnergyDCz_v_Angle_R%02d",j_hit.GetRing()),
			  4000,0,4000,s_energy_dc_z,
			  180,0,180,angle);
	*/
	
	
	TVector3 incBeam = TVector3(0.0,0.0,1.0);
	TVector3 reacPlane = j_pos.Cross(incBeam);
	TVector3 reacPlane_rot = j_pos_rot.Cross(incBeam);
	TVector3 detPlane = s_pos.Cross(incBeam);

	double reac_phi = reacPlane.Phi();
	if(reac_phi < 0)
	  {reac_phi += TMath::TwoPi();}

	double reac_phi_rot = reacPlane_rot.Phi();
	if(reac_phi_rot < 0)
	  {reac_phi_rot += TMath::TwoPi();}

	double det_phi = detPlane.Phi();
	if(det_phi < 0)
	  {det_phi += TMath::TwoPi();}

	double planeAng = reac_phi - det_phi;
	if(planeAng < 0)
	  {planeAng += TMath::TwoPi();}

	double planeAng_rot = reac_phi_rot - det_phi;
	if(planeAng_rot < 0)
	  {planeAng_rot += TMath::TwoPi();}

	//double thing1 = 33*TMath::Pi()/16.0;
	double thing1 = 65*TMath::Pi()/32.0;

	obj.FillHistogram(dirname,"GammaEnergy_v_PlaneAng",4000,0,4000,s_energy,32,0,thing1,planeAng);
	obj.FillHistogram(dirname,"GammaEnergyDCz_v_PlaneAng",8000,0,4000,s_energy_dc_z,32,0,thing1,planeAng);

	obj.FillHistogram(dirname,"GammaEnergy_v_PlaneAngRot",4000,0,4000,s_energy,32,0,thing1,planeAng_rot);
	obj.FillHistogram(dirname,"GammaEnergyDCz_v_PlaneAngRot",8000,0,4000,s_energy_dc_z_rot,32,0,thing1,planeAng_rot);

	/*
	obj.FillHistogram(Form("%sScattered%s_jRings_%s_PhiCorr",d_pref.c_str(),scatter_gate->GetName(),time_gate->GetName()),
			  Form("GammaEnergy_v_PlaneAngle_R%02d",j_hit.GetRing()),
			  4000,0,4000,s_energy,
			  32,0,thing1,planeAng);
	
	obj.FillHistogram(Form("%sScattered%s_jRings_%s_PhiCorr",d_pref.c_str(),scatter_gate->GetName(),time_gate->GetName()),
			  Form("GammaEnergyDCz_v_PlaneAngle_R%02d",j_hit.GetRing()),
			  4000,0,4000,s_energy_dc_z,
			  32,0,thing1,planeAng);

	
	obj.FillHistogram(Form("%sScattered%s_jRings_%s_PhiCorrRot",d_pref.c_str(),scatter_gate->GetName(),time_gate->GetName()),
			  Form("GammaEnergy_v_PlaneAngleRot_R%02d",j_hit.GetRing()),
			  4000,0,4000,s_energy,
			  32,0,thing1,planeAng_rot);
	
	obj.FillHistogram(Form("%sScattered%s_jRings_%s_PhiCorrRot",d_pref.c_str(),scatter_gate->GetName(),time_gate->GetName()),
			  Form("GammaEnergyDCz_v_PlaneAngleRot_R%02d",j_hit.GetRing()),
			  4000,0,4000,s_energy_dc_z_rot,
			  32,0,thing1,planeAng_rot);
	*/
	
	if(j_hit.Charge() > 4650 && j_hit.Charge() < 5500 && !gate_name.compare("CdS")) {

	  /*
	  if(janus.GetDetectorHits(detNum,lvl).size() < 3) {
	  obj.FillHistogram(jCh_dirname,
			    Form("GammaEnergy_v_Angle_jMult%d",(int)janus.GetDetectorHits(detNum,lvl).size()),
			    2000,0,4000,s_energy,
			    90,0,180,angle);
	  obj.FillHistogram(jCh_dirname,
			    Form("GammaEnergyDCz_v_Angle_jMult%d",(int)janus.GetDetectorHits(detNum,lvl).size()),
			    2000,0,4000,s_energy_dc_z,
			    180,0,180,angle);

	  obj.FillHistogram(jCh_dirname,
			    Form("GammaEnergy_v_PlaneAng_jMult%d",(int)janus.GetDetectorHits(detNum,lvl).size()),
			    2000,0,4000,s_energy,
			    32,0,thing1,planeAng);
	  obj.FillHistogram(jCh_dirname,
			    Form("GammaEnergyDCz_v_PlaneAng_jMult%d",(int)janus.GetDetectorHits(detNum,lvl).size()),
			    2000,0,4000,s_energy_dc_z,
			    32,0,thing1,planeAng);
	  }
	  else if(janus.GetDetectorHits(detNum,lvl).size() > 2) {
	    obj.FillHistogram(jCh_dirname,"GammaEnergy_v_Angle_jMult3+",
			      2000,0,4000,s_energy,
			      90,0,180,angle);
	    obj.FillHistogram(jCh_dirname,"GammaEnergyDCz_v_Angle_jMult3+",
			      2000,0,4000,s_energy_dc_z,
			      180,0,180,angle);

	    obj.FillHistogram(jCh_dirname,"GammaEnergy_v_PlaneAng_jMult3+",
			      2000,0,4000,s_energy,
			      32,0,thing1,planeAng);
	    obj.FillHistogram(jCh_dirname,"GammaEnergyDCz_v_PlaneAng_jMult3+",
			      2000,0,4000,s_energy_dc_z,
              		      32,0,thing1,planeAng);
	  }
	  */

	  /*
	  obj.FillHistogram(jCh_dirname,"GammaEnergy_v_Angle",
			    4000,0,4000,s_energy,
			    90,0,180,angle);
	  obj.FillHistogram(jCh_dirname,"GammaEnergyDCz_v_Angle",
			    4000,0,4000,s_energy_dc_z,
			    180,0,180,angle);
	  */
	  
	  obj.FillHistogram(jCh_dirname,"GammaEnergy_v_PlaneAng",
			    4000,0,4000,s_energy,
			    32,0,thing1,planeAng);
	  obj.FillHistogram(jCh_dirname,"GammaEnergyDCz_v_PlaneAng",
			    4000,0,4000,s_energy_dc_z,
			    32,0,thing1,planeAng);

	  obj.FillHistogram(jCh_dirname,"GammaEnergy_v_PlaneAngRot",
			    4000,0,4000,s_energy,
			    32,0,thing1,planeAng_rot);
	  obj.FillHistogram(jCh_dirname,"GammaEnergyDCz_v_PlaneAngRot",
			    4000,0,4000,s_energy_dc_z_rot,
			    32,0,thing1,planeAng_rot);

	  //obj.FillHistogram(Form("%sScattered%s_jRings_%s",d_pref.c_str(),scatter_gate->GetName(),time_gate->GetName()),
	  //		    Form("GammaEnergyDCz_v_PlaneAng_jRing%02d",j_hit.GetRing()),
	  //		    4000,0,4000,s_energy_dc_z,
	  //		    32,0,thing1,planeAng);

	}
	else if(j_hit.Charge() > 4500 && j_hit.Charge() < 5300 && (!gate_name.compare("TiS") || !gate_name.compare("TiSBig"))) {

	  /*
	  if(janus.GetDetectorHits(detNum,lvl).size() < 3) {
	  obj.FillHistogram(jCh_dirname,
			    Form("GammaEnergy_v_Angle_jMult%d",(int)janus.GetDetectorHits(detNum,lvl).size()),
			    2000,0,4000,s_energy,
			    90,0,180,angle);
	  obj.FillHistogram(jCh_dirname,
			    Form("GammaEnergyDCz_v_Angle_jMult%d",(int)janus.GetDetectorHits(detNum,lvl).size()),
			    2000,0,4000,s_energy_dc_z,
			    180,0,180,angle);

	  obj.FillHistogram(jCh_dirname,
			    Form("GammaEnergy_v_PlaneAng_jMult%d",(int)janus.GetDetectorHits(detNum,lvl).size()),
			    2000,0,4000,s_energy,
			    32,0,thing1,planeAng);
	  obj.FillHistogram(jCh_dirname,
			    Form("GammaEnergyDCz_v_PlaneAng_jMult%d",(int)janus.GetDetectorHits(detNum,lvl).size()),
			    2000,0,4000,s_energy_dc_z,
			    32,0,thing1,planeAng);
	  }
	  else if(janus.GetDetectorHits(detNum,lvl).size() > 2) {
	    obj.FillHistogram(jCh_dirname,"GammaEnergy_v_Angle_jMult3+",
			      2000,0,4000,s_energy,
			      90,0,180,angle);
	    obj.FillHistogram(jCh_dirname,"GammaEnergyDCz_v_Angle_jMult3+",
			      2000,0,4000,s_energy_dc_z,
			      180,0,180,angle);

	    obj.FillHistogram(jCh_dirname,"GammaEnergy_v_PlaneAng_jMult3+",
			      2000,0,4000,s_energy,
			      32,0,thing1,planeAng);
	    obj.FillHistogram(jCh_dirname,"GammaEnergyDCz_v_PlaneAng_jMult3+",
			      2000,0,4000,s_energy_dc_z,
              		      32,0,thing1,planeAng);
	  }
	  */

	  /*
	  obj.FillHistogram(jCh_dirname,"GammaEnergy_v_Angle",
			    4000,0,4000,s_energy,
			    90,0,180,angle);
	  obj.FillHistogram(jCh_dirname,"GammaEnergyDCz_v_Angle",
			    4000,0,4000,s_energy_dc_z,
			    180,0,180,angle);
	  */
	  
	  obj.FillHistogram(jCh_dirname,"GammaEnergy_v_PlaneAng",
			    4000,0,4000,s_energy,
			    32,0,thing1,planeAng);
	  obj.FillHistogram(jCh_dirname,"GammaEnergyDCz_v_PlaneAng",
			    4000,0,4000,s_energy_dc_z,
			    32,0,thing1,planeAng);

	  obj.FillHistogram(jCh_dirname,"GammaEnergy_v_PlaneAngRot",
			    4000,0,4000,s_energy,
			    32,0,thing1,planeAng_rot);
	  obj.FillHistogram(jCh_dirname,"GammaEnergyDCz_v_PlaneAngRot",
			    4000,0,4000,s_energy_dc_z_rot,
			    32,0,thing1,planeAng_rot);

	  //obj.FillHistogram(Form("%sScattered%s_jRings_%s",d_pref.c_str(),scatter_gate->GetName(),time_gate->GetName()),
	  //		    Form("GammaEnergyDCz_v_PlaneAng_jRing%02d",j_hit.GetRing()),
	  //		    4000,0,4000,s_energy_dc_z,
	  //		    32,0,thing1,planeAng);
	}
	else if(!gate_name.compare("CdS") || !gate_name.compare("TiS") || !gate_name.compare("TiSBig")) {

	  /*
	  obj.FillHistogram(jCh_dirname,"GammaEnergy_v_Angle_NotjChargeGate",
			    4000,0,4000,s_energy,
			    90,0,180,angle);
	  obj.FillHistogram(jCh_dirname,"GammaEnergyDCz_v_Angle_NotjChargeGate",
			    4000,0,4000,s_energy_dc_z,
			    180,0,180,angle);
	  */
	  
	  obj.FillHistogram(jCh_dirname,"GammaEnergy_v_PlaneAng_NotjChargeGate",
			    4000,0,4000,s_energy,
			    32,0,thing1,planeAng);
	  obj.FillHistogram(jCh_dirname,"GammaEnergyDCz_v_PlaneAng_NotjChargeGate",
			    4000,0,4000,s_energy_dc_z,
			    32,0,thing1,planeAng);

	  obj.FillHistogram(jCh_dirname,"GammaEnergy_v_PlaneAngRot_NotjChargeGate",
			    4000,0,4000,s_energy,
			    32,0,thing1,planeAng_rot);
	  obj.FillHistogram(jCh_dirname,"GammaEnergyDCz_v_PlaneAngRot_NotjChargeGate",
			    4000,0,4000,s_energy_dc_z_rot,
			    32,0,thing1,planeAng_rot);
	}
	  
	/*
	  obj.FillHistogram(dirname,"GammaEnergyDCz_allseg_summary",
	  512,1,513,32*(s_hit.GetDetnum() - 1) + s_hit.GetMainSegnum(),
	  35,600,670,s_energy_dc_z);

	  obj.FillHistogram(dirname,"GammaEnergy_SegMap_summary",
	  512,1,513,s_hit.GetMapSegnum(),
	  70,570,710,s_energy);

	  obj.FillHistogram(dirname,"GammaEnergy_PairMap_summary",
	  256,1,257,s_hit.GetMapPairnum(),
	  70,570,710,s_energy);

	  obj.FillHistogram(dirname,"GammaEnergy_SliceMap_summary",
	  128,1,129,s_hit.GetMapSlicenum(),
	  70,570,710,s_energy);

	  obj.FillHistogram(dirname,"GammaEnergyDCz_SegMap_summary",
	  512,1,513,s_hit.GetMapSegnum(),
	  35,600,670,s_energy_dc_z);

	  obj.FillHistogram(dirname,"GammaEnergyDCz_PairMap_summary",
	  256,1,257,s_hit.GetMapPairnum(),
	  35,600,670,s_energy_dc_z);

	  obj.FillHistogram(dirname,"GammaEnergyDCz_SliceMap_summary",
	  128,1,129,s_hit.GetMapSlicenum(),
	  35,600,670,s_energy_dc_z);
	*/
      }//end sega hit loop

      //obj.FillHistogram(dirname,"SeGA_Mult",25,0,25,s_count);
      
    }//end janus hit loop

    obj.FillHistogram(dirname,"Janus_HitMult",100,0,100,janus.GetDetectorHits(detNum,lvl).size());
    obj.FillHistogram(dirname,"Janus_GateMult",100,0,100,j_count);

    /*
    for(size_t i=0;i<good_jHits.size();i++) {
      TJanusDDASHit& j_hit = good_jHits.at(i);
      for(auto& s_hit : sega) {

	double s_energy = s_hit.GetEnergy();
	auto tdiff = (s_hit.Timestamp() - j_hit.Timestamp());

	obj.FillHistogram(dirname,Form("GammaEnergy_v_tdiff_NoTimeGate_jMult%d",j_count),
			  1000,-1000,3000,tdiff,
			  2500,0,5000,s_energy);
      }

      for(size_t j=i+1;j<good_jHits.size();j++) {
        TJanusDDASHit& j_hit1 = good_jHits.at(j);

	obj.FillHistogram(dirname,Form("All_MultiHit_TDiff_det%d",j_hit.GetDetnum()),
			  2750,-11000,11000,j_hit.Timestamp()-j_hit1.Timestamp());

	obj.FillHistogram(dirname,Form("All_MultiHit_TDiff_det%d_Mult%d",j_hit.GetDetnum(),(int)good_jHits.size()),
			  2750,-11000,11000,j_hit.Timestamp()-j_hit1.Timestamp());
	
      } //end second good janus hit loop
    } //end good janus hit loop
    */
    
  }//end janus level loop
}

void GammaGamma(TRuntimeObjects& obj, TSega& sega, TJanusDDAS& janus, GCutG* time_gate,
		GCutG* scatter_gate) {

  if(!scatter_gate || !time_gate) 
    {return;}

  std::string gate_name = scatter_gate->GetName();
  
  static TSRIM b_srim(system_map[GValue::Value("REAC_KEY")].at(2).c_str());
  static TSRIM t_srim(system_map[GValue::Value("REAC_KEY")].at(3).c_str());
  
  // Factors of 1e3 are because TNucleus and TReaction use MeV, while TSRIM uses keV.
  // (mg/cm^2) / (mg/cm^3) * (10^4 um/cm) distance energy loss calc
  static TReaction reac = TReaction(std::make_shared<TNucleus>(system_map[GValue::Value("REAC_KEY")].at(0).c_str()),
				    std::make_shared<TNucleus>(system_map[GValue::Value("REAC_KEY")].at(1).c_str()),
				    std::make_shared<TNucleus>(system_map[GValue::Value("REAC_KEY")].at(0).c_str()),
				    std::make_shared<TNucleus>(system_map[GValue::Value("REAC_KEY")].at(1).c_str()),
				    b_srim.GetAdjustedEnergy(GValue::Value("BEAM_ENERGY")*1000,
				    (GValue::Value("TARGET_THICKNESS")/GValue::Value("TARGET_DENSITY"))*10000.0/2.0)/1000.0
				   );
  int detNum = 1;
  
  std::string dirname = "Primary_GammaGamma_" + gate_name;
  if(time_gate) {
    dirname += Form("_%s",time_gate->GetName());
  }
    
  std::vector<int> times_used;
  times_used.resize(sega.Size());
  std::fill(times_used.begin(),times_used.end(),0);
  
  for(auto& j_hit : janus.GetHits(1)) {
    
    if(j_hit.GetDetnum() != detNum || !scatter_gate->IsInside(j_hit.GetRing(),j_hit.Charge())) {
      continue;
    }

    auto j_pos = j_hit.GetPosition(false,false);
    double beta = TJanusDDAS::GetBeta(j_pos.Theta(),reac,2,b_srim,false);

    auto recon_pos = j_hit.GetReconPosition(reac,3,2,false,false);
    double recon_beta = TJanusDDAS::GetBeta(recon_pos.Theta(),reac,2,b_srim,false);
    
    for(unsigned int i=0;i<sega.Size();i++) {
      auto& s_hit1 = sega.GetSegaHit(i);
	
      if(time_gate) {
	if(!time_gate->IsInside(s_hit1.Timestamp() - j_hit.Timestamp(),s_hit1.GetEnergy())) {
	  continue;
	}
      }

      if(times_used.at(i)) {
	continue;
      }
      times_used.at(i)++;

      double cen1 = s_hit1.GetEnergy();
      double dop1 = s_hit1.GetDoppler(beta,j_pos);
      double rec1 = s_hit1.GetDoppler(recon_beta,recon_pos);

      obj.FillHistogram(dirname,"All_CoreEn",8000,0,8000,cen1);
      obj.FillHistogram(dirname,"All_Dop",8000,0,4000,dop1);
      obj.FillHistogram(dirname,"All_Recon",8000,0,4000,rec1);
      
      for(unsigned int j=0;j<sega.Size();j++) {

	if(i==j) {
	  continue;
	}
	
	auto& s_hit2 = sega.GetSegaHit(j);
	
	if(time_gate) {
	  if(!time_gate->IsInside(s_hit2.Timestamp() - j_hit.Timestamp(),s_hit2.GetEnergy())) {
	    continue;
	  }
	}

	double cen2 = s_hit2.GetEnergy();
	double dop2 = s_hit2.GetDoppler(beta,j_pos);
	double rec2 = s_hit2.GetDoppler(recon_beta,recon_pos);

	obj.FillHistogram(dirname,"All_CoreEnGamma",
			  1000,0,2000,dop1,1000,0,2000,cen2);
	
	obj.FillHistogram(dirname,"All_GammaGamma",
			  1000,0,2000,dop1,1000,0,2000,dop2);

	obj.FillHistogram(dirname,"All_ReconCoreEnGamma",
			  1000,0,2000,rec1,1000,0,2000,cen2);
	
	obj.FillHistogram(dirname,"All_ReconGammaGamma",
			  1000,0,2000,rec1,1000,0,2000,rec2);

      } //End second SeGA loop	
	
    } //End first SeGA loop
      
  } //End janus loop

  std::fill(times_used.begin(),times_used.end(),0);
  std::vector<TJanusDDASHit> j_hits;
  std::vector< std::vector<TSegaHit> > s_hits;
  
  for(auto& j_hit : janus.GetHits(1)) {
    
    if(j_hit.GetDetnum() != detNum || !scatter_gate->IsInside(j_hit.GetRing(),j_hit.Charge())) {
      continue;
    }

    j_hits.push_back(j_hit);
    s_hits.emplace_back();
    
    for(unsigned int i=0;i<sega.Size();i++) {
      auto& s_hit1 = sega.GetSegaHit(i);
	
      if(time_gate) {
	if(!time_gate->IsInside(s_hit1.Timestamp() - j_hit.Timestamp(),s_hit1.GetEnergy())) {
	  continue;
	}
      }

      if(times_used.at(i)) {
	continue;
      }
      times_used.at(i)++;

      s_hits.back().push_back(s_hit1);
      
    }

  }

  for(unsigned int i=0;i<j_hits.size();i++) {

    auto& j_hit = j_hits.at(i);
    auto j_pos = j_hit.GetPosition(false,false);
    double beta = TJanusDDAS::GetBeta(j_pos.Theta(),reac,2,b_srim,false);

    auto recon_pos = j_hit.GetReconPosition(reac,3,2,false,false);
    double recon_beta = TJanusDDAS::GetBeta(recon_pos.Theta(),reac,2,b_srim,false);
    
    int s_mult = s_hits.at(i).size();
    for(int j=0;j<s_mult;j++) {
      auto& s_hit1 = s_hits.at(i).at(j);

      int det1 = s_hit1.GetDetnum();
      double cen1 = s_hit1.GetEnergy();
      double dop1 = s_hit1.GetDoppler(beta,j_pos);
      double rec1 = s_hit1.GetDoppler(recon_beta,recon_pos);

      obj.FillHistogram(dirname,"All_CoreEn_Check",8000,0,8000,cen1);

      for(int k=0;k<s_mult;k++) {

	if(j==k) {
	  continue;
	}
	auto& s_hit2 = s_hits.at(i).at(k);

	int det2 = s_hit2.GetDetnum();
	double cen2 = s_hit2.GetEnergy();
	double dop2 = s_hit2.GetDoppler(beta,j_pos);
	double rec2 = s_hit2.GetDoppler(recon_beta,recon_pos);

	obj.FillHistogram(dirname,"All_GammaGamma_Check",
			  1000,0,2000,dop1,1000,0,2000,dop2);

	if(s_mult == 2) {
	  if(nMap[det1].at(0) != det2 && nMap[det1].at(1) != det2) {
	    
	    obj.FillHistogram(dirname,"All_CoreEnGamma_M2_NN",
			      1000,0,2000,dop1,1000,0,2000,cen2);
	
	    obj.FillHistogram(dirname,"All_GammaGamma_M2_NN",
			      1000,0,2000,dop1,1000,0,2000,dop2);

	    obj.FillHistogram(dirname,"All_ReconCoreEnGamma_M2_NN",
			      1000,0,2000,rec1,1000,0,2000,cen2);
	
	    obj.FillHistogram(dirname,"All_ReconGammaGamma_M2_NN",
			      1000,0,2000,rec1,1000,0,2000,rec2);
	  }
	}
      
      } //End second SeGA loop
      
    } //End first SeGA loop
  } //End Janus loop

  return;
}

void MakeReconstructedSeGAJanus(TRuntimeObjects& obj, TSega& sega, TJanusDDAS& janus, GCutG* time_gate, GCutG* recon_gate) {

  if(!recon_gate || !time_gate) 
    {return;}
    
  static TSRIM b_srim(system_map[GValue::Value("REAC_KEY")].at(2).c_str());
  static TSRIM t_srim(system_map[GValue::Value("REAC_KEY")].at(3).c_str());
  
  // Factors of 1e3 are because TNucleus and TReaction use MeV, while TSRIM uses keV.
  // (mg/cm^2) / (mg/cm^3) * (10^4 um/cm) distance in um travelled through target, using SRIM to adjust reaction energy
  static TReaction reac = TReaction(std::make_shared<TNucleus>(system_map[GValue::Value("REAC_KEY")].at(0).c_str()),
				    std::make_shared<TNucleus>(system_map[GValue::Value("REAC_KEY")].at(1).c_str()),
				    std::make_shared<TNucleus>(system_map[GValue::Value("REAC_KEY")].at(0).c_str()),
				    std::make_shared<TNucleus>(system_map[GValue::Value("REAC_KEY")].at(1).c_str()),
				    b_srim.GetAdjustedEnergy(GValue::Value("BEAM_ENERGY")*1000,
				    (GValue::Value("TARGET_THICKNESS")/GValue::Value("TARGET_DENSITY"))*10000.0/2.0)/1000.0
				   );

  std::string gate_name = recon_gate->GetName();
  int det_part = -1;
  int recon_part = -1;
    
  if(!gate_name.compare("CdR") || !gate_name.compare("CdR_FW")) {
    det_part = 2;
    recon_part = 3;
  }
    
  else if(!gate_name.compare("TiR") || !gate_name.compare("TiRBig")  || !gate_name.compare("WR") || !gate_name.compare("Pb")
	  || !gate_name.compare("PbBig")) {
    det_part = 3;
    recon_part = 2;
  }
    
  else if(!gate_name.compare("UnKR")) {
    det_part = 3;
    recon_part = 2;
  }

  std::string d_pref = "";
  std::string dirname = "";
  std::string jCh_dirname = "";

  for(int lvl=min_level;lvl<max_level+1;lvl++) {

    if(lvl==0)
      {d_pref="Total_";}
    else if(lvl==1)
      {d_pref="Primary_";}
    else if(lvl==2)
      {d_pref="Secondary_";}
    else if(lvl==3)
      {d_pref="Tertiary_";}
    
    dirname = d_pref + "SeGAJanus_ReconstructingFrom" + recon_gate->GetName();
    jCh_dirname = d_pref + "SeGAJanus_jChargeGate_ReconstructingFrom" + recon_gate->GetName();
    if(time_gate){
      dirname += Form("_%s",time_gate->GetName());
      jCh_dirname += Form("_%s",time_gate->GetName());
    }

    std::vector<int> times_used;
    times_used.resize(sega.Size());
    std::fill(times_used.begin(),times_used.end(),0);
    
    int j_count = 0;
    for(auto& j_hit : janus.GetHits(lvl)) {
      if(j_hit.GetDetnum() != 1 || !recon_gate->IsInside(j_hit.GetRing(),j_hit.Charge())) 
        {continue;}
      j_count++;

      //BasicJanusSpectra(obj,j_hit,dirname,"All_");
 
      for(unsigned int i=0;i<sega.Size();i++) {
	auto& s_hit = sega.GetSegaHit(i);
	
        double s_energy = s_hit.GetEnergy();
        auto tdiff = (s_hit.Timestamp() - j_hit.Timestamp());

        if(time_gate) {
          if(!time_gate->IsInside(tdiff,s_energy)) 
            {continue;}
        }
	if(times_used.at(i))
	  {continue;}
	times_used.at(i)++;

	double num_segs = s_hit.GetNumSegments();

	if(num_segs == 0) {
	  obj.FillHistogram(dirname,"GammaEnergy_NoSegs",4000,0,4000,s_energy);
	}
	else if(num_segs > 0) {
	  obj.FillHistogram(dirname,"GammaEnergy_YesSegs",4000,0,4000,s_energy);
	}

	obj.FillHistogram(dirname,"GammaEnergy",4000,0,4000,s_energy);
	obj.FillHistogram(dirname,"GammaEnergy_summary",20,0,20,s_hit.GetDetnum(),4000,0,4000, s_energy);

	bool sol2 = false;
	if(!gate_name.compare("CdR") && j_hit.GetRing() < 3 && j_hit.Charge() < 3000)
	  {sol2=true;}
	
	//auto j_pos = j_hit.GetPosition(BEFORE);
	auto recon_pos = j_hit.GetReconPosition(reac,det_part,recon_part,sol2,BEFORE);

	//This will fix the doppler correction of the JANUS events with switched sectors
	//Histograms with the suffix jChargeGate show this nicely
	auto recon_pos_rot = j_hit.GetReconPosition(reac,det_part,recon_part,sol2,BEFORE);
	recon_pos_rot.SetPhi(recon_pos_rot.Phi() + TMath::Pi());
	
	//auto recon_pos_noS2 = j_hit.GetReconPosition(reac,det_part,recon_part,0,BEFORE);
        auto s_pos = s_hit.GetPosition();
	
        double recon_angle = s_pos.Angle(recon_pos)*TMath::RadToDeg();
	
	//double beta = get_beta(recon_pos.Theta(),recon_part,0);
	//double beta_noS2 = get_beta(recon_pos_noS2.Theta(),recon_part,0);

	double beta;
	//double beta_noS2;
	if(recon_part == 2) {
	  beta = TJanusDDAS::GetBeta(recon_pos.Theta(),reac,recon_part,b_srim,0);
          //beta_noS2 = TJanusDDAS::GetBeta(recon_pos_noS2.Theta(),reac,recon_part,b_srim,0);
	}
	else {
          beta = TJanusDDAS::GetBeta(recon_pos.Theta(),reac,recon_part,t_srim,0);
          //beta_noS2 = TJanusDDAS::GetBeta(recon_pos_noS2.Theta(),reac,recon_part,t_srim,0);
	}

	obj.FillHistogram(dirname,"Beta_v_Ring",26,0,26,j_hit.GetRing(),100,0,0.1,beta);
	obj.FillHistogram(dirname,"ReconAngle_v_Ring",26,0,26,j_hit.GetRing(),90,0,180,recon_angle);
	
        double s_energy_recon_dop = s_hit.GetDoppler(beta,recon_pos);
	double s_energy_recon_dop_rot = s_hit.GetDoppler(beta,recon_pos_rot);
	//double s_en_recon_dop_noS2 = s_hit.GetDoppler(beta_noS2,recon_pos_noS2);

	obj.FillHistogram(dirname,"GammaEnergyReconDC",12000,0,4000,s_energy_recon_dop);
	obj.FillHistogram(dirname,"GammaEnergyReconDC_summary",20,0,20,s_hit.GetDetnum(),8000,0,4000, s_energy_recon_dop);

	/*
	obj.FillHistogram(Form("%sRecoiling%s_jRings_%s",d_pref.c_str(),recon_gate->GetName(),time_gate->GetName()),
			  Form("GammaEnergy_R%02d",j_hit.GetRing()),
			  2000,0,4000,s_energy);

	obj.FillHistogram(Form("%sRecoiling%s_jRings_%s",d_pref.c_str(),recon_gate->GetName(),time_gate->GetName()),
			  Form("GammaEnergyReconDC_R%02d",j_hit.GetRing()),
			  2000,0,4000,s_energy_recon_dop);
	*/
	
	obj.FillHistogram(dirname,Form("ParticleCharge_v_GammaEnergy_jdet%d",j_hit.GetDetnum()),
			  2500,0,4000,s_energy,
			  2000,0,20000,j_hit.Charge());

	//obj.FillHistogram(dirname,Form("ParticleCharge_v_GammaEnergyReconDC_noS2_jdet%d",j_hit.GetDetnum()),
	//		  2500,0,4000,s_en_recon_dop_noS2,
	//		  2000,0,20000,j_hit.Charge());

	obj.FillHistogram(dirname,Form("ParticleCharge_v_GammaEnergyReconDC_jdet%d",j_hit.GetDetnum()),
			  2500,0,4000,s_energy_recon_dop,
			  2000,0,20000,j_hit.Charge());

	/*
	obj.FillHistogram(Form("%sRecoiling%s_jRings_%s",d_pref.c_str(),recon_gate->GetName(),time_gate->GetName()),
 			  Form("ParticleCharge_v_GammaEnergy_R%02d",j_hit.GetRing()),
 			  1000,0,2000,s_energy,
       	         	  3000,0,30000,j_hit.Charge());

 	obj.FillHistogram(Form("%sRecoiling%s_jRings_%s",d_pref.c_str(),recon_gate->GetName(),time_gate->GetName()),
 			  Form("ParticleCharge_v_GammaEnergyDCz_R%02d",j_hit.GetRing()),
			  1000,0,2000,s_energy_recon_dop,
			  2000,0,20000,j_hit.Charge());
	
	for(double b : BinCenters(180,0.0,0.09)) {
	  obj.FillHistogram(Form("%sRecoiling%s_jRings_BetaScan_%s",d_pref.c_str(),recon_gate->GetName(),time_gate->GetName()),
			    Form("GammaEnergyReconDC_BetaScan_R%02d",j_hit.GetRing()),
			    180,0,0.09,b,
			    2000,0,4000,s_hit.GetDoppler(b,recon_pos,TVector3(0,0,GValue::Value("Sega_Z_offset"))));
	}
	
	
	if(s_energy > 632.0 && s_energy < 634.0 && j_hit.GetRing() < 9) {

	  obj.FillHistogram(dirname,"GammaEnergy_EnGate_R01_R08",4000,0,4000,s_energy);
	  obj.FillHistogram(dirname,"GammaEnergyReconDC_EnGate_R01_R08",12000,0,4000,s_energy_recon_dop);
	  
	}
	*/
	
        obj.FillHistogram(dirname,"GammaEnergy_v_ReconAngle",4000,0,4000,s_energy,180,0,180,recon_angle);
	obj.FillHistogram(dirname,"GammaEnergyReconDC_v_ReconAngle",4000,0,4000,s_energy_recon_dop,180,0,180,recon_angle);

	obj.FillHistogram(dirname,"GammaEnergy_v_CosReconAngle",4000,0,4000,s_energy,
			  200,-1.1,1.1,TMath::Cos(recon_angle*TMath::DegToRad()));
	obj.FillHistogram(dirname,"GammaEnergyReconDC_v_CosReconAngle",4000,0,4000,s_energy_recon_dop,
			  200,-1.1,1.1,TMath::Cos(recon_angle*TMath::DegToRad()));
	
	/*
	obj.FillHistogram(Form("%sRecoiling%s_jRings_%s_ThetaCorr",d_pref.c_str(),recon_gate->GetName(),time_gate->GetName()),
			  Form("GammaEnergy_v_ReconAngle_R%02d",j_hit.GetRing()),
			  4000,0,4000,s_energy,
			  180,0,180,recon_angle);
	
	obj.FillHistogram(Form("%sRecoiling%s_jRings_%s_ThetaCorr",d_pref.c_str(),recon_gate->GetName(),time_gate->GetName()),
			  Form("GammaEnergyReconDC_v_ReconAngle_R%02d",j_hit.GetRing()),
			  4000,0,4000,s_energy_recon_dop,
			  180,0,180,recon_angle);
	*/
	
	TVector3 incBeam = TVector3(0.0,0.0,1.0);
        TVector3 reacPlane = recon_pos.Cross(incBeam);
	TVector3 reacPlane_rot = recon_pos_rot.Cross(incBeam);
        TVector3 detPlane = s_pos.Cross(incBeam);

	double reac_phi = reacPlane.Phi();
        if(reac_phi < 0)
	  {reac_phi += TMath::TwoPi();}

	double reac_phi_rot = reacPlane_rot.Phi();
        if(reac_phi_rot < 0)
	  {reac_phi_rot += TMath::TwoPi();}

        double det_phi = detPlane.Phi();
        if(det_phi < 0)
	  {det_phi += TMath::TwoPi();}

        double planeAng = reac_phi - det_phi;
        if(planeAng < 0)
	  {planeAng += TMath::TwoPi();}

	double planeAng_rot = reac_phi_rot - det_phi;
        if(planeAng_rot < 0)
	  {planeAng_rot += TMath::TwoPi();}

	//double thing1 = 33*TMath::Pi()/16.0;
	double thing1 = 65*TMath::Pi()/32.0;

	obj.FillHistogram(dirname,"GammaEnergy_v_PlaneAng",4000,0,4000,s_energy,32,0,thing1,planeAng);
	obj.FillHistogram(dirname,"GammaEnergyReconDC_v_PlaneAng",4000,0,4000,s_energy_recon_dop,32,0,thing1,planeAng);

	obj.FillHistogram(dirname,"GammaEnergy_v_PlaneAngRot",4000,0,4000,s_energy,32,0,thing1,planeAng_rot);
	obj.FillHistogram(dirname,"GammaEnergyReconDC_v_PlaneAngRot",4000,0,4000,s_energy_recon_dop_rot,32,0,thing1,planeAng_rot);

	/*
	obj.FillHistogram(Form("%sRecoiling%s_jRings_%s_PhiCorr",d_pref.c_str(),recon_gate->GetName(),time_gate->GetName()),
			  Form("GammaEnergy_v_PlaneAngle_R%02d",j_hit.GetRing()),
			  4000,0,4000,s_energy,
			  32,0,thing1,planeAng);
	
	obj.FillHistogram(Form("%sRecoiling%s_jRings_%s_PhiCorr",d_pref.c_str(),recon_gate->GetName(),time_gate->GetName()),
			  Form("GammaEnergyReconDC_v_PlaneAngle_R%02d",j_hit.GetRing()),
			  4000,0,4000,s_energy_recon_dop,
			  32,0,thing1,planeAng);

	obj.FillHistogram(Form("%sRecoiling%s_jRings_%s_PhiCorrRot",d_pref.c_str(),recon_gate->GetName(),time_gate->GetName()),
			  Form("GammaEnergy_v_PlaneAngleRot_R%02d",j_hit.GetRing()),
			  4000,0,4000,s_energy,
			  32,0,thing1,planeAng_rot);
	
	obj.FillHistogram(Form("%sRecoiling%s_jRings_%s_PhiCorrRot",d_pref.c_str(),recon_gate->GetName(),time_gate->GetName()),
			  Form("GammaEnergyReconDC_v_PlaneAngleRot_R%02d",j_hit.GetRing()),
			  4000,0,4000,s_energy_recon_dop_rot,
			  32,0,thing1,planeAng_rot);
	*/
	if(j_hit.Charge() > 4650 && j_hit.Charge() < 5500 && !gate_name.compare("CdR")) {

	  /*
	  if(janus.GetDetectorHits(1,lvl).size() < 3) {
	    obj.FillHistogram(jCh_dirname,
			      Form("GammaEnergy_v_ReconAngle_jMult%d",(int)janus.GetDetectorHits(1,lvl).size()),
			      2000,0,4000,s_energy,
			      90,0,180,recon_angle);
	    obj.FillHistogram(jCh_dirname,
			      Form("GammaEnergyReconDC_v_ReconAngle_jMult%d",
				   (int)janus.GetDetectorHits(1,lvl).size()),
			      2000,0,4000,s_energy_recon_dop,
			      180,0,180,recon_angle);

	    obj.FillHistogram(jCh_dirname,
			      Form("GammaEnergy_v_PlaneAng_jMult%d",(int)janus.GetDetectorHits(1,lvl).size()),
			      2000,0,4000,s_energy,
			      32,0,thing1,planeAng);
	    obj.FillHistogram(jCh_dirname,
			      Form("GammaEnergyReconDC_v_PlaneAng_jMult%d",(int)janus.GetDetectorHits(1,lvl).size()),
			      2000,0,4000,s_energy_recon_dop,
			      32,0,thing1,planeAng);
	  }
	  else if(janus.GetDetectorHits(1,lvl).size() > 2) {
	    obj.FillHistogram(jCh_dirname,"GammaEnergy_v_ReconAngle_jMult3+",
			      2000,0,4000,s_energy,
			      90,0,180,recon_angle);
	    obj.FillHistogram(jCh_dirname,"GammaEnergyReconDC_v_ReconAngle_jMult3+",
			      2000,0,4000,s_energy_recon_dop,
			      180,0,180,recon_angle);

	    obj.FillHistogram(jCh_dirname,"GammaEnergy_v_PlaneAng_jMult3+",
			      2000,0,4000,s_energy,
			      32,0,thing1,planeAng);
	    obj.FillHistogram(jCh_dirname,"GammaEnergyReconDC_v_PlaneAng_jMult3+",
			      2000,0,4000,s_energy_recon_dop,
			      32,0,thing1,planeAng);
	  }
	  */

	  /*
	  obj.FillHistogram(jCh_dirname,"GammaEnergy_v_ReconAngle",
			    4000,0,4000,s_energy,
			    90,0,180,recon_angle);
	  obj.FillHistogram(jCh_dirname,"GammaEnergyReconDC_v_ReconAngle",
			    4000,0,4000,s_energy_recon_dop,
			    180,0,180,recon_angle);
	  */

	  obj.FillHistogram(jCh_dirname,"GammaEnergy_v_PlaneAng",
			    4000,0,4000,s_energy,
			    32,0,thing1,planeAng);
	  obj.FillHistogram(jCh_dirname,"GammaEnergyReconDC_v_PlaneAng",
			    4000,0,4000,s_energy_recon_dop,
			    32,0,thing1,planeAng);

	  obj.FillHistogram(jCh_dirname,"GammaEnergy_v_PlaneAngRot",
			    4000,0,4000,s_energy,
			    32,0,thing1,planeAng_rot);
	  obj.FillHistogram(jCh_dirname,"GammaEnergyReconDC_v_PlaneAngRot",
			    4000,0,4000,s_energy_recon_dop_rot,
			    32,0,thing1,planeAng_rot);

	  //obj.FillHistogram(Form("%sRecoiling%s_jRings_%s",d_pref.c_str(),recon_gate->GetName(),time_gate->GetName()),
	  //		    Form("GammaEnergyReconDC_v_PlaneAng_jRing%02d",j_hit.GetRing()),
	  //		    4000,0,4000,s_energy_recon_dop,
	  //		    32,0,thing1,planeAng);
	}
	else if(j_hit.Charge() > 4500 && j_hit.Charge() < 5300 && (!gate_name.compare("TiR") || !gate_name.compare("TiRBig"))) {

	  /*
	  if(janus.GetDetectorHits(1,lvl).size() < 3) {
	    obj.FillHistogram(jCh_dirname,
			      Form("GammaEnergy_v_ReconAngle_jMult%d",(int)janus.GetDetectorHits(1,lvl).size()),
			      2000,0,4000,s_energy,
			      90,0,180,recon_angle);
	    obj.FillHistogram(jCh_dirname,
			      Form("GammaEnergyReconDC_v_ReconAngle_jMult%d",
				   (int)janus.GetDetectorHits(1,lvl).size()),
			      2000,0,4000,s_energy_recon_dop,
			      180,0,180,recon_angle);

	    obj.FillHistogram(jCh_dirname,
			      Form("GammaEnergy_v_PlaneAng_jMult%d",(int)janus.GetDetectorHits(1,lvl).size()),
			      2000,0,4000,s_energy,
			      32,0,thing1,planeAng);
	    obj.FillHistogram(jCh_dirname,
			      Form("GammaEnergyReconDC_v_PlaneAng_jMult%d",(int)janus.GetDetectorHits(1,lvl).size()),
			      2000,0,4000,s_energy_recon_dop,
			      32,0,thing1,planeAng);
	  }
	  else if(janus.GetDetectorHits(1,lvl).size() > 2) {
	    obj.FillHistogram(jCh_dirname,"GammaEnergy_v_ReconAngle_jMult3+",
			      2000,0,4000,s_energy,
			      90,0,180,recon_angle);
	    obj.FillHistogram(jCh_dirname,"GammaEnergyReconDC_v_ReconAngle_jMult3+",
			      2000,0,4000,s_energy_recon_dop,
			      180,0,180,recon_angle);

	    obj.FillHistogram(jCh_dirname,"GammaEnergy_v_PlaneAng_jMult3+",
			      2000,0,4000,s_energy,
			      32,0,thing1,planeAng);
	    obj.FillHistogram(jCh_dirname,"GammaEnergyReconDC_v_PlaneAng_jMult3+",
			      2000,0,4000,s_energy_recon_dop,
			      32,0,thing1,planeAng);
	  }
	  */

	  /*
	  obj.FillHistogram(jCh_dirname,"GammaEnergy_v_ReconAngle",
			    4000,0,4000,s_energy,
			    90,0,180,recon_angle);
	  obj.FillHistogram(jCh_dirname,"GammaEnergyReconDC_v_ReconAngle",
			    4000,0,4000,s_energy_recon_dop,
			    180,0,180,recon_angle);
	  */
	  
	  obj.FillHistogram(jCh_dirname,"GammaEnergy_v_PlaneAng",
			    4000,0,4000,s_energy,
			    32,0,thing1,planeAng);
	  obj.FillHistogram(jCh_dirname,"GammaEnergyReconDC_v_PlaneAng",
			    4000,0,4000,s_energy_recon_dop,
			    32,0,thing1,planeAng);

	  obj.FillHistogram(jCh_dirname,"GammaEnergy_v_PlaneAngRot",
			    4000,0,4000,s_energy,
			    32,0,thing1,planeAng_rot);
	  obj.FillHistogram(jCh_dirname,"GammaEnergyReconDC_v_PlaneAngRot",
			    4000,0,4000,s_energy_recon_dop_rot,
			    32,0,thing1,planeAng_rot);

	  //obj.FillHistogram(Form("%sRecoiling%s_jRings_%s",d_pref.c_str(),recon_gate->GetName(),time_gate->GetName()),
	  //		    Form("GammaEnergyReconDC_v_PlaneAng_jRing%02d",j_hit.GetRing()),
	  //		    4000,0,4000,s_energy_recon_dop,
	  //		    32,0,thing1,planeAng);
	}
	else if(!gate_name.compare("CdR") || !gate_name.compare("TiR") || !gate_name.compare("TiRBig")) {
	  
	  /*
	  obj.FillHistogram(jCh_dirname,"GammaEnergy_v_ReconAngle_NotjChargeGate",
			    4000,0,4000,s_energy,
			    90,0,180,recon_angle);
	  obj.FillHistogram(jCh_dirname,"GammaEnergyDCz_v_ReconAngle_NotjChargeGate",
			    4000,0,4000,s_energy_recon_dop,
			    180,0,180,recon_angle);
	  */
	  
	  obj.FillHistogram(jCh_dirname,"GammaEnergy_v_PlaneAng_NotjChargeGate",
			    4000,0,4000,s_energy,
			    32,0,thing1,planeAng);
	  obj.FillHistogram(jCh_dirname,"GammaEnergyReconDC_v_PlaneAng_NotjChargeGate",
			    4000,0,4000,s_energy_recon_dop,
			    32,0,thing1,planeAng);

	  obj.FillHistogram(jCh_dirname,"GammaEnergy_v_PlaneAngRot_NotjChargeGate",
			    4000,0,4000,s_energy,
			    32,0,thing1,planeAng_rot);
	  obj.FillHistogram(jCh_dirname,"GammaEnergyReconDC_v_PlaneAngRot_NotjChargeGate",
			    4000,0,4000,s_energy_recon_dop_rot,
			    32,0,thing1,planeAng_rot);
	}

	/*
	  obj.FillHistogram(dirname,"GammaEnergyReconDC_allseg_summary",
	  512,1,513,32*(s_hit.GetDetnum() - 1) + s_hit.GetMainSegnum(),
	  35,600,670,s_energy_recon_dop);

	  obj.FillHistogram(dirname,"GammaEnergy_SegMap_summary",
	  512,1,513,s_hit.GetMapSegnum(),
	  70,570,710,s_energy);

	  obj.FillHistogram(dirname,"GammaEnergy_PairMap_summary",
	  256,1,257,s_hit.GetMapPairnum(),
	  70,570,710,s_energy);

	  obj.FillHistogram(dirname,"GammaEnergy_SliceMap_summary",
	  128,1,129,s_hit.GetMapSlicenum(),
	  70,570,710,s_energy);

	  obj.FillHistogram(dirname,"GammaEnergyReconDC_SegMap_summary",
	  512,1,513,s_hit.GetMapSegnum(),
	  35,600,670,s_energy_recon_dop);

	  obj.FillHistogram(dirname,"GammaEnergyReconDC_PairMap_summary",
	  256,1,257,s_hit.GetMapPairnum(),
	  35,600,670,s_energy_recon_dop);

	  obj.FillHistogram(dirname,"GammaEnergyReconDC_SliceMap_summary",
	  128,1,129,s_hit.GetMapSlicenum(),
	  35,600,670,s_energy_recon_dop);
 
        */
      }//end sega hit loop

      //obj.FillHistogram(dirname,"SeGA_Mult",25,0,25,s_count);
      
    }//end janus hit loop

    obj.FillHistogram(dirname,"Janus_HitMult",100,0,100,janus.GetDetectorHits(1,lvl).size());
    obj.FillHistogram(dirname,"Janus_GateMult",100,0,100,j_count);

    //for(auto uses : times_used) {
    //obj.FillHistogram(dirname,"SeGA_Hit_Uses",10,0,10,uses);
    //}

  }//end janus level loop

}

int gates_loaded=0;
std::vector<GCutG*> scatter_gates;
std::vector<GCutG*> recon_gates;
std::vector<GCutG*> time_gates;
std::vector<GCutG*> janus_gates;

// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
  
  TSega* sega = obj.GetDetector<TSega>();
  TJanusDDAS* janus = obj.GetDetector<TJanusDDAS>();

  TList *gates = &(obj.GetGates());
  if(gates_loaded!=gates->GetSize()) {
    TIter iter(gates);
    while(TObject *obj = iter.Next()) {
      GCutG *gate = (GCutG*)obj;
      std::string tag = gate->GetTag();
      if(!tag.compare("Scatter")) {
        scatter_gates.push_back(gate);
	janus_gates.push_back(gate);
        std::cout << "Scatter Gate: << " << gate->GetName() << std::endl;
      }
      else if(!tag.compare("Reconstruct")) {
        recon_gates.push_back(gate);
	janus_gates.push_back(gate);
        std::cout << "Recon Gate: << " << gate->GetName() << std::endl;
      }
      else if(!tag.compare("time")) {
        time_gates.push_back(gate); 
        std::cout << "Time Gate: << " << gate->GetName() << std::endl;
      }
      gates_loaded++;
    }
  } 

  /*
  if(sega) {
    MakeSega(obj,*sega);
  }
  
  if(janus) {
    MakeJanus(obj,*janus);
    
    for(auto &gate : janus_gates){
      MakeGatedJanus(obj,*janus,gate);
    }

    if(janus_gates.size() > 1) {
      MakeDoubleGatedJanus(obj,*janus,janus_gates.at(0),janus_gates.at(1));
    }
  }
  */
  
  if(sega && janus) {
    //MakeSeGAJanus(obj,*sega,*janus);
      
    for(auto &time_gate : time_gates) {
      for(auto &scatter_gate : scatter_gates) {
	//MakeScatteredSeGAJanus(obj,*sega,*janus,time_gate,scatter_gate);
	GammaGamma(obj,*sega,*janus,time_gate,scatter_gate);
      }

      
      for(auto &recon_gate : recon_gates) {
	//MakeReconstructedSeGAJanus(obj,*sega,*janus,time_gate,recon_gate);
	GammaGamma(obj,*sega,*janus,time_gate,recon_gate);
      }
      
    }    
  }
  
  
}


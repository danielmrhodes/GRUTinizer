#include "TRuntimeObjects.h"
#include "TJanusDDAS.h"
#include "TSega.h"
#include "GCutG.h"

#include "TRandom.h"

//80Ge
/*
static const double beam_mass = 74441.6;
static const double beam_Z = 32;
static const double targ_mass = 182540.0;
static const double targ_Z = 78.0;
//static const double Ep = 268.97;
static const double beam_en = 281.6;
static const double Ex = 0.0;
*/

//80Kr

static const double beam_mass = 74442.6;
static const double beam_Z = 36;
static const double targ_mass = 182540.0;
static const double targ_Z = 78.0;
static const double beam_en = 281.6;
static const double Ex = 0.0;


/////////Kinematics/////////
double Theta_CM_FP(double ThetaLAB, double Ep = beam_en, bool sol2=false) {
  
  double tau = (beam_mass/targ_mass)/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  
  if(std::sin(ThetaLAB) > 1.0/tau) {
    ThetaLAB = std::asin(1.0/tau);

    if(ThetaLAB < 0) {
      ThetaLAB += TMath::Pi();
    }

    return std::asin(tau*std::sin(ThetaLAB)) + ThetaLAB;
  }

  if(!sol2) {
    return std::asin(tau*std::sin(ThetaLAB)) + ThetaLAB;
  }
  
  return std::asin(tau*std::sin(-ThetaLAB)) + ThetaLAB + TMath::Pi();

}

double Theta_CM_FR(double ThetaLAB, double Ep = beam_en, bool sol2=false) {
  
  double tau = 1.0/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  
  if(std::sin(ThetaLAB) > 1.0/tau) {
    ThetaLAB = std::asin(1.0/tau);

    if(ThetaLAB < 0) {
      ThetaLAB += TMath::Pi();
    }

    return std::asin(tau*std::sin(ThetaLAB)) + ThetaLAB;
  }

  if(!sol2) {
    return TMath::Pi() - (std::asin(tau*std::sin(ThetaLAB)) + ThetaLAB);
  }
  
  return -std::asin(tau*std::sin(-ThetaLAB)) - ThetaLAB;

}

double Theta_LAB(double thetaCM, double Ep = beam_en) {

  double tau = (beam_mass/targ_mass)/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  double tanTheta = std::sin(thetaCM)/(std::cos(thetaCM) + tau);

  if(tanTheta > 0) {
    return std::atan(tanTheta);
  }
  
  return std::atan(tanTheta) + TMath::Pi();
  
}

double Theta_LAB_Max(double Ep = beam_en) {

  double tau = (beam_mass/targ_mass)/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));

  if(tau < 1.0) {
    return TMath::Pi();
  }
  
  return std::asin(1.0/tau);
  
}

double Recoil_Theta_LAB(double thetaCM, double Ep = beam_en) {

  double tau = 1.0/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  double tanTheta = std::sin(TMath::Pi() - thetaCM)/(std::cos(TMath::Pi() - thetaCM) + tau);
  
  return std::atan(tanTheta);
  
}

double Recoil_Theta_LAB_Max(double Ep = beam_en) {

  double tau = 1.0/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  
  return std::asin(1.0/tau);
    
}

double KE_LAB(double thetaCM, double Ep = beam_en) {

  double tau = (beam_mass/targ_mass)/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));

  double term1 = std::pow(targ_mass/(beam_mass + targ_mass),2);
  double term2 = 1 + tau*tau + 2*tau*std::cos(thetaCM);
  double term3 = Ep - Ex*(1 + beam_mass/targ_mass);
  
  return term1*term2*term3;
}

double Recoil_KE_LAB(double thetaCM, double Ep = beam_en) {

  double tau = 1.0/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));

  double term1 = beam_mass*targ_mass/std::pow(beam_mass + targ_mass,2);
  double term2 = 1 + tau*tau + 2*tau*std::cos(TMath::Pi() - thetaCM);
  double term3 = Ep - Ex*(1 + beam_mass/targ_mass);
  
  return term1*term2*term3;
}

double Beta(double energy, double mass) {

  double gam = energy/mass + 1.0;
  
  return TMath::Sqrt(1.0 - 1.0/(gam*gam));
  
}
////////////////////////////

std::string Prefix(int level) {

  switch(level) {

  case 0:
    return "Total_";

  case 1:
    return "Primary_";

  case 2:
    return "Secondary_";

  case 3:
    return "Tertiary_";

  }

  return "_";
  
}

TRandom* rand_gen;
const std::vector<double> wT = {1.84304,1.78991,1.7357,1.68088,1.62586,1.57101,1.51664,1.46303,1.41041,
				1.35895,1.30881,1.26011,1.21292,1.16731,1.1233,1.08092,1.04017,1.00103,
				0.963486,0.927497,0.893027,0.860034,0.828469,0.798282};

const double wP = 360.0/32.0;

void BasicJanusSpectra(TRuntimeObjects& obj, TJanusDDASHit& hit, std::string dirname, std::string h_suf) {

  const int det = hit.GetDetnum();
  const int ring = hit.GetRing();
  const int sect = hit.GetSector();

  const int sCharge = hit.Charge();
  const int rCharge = hit.GetBackHit().Charge();
  
  const TVector3 pos = hit.GetPosition(false);
  const TVector3 pos_NoOff = hit.GetPosition(false,false);

  const double theta = pos.Theta()*TMath::RadToDeg();
  const double phi = pos.Phi()*TMath::RadToDeg();
  
  const double theta_r = rand_gen->Uniform(theta-(wT.at(ring-1)/2.0),theta+(wT.at(ring-1)/2.0));
  const double phi_r = rand_gen->Uniform(phi-(wP/2.0),phi+(wP/2.0));
  
  obj.FillHistogram(dirname,Form("%sSectorCharge_Hit_det%d",h_suf.c_str(),det),
		    3000,0,30000,sCharge);
  
  obj.FillHistogram(dirname,Form("%sRingCharge_Hit_det%d",h_suf.c_str(),det),
		    3000,0,30000,rCharge);

  /*
  obj.FillHistogram(dirname,Form("%sCharge_Summary_Hit",h_suf.c_str()),
                    128,0,128,hit.GetFrontChannel(),
                    3000,0,30000,sCharge);
  obj.FillHistogram(dirname,Form("%sCharge_Summary_Hit",h_suf.c_str()),
                    128,0,128,hit.GetBackChannel(),
                    3000,0,30000,rCharge);
  */
  
  obj.FillHistogram(dirname,Form("%sRing_det%d",h_suf.c_str(),det),26,0,26,ring);
  obj.FillHistogram(dirname,Form("%sSector_det%d",h_suf.c_str(),det),34,0,34,sect);

  obj.FillHistogram(dirname,Form("%sSectorCharge_v_Sector_det%d",h_suf.c_str(),det),
  		    34,0,34,sect,3000,0,30000,sCharge);

  obj.FillHistogram(dirname,Form("%sSectorCharge_v_Ring_det%d",h_suf.c_str(),det),
		    26,0,26,ring,3000,0,30000,sCharge);

  obj.FillHistogram(dirname,Form("%sSectorCharge_v_Theta_det%d",h_suf.c_str(),det),
		    200,-10,190,theta,3000,0,30000,sCharge);

  obj.FillHistogram(dirname,Form("%sRingCharge_v_Sector_det%d",h_suf.c_str(),det),
		    34,0,34,sect,3000,0,30000,rCharge);
  
  obj.FillHistogram(dirname,Form("%sRingCharge_v_Ring_det%d",h_suf.c_str(),det),
		    26,0,26,ring,3000,0,30000,rCharge);

  obj.FillHistogram(dirname,Form("%sRingCharge_v_Theta_det%d",h_suf.c_str(),det),
		    200,-10,190,theta,3000,0,30000,rCharge);

  obj.FillHistogram(dirname,Form("%sRing_v_Sector_det%d",h_suf.c_str(),det),
		    26,0,26,ring,34,0,34,sect);

  obj.FillHistogram(dirname,Form("%sPhi_v_Theta_det%d",h_suf.c_str(),det),
		    900,0,180,theta,
		    100,-200,200,phi);
  
  obj.FillHistogram(dirname,Form("%sPhi_v_Theta_Rand_det%d",h_suf.c_str(),det),
		    900,0,180,theta_r,
		    100,-200,200,phi_r);

  double shift = 0.5*TMath::TwoPi()/32.0;
  obj.FillHistogram(dirname,Form("%sPerp_v_Phi_NoOffset_det%d",h_suf.c_str(),det),
		    32,-TMath::Pi()+shift,TMath::Pi()+shift,pos_NoOff.Phi(),
		    24,1.1,3.5,pos_NoOff.Perp());

  obj.FillHistogram(dirname,Form("%sY_v_X_det%d",h_suf.c_str(),det),
		    200,-4,4,pos.X(),
		    200,-4,4,pos.Y());

  //obj.FillHistogram(dirname,Form("%sSectorCharge_v_RingCharge_Hit_det%d",h_suf.c_str(),det),
  //		    2000,0,30000,rCharge,
  //		    2000,0,30000,sCharge);

  /*
  auto tdiff = hit.GetBackHit().Timestamp() - hit.Timestamp();
   
  obj.FillHistogram(dirname,Form("%sFrontBack_Tdiff_Hit_det%d",h_suf.c_str(),det),
		    2750,-11000,11000,tdiff);

  auto cdiff = rCharge - sCharge;
  
  obj.FillHistogram(dirname,Form("%sCdiff_Hit_det%d",h_suf.c_str(),det),2000,0,20000,cdiff);

  obj.FillHistogram(dirname,Form("%sCdiff_v_Tdiff_Hit_det%d",h_suf.c_str(),det),
		    2500,-10000,10000,tdiff,
		    1000,0,30000,cdiff);

  obj.FillHistogram(dirname,Form("%sRingCharge_v_Tdiff_Hit_det%d",h_suf.c_str(),det),
		    2500,-10000,10000, tdiff,
		    1000,0,30000,rCharge);

  obj.FillHistogram(dirname,Form("%sRingCharge_v_Cdiff_Hit_det%d",h_suf.c_str(),det),
		    5000,0,20000,cdiff,
		    1000,0,30000,rCharge);

  obj.FillHistogram(dirname,Form("%sSectorCharge_v_Tdiff_Hit_det%d",h_suf.c_str(),det),
		    2500,-10000,10000,tdiff,
		    1000,0,30000,sCharge);

  obj.FillHistogram(dirname,Form("%sSectorCharge_v_Cdiff_Hit_det%d",h_suf.c_str(),det),
		    5000,0,20000,cdiff,
		    1000,0,30000,sCharge);
  */
  
  return;
 
}

void LeveledJanusSpectra(TRuntimeObjects& obj, TJanusDDAS& janus, int lvl) {

  std::string dirname = Prefix(lvl) + "Janus_Ungated";
  
  for(size_t i=0; i<janus.GetUnusedChannels(lvl).size(); i++) {
    
    auto &chan = janus.GetUnusedChannels(lvl).at(i);
    
    obj.FillHistogram(dirname,"Charge_Summary_UnusedChans",
		      128,0,128,chan.GetFrontChannel(),
		      3500,0,35000,chan.Charge());
    
    for(size_t j=i+1;j<janus.GetUnusedChannels(lvl).size();j++) {
      
      auto &chan1 = janus.GetUnusedChannels(lvl).at(j);
      
      if(chan1.GetDetnum() == chan.GetDetnum()) {
        obj.FillHistogram(dirname,Form("Tdiff_UnusedChans_det%d",chan.GetDetnum()),
	  		  2500,-11000,11000,chan.Timestamp()-chan1.Timestamp());
      }
        
    }//end second janus chan loop 
  }//end janus chan loop

  for(size_t i=0;i<janus.GetUnusedSectors(lvl).size();i++) {
    for(size_t j=0;j<janus.GetUnusedRings(lvl).size();j++) {
      
      auto &sect = janus.GetUnusedSectors(lvl).at(i);
      auto &ring = janus.GetUnusedRings(lvl).at(j);

      auto tdiff = ring.Timestamp() - sect.Timestamp();

      if(ring.GetDetnum() == sect.GetDetnum()) {

	obj.FillHistogram(dirname,Form("SectorCharge_v_TDiff_UnusedChans_det%d",sect.GetDetnum()),
		          2500,-11000,11000,tdiff,
		          2000,0,30000,sect.Charge());
	
        obj.FillHistogram(dirname,Form("SectorCharge_v_RingCharge_UnusedChans_det%d",sect.GetDetnum()),
		          2000,0,30000,ring.Charge(),
		          2000,0,30000,sect.Charge());

	obj.FillHistogram(dirname,Form("Sector_Ring_TDiff_UnusedChans_det%d",sect.GetDetnum()),
			  2500,-11000,11000,tdiff);
      }
    }
  }
  
  obj.FillHistogram(dirname,"All_HitMult_det1",50,0,50,janus.GetSpecificHits(1,lvl,0).size());   
  obj.FillHistogram(dirname,"All_HitMult_det0",50,0,50,janus.GetSpecificHits(0,lvl,0).size());
  for(auto& hit : janus.GetHits(lvl)) {
    BasicJanusSpectra(obj,hit,dirname,"All_"); 
  }
  
  obj.FillHistogram(dirname,"Single_HitMult_det1",50,0,50,janus.GetSpecificHits(1,lvl,1).size());   
  obj.FillHistogram(dirname,"Single_HitMult_det0",50,0,50,janus.GetSpecificHits(0,lvl,1).size());
  for(auto& hit : janus.GetSingleHits(lvl)) {
    BasicJanusSpectra(obj,hit,dirname,"Single_"); 
  }

  obj.FillHistogram(dirname,"Double_HitMult_det1",50,0,50,janus.GetSpecificHits(1,lvl,2).size());   
  obj.FillHistogram(dirname,"Double_HitMult_det0",50,0,50,janus.GetSpecificHits(0,lvl,2).size());
  for(auto& hit : janus.GetDoubleHits(lvl)) {
    BasicJanusSpectra(obj,hit,dirname,"Double_"); 
  }

  obj.FillHistogram(dirname,"Addback_HitMult_det1",50,0,50,janus.GetSpecificHits(1,lvl,3).size());   
  obj.FillHistogram(dirname,"Addback_HitMult_det0",50,0,50,janus.GetSpecificHits(0,lvl,3).size());
  for(auto& hit : janus.GetAddbackHits(lvl)) {
    BasicJanusSpectra(obj,hit,dirname,"Addback_"); 
  }
  
  return;

}

static long prev_ts = 0;
static long prev_chan_ts = 0;
void MakeJanus(TRuntimeObjects& obj, TJanusDDAS& janus) {
  
  std::string dirname = "Janus_Ungated";

  if(prev_ts > 0) {
    obj.FillHistogram(dirname,"Event_TDiff",100000,0,10000000000,janus.Timestamp()-prev_ts);
  }
  prev_ts = janus.Timestamp();
  
  std::vector<TJanusDDASHit> rings;
  std::vector<TJanusDDASHit> sectors;
  for(auto& chan : janus.GetAllChannels()) {

    obj.FillHistogram(dirname,"Charge_Summary_AllChans",
		      128,0,128,chan.GetFrontChannel(),
		      6000,0,70000,chan.Charge());

    if(chan.IsRing()) {
      rings.push_back(chan);
    }
    else {
      sectors.push_back(chan);
    }

  }//end all channel loop

  //std::string dir_tmp = "Janus_Ungated_Rings";
  //std::string dir_tmp1 = "Janus_Ungated_Sectors";
  for(auto& sect : sectors) {
    for(auto& ring : rings) {

      if(ring.GetDetnum() == sect.GetDetnum()) {
	
        obj.FillHistogram(dirname,Form("SectorCharge_v_RingCharge_AllChans_det%d",sect.GetDetnum()),
		          2000,0,36000,ring.Charge(),
		          2000,0,36000,sect.Charge());

	if(ring.GetSector() > 1) {

	  obj.FillHistogram(dirname,Form("SectorCharge_v_RingCharge_AllChans_NoR01_det%d",sect.GetDetnum()),
		          2000,0,36000,ring.Charge(),
		          2000,0,36000,sect.Charge());
	}

	/*
	obj.FillHistogram(dir_tmp,Form("SectorCharge_v_RingCharge_AllChans_det%d_R%02d",
				       sect.GetDetnum(),ring.GetSector()),
		          2000,0,36000,ring.Charge(),
		          2000,0,36000,sect.Charge());

	obj.FillHistogram(dir_tmp1,Form("SectorCharge_v_RingCharge_AllChans_det%d_S%02d",
				       sect.GetDetnum(),sect.GetSector()),
		          2000,0,36000,ring.Charge(),
		          2000,0,36000,sect.Charge());
	*/
      }
    } //end all ring loop
  } //end all sector loop
  
  for(size_t i=0; i<janus.GetChannels().size(); i++) {
    auto &chan = janus.GetChannels().at(i);

    if(prev_chan_ts > 0) {
      obj.FillHistogram(dirname,"Chan_TDiff",100000,0,10000000000,chan.Timestamp()-prev_chan_ts);
    }
    prev_chan_ts = chan.Timestamp();
    
    obj.FillHistogram(dirname,"Charge_Summary_UsableChans",
		      128,0,128,chan.GetFrontChannel(),
		      3500,0,35000,chan.Charge());

    for(size_t j=i+1;j<janus.GetChannels().size();j++) {
      auto &chan1 = janus.GetChannels().at(j);

      if(chan1.GetDetnum() == chan.GetDetnum()) {
        obj.FillHistogram(dirname,Form("Tdiff_UsableChans_det%d",chan.GetDetnum()),
	  		  2500,-11000,11000,chan.Timestamp()-chan1.Timestamp());
      }
      
    }//end second janus chan loop 
  }//end janus chan loop

  for(auto& sect : janus.GetSectors()) {
    for(auto& ring : janus.GetRings()) {
  
      auto tdiff = ring.Timestamp() - sect.Timestamp();
      
      if(ring.GetDetnum() == sect.GetDetnum()) {

	obj.FillHistogram(dirname,Form("SectorCharge_v_TDiff_UsableChans_det%d",sect.GetDetnum()),
		          2500,-11000,11000,tdiff,
		          2000,0,30000,sect.Charge());
	
        obj.FillHistogram(dirname,Form("SectorCharge_v_RingCharge_UsableChans_det%d",sect.GetDetnum()),
		          2000,0,30000,ring.Charge(),
		          2000,0,30000,sect.Charge());

	obj.FillHistogram(dirname,Form("Sector_Ring_TDiff_UsableChans_det%d",sect.GetDetnum()),
			  6000,-11000,11000,tdiff);
	
	if(ring.GetSector() > 1) {
	  obj.FillHistogram(dirname,Form("SectorCharge_v_RingCharge_UsableChans_NoR01_det%d",
					 sect.GetDetnum()),
			    2000,0,30000,ring.Charge(),
			    2000,0,30000,sect.Charge());

	  obj.FillHistogram(dirname,Form("SectorCharge_v_TDiff_UsableChans_NoR01_det%d",sect.GetDetnum()),
		          2500,-11000,11000,tdiff,
		          2000,0,30000,sect.Charge());
	}
	
      } 
    }//end ring loop
  }//end sector loop
  
  LeveledJanusSpectra(obj,janus,1);

  return;
  
}

void MakeGatedJanus(TRuntimeObjects& obj, TJanusDDAS& janus, GCutG* gate) {

  if(!gate)
    {return;}

  int lvl = 1;
  std::string dirname = Prefix(lvl) + Form("Janus_%s",gate->GetName());

  int detNum = 1;
  std::string tag = gate->GetTag();
  if(!tag.compare("D0")) {
    detNum = 0;
  }

  int mult = 0;
  for(auto& hit : janus.GetHits(lvl)) {
      
    if(hit.GetDetnum() != detNum || !gate->IsInside(hit.GetRing(),hit.Charge())) {
      continue;
    }

    mult++;
    BasicJanusSpectra(obj,hit,dirname,"All_");

  }
  obj.FillHistogram(dirname,Form("All_GateMult_det%d",detNum),20,0,20,mult);

  int multS = 0;
  for(auto& hit : janus.GetSingleHits(lvl)) {
      
    if(hit.GetDetnum() != detNum || !gate->IsInside(hit.GetRing(),hit.Charge())) {
      continue;
    }

    multS++;
    BasicJanusSpectra(obj,hit,dirname,"Single_");

  }
  obj.FillHistogram(dirname,Form("Single_GateMult_det%d",detNum),20,0,20,multS);

  /*
  int multD = 0;
  for(auto& hit : janus.GetDoubleHits(lvl)) {
      
    if(hit.GetDetnum() != detNum || !gate->IsInside(hit.GetRing(),hit.Charge())) {
      continue;
    }

    multD++;
    BasicJanusSpectra(obj,hit,dirname,"Double_");

  }
  obj.FillHistogram(dirname,Form("Double_GateMult_det%d",detNum),20,0,20,multD);
  */

  int multA = 0;
  for(auto& hit : janus.GetAddbackHits(lvl)) {
      
    if(hit.GetDetnum() != detNum || !gate->IsInside(hit.GetRing(),hit.Charge())) {
      continue;
    }

    multA++;
    BasicJanusSpectra(obj,hit,dirname,"Addback_");

  }
  obj.FillHistogram(dirname,Form("Addback_GateMult_det%d",detNum),20,0,20,multA);

  return;

}

void MakeSega(TRuntimeObjects& obj, TSega& sega) {

  std::string dirname = "SeGA";

  obj.FillHistogram(dirname,"Mult",20,0,20,sega.Size());

  //int index = -1;
  //for(auto hit : sega) {
  for(unsigned int i=0;i<sega.Size();i++) {

    TSegaHit hit = sega.GetSegaHit(i);
    //int num_segs = hit.GetNumSegments();

    //if(hit.GetDetnum() == 1) {
    //index = i;
    //}

    obj.FillHistogram(dirname,"CoreEnergy",4000,0,4000,hit.GetEnergy());
    obj.FillHistogram(dirname,"CoreCharge",10000,0,33000,hit.Charge());
    
    obj.FillHistogram("SeGA_Detectors",Form("CoreEnergy_det%02d",hit.GetDetnum()),
		      4000,0,4000,hit.GetEnergy());
    
    obj.FillHistogram("SeGA_Detectors",Form("CoreCharge_det%02d",hit.GetDetnum()),
		      9000,0,36000,hit.Charge());

    //obj.FillHistogram("SeGA_Detectors_Segs",Form("SegSummary_det%02d",hit.GetDetnum()),
    //		      34,0,34,hit.GetMainSegnum(),3000,0,36000,hit.Charge());

    obj.FillHistogram(dirname,"CoreEnergy_Summary",18,0,18,hit.GetDetnum(),
		      4000,0,4000,hit.GetEnergy());
    
    obj.FillHistogram(dirname,"CoreCharge_Summary",18,0,18,hit.GetDetnum(),5000,0,33000,hit.Charge());
    
    //obj.FillHistogram("SeGA_Detectors_Segs",Form("SegMult_det%02d",hit.GetDetnum()),33,0,33,num_segs);

    //obj.FillHistogram("SeGA_Detectors",Form("CoreEnergy_v_SegMult_det%02d",hit.GetDetnum()),
    //		      33,0,33,num_segs,
    //		      2000,0,4000,hit.GetEnergy());

    /*
    if(!num_segs) {
      obj.FillHistogram(dirname,"CoreEnergy_NoSegs",4000,0,4000,hit.GetEnergy());
      obj.FillHistogram("SeGA_Detectors_NoSegs",Form("GammaEnergy_det%02d",hit.GetDetnum()),
			4000,0,4000,hit.GetEnergy());
    }
    else {
      obj.FillHistogram(dirname,"CoreEnergy_YesSegs",4000,0,4000,hit.GetEnergy());
      obj.FillHistogram("SeGA_Detectors_YesSegs",Form("GammaEnergy_det%02d",hit.GetDetnum()),
			4000,0,4000,hit.GetEnergy());
    }
    */
    
    long cc_ts = hit.Timestamp();
    //double seg_sum_charge = 0;

    obj.FillHistogram("SeGA_Detectors_Time",Form("CoreEnergy_v_TimeStamp_det%02d",hit.GetDetnum()),
		      500,0,4000000000000,cc_ts,
		      2000,0,4000,hit.GetEnergy());

    obj.FillHistogram("SeGA_Detectors_Time",Form("CoreCharge_v_TimeStamp_det%02d",hit.GetDetnum()),
		      500,0,4000000000000,cc_ts,
		      1500,0,30000,hit.Charge());
    
    
    /*
    for(int j=0;j<num_segs;j++) {

      TSegaSegmentHit seg = hit.GetSegment(j);
      long seg_ts = seg.Timestamp();
      seg_sum_charge += seg.Charge();

      obj.FillHistogram(dirname,"Core_Seg_tdiff",2000,-20000,20000,seg_ts-cc_ts);

      obj.FillHistogram(dirname,"CoreEn_v_Core_Seg_tdiff",2000,-20000,20000,cc_ts-seg_ts,
			2000,0,4000,hit.GetEnergy());

    }//end segment loop
    
    if(hit.Charge() <= 0) {

      obj.FillHistogram("SeGA_Detectors_NoCoreCharge",Form("SegMult_det%02d",hit.GetDetnum()),
			33,0,33,num_segs);
      
      obj.FillHistogram("SeGA_Detectors_NoCoreCharge",Form("SegSumCharge_det%02d",hit.GetDetnum()),
			33000,0,33000,seg_sum_charge);
      
    }
    obj.FillHistogram("SeGA_Core_Seg",Form("SegSum_v_CoreCharge_det%02d",hit.GetDetnum()),
		      2048,0,32768,hit.Charge(),2048,0,32768,seg_sum_charge);
    */
   
  }//end hit loop

  /*
  if(index > -1 && sega.Size() > 1) {

    TSegaHit ref_hit = sega.GetSegaHit(index);
    auto tref = ref_hit.Timestamp();
    
    for(unsigned int i=0;i<sega.Size();i++) {

      if(i == index) {
	continue;
      }

      TSegaHit hit = sega.GetSegaHit(i);

      auto tdiff = hit.Timestamp() - tref;

      obj.FillHistogram(dirname,"CoreEnergy_v_tdiff",
			300,-3000,3000,tdiff,
			2500,0,5000,hit.GetEnergy());

    }
  }
  */

  return;
}

void MakeSeGAJanus(TRuntimeObjects& obj, TSega& sega, TJanusDDAS& janus) {

  int lvl = 1;
  std::string dirname = Prefix(lvl) + "SeGAJanus";
 
  obj.FillHistogram(dirname,"SegaMult",20,0,20,sega.Size());
  
  obj.FillHistogram(dirname,"All_jHitMult_det1",40,0,40,janus.GetSpecificHits(1,lvl,0).size());
  obj.FillHistogram(dirname,"All_jHitMult_det0",40,0,40,janus.GetSpecificHits(0,lvl,0).size());

  obj.FillHistogram(dirname,"Single_jHitMult_det1",50,0,50,janus.GetSpecificHits(1,lvl,1).size());   
  obj.FillHistogram(dirname,"Single_jHitMult_det0",50,0,50,janus.GetSpecificHits(0,lvl,1).size());

  //obj.FillHistogram(dirname,"Double_jHitMult_det1",50,0,50,janus.GetSpecificHits(1,lvl,2).size());   
  //obj.FillHistogram(dirname,"Double_jHitMult_det0",50,0,50,janus.GetSpecificHits(0,lvl,2).size());

  obj.FillHistogram(dirname,"Addback_jHitMult_det1",50,0,50,janus.GetSpecificHits(1,lvl,3).size());   
  obj.FillHistogram(dirname,"Addback_jHitMult_det0",50,0,50,janus.GetSpecificHits(0,lvl,3).size());
  
  for(auto& j_hit : janus.GetHits(lvl)) {

    BasicJanusSpectra(obj,j_hit,dirname,"All_");
     
    for(auto& s_hit : sega) {

      double s_energy = s_hit.GetEnergy();
      
      obj.FillHistogram(dirname,"Core_Energy",5000,0,5000,s_energy);
      obj.FillHistogram(dirname,"Core_Energy_Summary",18,0,18,s_hit.GetDetnum(),
			25000,0,5000,s_energy);

      auto tdiff = (s_hit.Timestamp() - j_hit.Timestamp());
      obj.FillHistogram(dirname,Form("Particle_Gamma_Tdiff_det%d",j_hit.GetDetnum()),
			5000,-11000,11000,tdiff);
      
      obj.FillHistogram(dirname,Form("CoreEnergy_v_tdiff_det%d",j_hit.GetDetnum()),
			1000,-1000,3000,tdiff,
			2500,0,5000,s_energy);
       
    }//end sega hit loop
  }//end janus hit loop

  return;
}

void MakeGatedSeGAJanus(TRuntimeObjects& obj, TSega& sega, TJanusDDAS& janus, GCutG* janus_gate,
			GCutG* time_gate) {

  //static TSRIM srimB("80Ge_on_196Pt");
  static TSRIM srimB("80Kr_on_196Pt");
  static TSRIM srimT("196Pt_on_196Pt");

  if(!janus_gate){
    return;
  }

  int lvl = 1;
  std::string name = janus_gate->GetName();
  std::string dirname = Prefix(lvl) + Form("SeGAJanus_%s",name.c_str());
  
  if(time_gate) {
    dirname += Form("_%s",time_gate->GetName());
  }

  std::string dirname1 = dirname + "_Sectors";
  std::string dirname2 = dirname + "_Rings";
  std::string dirname3 = dirname + "_sDets";

  std::string dirnameAS = dirname + "_AScan";
  std::string dirnameXS = dirname + "_XScan";
  std::string dirnameYS = dirname + "_YScan";
  std::string dirnameZS = dirname + "_ZScan";
  
  std::string dirnameY = Prefix(lvl) + Form("Yields_%s",name.c_str());
  
  bool recoil = false;
  if(!name.compare("196Pt")) {
    recoil = true;
  }

  int detNum = 1;
  std::string tag = janus_gate->GetTag();
  if(!tag.compare("D0")) {
    detNum = 0;
  }
  
  TVector3 incBeam = TVector3(0.0,0.0,1.0);
  int jMult = 0;
  
  for(auto& j_hit : janus.GetHits(lvl)) {

    if(j_hit.GetDetnum() != detNum || !janus_gate->IsInside(j_hit.GetRing(),j_hit.Charge())) {
      continue;
    }

    jMult++;
    //BasicJanusSpectra(obj,j_hit,Prefix(lvl) + Form("SeGAJanus_%s_JanSpec",name.c_str()),"All_");
    
    TVector3 j_pos = j_hit.GetPosition(false);
    double beta;
    TVector3 recon_pos(0,0,1);
    double recon_beta;
    if(recoil) {

      double theta_LAB = j_pos.Theta();

      double width = 0.369; //196Pt
      double reac_en = 0.001*srimB.GetAdjustedEnergy(1000.0*beam_en,width,0.05);
      
      double theta_CM = Theta_CM_FR(theta_LAB,reac_en);
      double En = Recoil_KE_LAB(theta_CM,reac_en);
      
      double distance = width/std::abs(std::cos(theta_LAB));
      
      En = 0.001*srimT.GetAdjustedEnergy(1000.0*En,distance,0.05);
      beta = Beta(En,targ_mass);

      double Recon_Theta_LAB = Theta_LAB(theta_CM,reac_en);
      double Recon_En = KE_LAB(theta_CM,reac_en);
      recon_pos.SetTheta(Recon_Theta_LAB);
      
      double recon_distance = width/std::abs(std::cos(Recon_Theta_LAB));
      Recon_En = 0.001*srimB.GetAdjustedEnergy(1000.0*Recon_En,recon_distance,0.05);
      
      recon_beta = Beta(Recon_En,beam_mass);  
	
    }
    else {

      double Theta_LAB = j_pos.Theta();

      double width = 0.369; //196Pt
      double reac_en = 0.001*srimB.GetAdjustedEnergy(1000.0*beam_en,width,0.05);
      
      double theta_CM = Theta_CM_FP(Theta_LAB,reac_en);
      double En = KE_LAB(theta_CM,reac_en);

      double distance = width/std::abs(std::cos(Theta_LAB));
	
      En = 0.001*srimB.GetAdjustedEnergy(1000.0*En,distance,0.05);
      beta = Beta(En,beam_mass);

      double Recon_Theta_LAB = Recoil_Theta_LAB(theta_CM,reac_en);
      double Recon_En = Recoil_KE_LAB(theta_CM,reac_en);
      recon_pos.SetTheta(Recon_Theta_LAB);
      
      double recon_distance = width/std::abs(std::cos(Recon_Theta_LAB));
      Recon_En = 0.001*srimT.GetAdjustedEnergy(1000.0*Recon_En,recon_distance,0.05);
      
      recon_beta = Beta(Recon_En,targ_mass);
	
    }
    recon_pos.SetPhi(j_pos.Phi() + TMath::Pi());

    const double theta = j_pos.Theta()*TMath::RadToDeg();
    const double theta_r = rand_gen->Uniform(theta-(wT.at(j_hit.GetRing()-1)/2.0),
					     theta+(wT.at(j_hit.GetRing()-1)/2.0));
    
    for(auto& s_hit : sega) {
      
      auto tdiff = (s_hit.Timestamp() - j_hit.Timestamp());
      double s_energy = s_hit.GetEnergy();

      if(time_gate) {
	if(!time_gate->IsInside(tdiff,s_energy)) {
	  continue;
	}
      }

      TVector3 s_pos = s_hit.GetPosition2(true,true);
      double dop_en = s_hit.GetDoppler2(beta,j_pos);

      int slice = s_hit.GetMainSlicenum();
      int quad = s_hit.GetMainQuadnum();

      /*
      obj.FillHistogram(dirname,Form("CoreEnergy_v_tdiff_jdet%d",detNum),
			1000,-1000,3000,tdiff,
			2500,0,5000,s_energy);
      */
      
      obj.FillHistogram(dirname,Form("ParticleCharge_v_CoreEnergy_jdet%d",detNum),
			2000,0,2000,s_energy,
			800,5000,20000,j_hit.Charge());
      
      obj.FillHistogram(dirname,Form("ParticleCharge_v_DopplerEnergy_jdet%d",detNum),
			2000,0,2000,dop_en,
			800,5000,20000,j_hit.Charge());

      obj.FillHistogram(dirname,Form("ParticleTheta_v_CoreEnergy_jdet%d",detNum),
			2000,0,2000,s_energy,
			120,10,70,theta_r);

      obj.FillHistogram(dirname,Form("ParticleTheta_v_DopplerEnergy_jdet%d",detNum),
			2000,0,2000,dop_en,
			120,10,70,theta_r);

      obj.FillHistogram(dirname,Form("CoreEnergy_v_Ring_jdet%d",detNum),
			24,1,25,j_hit.GetRing(),
			2000,0,2000,s_energy);
      
      obj.FillHistogram(dirname,Form("DopplerEnergy_v_Ring_jdet%d",detNum),
			24,1,25,j_hit.GetRing(),
			2000,0,2000,dop_en);
      
      obj.FillHistogram(dirname,"Core_Energy",4000,0,4000,s_energy);
      obj.FillHistogram(dirname,"Core_Energy_Summary",18,0,18,s_hit.GetDetnum(),
			25000,0,5000,s_energy);

      //obj.FillHistogram(dirname1,Form("Core_Energy_Sec%02d",j_hit.GetSector()),4000,0,4000,s_energy);
      obj.FillHistogram(dirname2,Form("Core_Energy_Ring%02d",j_hit.GetRing()),4000,0,4000,s_energy);
      obj.FillHistogram(dirname3,Form("Core_Energy_sDets%02d",s_hit.GetDetnum()),4000,0,4000,s_energy);
      
      obj.FillHistogram(dirname,"Doppler_Energy",4000,0,4000,dop_en);
      obj.FillHistogram(dirname,"Doppler_Energy_Summary",18,0,18,s_hit.GetDetnum(),4000,0,4000,dop_en);

      obj.FillHistogram(dirname,"Doppler_Energy_SegSum",516,0,516,s_hit.GetMapSegnum(),2000,0,4000,dop_en);

      std::string tmp_name = dirname + Form("_sDet%02d",s_hit.GetDetnum());
      if(slice >= 0) {
	obj.FillHistogram(tmp_name,Form("Doppler_Energy_Slice%02d",slice),4000,0,4000,dop_en);
      }
      if(quad >= 0) {
	obj.FillHistogram(tmp_name,Form("Doppler_Energy_Quad%02d",quad),4000,0,4000,dop_en);
      }
      
      obj.FillHistogram(dirname1,Form("Doppler_Energy_Sec%02d",j_hit.GetSector()),4000,0,4000,dop_en);
      obj.FillHistogram(dirname2,Form("Doppler_Energy_Ring%02d",j_hit.GetRing()),4000,0,4000,dop_en);
      obj.FillHistogram(dirname3,Form("Doppler_Energy_sDet%02d",s_hit.GetDetnum()),4000,0,4000,dop_en);

      //obj.FillHistogram(dirname,Form("Doppler_Energy_v_tdiff"),1000,-1000,3000,tdiff,2500,0,5000,dop_en);
      obj.FillHistogram(dirname,Form("Beta_v_Ring_jdet%d",detNum),26,0,26,j_hit.GetRing(),100,0,0.1,beta);
      
      double angle = j_pos.Angle(s_pos)*TMath::RadToDeg();

      obj.FillHistogram(dirname,"Theta_Correlation",4000,0,4000,s_energy,180,0,180,angle);
      obj.FillHistogram(dirname,"Theta_Correction",4000,0,4000,dop_en,180,0,180,angle);

      
      obj.FillHistogram(dirname2,Form("Theta_Correlation_R%02d",j_hit.GetRing()),
			4000,0,4000,s_energy,180,0,180,angle);
      obj.FillHistogram(dirname2,Form("Theta_Correction_R%02d",j_hit.GetRing()),
			4000,0,4000,dop_en,180,0,180,angle);

      obj.FillHistogram(dirname1,Form("Theta_Correlation_S%02d",j_hit.GetSector()),
			4000,0,4000,s_energy,180,0,180,angle);
      obj.FillHistogram(dirname1,Form("Theta_Correction_S%02d",j_hit.GetSector()),
			4000,0,4000,dop_en,180,0,180,angle);

      obj.FillHistogram(dirname3,Form("Theta_Correlation_sDet%02d",s_hit.GetDetnum()),
			4000,0,4000,s_energy,
			180,0,180,angle);
      obj.FillHistogram(dirname3,Form("Theta_Correction_sDet%02d",s_hit.GetDetnum()),
			4000,0,4000,dop_en,
			180,0,180,angle);

      if(slice >= 0) {
	obj.FillHistogram(tmp_name,Form("Theta_Correlation_Slice%02d",slice),
			  4000,0,4000,s_energy,
			  180,0,180,angle);
	obj.FillHistogram(tmp_name,Form("Theta_Correction_Slice%02d",slice),
			  4000,0,4000,dop_en,
			  180,0,180,angle);
      }

      if(quad >= 0) {
	obj.FillHistogram(tmp_name,Form("Theta_Correlation_Quad%02d",quad),
			  4000,0,4000,s_energy,
			  180,0,180,angle);
	obj.FillHistogram(tmp_name,Form("Theta_Correction_Quad%02d",quad),
			  4000,0,4000,dop_en,
			  180,0,180,angle);
      }
      
      TVector3 reacPlane = j_pos.Cross(incBeam);
      TVector3 detPlane = s_pos.Cross(incBeam);

      double reac_phi = reacPlane.Phi();
      double det_phi = detPlane.Phi();

      double planeAng = reac_phi - det_phi;
      if(planeAng < 0) {
	planeAng += TMath::TwoPi();
      }

      double thing = 65*TMath::Pi()/32.0;
      obj.FillHistogram(dirname,"Phi_Correlation",4000,0,4000,s_energy,32,0,thing,planeAng);
      obj.FillHistogram(dirname,"Phi_Correction",4000,0,4000,dop_en,32,0,thing,planeAng);

      
      obj.FillHistogram(dirname2,Form("Phi_Correlation_R%02d",j_hit.GetRing()),
			4000,0,4000,s_energy,32,0,thing,planeAng);
      obj.FillHistogram(dirname2,Form("Phi_Correction_R%02d",j_hit.GetRing()),
			4000,0,4000,dop_en,32,0,thing,planeAng);

      obj.FillHistogram(dirname1,Form("Phi_Correlation_S%02d",j_hit.GetSector()),
			4000,0,4000,s_energy,8,0,thing,planeAng);
      obj.FillHistogram(dirname1,Form("Phi_Correction_S%02d",j_hit.GetSector()),
			4000,0,4000,dop_en,8,0,thing,planeAng);
      
      obj.FillHistogram(dirname3,Form("Phi_Correlation_sDet%02d",s_hit.GetDetnum()),
			4000,0,4000,s_energy,
			32,0,thing,planeAng);
      obj.FillHistogram(dirname3,Form("Phi_Correction_sDet%02d",s_hit.GetDetnum()),
			4000,0,4000,dop_en,
			32,0,thing,planeAng);

      if(slice >= 0) {
	obj.FillHistogram(tmp_name,Form("Phi_Correlation_Slice%02d",slice),
			  4000,0,4000,s_energy,
			  32,0,thing,planeAng);
	obj.FillHistogram(tmp_name,Form("Phi_Correction_Slice%02d",slice),
			  4000,0,4000,dop_en,
			  32,0,thing,planeAng);
      }

      if(quad >= 0) {
	obj.FillHistogram(tmp_name,Form("Phi_Correlation_Quad%02d",quad),
			  4000,0,4000,s_energy,
			  32,0,thing,planeAng);
	obj.FillHistogram(tmp_name,Form("Phi_Correction_Quad%02d",quad),
			  4000,0,4000,dop_en,
			  32,0,thing,planeAng);    
      }

      /*
      for(int i=0;i<100;i++) {
	double tmp_xoff = -5.0 + 0.1*(i+0.5);
	obj.FillHistogram(dirnameXS,Form("sX_Scan_%02d",s_hit.GetDetnum()),100,-5.0,5.0,tmp_xoff,
			  2000,0,4000,s_hit.GetDoppler2(beta,j_pos,TVector3(tmp_xoff,0,0)));
      }

      for(int i=0;i<100;i++) {
	double tmp_yoff = -5.0 + 0.1*(i+0.5);
	obj.FillHistogram(dirnameYS,Form("sY_Scan_%02d",s_hit.GetDetnum()),100,-5.0,5.0,tmp_yoff,
			  2000,0,4000,s_hit.GetDoppler2(beta,j_pos,TVector3(0,tmp_yoff,0)));
      }

      for(int i=0;i<100;i++) {
	double tmp_zoff = -5.0 + 0.1*(i+0.5);
	obj.FillHistogram(dirnameZS,Form("sZ_Scan_%02d",s_hit.GetDetnum()),100,-5.0,5.0,tmp_zoff,
			  2000,0,4000,s_hit.GetDoppler2(beta,j_pos,TVector3(0,0,tmp_zoff)));
      } 

      for(int i=0;i<40;i++) {
	double tmp_ang = -20.0 + i + 0.5;
	obj.FillHistogram(dirnameAS,Form("sA_Scan_%02d",s_hit.GetDetnum()),40,-20.0,20.0,tmp_ang,
			  2000,0,4000,s_hit.GetDoppler2(beta,j_pos,TVector3(0,0,0),tmp_ang));
      }
      */
      
      /*
      for(int i=1;i<101;i++) {
	double tmp_xoff = -0.5 + 0.01*i;
	TVector3 tmp_pos = j_hit.GetPosition(false,true,TVector3(tmp_xoff,0,0));
	
	double theta_CM = Theta_CM_FP(tmp_pos.Theta());
	double En = KE_LAB(theta_CM);
	
	double distance = 0.369/std::abs(std::cos(tmp_pos.Theta()));
	En += 0.001*srimB.GetEnergyChange(1000.*En,distance);
      
	double tmp_beta = Beta(En,beam_mass);
	
	obj.FillHistogram(dirname,"jX_Scan",100,-0.5,0.51,tmp_xoff,
			  2000,0,4000,s_hit.GetDoppler(tmp_beta,tmp_pos));	
      }

      for(int i=1;i<101;i++) {
	double tmp_yoff = -0.5 + 0.01*i;
	TVector3 tmp_pos = j_hit.GetPosition(false,true,TVector3(0,tmp_yoff,0));
	
	double theta_CM = Theta_CM_FP(tmp_pos.Theta());
	double En = KE_LAB(theta_CM);
	
	double distance = 0.369/std::abs(std::cos(tmp_pos.Theta()));
	En += 0.001*srimB.GetEnergyChange(1000.*En,distance);
      
	double tmp_beta = Beta(En,beam_mass);
	
	obj.FillHistogram(dirname,"jY_Scan",100,-0.5,0.51,tmp_yoff,
			  2000,0,4000,s_hit.GetDoppler(tmp_beta,tmp_pos));	
      }
      
      for(int i=1;i<101;i++) {
	double tmp_zoff = -0.5 + 0.01*i;
	TVector3 tmp_pos = j_hit.GetPosition(false,true,TVector3(0,0,tmp_zoff));
	
	double theta_CM = Theta_CM_FP(tmp_pos.Theta());
	double En = KE_LAB(theta_CM);
	
	double distance = 0.369/std::abs(std::cos(tmp_pos.Theta()));
	En += 0.001*srimB.GetEnergyChange(1000.*En,distance);
      
	double tmp_beta = Beta(En,beam_mass);
	
	obj.FillHistogram(dirname,"jZ_Scan",100,-0.5,0.51,tmp_zoff,
			  2000,0,4000,s_hit.GetDoppler(tmp_beta,tmp_pos));	
      }

      for(int i=0;i<100;i++) {

	double tmp_Ep = 220.0 + i + 0.5;
	double theta_CM = Theta_CM_FP(j_pos.Theta(),tmp_Ep);
	double En = KE_LAB(theta_CM,tmp_Ep);
	
	double distance = 0.369/std::abs(std::cos(j_pos.Theta()));
	En += 0.001*srimB.GetEnergyChange(1000.*En,distance);
	
	double tmp_beta = Beta(En,beam_mass);
	
	obj.FillHistogram(dirname2,Form("Ep_Scan_R%02d",j_hit.GetRing()),100,220,320,tmp_Ep,
			  2000,0,2000,s_hit.GetDoppler(tmp_beta,j_pos));	
      }

      for(int i=0;i<30;i++) {
	
	double tmp_beta = 0.09 - 0.001*(i+0.5);
	
	obj.FillHistogram(dirname2,Form("Beta_Scan_R%02d",j_hit.GetRing()),30,0.06,0.09,tmp_beta,
			  2000,0,2000,s_hit.GetDoppler2(tmp_beta,j_pos));	
      }
      */
      
      obj.FillHistogram(dirname,Form("ReconBeta_v_Ring_jdet%d",detNum),26,0,26,j_hit.GetRing(),
      			100,0,0.1,recon_beta);
      
      double recon_en = s_hit.GetDoppler(recon_beta,recon_pos);
      obj.FillHistogram(dirname,"Recon_Energy",4000,0,4000,recon_en);
      obj.FillHistogram(dirname,"Recon_Energy_Summary",18,0,18,s_hit.GetDetnum(),4000,0,4000,recon_en);

      obj.FillHistogram(dirname1,Form("Recon_Energy_Sec%02d",j_hit.GetSector()),4000,0,4000,recon_en);
      obj.FillHistogram(dirname2,Form("Recon_Energy_Ring%02d",j_hit.GetRing()),4000,0,4000,recon_en);
      obj.FillHistogram(dirname3,Form("Recon_Energy_sDet%02d",s_hit.GetDetnum()),4000,0,4000,recon_en);

      if(slice >= 0) {
	obj.FillHistogram(tmp_name,Form("Recon_Energy_Slice%02d",slice),4000,0,4000,recon_en);
      }
      if(quad >= 0) {
	obj.FillHistogram(tmp_name,Form("Recon_Energy_Quad%02d",quad),4000,0,4000,recon_en);
      }
      
      obj.FillHistogram(dirname,Form("ParticleCharge_v_ReconEnergy_jdet%d",detNum),
			2000,0,2000,recon_en,
			900,0,20000,j_hit.Charge());

      obj.FillHistogram(dirname,Form("ParticleTheta_v_ReconEnergy_jdet%d",detNum),
			2000,0,2000,recon_en,
			120,10,70,theta_r);
      
      obj.FillHistogram(dirname,Form("ReconEnergy_v_Ring_jdet%d",detNum),
			24,1,25,j_hit.GetRing(),
			2000,0,2000,recon_en);

      double recon_angle = recon_pos.Angle(s_pos)*TMath::RadToDeg();
      obj.FillHistogram(dirname,"ReconTheta_Correlation",4000,0,4000,s_energy,180,0,180,recon_angle);
      obj.FillHistogram(dirname,"ReconTheta_Correction",4000,0,4000,recon_en,180,0,180,recon_angle);

      obj.FillHistogram(dirname2,Form("ReconTheta_Correlation_R%02d",j_hit.GetRing()),
			4000,0,4000,s_energy,180,0,180,recon_angle);
      obj.FillHistogram(dirname2,Form("ReconTheta_Correction_R%02d",j_hit.GetRing()),
			4000,0,4000,recon_en,180,0,180,recon_angle);

      obj.FillHistogram(dirname1,Form("ReconTheta_Correlation_S%02d",j_hit.GetSector()),
			4000,0,4000,s_energy,180,0,180,recon_angle);
      obj.FillHistogram(dirname1,Form("ReconTheta_Correction_S%02d",j_hit.GetSector()),
			4000,0,4000,recon_en,180,0,180,recon_angle);

      if(slice >= 0) {
	obj.FillHistogram(tmp_name,Form("ReconTheta_Correlation_Slice%02d",slice),
			  4000,0,4000,s_energy,
			  180,0,180,recon_angle);
	obj.FillHistogram(tmp_name,Form("ReconTheta_Correction_Slice%02d",slice),
			  4000,0,4000,recon_en,
			  180,0,180,recon_angle);
      }

      if(quad >= 0) {
	obj.FillHistogram(tmp_name,Form("ReconTheta_Correlation_Quad%02d",quad),
			  4000,0,4000,s_energy,
			  180,0,180,recon_angle);
	obj.FillHistogram(tmp_name,Form("ReconTheta_Correction_Quad%02d",quad),
			  4000,0,4000,recon_en,
			  180,0,180,recon_angle);
      }

      TVector3 reacPlane_rec = recon_pos.Cross(incBeam);
      double reac_phi_rec = reacPlane_rec.Phi();

      double planeAng_rec = reac_phi_rec - det_phi;
      if(planeAng_rec < 0) {
	planeAng_rec += TMath::TwoPi();
      }

      obj.FillHistogram(dirname,"ReconPhi_Correlation",4000,0,4000,s_energy,32,0,thing,planeAng_rec);
      obj.FillHistogram(dirname,"ReconPhi_Correction",4000,0,4000,recon_en,32,0,thing,planeAng_rec);

      obj.FillHistogram(dirname2,Form("ReconPhi_Correlation_R%02d",j_hit.GetRing()),
			4000,0,4000,s_energy,8,0,thing,planeAng_rec);
      obj.FillHistogram(dirname2,Form("ReconPhi_Correction_R%02d",j_hit.GetRing()),
			4000,0,4000,recon_en,8,0,thing,planeAng_rec);

      obj.FillHistogram(dirname1,Form("ReconPhi_Correlation_S%02d",j_hit.GetSector()),
			4000,0,4000,s_energy,8,0,thing,planeAng_rec);
      obj.FillHistogram(dirname1,Form("ReconPhi_Correction_S%02d",j_hit.GetSector()),
			4000,0,4000,recon_en,8,0,thing,planeAng_rec);

      if(slice >= 0) {
	obj.FillHistogram(tmp_name,Form("ReconPhi_Correlation_Slice%02d",slice),
			  4000,0,4000,s_energy,
			  32,0,thing,planeAng_rec);
	obj.FillHistogram(tmp_name,Form("ReconPhi_Correction_Slice%02d",slice),
			  4000,0,4000,recon_en,
			  32,0,thing,planeAng_rec);
      }

      if(quad >= 0) {
	obj.FillHistogram(tmp_name,Form("ReconPhi_Correlation_Quad%02d",quad),
			  4000,0,4000,s_energy,
			  32,0,thing,planeAng_rec);
	obj.FillHistogram(tmp_name,Form("ReconPhi_Correction_Quad%02d",quad),
			  4000,0,4000,recon_en,
			  32,0,thing,planeAng_rec);    
      }
      
      if(recoil) {
	if(theta_r < 43.0) {
	  obj.FillHistogram(dirnameY,"Doppler_Energy_Y1",4000,0,4000,dop_en);
	  obj.FillHistogram(dirnameY,"Recon_Energy_Y1",4000,0,4000,recon_en);
	}
	else if(theta_r < 60.0) {
	  obj.FillHistogram(dirnameY,"Doppler_Energy_Y2",4000,0,4000,dop_en);
	  obj.FillHistogram(dirnameY,"Recon_Energy_Y2",4000,0,4000,recon_en);
	}
      }
      else {
	if(theta_r < 34.0) {
	  obj.FillHistogram(dirnameY,"Doppler_Energy_Y1",4000,0,4000,dop_en);
	  obj.FillHistogram(dirnameY,"Recon_Energy_Y1",4000,0,4000,recon_en);
	}
	else if(theta_r < 42.0) {
	  obj.FillHistogram(dirnameY,"Doppler_Energy_Y2",4000,0,4000,dop_en);
	  obj.FillHistogram(dirnameY,"Recon_Energy_Y2",4000,0,4000,recon_en);
	}
	else if(theta_r < 60.0) {
	  obj.FillHistogram(dirnameY,"Doppler_Energy_Y3",4000,0,4000,dop_en);
	  obj.FillHistogram(dirnameY,"Recon_Energy_Y3",4000,0,4000,recon_en);
	}
	
      }
      
    }//end sega hit loop
    
  }//end janus hit loop

  obj.FillHistogram(Prefix(lvl) + Form("SeGAJanus_%s_JanSpec",name.c_str()),
		    Form("All_Gated_jMult_det%d",detNum),20,0,20,jMult);
  
  return;
}

int gates_loaded=0;
GCutG* time_gate=0;
std::vector<GCutG*> janus_gates;

// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
  
  TSega* sega = obj.GetDetector<TSega>();
  TJanusDDAS* janus = obj.GetDetector<TJanusDDAS>();

  if(!rand_gen) {
    rand_gen = new TRandom(50747227);
  }
  
  TList *gates = &(obj.GetGates());
  if(gates_loaded!=gates->GetSize()) {
    TIter iter(gates);
    while(TObject *obj = iter.Next()) {
      GCutG *gate = (GCutG*)obj;
      std::string tag = gate->GetTag();
      std::string name = gate->GetName();
      if(!tag.compare("time")) {
        time_gate = gate;
	gates_loaded++;
        std::cout << "Time Gate: " << name << std::endl;
      }
      else if(!tag.compare("D0") || !tag.compare("D1")) {
        janus_gates.push_back(gate);
	gates_loaded++;
        std::cout << "Janus Gate: " << tag << " " << name << std::endl;
      }
      else {
	std::cout << "Unknown Gate: Name = " << name << ", Tag = " << tag << std::endl;
      }
    }
  }

  /*  
  if(sega) {
    MakeSega(obj,*sega);
  }
  
  if(janus) {
    MakeJanus(obj,*janus);

    for(auto &gate : janus_gates) {
      MakeGatedJanus(obj,*janus,gate);
    }
  }
  */
  
  
  if(sega && janus) {

    //MakeSeGAJanus(obj,*sega,*janus);

    for(auto &gate : janus_gates) {
      MakeGatedSeGAJanus(obj,*sega,*janus,gate,time_gate);
    }
  }
  
  
  return;
}


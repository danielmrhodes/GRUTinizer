#include "TSegaHit.h"

#include <algorithm>
#include <iostream>

#include "TString.h"

#include "GCanvas.h"
#include "GValue.h"
#include "TSega.h"

#include <GH1D.h>

TSegaHit::TSegaHit() {
  Clear();
}

void TSegaHit::Copy(TObject& obj) const{
  TDetectorHit::Copy(obj);

  TSegaHit& sega = (TSegaHit&)obj;
  sega.fTrace = fTrace;
}

void TSegaHit::Draw(Option_t* opt) {
  TString option = opt;
  if(!gPad || option.Contains("new", TString::ECaseCompare::kIgnoreCase)) {
    new GCanvas;
  } else {
    gPad->Clear();
  }

  TVirtualPad* base_pad = gPad;

  bool draw_all = option.Contains("all", TString::ECaseCompare::kIgnoreCase);
  if(draw_all){
    base_pad->Divide(1, 2, 0.005, 0.005);
    base_pad->cd(2)->Divide(4,8,0,0);
    base_pad->cd(1);
  }

  DrawTrace(0);

  if(draw_all){
    for(int i=1; i<33; i++) {
      base_pad->cd(2)->cd(i);
      DrawTrace(i);
    }
  }
  base_pad->cd();
}

void TSegaHit::DrawTrace(int segnum) {
  std::vector<unsigned short>* trace = GetTrace(segnum);
  if(!trace){
    std::cout << "No segment trace found for segment " << segnum << std::endl;
    return;
  }

  GH1D hist("sega_wave", "", trace->size(), 0, 10*trace->size());
  hist.SetStats(false);

  if(segnum==0){
    hist.SetTitle(Form("SeGA Detector %d at %ld ns", GetDetnum(), Timestamp()));
    hist.GetXaxis()->SetTitle("Time (ns)");
    hist.GetYaxis()->SetTitle("ADC units");
  }

  for(size_t i=0; i<trace->size(); i++) {
    hist.SetBinContent(i+1,(*trace)[i]);
  }
  hist.DrawCopy();
}

void TSegaHit::Clear(Option_t *opt) {
  TDetectorHit::Clear(opt);
  fTrace.clear();
}

void TSegaHit::Print(Option_t *opt) const {
  std::cout << "TSegaHit:\n"
            << "\tChannel: " << GetChannel() << "\n"
            << "\tCharge: " << Charge() << "\n"
            << std::flush;
}

void TSegaHit::SetTrace(unsigned int trace_length, const unsigned short* trace) {
  if(!trace){
    fTrace.clear();
    return;
  }

  fTrace.clear();
  return;
  fTrace.reserve(trace_length);
  for(unsigned int i=0; i<trace_length; i++){
    fTrace.push_back(trace[i]);
  }
}

bool TSegaHit::HasCore() const {
  return fCharge != -1;
}

int TSegaHit::GetDetnum() const {
  TChannel* chan = TChannel::GetChannel(fAddress);
  int output = -1;
  if(chan && fAddress!=-1){
    output = chan->GetArrayPosition();
  } else if(fSegments.size()) {
    output = fSegments[0].GetDetnum();
  } else {
    // std::cout << "Unknown address: " << std::hex << fAddress << std::dec
    //           << std::endl;
    output = -1;
  }

  if(output == -1 && chan){
    // std::cout << "Chan with det=-1: " << chan->GetName() << std::endl;
    // std::cout << "address: " << fAddress << std::endl;
  }

  return output;
}

int TSegaHit::GetCrate() const {
  return (fAddress&0x00ff0000)>>16;
}

int TSegaHit::GetSlot() const {
  return (fAddress&0x0000ff00)>>8;
}

int TSegaHit::GetChannel() const {
  return (fAddress&0x000000ff)>>0;
}

TSegaSegmentHit& TSegaHit::MakeSegmentByAddress(unsigned int address){
  // for(auto& segment : fSegments){
  //   if(segment.Address() == address){
  //     return segment;
  //   }
  // }

  fSegments.emplace_back();
  TSegaSegmentHit& output = fSegments.back();
  output.SetAddress(address);
  return output;
}

std::vector<unsigned short>* TSegaHit::GetTrace(int segnum) {
  if(segnum == 0){
    return &fTrace;
  }
  for(auto& seg : fSegments) {
    if(seg.GetSegnum() == segnum) {
      return &seg.GetTrace();
    }
  }
  return NULL;
}

int TSegaHit::GetMainSegnum() const {
  int output = 0;
  double max_energy = -9e99;
  for(auto& segment : fSegments){
    if(segment.GetEnergy() > max_energy){
      output = segment.GetSegnum();
      max_energy = segment.GetEnergy();
    }
  }
  return output;
}

int TSegaHit::GetMainSlicenum() const {

  if(GetNumSegments()<1) {
    return std::sqrt(-1);
  }
  
  return int(GetMainSegnum() - 1)/int(4);
  
}

int TSegaHit::GetMainQuadnum() const {

  if(GetNumSegments()<1) {
    return std::sqrt(-1);
  }
  
  int segnum = GetMainSegnum();
  if(segnum < 5) {
    return segnum-1;
  }
  else if(segnum < 9) {
    return segnum-5;
  }
  else if(segnum < 13) {
    return segnum-9;
  }
  else if(segnum < 17) {
    return segnum-13;
  }
  else if(segnum < 21) {
    return segnum-17;
  }
  else if(segnum < 25) {
    return segnum-21;
  }
  else if(segnum < 29) {
    return segnum-25;
  }
  else {
    return segnum-29;
  }
  
}

//Mapped Numbers
int TSegaHit::GetMapSegnum() const {return TSega::MappedSegnum(GetDetnum(),GetMainSegnum());}
int TSegaHit::GetMapPairnum() const {return TSega::MappedPairnum(GetDetnum(),GetMainSegnum());}
int TSegaHit::GetMapSlicenum() const {return TSega::MappedSlicenum(GetDetnum(),GetMainSegnum());}

TVector3 TSegaHit::GetPosition(bool apply_array_offset, bool apply_det_offset, TVector3 offset) const {
  
  TVector3 array_pos = TSega::GetSegmentPosition(GetDetnum(),GetMainSegnum());
  if(!std::isnan(offset.X()) && !std::isnan(offset.Y()) && !std::isnan(offset.Z())) {
      array_pos += offset;
  }

  if(apply_array_offset) {

    double X_off = GValue::Value("Sega_X_offset");
    if(!std::isnan(X_off)) {
      array_pos.SetX(array_pos.X() + X_off);
    }

    double Y_off = GValue::Value("Sega_Y_offset");
    if(!std::isnan(Y_off)) {
      array_pos.SetY(array_pos.Y() + Y_off);
    }
    
    double Z_off = GValue::Value("Sega_Z_offset");
    if(!std::isnan(Z_off)) {
      array_pos.SetZ(array_pos.Z() + Z_off);
    }

  }

  if(apply_det_offset) {
    
    double X_off = GValue::Value(Form("Sega_X_offset_%02d",GetDetnum()));
    if(!std::isnan(X_off)) {
      array_pos.SetX(array_pos.X() + X_off);
    }

    double Y_off = GValue::Value(Form("Sega_Y_offset_%02d",GetDetnum()));
    if(!std::isnan(Y_off)) {
      array_pos.SetY(array_pos.Y() + Y_off);
    }
    
    double Z_off = GValue::Value(Form("Sega_Z_offset_%02d",GetDetnum()));
    if(!std::isnan(Z_off)) {
      array_pos.SetZ(array_pos.Z() + Z_off);
    }
    
  }
  
  return array_pos;
}

TVector3 TSegaHit::GetPosition2(bool apply_array_offset, bool apply_det_offset, const TVector3 offset,
				const double angle) const {

  double ang = 0.0; 
  if(apply_det_offset) {

    ang = GValue::Value(Form("Sega_Rotation_%02d",GetDetnum()));
    if(std::isnan(ang)) {
      ang = 0.0;
    }
    
  }
  
  TVector3 array_pos = TSega::GetSegmentPosition2(GetDetnum(),GetMainSegnum(),ang);
  if(!std::isnan(offset.X()) && !std::isnan(offset.Y()) && !std::isnan(offset.Z())) {
      array_pos += offset;
  }

  if(apply_array_offset) {

    double X_off = GValue::Value("Sega_X_offset");
    if(!std::isnan(X_off)) {
      array_pos.SetX(array_pos.X() + X_off);
    }

    double Y_off = GValue::Value("Sega_Y_offset");
    if(!std::isnan(Y_off)) {
      array_pos.SetY(array_pos.Y() + Y_off);
    }
    
    double Z_off = GValue::Value("Sega_Z_offset");
    if(!std::isnan(Z_off)) {
      array_pos.SetZ(array_pos.Z() + Z_off);
    }

  }
  
  if(apply_det_offset) {
    
    double X_off = GValue::Value(Form("Sega_X_offset_%02d",GetDetnum()));
    if(!std::isnan(X_off)) {
      array_pos.SetX(array_pos.X() + X_off);
    }

    double Y_off = GValue::Value(Form("Sega_Y_offset_%02d",GetDetnum()));
    if(!std::isnan(Y_off)) {
      array_pos.SetY(array_pos.Y() + Y_off);
    }
    
    double Z_off = GValue::Value(Form("Sega_Z_offset_%02d",GetDetnum()));
    if(!std::isnan(Z_off)) {
      array_pos.SetZ(array_pos.Z() + Z_off);
    }
    
  }

  array_pos.RotateY(angle*TMath::DegToRad());
  
  return array_pos;
}

Int_t TSegaHit::Charge() const {
  if(fCharge > 30000) {
    return fCharge - 32768;
  } else {
    return fCharge;
  }
}

double TSegaHit::GetDoppler(double beta,const TVector3& particle_vec, const TVector3& sega_offset) const {
  if(GetNumSegments()<1) {
    return std::sqrt(-1);
  }

  double gamma = 1/(sqrt(1-pow(beta,2)));
  TVector3 pos = GetPosition(true,true,sega_offset);
  double cos_angle = TMath::Cos(pos.Angle(particle_vec));
  double dc_en = GetEnergy()*gamma *(1 - beta*cos_angle);
  return dc_en;
}

double TSegaHit::GetDoppler2(double beta,const TVector3& particle_vec, const TVector3& sega_offset,
			     const double angle) const {
  if(GetNumSegments()<1) {
    return std::sqrt(-1);
  }

  double gamma = 1/(sqrt(1-pow(beta,2)));
  TVector3 pos = GetPosition2(true,true,sega_offset,angle);
  double cos_angle = TMath::Cos(pos.Angle(particle_vec));
  double dc_en = GetEnergy()*gamma *(1 - beta*cos_angle);
  return dc_en;
}

double TSegaHit::GetTraceHeight() const {
  if(fTrace.size() < 20){
    return std::sqrt(-1);
  }

  double low = 0;
  double high = 0;
  for(unsigned int i=0; i<10; i++){
    low += fTrace[i];
    high += fTrace[fTrace.size()-i-1];
  }

  return (high-low)/10;
}

double TSegaHit::GetTraceHeightDoppler(double beta,const TVector3& vec) const {
  if(GetNumSegments()<1) {
    return std::sqrt(-1);
  }

  double gamma = 1/(sqrt(1-pow(beta,2)));
  TVector3 pos = GetPosition();
  double cos_angle = TMath::Cos(pos.Angle(vec));
  double dc_en = GetTraceHeight()*gamma *(1 - beta*cos_angle);
  return dc_en;
}

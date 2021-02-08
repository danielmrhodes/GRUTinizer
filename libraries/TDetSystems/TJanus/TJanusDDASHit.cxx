#include "TJanusDDASHit.h"

#include "GValue.h"
#include "TJanusDDAS.h"
#include "TReaction.h"
#include "TSRIM.h"

#include "TRandom.h"

TJanusDDASHit::TJanusDDASHit(const TJanusDDASHit& hit) {
  hit.Copy(*this);
}

TJanusDDASHit& TJanusDDASHit::operator=(const TJanusDDASHit& hit) {
  hit.Copy(*this);
  return *this;
}

bool TJanusDDASHit::operator==(const TJanusDDASHit& rhs) {

  return Address()==rhs.Address() && Charge()==rhs.Charge() && Timestamp()==rhs.Timestamp();
}

void TJanusDDASHit::Copy(TObject& obj) const {
  TDetectorHit::Copy(obj);

  TJanusDDASHit& hit = (TJanusDDASHit&)obj;
  back_hit.Copy(hit.back_hit);
}

void TJanusDDASHit::Clear(Option_t* opt) {
  TDetectorHit::Clear(opt);
}

float TJanusDDASHit::RawCharge() const {
  if(fCharge > 30000) {
    return fCharge - 32768;
  } else {
    return fCharge;
  }
}

int TJanusDDASHit::GetFrontChannel() const {
  int address = Address();

  int crate_id = (address & 0x00ff0000) >> 16;
  int slot_id  = (address & 0x0000ff00) >> 8;
  int chan_id  = (address & 0x000000ff) >> 0;

  return (slot_id-2)*16 + chan_id;
}

int TJanusDDASHit::GetBackChannel() const {
  int address = back_hit.Address();

  int crate_id = (address & 0x00ff0000) >> 16;
  int slot_id  = (address & 0x0000ff00) >> 8;
  int chan_id  = (address & 0x000000ff) >> 0;

  return (slot_id-2)*16 + chan_id;
}

int TJanusDDASHit::GetDetnum() const {
  TChannel* chan = TChannel::GetChannel(fAddress);
  if(chan){
    return chan->GetArrayPosition();
  } else {
    return -1;
  }
}

int TJanusDDASHit::GetRing() const {
  TChannel* chan = TChannel::GetChannel(back_hit.Address());
  if(chan){
    return chan->GetSegment();
  } else {
    return 0;
  }
}

int TJanusDDASHit::GetSector() const {
  TChannel* chan = TChannel::GetChannel(fAddress);
  if(chan){
    return chan->GetSegment();
  } else {
    return 0;
  }
}

void TJanusDDASHit::Print(Option_t *opt) const {
  printf("Det %i JanusDDASHit[%03i,%03i]  Ring: %02i  Sector: %02i  Charge: %i   Energy: %.02f\n",
          GetDetnum(),GetFrontChannel(),GetBackChannel(),GetRing(),GetSector(),Charge(),GetEnergy());
}

bool TJanusDDASHit::IsRing() const {
  TChannel* chan = TChannel::GetChannel(Address());
  if(chan){
    return (*chan->GetArraySubposition() == 'F');
  } else {
    return false;
  }
}

bool TJanusDDASHit::IsSector() const {
  TChannel* chan = TChannel::GetChannel(Address());
  if(chan){
    return (*chan->GetArraySubposition() != 'F');
  } else {
    return false;
  }
}

TVector3 TJanusDDASHit::GetPosition(bool before, bool apply_array_offset, TVector3 offset) const {
  
  TVector3 output = TJanusDDAS::GetPosition(GetDetnum(), GetRing(), GetSector(), before);

  TVector3 off(0,0,0);
  if(!std::isnan(offset.X()) && !std::isnan(offset.Y()) && !std::isnan(offset.Z())) {
      off += offset;
  }
  
  if(apply_array_offset) {

    if(GetDetnum()) {
      double xoff = GValue::Value("JanusDS_X_offset");
      if(!std::isnan(xoff)) {
	off.SetX(off.X() + xoff);
      }

      double yoff = GValue::Value("JanusDS_Y_offset");
      if(!std::isnan(yoff)) {
	off.SetY(off.Y() + yoff);
      }

      double zoff = GValue::Value("JanusDS_Z_offset");
      if(!std::isnan(zoff)) {
	off.SetZ(off.Z() + zoff);
      }
    }
    else {
      double xoff = GValue::Value("JanusUS_X_offset");
      if(!std::isnan(xoff)) {
	off.SetX(off.X() + xoff);
      }

      double yoff = GValue::Value("JanusUS_Y_offset");
      if(!std::isnan(yoff)) {
	off.SetY(off.Y() + yoff);
      }

      double zoff = GValue::Value("JanusUS_Z_offset");
      if(!std::isnan(zoff)) {
	off.SetZ(off.Z() + zoff);
      }
    }
    
  }
    
  return output + off;
}

TVector3 TJanusDDASHit::GetReconPosition(TReaction& reac, int d_p, int r_p, bool s2, bool before,
					 bool apply_offset) const {

  TVector3 det_pos = GetPosition(before,apply_offset);

  return TJanusDDAS::GetReconPosition(det_pos.Theta(),det_pos.Phi(),reac,d_p,r_p,s2);

}

double TJanusDDASHit::GetLabSolidAngle() const { return TJanusDDAS::LabSolidAngle(GetDetnum(),GetRing()); }

double TJanusDDASHit::GetCmSolidAngle(TReaction& reac, int part, bool before) const {
  return TJanusDDAS::CmSolidAngle(GetDetnum(),GetRing(),reac,part,before);
}





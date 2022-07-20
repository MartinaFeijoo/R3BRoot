/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

// ----------------------------------------------------------------
// -----       R3BAnalysisTrackerFragment source file         -----
// -----       Created 13/06/2022 by M. Feijoo Fontan         -----
// ----------------------------------------------------------------

/*
 * This task should make the analysis of the fragments after the target
 */

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

#include "R3BAnalysisTrackerFragment.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BLosHitData.h"
#include "R3BMusicHitData.h"
#include "R3BTofdHitData.h"
#include "R3BFiberMAPMTHitData.h"
#include "R3BSofGladFieldPar.h"
#include "R3BSofTrackingData.h"
#include "R3BTGeoPar.h"

#include "TClonesArray.h"
#include "TMath.h"

const Double_t c = 29.9792458;

R3BAnalysisTrackerFragment::R3BAnalysisTrackerFragment()
    : R3BAnalysisTrackerFragment("AnalysisTrackerFragment", 1)
{
}

R3BAnalysisTrackerFragment::R3BAnalysisTrackerFragment(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fHeader(NULL)
    , fHitItemsMus(NULL)
    , fHitItemsLos(NULL)
    , fHitItemsFib10(NULL)
    , fHitItemsFib11(NULL)
    , fHitItemsFib12(NULL)
    , fHitItemsTofd(NULL)
    , fFieldCentre(0)
    , fEffLength(0.)
    , fBfield_Glad(0.)
    , fOnline(kFALSE)
{
}

R3BAnalysisTrackerFragment::~R3BAnalysisTrackerFragment()
{
    R3BLOG(DEBUG1, "");
    if (fTrackingDataCA)
        delete fTrackingDataCA;
}

void R3BAnalysisTrackerFragment::SetParContainers()
{
    R3BLOG(INFO, "");
    // Reading TrackerFragmentPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(FATAL, !rtdb, "FairRuntimeDb not found");

    fGladPar = (R3BSofGladFieldPar*)rtdb->getContainer("GladFieldPar");
    if (!fGladPar)
    {
        LOG(ERROR) << "R3BAnalysisTrackerFragment::SetParContainers() Couldn't get handle on GladFieldPar container";
    }

    fTargetGeoPar = (R3BTGeoPar*)rtdb->getContainer("TargetGeoPar");
    if (!fTargetGeoPar)
    {
        LOG(WARNING) << "R3BAnalysisTrackerFragment::SetParContainers() : Could not get access to TargetGeoPar container.";
        return;
    }
    else
        LOG(INFO) << "R3BAnalysisTrackerFragment::SetParContainers() : Container TargetGeoPar found.";

    fMusicGeoPar = (R3BTGeoPar*)rtdb->getContainer("MusicGeoPar");
    if (!fMusicGeoPar)
    {
        LOG(ERROR) << "R3BAnalysisTrackerFragment::SetParContainers() : Could not get access to MusicGeoPar container.";
        return;
    }
    else
        LOG(INFO) << "R3BAnalysisTrackerFragment::SetParContainers() : Container MusicGeoPar found.";

    fFib10GeoPar = (R3BTGeoPar*)rtdb->getContainer("Fi10GeoPar");
    if (!fFib10GeoPar)
    {
        LOG(ERROR) << "R3BAnalysisTrackerFragment::SetParContainers() : Could not get access to Fi10GeoPar container.";
        return;
    }
    else
        LOG(INFO) << "R3BAnalysisTrackerFragment::SetParContainers() : Container Fi10GeoPar found.";

    fFib11GeoPar = (R3BTGeoPar*)rtdb->getContainer("Fi11GeoPar");
    if (!fFib11GeoPar)
    {
        LOG(ERROR) << "R3BAnalysisTrackerFragment::SetParContainers() : Could not get access to Fi11GeoPar container.";
        return;
    }
    else
        LOG(INFO) << "R3BAnalysisTrackerFragment::SetParContainers() : Container Fi11GeoPar found.";

    fFib12GeoPar = (R3BTGeoPar*)rtdb->getContainer("Fi12GeoPar");
    if (!fFib12GeoPar)
    {
        LOG(ERROR) << "R3BAnalysisTrackerFragment::SetParContainers() : Could not get access to Fi12GeoPar container.";
        return;
    }
    else
        LOG(INFO) << "R3BAnalysisTrackerFragment::SetParContainers() : Container Fi12GeoPar found.";

    fTofDGeoPar = (R3BTGeoPar*)rtdb->getContainer("tofdGeoPar");
    if (!fTofDGeoPar)
    {
        LOG(ERROR) << "R3BAnalysisTrackerFragment::SetParContainers() : Could not get access to tofdGeoPar container.";
        return;
    }
    else
        LOG(INFO) << "R3BAnalysisTrackerFragment::SetParContainers() : Container tofdGeoPar found.";

    fMwpc0GeoPar = (R3BTGeoPar*)rtdb->getContainer("Mwpc0GeoPar");
    if (!fMwpc0GeoPar)
    {
        LOG(ERROR) << "R3BAnalysisTrackerFragment::SetParContainers() : Could not get access to Mwpc0GeoPar container.";
        return;
    }
    else
        LOG(INFO) << "R3BAnalysisTrackerFragment::SetParContainers() : Container Mwpc0GeoPar found.";

    return;
}

void R3BAnalysisTrackerFragment::SetParameter()
{
    //--- Parameter Container ---
    fFieldCentre = fGladPar->GetFieldCentre();
    fEffLength = fGladPar->GetEffectiveLength();
    fBfield_Glad = fGladPar->GetMagneticField();
    return;
}

InitStatus R3BAnalysisTrackerFragment::Init()
{
    R3BLOG(INFO, "");

    FairRootManager* mgr = FairRootManager::Instance();
    R3BLOG_IF(FATAL, NULL == mgr, "FairRootManager not found");

    fHeader = (R3BEventHeader*)mgr->GetObject("EventHeader.");

    // Get access to hit data of the MUSIC
    fHitItemsMus = (TClonesArray*)mgr->GetObject("MusicHitData");
    R3BLOG_IF(WARNING, !fHitItemsMus, "MusicHitData not found");

    // // Get access to hit data of the LOS
    // fHitItemsLos = (TClonesArray*)mgr->GetObject("LosHit");
    // R3BLOG_IF(WARNING, !fHitItemsLos, "LosHit not found");

    // Get access to hit data of Fiber10
    fHitItemsFib10 = (TClonesArray*)mgr->GetObject("Fi10Hit");
    R3BLOG_IF(WARNING, !fHitItemsFib10, "Fi10Hit not found");

    // Get access to hit data of Fiber11
    fHitItemsFib11 = (TClonesArray*)mgr->GetObject("Fi11Hit");
    R3BLOG_IF(WARNING, !fHitItemsFib11, "Fi11Hit not found");

    // Get access to hit data of Fiber12
    fHitItemsFib12 = (TClonesArray*)mgr->GetObject("Fi12Hit");
    R3BLOG_IF(WARNING, !fHitItemsFib12, "Fi12Hit not found");

    // Get access to hit data of the ToFD
    fHitItemsTofd = (TClonesArray*)mgr->GetObject("TofdHit");
    R3BLOG_IF(WARNING, !fHitItemsTofd, "TofdHit not found");

    // Output data
    fTrackingDataCA = (TClonesArray*)mgr->GetObject("TrackingData");
    if (fTrackingDataCA == NULL)
    {
        fTrackingDataCA = new TClonesArray("R3BSofTrackingData");
        mgr->Register("TrackingData", "Analysis Tracking", fTrackingDataCA, !fOnline);
    }

    SetParameter();

    return kSUCCESS;
}

InitStatus R3BAnalysisTrackerFragment::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

void R3BAnalysisTrackerFragment::Exec(Option_t* option)
{
  // Reset entries in output arrays, local arrays
  std::cout << "exec" << std::endl;
  Reset();
  Int_t nHitMusic = fHitItemsMus->GetEntries();
  Int_t nHitFib10 = fHitItemsFib10->GetEntries();
  Int_t nHitFib11 = fHitItemsFib11->GetEntries();
  Int_t nHitFib12 = fHitItemsFib12->GetEntries();
  Int_t nHitTofD = fHitItemsTofd->GetEntries();

  if (nHitMusic == 0 || nHitTofD == 0 || nHitFib10 == 0 || nHitFib11 == 0 || nHitFib12 == 0)
      return;

  TVector3 pos1[2];
  pos1[0].SetXYZ(-100., 0., 0.);
  pos1[1].SetXYZ(-100., 0., 0.);
  TVector3 pos2[2];
  pos2[0].SetXYZ(-100., 0., 0.);
  pos2[1].SetXYZ(-100., 0., 0.);
  TVector3 pos3[2];
  pos3[0].SetXYZ(-1000., 0., 0.);
  pos3[1].SetXYZ(-1000., 0., 0.);

  Double_t Length = 0.; Double_t Brho = 0.;

  Float_t tofd_x = 0.; Float_t tofd_y = 0.;
  Float_t tof = 0.;
  Float_t fib10_x = 0.; Float_t fib10_y = 0.;
  Float_t fib11_x = 0.; Float_t fib11_y = 0.;
  Float_t fib12_x = 0.; Float_t fib12_y = 0.;


  for (Int_t ihit = 0; ihit < nHitMusic; ihit++)
  {
      auto hit = (R3BMusicHitData*)fHitItemsMus->At(ihit);
      if (!hit)
          continue;
      music_ang = hit->GetTheta(); // mrad
      music_z = hit->GetZcharge();
  }

  for (Int_t ihit = 0; ihit < nHitTofD; ihit++)
  {
      auto hit = (R3BTofdHitData*)fHitItemsTofd->At(ihit);
      if (!hit)
          continue;
      tofd_x = hit->GetX(); tofd_y = hit->GetY();
      tof = hit->GetTof();
  }

  for (Int_t ihit = 0; ihit < nHitFib10; ihit++)
  {
      auto hit = (R3BFiberMAPMTHitData*)fHitItemsFib10->At(ihit);
      if (!hit)
          continue;
      fib10_x = hit->GetX(); fib10_y = hit->GetY();
  }
  for (Int_t ihit = 0; ihit < nHitFib11; ihit++)
  {
      auto hit = (R3BFiberMAPMTHitData*)fHitItemsFib11->At(ihit);
      if (!hit)
          continue;
      fib11_x = hit->GetX(); fib11_y = hit->GetY();
  }
  for (Int_t ihit = 0; ihit < nHitFib12; ihit++)
  {
      auto hit = (R3BFiberMAPMTHitData*)fHitItemsFib12->At(ihit);
      if (!hit)
          continue;
      fib12_x = hit->GetX(); fib12_y = hit->GetY();
  }

  Length = GetLength(pos1[1].X() / 10., pos2[1].X() / 10., pos3[1].X() / 10.);
  Brho = GetBrho(pos1[1].X() / 10., pos2[1].X() / 10., pos3[1].X() / 10.);

  Double_t v = Length / tof / c;
  Double_t gamma = 1. / sqrt(1. - v * v);

  AddData(music_z, Brho / v / gamma / 3.107, v, Length, Brho, 0);
}


Double_t R3BAnalysisTrackerFragment::GetBrho(Double_t position1, Double_t position2, Double_t position3)
{
    Double_t brho = 0.;
    Double_t L = fEffLength; // cm
    Double_t alpha = -14. * TMath::DegToRad();
    Double_t beta = -18. * TMath::DegToRad();
    Double_t theta = music_ang;

    Double_t xc = 1. / (1. + tan(alpha) * tan(theta)) * (position1 + (fFieldCentre - fMusicGeoPar->GetPosZ()) * tan(theta));
    Double_t zc = fFieldCentre - xc * tan(alpha);

    TVector3 posc = { xc, 0., zc };

    // std::cout<<"c: "<<xc<<" "<<zc<<std::endl;

    Double_t xb = xc - L / 2. * sin(theta) / cos(theta - alpha);
    Double_t zb = zc - L / 2. * cos(theta) / cos(theta - alpha);

    TVector3 posb = { xb, 0., zb };

    // std::cout<<"b: "<<xb<<" "<<zb<<std::endl;

    TVector3 ftrans = { fFib10GeoPar->GetPosX(), fFib10GeoPar->GetPosY(), fFib10GeoPar->GetPosZ() };
    TRotation frot;
    frot.RotateX(-1. * fFib10GeoPar->GetRotX() * TMath::DegToRad());
    frot.RotateY(-1. * fFib10GeoPar->GetRotY() * TMath::DegToRad());
    frot.RotateZ(-1. * fFib10GeoPar->GetRotZ() * TMath::DegToRad());

    TVector3 pos;
    pos.SetXYZ(position3, 0., 0.);

    auto pos2 = frot * pos + ftrans;

    // std::cout<<"E: "<<pos2.X()<<" "<<pos2.Y()<<" "<<pos2.Z()<<std::endl;

    auto pos3 = pos2 - posc;

    Double_t eta = atan(pos3.X() / pos3.Z());

    // std::cout<<eta<<std::endl;

    Double_t xd = xc + L / 2. * sin(eta) / cos(eta + alpha);
    Double_t zd = zc + L / 2. * cos(eta) / cos(eta + alpha);

    TVector3 posd = { xd, 0., zd };

    // std::cout<<"d: "<<xd<<" "<<zd<<std::endl;

    Double_t xf = pos2.X() + (fTofDGeoPar->GetPosZ() - fFib10GeoPar->GetPosZ()) * sin(eta) / cos(eta - beta);
    Double_t zf = pos2.Z() + (fTofDGeoPar->GetPosZ() - fFib10GeoPar->GetPosZ()) * cos(eta) / cos(eta - beta);

    TVector3 posf = { xf, 0., zf };

    // std::cout<<"f: "<<xf<<" "<<zf<<std::endl;

    Double_t rho = 0.5 * L / sin((theta - eta) / 2.0);
    brho = rho * fBfield_Glad;

    return brho / 100.; //[mT]
}

Double_t R3BAnalysisTrackerFragment::GetLength(Double_t position1, Double_t position2, Double_t position3)
{
    Double_t length = 0.;
    Double_t L = fEffLength; // cm
    Double_t alpha = -14. * TMath::DegToRad();
    Double_t beta = -20. * TMath::DegToRad();

    //Double_t theta = atan((pos2 - pos1) / (fMw2GeoPar->GetPosZ() - fMwpc0GeoPar->GetPosZ()));
    Double_t theta = music_ang;

    // std::cout<<mw1<<" "<<mw2<<" "<<fMwpc0GeoPar->GetPosZ()<<" "<<fMw2GeoPar->GetPosZ()<<std::endl;
    // std::cout<<alpha<<" "<<theta<<std::endl;

    Double_t xc = 1. / (1. + tan(alpha) * tan(theta)) * (position1 + (fFieldCentre - fMusicGeoPar->GetPosZ()) * tan(theta));
    Double_t zc = fFieldCentre - xc * tan(alpha);

    TVector3 posc = { xc, 0., zc };

    // std::cout<<"c: "<<xc<<" "<<zc<<std::endl;

    Double_t xb = xc - L / 2. * sin(theta) / cos(theta - alpha);
    Double_t zb = zc - L / 2. * cos(theta) / cos(theta - alpha);

    TVector3 posb = { xb, 0., zb };

    // std::cout<<"b: "<<xb<<" "<<zb<<std::endl;

    TVector3 ftrans = { fFib10GeoPar->GetPosX(), fFib10GeoPar->GetPosY(), fFib10GeoPar->GetPosZ() };
    TRotation frot;
    frot.RotateX(-1. * fFib10GeoPar->GetRotX() * TMath::DegToRad());
    frot.RotateY(-1. * fFib10GeoPar->GetRotY() * TMath::DegToRad());
    frot.RotateZ(-1. * fFib10GeoPar->GetRotZ() * TMath::DegToRad());
    TVector3 pos;
    pos.SetXYZ(position3, 0., 0.);

    auto pos2 = frot * pos + ftrans;

    // std::cout<<"E: "<<pos2.X()<<" "<<pos2.Y()<<" "<<pos2.Z()<<std::endl;

    auto pos3 = pos2 - posc;

    Double_t eta = atan(pos3.X() / pos3.Z());

    // std::cout<<eta<<std::endl;

    Double_t xd = xc + L / 2. * sin(eta) / cos(eta + alpha);
    Double_t zd = zc + L / 2. * cos(eta) / cos(eta + alpha);

    TVector3 posd = { xd, 0., zd };

    // std::cout<<"d: "<<xd<<" "<<zd<<std::endl;

    Double_t xf = pos2.X() + (fTofDGeoPar->GetPosZ() - fFib10GeoPar->GetPosZ()) * sin(eta) / cos(eta - beta);
    Double_t zf = pos2.Z() + (fTofDGeoPar->GetPosZ() - fFib10GeoPar->GetPosZ()) * cos(eta) / cos(eta - beta);

    TVector3 posf = { xf, 0., zf };

    // std::cout<<"f: "<<xf<<" "<<zf<<std::endl;

    // double rho = -0.5*L/sin((eta-theta)/2.0)/cos(alpha-(theta+eta)/2.0);

    Double_t rho = 0.5 * L / sin((theta - eta) / 2.0);

    // double omega = 2.0/sin(sqrt( (posd.Z()-posb.Z())*(posd.Z()-posb.Z()) + (posd.X()-posb.X())*(posd.X()-posb.X())
    // )/2.0/rho);

    Double_t omega =
        2.0 * asin(sqrt((posd.Z() - posb.Z()) * (posd.Z() - posb.Z()) + (posd.X() - posb.X()) * (posd.X() - posb.X())) /
                   2.0 / rho);

    length = (zb - fTargetGeoPar->GetPosZ()) / cos(theta) + omega * rho + (posf.Z() - posd.Z()) / cos(-eta) + 98.75;

    // std::cout<<mw2<<" "<<mw3<<" "<< rho <<" "<<omega <<" "<<omega*rho <<" "<<(posd-posb).Mag()<<" "<<
    // fTargetGeoPar->GetPosZ()<<" "<<length<<std::endl;

    return length; // [cm]
}

void R3BAnalysisTrackerFragment::Finish() {}
void R3BAnalysisTrackerFragment::FinishEvent() {}


void R3BAnalysisTrackerFragment::Reset()
{
    LOG(DEBUG) << "Clearing AnalysisTrackerFragments Structures";
    if (fHitItemsLos)
        fHitItemsLos->Clear();
    if (fHitItemsMus)
        fHitItemsMus->Clear();
    if (fHitItemsFib10)
        fHitItemsFib10->Clear();
    if (fHitItemsFib11)
        fHitItemsFib11->Clear();
    if (fHitItemsFib12)
        fHitItemsFib12->Clear();
    if (fHitItemsTofd)
        fHitItemsTofd->Clear();
    if (fTrackingDataCA)
        fTrackingDataCA->Clear();
}

R3BSofTrackingData* R3BAnalysisTrackerFragment::AddData(Double_t z,
                                                   Double_t aq,
                                                   Double_t beta,
                                                   Double_t length,
                                                   Double_t brho,
                                                   Int_t paddle)
{
    // It fills the R3BSofTrackingData
    TClonesArray& clref = *fTrackingDataCA;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BSofTrackingData(z, aq, beta, length, brho, paddle);
}

ClassImp(R3BAnalysisTrackerFragment);

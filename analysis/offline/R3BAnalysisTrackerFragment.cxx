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
#include "R3BFragmentData.h"
#include "R3BTGeoPar.h"
#include "R3BMusicPoint.h"

#include "TClonesArray.h"
#include "TMath.h"

const Double_t c = 29.9792458;

TF1 *a, *b;
Double_t finter(double *x, double*par) {
   return TMath::Abs(a->EvalPar(x,par) - b->EvalPar(x,par));}


R3BAnalysisTrackerFragment::R3BAnalysisTrackerFragment()
    : R3BAnalysisTrackerFragment("AnalysisTrackerFragment", 1)
{
}

R3BAnalysisTrackerFragment::R3BAnalysisTrackerFragment(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fHeader(NULL)
    , fHitItemsMus(NULL)
    , fPointItemsMus(NULL)
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
    }
    else
        LOG(INFO) << "R3BAnalysisTrackerFragment::SetParContainers() : Container TargetGeoPar found.";

    fMusicGeoPar = (R3BTGeoPar*)rtdb->getContainer("MusicGeoPar");
    if (!fMusicGeoPar)
    {
        LOG(ERROR) << "R3BAnalysisTrackerFragment::SetParContainers() : Could not get access to MusicGeoPar container.";
    }
    else
        LOG(INFO) << "R3BAnalysisTrackerFragment::SetParContainers() : Container MusicGeoPar found.";

    fFib10GeoPar = (R3BTGeoPar*)rtdb->getContainer("Fi10GeoPar");
    if (!fFib10GeoPar)
    {
        LOG(ERROR) << "R3BAnalysisTrackerFragment::SetParContainers() : Could not get access to Fi10GeoPar container.";
    }
    else
        LOG(INFO) << "R3BAnalysisTrackerFragment::SetParContainers() : Container Fi10GeoPar found.";

    fFib11GeoPar = (R3BTGeoPar*)rtdb->getContainer("Fi11GeoPar");
    if (!fFib11GeoPar)
    {
        LOG(ERROR) << "R3BAnalysisTrackerFragment::SetParContainers() : Could not get access to Fi11GeoPar container.";
    }
    else
        LOG(INFO) << "R3BAnalysisTrackerFragment::SetParContainers() : Container Fi11GeoPar found.";

    fFib12GeoPar = (R3BTGeoPar*)rtdb->getContainer("Fi12GeoPar");
    if (!fFib12GeoPar)
    {
        LOG(ERROR) << "R3BAnalysisTrackerFragment::SetParContainers() : Could not get access to Fi12GeoPar container.";
    }
    else
        LOG(INFO) << "R3BAnalysisTrackerFragment::SetParContainers() : Container Fi12GeoPar found.";

    fTofDGeoPar = (R3BTGeoPar*)rtdb->getContainer("tofdGeoPar");
    if (!fTofDGeoPar)
    {
        LOG(ERROR) << "R3BAnalysisTrackerFragment::SetParContainers() : Could not get access to tofdGeoPar container.";
    }
    else
        LOG(INFO) << "R3BAnalysisTrackerFragment::SetParContainers() : Container tofdGeoPar found.";

    // fMwpc0GeoPar = (R3BTGeoPar*)rtdb->getContainer("Mwpc0GeoPar");
    // if (!fMwpc0GeoPar)
    // {
    //     LOG(ERROR) << "R3BAnalysisTrackerFragment::SetParContainers() : Could not get access to Mwpc0GeoPar container.";
    // }
    // else
    //     LOG(INFO) << "R3BAnalysisTrackerFragment::SetParContainers() : Container Mwpc0GeoPar found.";

    return;
}

void R3BAnalysisTrackerFragment::SetParameter()
{
    //--- Parameter Container ---
    // fFieldCentre = fGladPar->GetFieldCentre();
    // fEffLength = fGladPar->GetEffectiveLength();
    // fBfield_Glad = fGladPar->GetMagneticField();

    fFieldCentre = 271.4;
    fEffLength =  230.;
    fBfield_Glad = 2.155;
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

    // fPointItemsMus = (TClonesArray*)mgr->GetObject("MusicPoint");
    // R3BLOG_IF(WARNING, !fPointItemsMus, "MusicPoint not found");

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
        fTrackingDataCA = new TClonesArray("R3BFragmentData");
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
  Reset();
  Int_t nHitMusic = fHitItemsMus->GetEntriesFast();
  Int_t nHitFib10 = fHitItemsFib10->GetEntriesFast();
  Int_t nHitFib11 = fHitItemsFib11->GetEntriesFast();
  Int_t nHitFib12 = fHitItemsFib12->GetEntriesFast();
  Int_t nHitTofD = fHitItemsTofd->GetEntriesFast();

  if (nHitMusic == 0 || nHitTofD == 0 || nHitFib10 == 0 || nHitFib11 == 0 || nHitFib12 == 0)
      return;

  Double_t Length = 0.; Double_t Brho = 0.;

  Float_t tofd_x = 0.; Float_t tofd_y = 0.;
  Float_t tof = 0.;
  Float_t fib10_x = 0.; Float_t fib10_y = 0.;
  Float_t fib11_x = 0.; Float_t fib11_y = 0.;
  Float_t fib12_x = 0.; Float_t fib12_y = 0.;
  Float_t x_mus = 0.;

  for (Int_t ihit = 0; ihit < nHitMusic; ihit++)
  {
      auto hit = (R3BMusicHitData*)fHitItemsMus->At(ihit);
      if (!hit)
          continue;
      music_ang = hit->GetTheta(); // mrad
      music_z = hit->GetZcharge();
      x_mus = hit->GetGoodDt();
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

  Length = GetLength(x_mus/10., fib10_x/10., fib11_x/10., music_ang);
  Brho = GetBrho(x_mus/10., fib10_x/10., fib11_x/10., music_ang);

  Double_t v = Length / tof / c;
  Double_t gamma = 1. / sqrt(1. - v*v);

  AddData(music_z, Brho / v / gamma / 3.107, v, Length, Brho);
}


Double_t R3BAnalysisTrackerFragment::GetBrho(Double_t position1, Double_t position2, Double_t position3, Double_t theta)
{
    Double_t brho = 0.;
    Double_t L = fEffLength; // cm
    Double_t alpha = -14. * TMath::DegToRad();
    Double_t beta = -18. * TMath::DegToRad();

    Double_t eta = atan((position3-position2)/9662.0) + beta;

    Double_t xc = 1. / (1.+tan(alpha)*tan(theta)) * (position1 + (fFieldCentre - 63.2) * tan(theta));
    Double_t zc = fFieldCentre - xc * tan(alpha);
    TVector3 posc = { xc, 0., zc };

    Double_t xb = xc - L/2. * sin(theta) / cos(theta - alpha);
    Double_t zb = zc - L/2. * cos(theta) / cos(theta - alpha);
    TVector3 posb = { xb, 0., zb };

    Double_t xd = xc + L/2. * sin(eta) / cos(eta - alpha);
    Double_t zd = zc + L/2. * cos(eta) / cos(eta - alpha);
    TVector3 posd = { xd, 0., zd };

    a=new TF1("a","pol1",-500.,1650.);
    b=new TF1("b","pol1",-500.,1650.);
    b->SetParameter(0,xd-tan(eta-TMath::Pi()/2)*zd); b->SetParameter(1,tan(eta-TMath::Pi()/2));
    Double_t zminimum_rho = 0.;
    if (theta != 0) {
      a->SetParameter(0,xb-tan(theta-TMath::Pi()/2)*zb); a->SetParameter(1,tan(theta-TMath::Pi()/2));
      TF1 *rho_inter = new TF1("rho_inter",finter,-300.,1650.);
      zminimum_rho = rho_inter->GetMinimumX();
      }
    if (theta == 0) {
      zminimum_rho = zb;
    }
    Double_t x_rho = b->Eval(zminimum_rho);

    //Double_t rho = abs(0.5*L/sin((eta-theta)/2.0)/cos(alpha-(theta+eta)/2.0));
    //Double_t rho = abs(0.5*L*sin((eta-theta)/2.0)/cos((eta-theta)/2.0))/(1-cos(alpha));
    Double_t rho = 0.5 * (sqrt(pow(x_rho-xb,2)+pow(zminimum_rho-zb,2)) + sqrt(pow(x_rho-xd,2)+pow(zminimum_rho-zd,2)));

    brho = rho * fBfield_Glad;

    return brho / 100.; //[mT]
}

//Double_t R3BAnalysisTrackerFragment::finter()

Double_t R3BAnalysisTrackerFragment::GetLength(Double_t position1, Double_t position2, Double_t position3, Double_t theta)
{
    Double_t length = 0.;
    Double_t L = fEffLength; // cm
    Double_t alpha = -14. * TMath::DegToRad();
    Double_t beta = -18. * TMath::DegToRad();

    Double_t xc = 1. / (1. + tan(alpha) * tan(theta)) * (position1 + (fFieldCentre - fMusicGeoPar->GetPosZ()) * tan(theta));
    Double_t zc = fFieldCentre - xc * tan(alpha);
    TVector3 posc = { xc, 0., zc };

    Double_t xb = xc - L / 2. * sin(theta) / cos(theta - alpha);
    Double_t zb = zc - L / 2. * cos(theta) / cos(theta - alpha);
    TVector3 posb = { xb, 0., zb };

    Double_t eta = atan((position3-position2)/9662.0) + beta;
    Double_t xd = xc + L / 2. * sin(eta) / cos(eta - alpha);
    Double_t zd = zc + L / 2. * cos(eta) / cos(eta - alpha);
    TVector3 posd = { xd, 0., zd };

    a=new TF1("a","pol1",-500.,1650.);
    b=new TF1("b","pol1",-500.,1650.);
    // TF1 *a=new TF1("a","pol1",-200.,600.);
    // TF1 *b=new TF1("b","pol1",-200.,600.);
    b->SetParameter(0,xd-tan(eta-TMath::Pi()/2)*zd); b->SetParameter(1,tan(eta-TMath::Pi()/2));
    Double_t zminimum_rho = 0.;
    if (theta != 0) {
      a->SetParameter(0,xb-tan(theta-TMath::Pi()/2)*zb); a->SetParameter(1,tan(theta-TMath::Pi()/2));
      TF1 *rho_inter = new TF1("rho_inter",finter,-300.,1650.);
      zminimum_rho = rho_inter->GetMinimumX();
      }
    if (theta == 0) {
      zminimum_rho = zb;
    }
    Double_t x_rho = b->Eval(zminimum_rho);


    TVector3 ftrans_fi11 = {-416.31, 0.0, 1552.66 };
    TVector3 ftrans_fi10 = {-117.74, 0.0, 633.75 };
    TRotation frot;
    frot.RotateY(-1.*18*TMath::DegToRad());
    TVector3 pos_fi11; TVector3 pos_fi10;
    pos_fi11.SetXYZ(position3/10., 0., 0.); pos_fi10.SetXYZ(position2/10., 0., 0.);
    auto pos_fi11lab = frot*pos_fi11+ftrans_fi11;
    auto pos_fi10lab = frot*pos_fi10+ftrans_fi10;

    TGraph *g_fi = new TGraph();
    g_fi->SetMarkerColor(kBlack); g_fi->SetMarkerStyle(8);
    g_fi->SetPoint(0,pos_fi10lab.Z(),pos_fi10lab.X());
    g_fi->SetPoint(1,pos_fi11lab.Z(),pos_fi11lab.X());

    a=new TF1("a","pol1",100.,1750.);
    b=new TF1("b","pol1",100.,1750.);
    g_fi->Fit("a","q+r0");
    b->SetParameter(0, fTofDGeoPar->GetPosX()-tan(TMath::Pi()/2+beta)*fTofDGeoPar->GetPosZ());
    b->SetParameter(1,tan(TMath::Pi()/2+beta));
    TF1 *f_intersection = new TF1("c",finter,200.,1650.);
    Double_t zf = f_intersection->GetMinimumX();
    Double_t xf = a->Eval(zf);
    TVector3 posf = { xf, 0., zf };

    //Double_t rho = abs(0.5*L/sin((eta-theta)/2.0)/cos(alpha-(theta+eta)/2.0));
    Double_t rho = 0.5 * (sqrt(pow(x_rho-xb,2)+pow(zminimum_rho-zb,2)) + sqrt(pow(x_rho-xd,2)+pow(zminimum_rho-zd,2)));

    Double_t omega =
        2.0*asin(sqrt((posd.Z()-posb.Z()) * (posd.Z()-posb.Z()) + (posd.X()-posb.X()) * (posd.X()-posb.X()))/2.0/rho);

    length = (zb - fTargetGeoPar->GetPosZ()) / cos(theta) + omega*rho + (posf.Z() - posd.Z()) / cos(eta);

    return length; // [cm]
}


void R3BAnalysisTrackerFragment::Finish() {}
void R3BAnalysisTrackerFragment::FinishEvent() {}


void R3BAnalysisTrackerFragment::Reset()
{
    LOG(DEBUG) << "Clearing AnalysisTrackerFragments Structures";
    // if (fHitItemsMus)
    //     fHitItemsMus->Clear();
    // if (fHitItemsFib10)
    //     fHitItemsFib10->Clear();
    // if (fHitItemsFib11)
    //     fHitItemsFib11->Clear();
    // if (fHitItemsFib12)
    //     fHitItemsFib12->Clear();
    // if (fHitItemsTofd)
    //     fHitItemsTofd->Clear();
    if (fTrackingDataCA)
        fTrackingDataCA->Clear();
}

R3BFragmentData* R3BAnalysisTrackerFragment::AddData(Double_t z,
                                                   Double_t aq,
                                                   Double_t beta,
                                                   Double_t length,
                                                   Double_t brho)
{
    // It fills the R3BFragmentData
    TClonesArray& clref = *fTrackingDataCA;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BFragmentData(z, aq, beta, length, brho);
}

ClassImp(R3BAnalysisTrackerFragment);

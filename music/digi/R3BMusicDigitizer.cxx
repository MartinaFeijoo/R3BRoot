/******************************************************************************
 *   Copyright (C) 2021 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2021-2023 Members of R3B Collaboration                     *
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
// -----            R3BMusicDigitizer source file             -----
// -----          Created 18/10/21 by JL Rodriguez            -----
// ----------------------------------------------------------------

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "TClonesArray.h"

#include "TMath.h"
#include "TRandom.h"
#include <iostream>
#include <string>

#include "R3BLogger.h"
#include "R3BMCTrack.h"
#include "R3BMusicHitData.h"
#include "R3BMusicDigitizer.h"
#include "R3BMusicPoint.h"

// R3BMusicDigitizer: Default Constructor --------------------------
R3BMusicDigitizer::R3BMusicDigitizer()
    : FairTask("R3BMusic Digitizer", 1)
    , fName("Music")
    , fMCTrack(NULL)
    , fMusicPoints(NULL)
    , fMusicHit(NULL)
    , fsigma_x(0.003) // sigma=0.03mm
{
}

// R3BMusicDigitizer: Standard Constructor --------------------------
R3BMusicDigitizer::R3BMusicDigitizer(const TString& name, Int_t iVerbose)
    : FairTask("R3B" + name + "Digitizer", iVerbose)
    , fName(name)
    , fMCTrack(NULL)
    , fMusicPoints(NULL)
    , fMusicHit(NULL)
    , fsigma_x(0.003)
{
}

// Virtual R3BMusicDigitizer: Destructor ----------------------------
R3BMusicDigitizer::~R3BMusicDigitizer()
{
    LOG(info) << "R3B" + fName + "Digitizer: Delete instance";
    if (fMusicPoints)
        delete fMusicPoints;
    if (fMusicHit)
        delete fMusicHit;
}

// ----   Public method Init  -----------------------------------------
InitStatus R3BMusicDigitizer::Init()
{
    R3BLOG(info, "for " << fName);

    // Get input array
    FairRootManager* ioman = FairRootManager::Instance();
    R3BLOG_IF(fatal, !ioman, "FairRootManager not found.");

    fMCTrack = dynamic_cast<TClonesArray*>(ioman->GetObject("MCTrack"));
    R3BLOG_IF(fatal, !fMCTrack, "MCTrack not found.");

    fMusicPoints = dynamic_cast<TClonesArray*>(ioman->GetObject(fName + "Point"));
    R3BLOG_IF(fatal, !fMusicPoints, fName << "Point not found.");

    // Register output array fMusicHit
    fMusicHit = new TClonesArray("R3BMusicHitData", 10);
    ioman->Register(fName + "HitData", "Digital response in " + fName, fMusicHit, kTRUE);

    return kSUCCESS;
}

// -----   Public method Execution   --------------------------------------------
void R3BMusicDigitizer::Exec(Option_t*)
{
    Reset();
    // Reading the Input -- Point Data --
    int nHits = fMusicPoints->GetEntriesFast();
    if (nHits == 0)
    {
        return;
    // Data from Point level
    R3BMusicPoint** pointData;
    pointData = new R3BMusicPoint*[nHits];
    Int_t TrackId = 0, PID = 0, anodeId = 0;
    Double_t x[nHits], fZ_out[nHits];
    Double_t Eave = 0.; Double_t z = 0.; Double_t theta = 0.;
    Double_t x_save = 0.;
    for (Int_t i = 0; i < nHits; i++)
    {
        x[i] = 0.; fZ_out[i] = 0.;
    }

    for (int i = 0; i < nHits; i++)
    {
        // Data from Point
        auto pointData = dynamic_cast<R3BMusicPoint*>(fMusicPoints->At(i));
        TrackId = pointData->GetTrackID();

        auto Track = dynamic_cast<R3BMCTrack*>(fMCTrack->At(TrackId));
        PID = Track->GetPdgCode();
        anodeId = pointData[i]->GetDetCopyID();

        Double_t fX_in = pointData[i]->GetXIn();
        Double_t fX_out = pointData[i]->GetXOut();
        x[i] = (fX_out + fX_in) / 2. ;
        if (x_save == 0 && anodeId == 3)
          x_save = x[anodeId];
        z += pointData[i]->GetZFF();
        Eave += pointData[i]->GetEnergyLoss();
        fZ_out[i] = pointData[i]->GetZOut();

        //std::cout << i << " " << x[i] << " " << x_save << " z=" << fZ_out[i] << std::endl;
    }

    if (PID > 1000200400) // Z=20 and A=40
    {
      theta = atan((x[nHits-1]-x[0])/(fZ_out[nHits-1]-fZ_out[0]));
      z = z/nHits;

      //std::cout << theta << " " << z << " " << Eave << " " << x[3] << std::endl;
      if (Eave > 0.01)
          AddHitData(theta, z, Eave, x_save);
    }

    if (pointData)
        delete pointData;

    return;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BMusicDigitizer::ReInit() { return kSUCCESS; }

// -----   Public method Reset   -----------------------------------------------
void R3BMusicDigitizer::Reset()
{
    LOG(DEBUG) << "Clearing R3B" + fName + "Digitizer Structure";
    if (fMusicHit)
        fMusicHit->Clear();
}

// -----   Private method AddR3BHitData  -------------------------------------------
R3BMusicHitData* R3BMusicDigitizer::AddHitData(Double_t theta, Double_t z, Double_t ene, Double_t good_dt)
{
    // It fills the R3BMusicHitData
    TClonesArray& clref = *fMusicHit;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BMusicHitData(theta, z, ene, good_dt);
}

ClassImp(R3BMusicDigitizer);

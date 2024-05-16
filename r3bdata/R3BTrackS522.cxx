/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2023 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

// -----------------------------------------------------------------------------
// -----                              R3BTrack                             -----
// -----                 Created on 12.03.2023 by A.Lagni for S522  -----
// -----------------------------------------------------------------------------

#include "R3BTrackS522.h"

R3BTrackS522::R3BTrackS522()
    : fX(0.)
    , fY(0.)
    , fZ(0.)
    , fPx(0.)
    , fPy(0.)
    , fPz(0.)
    , fQ(0.)
    , fQ_out(0.)
    , fAoZ(0.)
    , fChix(0.)
    , fToFm(0.)
    , fFlightPathm(0.)
    , fAng(0.,0.,0.)
    , fOutAng(0.,0.,0.)
    , fFibId(0)
    , fChiy(0.)
    , fQuality(0.)
    , fFootAngle(0.,0.,0.)
    , fFoot1Pos(0.,0.,0.)
    , fFoot2Pos(0.,0.,0.)
    , fFoot15Pos(0.,0.,0.)
    , fFoot16Pos(0.,0.,0.)
{
}

R3BTrackS522::R3BTrackS522(Double_t x,
                   Double_t y,
                   Double_t z,
                   Double_t px,
                   Double_t py,
                   Double_t pz,
                   Double_t q,
                   Double_t q_out,
                   Double_t AoZ,
                   Double_t ToFm,
                   Double_t FlightPathm,
                   TVector3 Ang,
                   TVector3 OutAng,
                   Int_t FibId,
                   Double_t chix,
                   Double_t chiy,
                   Int_t quality,
                   TVector3 FootAngle,
                   TVector3 Foot1Pos,
                   TVector3 Foot2Pos,
                   TVector3 Foot15Pos,
                   TVector3 Foot16Pos
		   )
    : fX(x)
    , fY(y)
    , fZ(z)
    , fPx(px)
    , fPy(py)
    , fPz(pz)
    , fQ(q)
    , fQ_out(q_out)
    , fAoZ(AoZ)
    , fToFm(ToFm)
    , fFlightPathm(FlightPathm)
    , fAng(Ang)
    , fOutAng(OutAng)
    , fFibId(FibId)
    , fChix(chix)
    , fChiy(chiy)
    , fQuality(quality)
    , fFootAngle(FootAngle)
    , fFoot1Pos(Foot1Pos)
    , fFoot2Pos(Foot2Pos)
    , fFoot15Pos(Foot15Pos)
    , fFoot16Pos(Foot16Pos)
{
}

R3BTrackS522::~R3BTrackS522() {}

ClassImp(R3BTrackS522)

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
// -----                              R3BTrackS522                             -----
// -----                 Created on 22.03.2023 by A.Lagni for S522     -----
// -----------------------------------------------------------------------------

#ifndef R3BTrackS522_H
#define R3BTrackS522_H

#include "TObject.h"
#include "TVector3.h"

class R3BTrackS522 : public TObject
{
  public:
    R3BTrackS522();
    R3BTrackS522(Double_t x,
             Double_t y,
             Double_t z,
             Double_t px,
             Double_t py,
             Double_t pz,
             Double_t q,
             Double_t q_out,
             Double_t AoZ,
             Double_t ToFm,
             Double_t Flight_Pathm,
             TVector3 Ang,
             TVector3 OutAng,
             Int_t FibId,
	     Double_t chix,
             Double_t chiy,
             Int_t quality);
    virtual ~R3BTrackS522();

    inline const Double_t& GetX() const { return fX; }
    inline const Double_t& GetY() const { return fY; }
    inline const Double_t& GetZ() const { return fZ; }
    inline const Double_t& GetPx() const { return fPx; }
    inline const Double_t& GetPy() const { return fPy; }
    inline const Double_t& GetPz() const { return fPz; }
    inline const Double_t& GetToF() const { return fToFm; }
    inline const Double_t& GetFlightPath() const { return fFlightPathm; }
    inline const Double_t& GetQ() const { return fQ; }
    inline const Double_t& GetQout() const { return fQ_out; }
    inline const Double_t& GetAoZ() const { return fAoZ; }
    inline const TVector3 GetAngles() const { return fAng; }
    inline const TVector3 GetOutAngles() const { return fOutAng; }
    inline const Int_t& GetFibId() const { return fFibId; }
    inline const Double_t& GetChix() const { return fChix; }
    inline const Double_t& GetChiy() const { return fChiy; }
    inline const Int_t& GetQuality() const { return fQuality; }

  protected:
    Double_t fX;
    Double_t fY;
    Double_t fZ;
    Double_t fPx;
    Double_t fPy;
    Double_t fPz;
    Double_t fQ;
    Double_t fQ_out;
    Double_t fAoZ;
    Double_t fToFm;
    Double_t fFlightPathm;
    TVector3 fAng;
    TVector3 fOutAng;
    Int_t fFibId;
    Double_t fChix;
    Double_t fChiy;
    Int_t fQuality;

  public:
    ClassDef(R3BTrackS522, 1)
};

#endif

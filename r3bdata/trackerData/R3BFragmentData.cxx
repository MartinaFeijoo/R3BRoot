// ---------------------------------------------------------------------------
// -----                                                                 -----
// -----                        R3BFragmentData                          -----
// -----                Created 13/10/2022 by M. Feijoo Fontán           -----
// -----                                                                 -----
// ---------------------------------------------------------------------------

#include "R3BFragmentData.h"

R3BFragmentData::R3BFragmentData()
    : fZ(0.)
    , fAq(0.)
    , fBeta(0.)
    , fLength(0.)
    , fBrho(0.)
{
}

//------------------------------

R3BFragmentData::R3BFragmentData(Double_t z,
                                       Double_t aq,
                                       Double_t beta,
                                       Double_t length,
                                       Double_t brho)
    : fZ(z)
    , fAq(aq)
    , fBeta(beta)
    , fLength(length)
    , fBrho(brho)
{
}

ClassImp(R3BFragmentData)

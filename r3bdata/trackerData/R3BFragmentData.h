// ---------------------------------------------------------------------------
// -----                                                                 -----
// -----                        R3BFragmentData                          -----
// -----                Created 13/10/2022 by M. Feijoo Fontán           -----
// -----                                                                 -----
// ---------------------------------------------------------------------------

#ifndef R3BFragmentData_H
#define R3BFragmentData_H
#include "TObject.h"

class R3BFragmentData : public TObject
{

  public:
    // Default Constructor
    R3BFragmentData();

    /** Standard Constructor
     *@param fZ      Z of fragments
     *@param fAq     A/q of fragments
     *@param fBeta   Beta of fragments
     *@param fLength Path length of fragments
     *@param fBrho   Brho of fragments
     **/

    R3BFragmentData(Double_t z, Double_t aq, Double_t beta, Double_t length, Double_t brho);

    // Destructor
    virtual ~R3BFragmentData() {}

    // Getters
    inline const Double_t GetZ() const { return fZ; }
    inline const Double_t GetAq() const { return fAq; }
    inline const Double_t GetBeta() const { return fBeta; }
    inline const Double_t GetBrho() const { return fBrho; }
    inline const Double_t GetLength() const { return fLength; }

  protected:
    Double_t fZ, fAq; // ID
    Double_t fBeta, fBrho, fLength;

  public:
    ClassDef(R3BFragmentData, 1)
};

#endif

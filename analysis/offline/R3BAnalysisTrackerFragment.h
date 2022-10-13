// ----------------------------------------------------------------------
// -----          R3BAnalysisTrackerFragment source file            -----
// -----          Created 13/06/2022 by M. Feijoo Fontan            -----
// ----------------------------------------------------------------------

#ifndef R3BAnalysisTrackerFragment_H
#define R3BAnalysisTrackerFragment_H

// ROOT headers
#include "TArrayF.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TH1F.h"
#include "TMath.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TRotation.h"
#include "TVector3.h"
#include <TRandom.h>
#include <iomanip>

// Fair headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "FairTask.h"

#include "R3BFragmentData.h"


class TClonesArray;
class R3BTGeoPar;
class R3BSofGladFieldPar;
class R3BEventHeader;

class R3BAnalysisTrackerFragment : public FairTask
{

  public:
    /** Default constructor **/
    R3BAnalysisTrackerFragment();

    /** Standard constructor **/
    R3BAnalysisTrackerFragment(const TString& name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BAnalysisTrackerFragment();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* option);

    /** Virtual method Reset **/
    virtual void Reset();
    /**
     * A method for finish of processing of an event.
     * Is called by the framework for each event after executing
     * the tasks.
     */
    virtual void FinishEvent();

    virtual void SetParContainers();
    //virtual void SetParamater();

    // Fair specific
    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /** Virtual method Finish **/
    virtual void Finish();

    /** Accessor to select online mode **/
    void SetOnline(Bool_t option) { fOnline = option; }

  private:
    void SetParameter();
    Double_t GetLength(Double_t pos1, Double_t pos2, Double_t pos3, Double_t theta);
    Double_t GetBrho(Double_t pos1, Double_t pos2, Double_t pos3, Double_t theta);
    Double_t GetVelocity(Double_t len, Double_t tof);
    Double_t GetAoverq(Double_t brho, Double_t vel);



    Bool_t fOnline; // Don't store data for online

    // Parameters set with accessor functions
    Double_t fFieldCentre, fEffLength, fBfield_Glad;

    R3BSofGladFieldPar* fGladPar;
    R3BTGeoPar* fTargetGeoPar;
    R3BTGeoPar* fTofDGeoPar;
    R3BTGeoPar* fMusicGeoPar;
    R3BTGeoPar* fFib10GeoPar;
    R3BTGeoPar* fFib11GeoPar;
    R3BTGeoPar* fFib12GeoPar;
    R3BTGeoPar* fpospc0GeoPar;

    TClonesArray* fHitItemsMus;
    TClonesArray* fPointItemsMus;
    TClonesArray* fHitItemsFib10;
    TClonesArray* fHitItemsFib11;
    TClonesArray* fHitItemsFib12;
    TClonesArray* fHitItemsLos;
    TClonesArray* fHitItemsTofd;  /**< Array with ToF Hit-input data. >*/
    TClonesArray* fTrackingDataCA; /**< Array with Tracking-output data. >*/
    R3BEventHeader* fHeader; // Event header

    Float_t music_ang = 0.; Float_t music_z = 0.;
    // Private method TrackingData
    R3BFragmentData* AddData(Double_t z, Double_t aq, Double_t beta, Double_t length, Double_t brho);
    



  public:
    // Class definition
    ClassDef(R3BAnalysisTrackerFragment, 1)
};

#endif /* R3BAnalysisTrackerFragment_H */

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

// ------------------------------------------------------------------
// -----                  R3BLosCal2Hit                         -----
// -----            Created March 10th 2016 by R.Plag           -----
// ----- Convert time calibrated data to hit level (single time) ----
// ------------------------------------------------------------------

#ifndef R3BLOSCAL2HIT
#define R3BLOSCAL2HIT

#include <map>

#include "FairTask.h"

class TClonesArray;
class TH1F;
class TH2F;
class R3BEventHeader;

/**
 * TODO: This explanation is humbug.
 * An analysis task to apply TCAL calibration for NeuLAND.
 * This class reads NeuLAND mapped items with TDC values and
 * produces time items with time in [ns]. It requires TCAL
 * calibration parameters, which are produced in a separate
 * analysis run containing R3BLosCal2HitFill task.
 */
class R3BLosCal2Hit : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BLosCal2Hit();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BLosCal2Hit(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BLosCal2Hit();

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

    /**
     * Method for re-initialization of parameter containers
     * in case the Run ID has changed.
     */
    virtual InitStatus ReInit();

    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    virtual void Exec(Option_t* option);

    /**
     * A method for finish of processing of an event.
     * Is called by the framework for each event after executing
     * the tasks.
     */
    virtual void FinishEvent();

    /**
     * A method to create histograms
     **/
    void CreateHisto();

    /**
     * A method to fill histograms
     **/
    void FillHisto();

    /**
     * Methods for setting position offset and effective velocity of light
     */
    inline void SetLosParamMCFD(Double_t offsetX, Double_t offsetY, Double_t veffX, Double_t veffY)
    {
        flosOffsetX = offsetX;
        flosOffsetY = offsetY;
        flosVeffX = veffX;
        flosVeffY = veffY;
    }

    inline void SetLosParamToT(Double_t offsetXQ, Double_t offsetYQ, Double_t veffXQ, Double_t veffYQ)
    {
        flosOffsetXQ = offsetXQ;
        flosOffsetYQ = offsetYQ;
        flosVeffXQ = veffXQ;
        flosVeffYQ = veffYQ;
    }

    inline void SetLosParamTAMEX(Double_t offsetXT, Double_t offsetYT, Double_t veffXT, Double_t veffYT)
    {
        flosOffsetXT = offsetXT;
        flosOffsetYT = offsetYT;
        flosVeffXT = veffXT;
        flosVeffYT = veffYT;
    }

    /**
     * Methods for setting input files
     */
    inline void SetLosInput(Int_t iOptHisto, std::string const& walk_param_file, std::string const& tot_param_file)
    {
        OptHisto = iOptHisto;
        fwalk_param_file = walk_param_file;
        ftot_param_file = tot_param_file;
    }

    /**
     * Method for setting the trigger value.
     * @param trigger 1 - physics, 2 - offspill, -1 - all events.
     */
    inline void SetTrigger(Int_t trigger) { fTrigger = trigger; }
    inline void SetTpat(Int_t tpat) { fTpat = tpat; }

    /**
     * Method for walk calculation.
     */
    virtual Double_t walk(Int_t inum, Double_t tot);

    /**
     * Method for saturation correction.
     */
    virtual Double_t satu(Int_t inum, Double_t tot, Double_t dt);

    /**
     * Method for finish of the task execution.
     * Is called by the framework after processing the event loop.
     */
    virtual void FinishTask();

  private:
    TClonesArray* fCalItems; /**< Array with Cal items - input data. */
    TClonesArray* fHitItems; /**< Array with Hit items - output data. */

    // check for trigger
    R3BEventHeader* header; /**< Event header. */
    Int_t fTrigger;         /**< Trigger value. */
    Int_t fTpat;

    UInt_t fNofHitItems; /**< Number of hit items for cur event. */
    Double_t fClockFreq; /**< Clock cycle in [ns]. */
    Double_t flosVeffX;
    Double_t flosVeffY;
    Double_t flosOffsetX;
    Double_t flosOffsetY;
    Double_t flosVeffXQ;
    Double_t flosVeffYQ;
    Double_t flosOffsetXQ;
    Double_t flosOffsetYQ;
    Double_t flosVeffXT;
    Double_t flosVeffYT;
    Double_t flosOffsetXT;
    Double_t flosOffsetYT;
    Double_t walk_par[16][11]{}; // Array containing walk parameters: x=PM, y=min,max,p0...p9; MCFD and TAMEX considered
    Double_t tot_par[8][4]{};    // Array containing walk parameters: x=PM, y=p0...p3;
    Int_t OptHisto;
    std::string fwalk_param_file;
    std::string ftot_param_file;

    TClonesArray* fMapped; /**< Array with mapped data - input data. */



    Int_t Icount = 0;

  public:
    ClassDef(R3BLosCal2Hit, 1)
};

#endif

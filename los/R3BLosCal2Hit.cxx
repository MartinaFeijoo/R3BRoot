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

// ------------------------------------------------------------
// -----                  R3BLosCal2Hit                   -----
// -----          Created Mar 10th 2016 by R.Plag         -----
// ------------------------------------------------------------

#include "R3BLosCal2Hit.h"
#include "FairLogger.h"
#include "R3BEventHeader.h"
#include "R3BLosCalData.h"
#include "R3BLosHitData.h"
#include "R3BLosMapped2Cal.h"
#include "R3BLosMappedData.h"
#include "R3BTCalEngine.h"
#include "R3BTCalPar.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include <fstream>
#include <iomanip>
#include <iostream>

#include <algorithm>
#include <cstdlib>
#include <iterator>

using namespace std;
#define IS_NAN(x) TMath::IsNaN(x)

R3BLosCal2Hit::R3BLosCal2Hit()
    : FairTask("LosCal2Hit", 1)
    , fCalItems(NULL)
    , fHitItems(new TClonesArray("R3BLosHitData"))
    , fNofHitItems(0)
    , fTrigger(-1)
    , fTpat(-1)
    , flosVeffX(1.)
    , flosVeffY(1.)
    , flosOffsetX(0.)
    , flosOffsetY(0.)
    , flosVeffXT(1.)
    , flosVeffYT(1.)
    , flosOffsetXT(0.)
    , flosOffsetYT(0.)
    , flosVeffXQ(1.)
    , flosVeffYQ(1.)
    , flosOffsetXQ(0.)
    , flosOffsetYQ(0.)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
{}

R3BLosCal2Hit::R3BLosCal2Hit(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fCalItems(NULL)
    , fHitItems(new TClonesArray("R3BLosHitData"))
    , fNofHitItems(0)
    , fTrigger(-1)
    , fTpat(-1)
    , flosVeffX(1.)
    , flosVeffY(1.)
    , flosOffsetX(0.)
    , flosOffsetY(0.)
    , flosVeffXT(1.)
    , flosVeffYT(1.)
    , flosOffsetXT(0.)
    , flosOffsetYT(0.)
    , flosVeffXQ(1.)
    , flosVeffYQ(1.)
    , flosOffsetXQ(0.)
    , flosOffsetYQ(0.)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
{}

R3BLosCal2Hit::~R3BLosCal2Hit()
{
    if (fHitItems)
    {
        delete fHitItems;
        fHitItems = NULL;
    }
}

InitStatus R3BLosCal2Hit::Init()
{
    // get access to Cal data
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(ERROR) << "FairRootManager not found";

    header = (R3BEventHeader*)mgr->GetObject("R3BEventHeader");

    fCalItems = (TClonesArray*)mgr->GetObject("LosCal");
    if (NULL == fCalItems)
        LOG(ERROR) << "Branch LosCal not found";

    // request storage of Hit data in output tree
    mgr->Register("LosHit", "Land", fHitItems, kTRUE);

    Icount = 0;

    // file with walk-correction parameters
    ifstream infile(fwalk_param_file.c_str());
    if (infile.is_open())
    {
        for (Int_t ivec = 0; ivec < 16; ivec++)
        {
            /*
             ************* Parameters 0-7 MCFD  **********************
             ************* Parameters 8-15 TAMEX *********************
             */

            infile >> walk_par[ivec][0] >> walk_par[ivec][1] >> walk_par[ivec][2] >> walk_par[ivec][3] >>
                walk_par[ivec][4] >> walk_par[ivec][5] >> walk_par[ivec][6] >> walk_par[ivec][7] >> walk_par[ivec][8] >>
                walk_par[ivec][9] >> walk_par[ivec][10];

            /*
                 cout<<setprecision(10)<< ivec<<", "<< walk_par[ivec][0] <<", "<< walk_par[ivec][1] <<", "<<
               walk_par[ivec][2] <<", "<< walk_par[ivec][3] <<", "<< walk_par[ivec][4] <<", "<< walk_par[ivec][5] <<",
               "<< walk_par[ivec][6] <<", "<< walk_par[ivec][7] <<", "<< walk_par[ivec][8] <<", "<< walk_par[ivec][9]
               <<", "<< walk_par[ivec][10] <<endl;
             */
        }
    }
    else // cerr << "Unable to open file \""<<fwalk_param_file<<"\" with walk parameters!"<<endl;
    {
        cout << "*****************************************************************" << endl;
        cout << "UNABLE TO OPEN FILE WITH WALK PARAMETERS! Parameters set to zero!" << endl;
        cout << "*****************************************************************" << endl;
        for (Int_t ivec = 0; ivec < 16; ivec++)
        {
            walk_par[ivec][0] = 10.;
            walk_par[ivec][1] = 1000.;
            for (Int_t iv = 2; iv < 11; iv++)
            {
                walk_par[ivec][iv] = 0.;
                ;
            }
        }
    }
    // file with tot-correction parameters
    ifstream infile1(ftot_param_file.c_str());

    if (infile1.is_open())
    {
        for (Int_t ivec = 0; ivec < 8; ivec++)
        {
            infile1 >> tot_par[ivec][0] >> tot_par[ivec][1] >> tot_par[ivec][2] >> tot_par[ivec][3];
        }
    }
    else // cerr << "Unable to open file \""<<ftot_param_file<<"\" with tot parameters!"<<endl;
    {
        cout << "****************************************************************" << endl;
        cout << "UNABLE TO OPEN FILE WITH ToT PARAMETERS! Parameters set to zero!" << endl;
        cout << "****************************************************************" << endl;
        for (Int_t ivec = 0; ivec < 8; ivec++)
        {
            tot_par[ivec][0] = 0.;
            tot_par[ivec][1] = 0.;
            tot_par[ivec][2] = 0.;
            tot_par[ivec][3] = 1.; // Normalization factor
        }
    }

    // cout << "R3BLosCal2Hit::Init END" << endl;
    return kSUCCESS;
}

InitStatus R3BLosCal2Hit::ReInit() { return kSUCCESS; }

/* Calculate a single hit time for each LOS detector
 *
 * Remember: The times of individual channels depend on the position of
 * the particle on the scintillator. To obtain a precise time of the
 * particle, we need to average either over all four signals (right, top,
 * left, bottom) or over two opposite signals (left+right or top+bottom).
 */
void R3BLosCal2Hit::Exec(Option_t* option)
{
    // cout << "R3BLosCal2Hit::Exec BEGIN: " << Icount << endl;

    // ofstream myFile("data_s473_run197.dat",ios_base::out|ios_base::app);

    // check for requested trigger (Todo: should be done globablly / somewhere else)
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;

    // fTpat = 1-16; fTpat_bit = 0-15
    Int_t fTpat_bit = fTpat - 1;
    if (fTpat_bit >= 0)
    {
        Int_t itpat = header->GetTpat();
        Int_t tpatvalue = (itpat & (1 << fTpat_bit)) >> fTpat_bit;
        if ((header) && (tpatvalue == 0))
            return;
    }

    Int_t nHits = fCalItems->GetEntries();

    if (nHits < 1)
        return;

    // missing times are NAN, hence other times will also
    // be NAN if one time is missing.
    //  cout<<nHits<<"; "<<narr<<endl;

    Double_t time_V[nHits][8]; // [multihit][pm]
    Double_t time_V_temp[nHits][8];
    Double_t time_L[nHits][8];
    Double_t time_T[nHits][8];
    Double_t time_V_corr[nHits][8]; // [multihit][pm]
    Double_t time_L_corr[nHits][8];
    Double_t timeLosM[nHits];
    Double_t LosTresM[nHits];
    Double_t timeLosT[nHits];
    Double_t LosTresT[nHits];
    Double_t timeLos[nHits];
    Double_t timeLosM_corr[nHits];
    Double_t LosTresM_corr[nHits];
    Double_t timeLosT_corr[nHits];
    Double_t LosTresT_corr[nHits];
    Double_t timeLos_corr[nHits];
    Double_t totsum[nHits];
    Double_t tot[nHits][8];
    Double_t totsum_corr[nHits];
    Double_t tot_corr[nHits][8];
    Double_t xT_cm[nHits];
    Double_t yT_cm[nHits];
    Double_t xV_cm[nHits];
    Double_t yV_cm[nHits];
    Double_t xToT_cm[nHits];
    Double_t yToT_cm[nHits];
    Double_t x_cm[nHits];
    Double_t y_cm[nHits];
    Double_t Z[nHits];
    Double_t t_hit[nHits]; // NAN

    memset(time_V, 0. / 0., nHits * 8 * sizeof(Double_t));
    memset(time_V_temp, 0. / 0., nHits * 8 * sizeof(Double_t));
    memset(time_L, 0. / 0., nHits * 8 * sizeof(Double_t));
    memset(time_T, 0. / 0., nHits * 8 * sizeof(Double_t));
    memset(time_V_corr, 0. / 0., nHits * 8 * sizeof(Double_t));
    memset(time_L_corr, 0. / 0., nHits * 8 * sizeof(Double_t));
    memset(timeLosM, 0., nHits * sizeof(Double_t));
    memset(LosTresM, 0. / 0., nHits * sizeof(Double_t));
    memset(timeLosT, 0., nHits * sizeof(Double_t));
    memset(LosTresT, 0. / 0., nHits * sizeof(Double_t));
    memset(timeLos, 0., nHits * sizeof(Double_t));
    memset(timeLosM_corr, 0., nHits * sizeof(Double_t));
    memset(LosTresM_corr, 0. / 0., nHits * sizeof(Double_t));
    memset(timeLosT_corr, 0., nHits * sizeof(Double_t));
    memset(LosTresT_corr, 0. / 0., nHits * sizeof(Double_t));
    memset(timeLos_corr, 0., nHits * sizeof(Double_t));
    memset(totsum, 0., nHits * sizeof(Double_t));
    memset(tot, 0. / 0., nHits * 8 * sizeof(Double_t));
    memset(totsum_corr, 0., nHits * sizeof(Double_t));
    memset(tot_corr, 0. / 0., nHits * 8 * sizeof(Double_t));
    memset(xT_cm, 0. / 0., nHits * sizeof(Double_t));
    memset(yT_cm, 0. / 0., nHits * sizeof(Double_t));
    memset(xV_cm, 0. / 0., nHits * sizeof(Double_t));
    memset(yV_cm, 0. / 0., nHits * sizeof(Double_t));
    memset(xToT_cm, 0. / 0., nHits * sizeof(Double_t));
    memset(yToT_cm, 0. / 0., nHits * sizeof(Double_t));
    memset(x_cm, 0. / 0., nHits * sizeof(Double_t));
    memset(y_cm, 0. / 0., nHits * sizeof(Double_t));
    memset(Z, 0., nHits * sizeof(Double_t));
    memset(t_hit, 0. / 0., nHits * sizeof(Double_t));

    Int_t nDet = 0;


    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {

        R3BLosCalData* calItem = (R3BLosCalData*)fCalItems->At(ihit);

        nDet = calItem->GetDetector();

        // lt=0, l=1,lb=2,b=3,rb=4,r=5,rt=6,t=7
        for (Int_t iCha = 0; iCha < 8; iCha++)
        {
            time_V[ihit][iCha] = 0. / 0.;
            if (!(IS_NAN(calItem->GetTimeV_ns(iCha))))
            { // VFTX
                time_V[ihit][iCha] = calItem->GetTimeV_ns(iCha);
            }
            time_L[ihit][iCha] = 0. / 0.;
            if (!(IS_NAN(calItem->GetTimeL_ns(iCha))))
            { // TAMEX leading
                time_L[ihit][iCha] = calItem->GetTimeL_ns(iCha);
            }
            time_T[ihit][iCha] = 0. / 0.;
            if (!(IS_NAN(calItem->GetTimeT_ns(iCha))))
            { // TAMEX trailing
                time_T[ihit][iCha] = calItem->GetTimeT_ns(iCha);
            }
        }
        if (!calItem)
        {
            cout << " !calItem" << endl;
            continue; // can this happen?
        }
    }

    // Sorting VFTX data:

    std::qsort(time_V, nHits, sizeof(*time_V), [](const void* arg1, const void* arg2) -> int {
        double const* lhs = static_cast<double const*>(arg1);
        double const* rhs = static_cast<double const*>(arg2);

        return (lhs[0] < rhs[0]) ? -1 : ((rhs[0] < lhs[0]) ? 1 : 0);
    });
    // End sorting

    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {
        Bool_t iLOSTypeMCFD = false;
        Bool_t iLOSTypeTAMEX = false;
        Bool_t iLOSType = false;
        UInt_t Igood_event = 0;

        if (time_V[ihit][0] > 0. && !(IS_NAN(time_V[ihit][0])) && time_V[ihit][1] > 0. && !(IS_NAN(time_V[ihit][1])) &&
            time_V[ihit][2] > 0. && !(IS_NAN(time_V[ihit][2])) && time_V[ihit][3] > 0. && !(IS_NAN(time_V[ihit][3])) &&
            time_V[ihit][4] > 0. && !(IS_NAN(time_V[ihit][4])) && time_V[ihit][5] > 0. && !(IS_NAN(time_V[ihit][5])) &&
            time_V[ihit][6] > 0. && !(IS_NAN(time_V[ihit][6])) && time_V[ihit][7] > 0. && !(IS_NAN(time_V[ihit][7])))
        {
            iLOSTypeMCFD = true; // all 8 MCFD times
        }

        if (time_L[ihit][0] > 0. && !(IS_NAN(time_L[ihit][0])) && time_L[ihit][1] > 0. && !(IS_NAN(time_L[ihit][1])) &&
            time_L[ihit][2] > 0. && !(IS_NAN(time_L[ihit][2])) && time_L[ihit][3] > 0. && !(IS_NAN(time_L[ihit][3])) &&
            time_L[ihit][4] > 0. && !(IS_NAN(time_L[ihit][4])) && time_L[ihit][5] > 0. && !(IS_NAN(time_L[ihit][5])) &&
            time_L[ihit][6] > 0. && !(IS_NAN(time_L[ihit][6])) && time_L[ihit][7] > 0. && !(IS_NAN(time_L[ihit][7])) &&

            time_T[ihit][0] > 0. && !(IS_NAN(time_T[ihit][0])) && time_T[ihit][1] > 0. && !(IS_NAN(time_T[ihit][1])) &&
            time_T[ihit][2] > 0. && !(IS_NAN(time_T[ihit][2])) && time_T[ihit][3] > 0. && !(IS_NAN(time_T[ihit][3])) &&
            time_T[ihit][4] > 0. && !(IS_NAN(time_T[ihit][4])) && time_T[ihit][5] > 0. && !(IS_NAN(time_T[ihit][5])) &&
            time_T[ihit][6] > 0. && !(IS_NAN(time_T[ihit][6])) && time_T[ihit][7] > 0. && !(IS_NAN(time_T[ihit][7])))
        {
            iLOSTypeTAMEX = true; // all 8 leading and trailing times
        }

        // We will consider only events in which booth MCFD and TAMEX see same number of channels:
        if (iLOSTypeTAMEX && iLOSTypeMCFD)
            iLOSType = true;
        //   if(iLOSTypeMCFD ) LOSType = true;

        if (iLOSType)
        {
            // calculate time over threshold and check if clock counter went out of range

            Igood_event = 1;

            int nPMT = 0;
            int nPMV = 0;
            for (int ipm = 0; ipm < 8; ipm++)
            {
                tot[ihit][ipm] = 0. / 0.;
                if (time_T[ihit][ipm] > 0. && time_L[ihit][ipm] > 0. && !(IS_NAN(time_T[ihit][ipm])) &&
                    !(IS_NAN(time_L[ihit][ipm])))
                {
                    while (time_T[ihit][ipm] - time_L[ihit][ipm] <= 0.)
                    {
                        time_T[ihit][ipm] = time_T[ihit][ipm] + 2048. * fClockFreq;
                    }

                    nPMT = nPMT + 1;
                    tot[ihit][ipm] = time_T[ihit][ipm] - time_L[ihit][ipm];
                }

                if (tot[ihit][ipm] != 0. && !(IS_NAN(tot[ihit][ipm])))
                    totsum[ihit] += tot[ihit][ipm];

                // Corrected for saturation and absorption
                if (ihit > 0)
                {
                    Double_t dthit = time_V[ihit][ipm] - time_V[ihit - 1][ipm];
                    tot_corr[ihit][ipm] = satu(ipm, tot[ihit][ipm], dthit);
                }
                else
                    tot_corr[ihit][ipm] = tot[ihit][ipm];

                totsum_corr[ihit] += tot_corr[ihit][ipm];

                if (time_L[ihit][ipm] > 0. && !(IS_NAN(time_L[ihit][ipm])))
                    timeLosT[ihit] += time_L[ihit][ipm];

                if (time_V[ihit][ipm] > 0. && !(IS_NAN(time_V[ihit][ipm])))
                {
                    timeLosM[ihit] += time_V[ihit][ipm];
                    nPMV = nPMV + 1;
                }
            }

            timeLosM[ihit] = timeLosM[ihit] / nPMV;
            timeLosT[ihit] = timeLosT[ihit] / nPMT;

            totsum[ihit] = totsum[ihit] / nPMT;
            totsum_corr[ihit] = totsum_corr[ihit] / nPMT;

            // Time resolution TAMEX
            LosTresT[ihit] = ((time_L[ihit][0] + time_L[ihit][2] + time_L[ihit][4] + time_L[ihit][6]) -
                              (time_L[ihit][1] + time_L[ihit][3] + time_L[ihit][5] + time_L[ihit][7])) /
                             4.;

            // Time resolution MCFD
            LosTresM[ihit] = ((time_V[ihit][0] + time_V[ihit][2] + time_V[ihit][4] + time_V[ihit][6]) -
                              (time_V[ihit][1] + time_V[ihit][3] + time_V[ihit][5] + time_V[ihit][7])) /
                             4.;

            // Position TAMEX:
            xT_cm[ihit] = ((time_L[ihit][5] + time_L[ihit][6]) / 2. - (time_L[ihit][1] + time_L[ihit][2]) / 2.) * (-1.);
            yT_cm[ihit] = ((time_L[ihit][7] + time_L[ihit][0]) / 2. - (time_L[ihit][3] + time_L[ihit][4]) / 2.) * (-1.);
            xT_cm[ihit] = (xT_cm[ihit] - flosOffsetXT) * flosVeffXT;
            yT_cm[ihit] = (yT_cm[ihit] - flosOffsetYT) * flosVeffYT;

            // Position MCFD:
            xV_cm[ihit] = ((time_V[ihit][5] + time_V[ihit][6]) / 2. - (time_V[ihit][1] + time_V[ihit][2]) / 2.) * (-1.);
            yV_cm[ihit] = ((time_V[ihit][7] + time_V[ihit][0]) / 2. - (time_V[ihit][3] + time_V[ihit][4]) / 2.) * (-1.);
            xV_cm[ihit] = (xV_cm[ihit] - flosOffsetX) * flosVeffX;
            yV_cm[ihit] = (yV_cm[ihit] - flosOffsetY) * flosVeffY;

            /*
                    for(int ipm=0; ipm<8; ipm++)
                    {
                              if(Icount < 500000)
                     {
                        myFile<<setprecision(10)<<ihit<<" "<<ipm<<" "<<time_V[ihit][ipm]<<" "<<time_L[ihit][ipm]<<"
               "<<tot[ihit][ipm]<<endl;
                      }
                    }
            */

            // Walk correction for MCFD and TAMEX
            for (int ipm = 0; ipm < 8; ipm++)
            {
                time_V_corr[ihit][ipm] = time_V[ihit][ipm] - walk(ipm, tot[ihit][ipm]);
                time_L_corr[ihit][ipm] = time_L[ihit][ipm] - walk(ipm + 8, tot[ihit][ipm]);
                timeLosM_corr[ihit] += time_V_corr[ihit][ipm];
                timeLosT_corr[ihit] += time_L_corr[ihit][ipm];
            }

            // Walk-corrected time-properties

            // TAMEX:
            LosTresT_corr[ihit] =
                ((time_L_corr[ihit][0] + time_L_corr[ihit][2] + time_L_corr[ihit][4] + time_L_corr[ihit][6]) -
                 (time_L_corr[ihit][1] + time_L_corr[ihit][3] + time_L_corr[ihit][5] + time_L_corr[ihit][7])) /
                4.;
            timeLosT_corr[ihit] = timeLosT_corr[ihit] / nPMT;

            // MCFD
            LosTresM_corr[ihit] = ((time_V_corr[ihit][0] + time_V[ihit][2] + time_V[ihit][4] + time_V[ihit][6]) -
                                   (time_V_corr[ihit][1] + time_V[ihit][3] + time_V[ihit][5] + time_V[ihit][7])) /
                                  4.;
            timeLosM_corr[ihit] = timeLosM_corr[ihit] / nPMV;

            // Position from ToT:
            xToT_cm[ihit] = (((tot[ihit][5] + tot[ihit][6]) / 2. - (tot[ihit][1] + tot[ihit][2]) / 2.) /
                             ((tot[ihit][1] + tot[ihit][2] + tot[ihit][5] + tot[ihit][6]) / 4.));

            yToT_cm[ihit] = (((tot[ihit][0] + tot[ihit][7]) / 2. - (tot[ihit][3] + tot[ihit][4]) / 2.) /
                             ((tot[ihit][7] + tot[ihit][0] + tot[ihit][3] + tot[ihit][4]) / 4.));

            xToT_cm[ihit] = (xToT_cm[ihit] - flosOffsetXQ) * flosVeffXQ;
            yToT_cm[ihit] = (yToT_cm[ihit] - flosOffsetYQ) * flosVeffYQ;

            x_cm[ihit] = xV_cm[ihit];
            y_cm[ihit] = yV_cm[ihit];
            Z[ihit] = totsum_corr[ihit];
            t_hit[ihit] = timeLosM_corr[ihit];



        } // end of LOSTYPE = true

        new ((*fHitItems)[fNofHitItems]) R3BLosHitData(nDet, t_hit[ihit], x_cm[ihit], y_cm[ihit], Z[ihit]);
        fNofHitItems += 1;

        Icount++;
     }

    // cout << "R3BLosCal2Hit::Exec END: " << Icount << endl;
}

void R3BLosCal2Hit::FinishEvent()
{

    if (fHitItems)
    {
        fHitItems->Clear();
        fNofHitItems = 0;
    }
}

void R3BLosCal2Hit::FinishTask() {}

Double_t R3BLosCal2Hit::walk(Int_t inum, Double_t tot)
{

    Double_t y = 0. / 0., ysc = 0. / 0., term[8] = { 0. };
    Double_t x;

    x = tot;
    term[0] = x;
    for (Int_t i = 0; i < 7; i++)
    {
        term[i + 1] = term[i] * x;
    }

    ysc = walk_par[inum][2] + walk_par[inum][3] * term[0] + walk_par[inum][4] * term[1] + walk_par[inum][5] * term[2] +
          walk_par[inum][6] * term[3] + walk_par[inum][7] * term[4] + walk_par[inum][8] * term[5] +
          walk_par[inum][9] * term[6] + walk_par[inum][10] * term[7];

    if (tot < walk_par[inum][0] || tot > walk_par[inum][1])
        ysc = 0.0 / 0.0;

    /*
       cout<< inum<<", "<< walk_par[inum][0] <<", "<< walk_par[inum][1] <<", "<< walk_par[inum][2] <<", "<<
                           walk_par[inum][3] <<", "<< walk_par[inum][4] <<", "<< walk_par[inum][5] <<", "<<
                           walk_par[inum][6] <<", "<< walk_par[inum][7] <<", "<< walk_par[inum][8] <<", "<<
                           walk_par[inum][9] <<", "<< walk_par[inum][10] <<endl;
    */
    return ysc;
}

Double_t R3BLosCal2Hit::satu(Int_t inum, Double_t tot, Double_t dt)
{

    Double_t ysc = 0. / 0.;

    // if(tot_par[inum][0] > 0.)
    // ysc  = (tot_par[inum][0]*tot+tot_par[inum][1])/(tot_par[inum][2]-tot)*tot_par[inum][3] ;

    ysc = tot_par[inum][0] + tot_par[inum][1] * (1. - 1. / (exp((dt - tot_par[inum][2]) / tot_par[inum][3]) + 1.));
    ysc = tot / ysc;
    ysc = ysc * (tot_par[inum][0] + tot_par[inum][1]);

    /*
     cout<< inum<<", "<< tot_par[inum][0] <<", "<< tot_par[inum][1] <<", "<< tot_par[inum][2]<<", "<< tot_par[inum][3]
     <<endl; cout<<tot<<", "<<ysc<<endl;
  */
    return ysc;
}

ClassImp(R3BLosCal2Hit)

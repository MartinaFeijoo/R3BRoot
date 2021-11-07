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
// -----             R3BCalifaJulichOnlineSpectra                 -----
// -----    Created 16/07/21  by J.L. Rodriguez-Sanchez   -----
// -----          Fill CalifaJulich online histograms             -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with CalifaJulich online data
 */

#include "R3BCalifaJulichOnlineSpectra.h"
#include "R3BAmsMappedData.h"
#include "R3BAmsStripCalData.h"
#include "R3BAmsHitData.h"
#include "R3BCalifaMappedData.h"
#include "R3BCalifaCrystalCalData.h"
#include "R3BCalifaHitData.h"

#include "R3BEventHeader.h"
#include "THttpServer.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"
#include "TCanvas.h"
#include "TFolder.h"
#include "TH1F.h"
#include "TH2F.h"

#include "TClonesArray.h"
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>

R3BCalifaJulichOnlineSpectra::R3BCalifaJulichOnlineSpectra()
    : FairTask("CalifaJulichOnlineSpectra", 1)
    , fMappedItemsCalifa(NULL)
    , fCalItemsCalifa(NULL)
    , fHitItemsCalifa(NULL)
    , fMappedItemsSi(NULL)
    , fCalItemsSi(NULL)
    , fHitItemsSi(NULL)
    , fTrigger(-1)
    , fNEvents(0)
    , fNbDet(1)
    , fNbCrystals(8)
{
}

R3BCalifaJulichOnlineSpectra::R3BCalifaJulichOnlineSpectra(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMappedItemsCalifa(NULL)
    , fCalItemsCalifa(NULL)
    , fHitItemsCalifa(NULL)
    , fMappedItemsSi(NULL)
    , fCalItemsSi(NULL)
    , fHitItemsSi(NULL)
    , fTrigger(-1)
    , fNEvents(0)
    , fNbDet(1)
    , fNbCrystals(8)
{
}

R3BCalifaJulichOnlineSpectra::~R3BCalifaJulichOnlineSpectra()
{
    LOG(DEBUG) << "R3BCalifaJulichOnlineSpectra::Delete instance";
    if (fMappedItemsCalifa)
        delete fMappedItemsCalifa;
    if (fMappedItemsSi)
        delete fMappedItemsSi;
    if (fCalItemsCalifa)
        delete fCalItemsCalifa;
    if (fCalItemsSi)
        delete fCalItemsSi;
    if (fHitItemsSi)
        delete fHitItemsSi;
}

InitStatus R3BCalifaJulichOnlineSpectra::Init()
{
    LOG(INFO) << "R3BCalifaJulichOnlineSpectra::Init()";

    // Looking for FairRootManager
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(FATAL) << "R3BCalifaJulichOnlineSpectra::FairRootManager not found";

    // Create histograms for all the detectors

    // Get access to Mapped data
    fMappedItemsSi = (TClonesArray*)mgr->GetObject("AmsMappedData");
    if (!fMappedItemsSi)
    {
        LOG(FATAL) << "R3BCalifaJulichOnlineSpectra::AmsMappedData not found";
        return kFATAL;
    }

    // Get access to Mapped data
    fMappedItemsCalifa = (TClonesArray*)mgr->GetObject("CalifaMappedData");
    if (!fMappedItemsCalifa)
    {
        LOG(FATAL) << "R3BCalifaJulichOnlineSpectra::CalifaMappedData not found";
        return kFATAL;
    }

    // Get access to Cal data
    fCalItemsCalifa = (TClonesArray*)mgr->GetObject("CalifaCrystalCalData");
    if (!fCalItemsCalifa)
        LOG(WARNING) << "R3BCalifaJulichOnlineSpectra::CalifaCrystalCalData not found";

    fCalItemsSi = (TClonesArray*)mgr->GetObject("AmsStripCalData");
    if (!fCalItemsSi)
        LOG(WARNING) << "R3BCalifaJulichOnlineSpectra::AmsStripCalData not found";

      fHitItemsSi = (TClonesArray*)mgr->GetObject("AmsHitData");
      if (!fHitItemsSi)
          LOG(WARNING) << "R3BCalifaJulichOnlineSpectra::AmsHitData not found";

    // // Get access to Hit data
    // fHitItems = (TClonesArray*)mgr->GetObject("CalifaJulichSiHitData");
    // if (!fHitItems)
    //     LOG(WARNING) << "R3BCalifaJulichOnlineSpectra::CalifaJulichSiHitData not found";

    // Energy range for strips
    Double_t binsE = 200;
    Double_t minE = 0;
    Double_t maxE = 3500;
    char Name1[255];
    char Name2[255];

    // MAIN FOLDERs
    TFolder* mainfolCalifa = new TFolder("Califa", "Califa info");
    TFolder* mainfolSi = new TFolder("Si", "Si info");
    // Folders for mapped data
    TFolder* mapfolCalifa = new TFolder("Mapped", "Califa Mapped");
    TFolder* mapfolSi = new TFolder("MappedSi", "Si Mapped");
    // Folders for cal data
    TFolder* calfolCalifa = new TFolder("Cal", "Cal Califa info");
    TFolder* calfolSi = new TFolder("CalSi", "Cal Si info");
    // Folders for hit data and correlations
    TFolder* hitfolCalifa = new TFolder("Correlations", "Hit Califa info");
    TFolder* hitfolSi = new TFolder("HitSi", "Hit Si info");

    // Mapped data Si
    fh2_EnergyVsStrip.resize(fNbDet);
    for (Int_t i = 0; i < fNbDet; i++)
    { // one histo per detector
        sprintf(Name1, "fh2_energy_vs_strip_det_%d", i + 1);
        sprintf(Name2, "Energy vs strip number for Si Det: %d", i + 1);
        fh2_EnergyVsStrip[i] = new TH2F(Name1, Name2, 64, 0, 64, binsE, minE, maxE);
        fh2_EnergyVsStrip[i]->GetXaxis()->SetTitle("Strip number");
        fh2_EnergyVsStrip[i]->GetYaxis()->SetTitle("Energy [channels]");
        fh2_EnergyVsStrip[i]->GetYaxis()->SetTitleOffset(1.4);
        fh2_EnergyVsStrip[i]->GetXaxis()->CenterTitle(true);
        fh2_EnergyVsStrip[i]->GetYaxis()->CenterTitle(true);
        fh2_EnergyVsStrip[i]->Draw("col");
        mapfolSi->Add(fh2_EnergyVsStrip[i]);
    }
    mainfolSi->Add(mapfolSi);

    // Mapped data Califa
    fh1_EnergyCalifaCrystals.resize(fNbCrystals);
    for (Int_t i=0; i<fNbCrystals; i++)
    {
      sprintf(Name1, "fh1_energy_CalifaCrystal_%d", i + 1);
      sprintf(Name2, "Energy in Califa crystal: %d", i + 1);
      fh1_EnergyCalifaCrystals[i] = new TH1F(Name1, Name2, 3000, 0, 3000);
      fh1_EnergyCalifaCrystals[i]->GetXaxis()->SetTitle("Energy [channels]");
      fh1_EnergyCalifaCrystals[i]->GetYaxis()->SetTitle("counts");
      fh1_EnergyCalifaCrystals[i]->GetYaxis()->SetTitleOffset(1.4);
      fh1_EnergyCalifaCrystals[i]->GetXaxis()->CenterTitle(true);
      fh1_EnergyCalifaCrystals[i]->GetYaxis()->CenterTitle(true);
      fh1_EnergyCalifaCrystals[i]->Draw("col");
      mapfolCalifa->Add(fh1_EnergyCalifaCrystals[i]);
    }
    mainfolCalifa->Add(mapfolCalifa);

    // Cal data Si
    fh2_EnergyCalVsStrip.resize(fNbDet);
    for (Int_t i = 0; i < fNbDet; i++)
    { // one histo per detector
        sprintf(Name1, "fh2_energyCal_vs_strip_det_%d", i + 1);
        sprintf(Name2, "Energy calibrated vs strip number for Si Det: %d", i + 1);
        fh2_EnergyCalVsStrip[i] = new TH2F(Name1, Name2, 64, 0, 64, binsE, minE, maxE);
        fh2_EnergyCalVsStrip[i]->GetXaxis()->SetTitle("Strip number");
        fh2_EnergyCalVsStrip[i]->GetYaxis()->SetTitle("Energy [keV]");
        fh2_EnergyCalVsStrip[i]->GetYaxis()->SetTitleOffset(1.4);
        fh2_EnergyCalVsStrip[i]->GetXaxis()->CenterTitle(true);
        fh2_EnergyCalVsStrip[i]->GetYaxis()->CenterTitle(true);
        fh2_EnergyCalVsStrip[i]->Draw("col");
        calfolSi->Add(fh2_EnergyCalVsStrip[i]);
    }
    mainfolSi->Add(calfolSi);

    // Cal data CALIFA
    fh1_EnergyCalCalifaCrystals.resize(fNbCrystals);
    for (Int_t i=0; i<fNbCrystals; i++)
    {
      sprintf(Name1, "fh1_energyCal_CalifaCrystal_%d", i + 1);
      sprintf(Name2, "Energy calibrated in Califa crystal: %d", i + 1);
      fh1_EnergyCalCalifaCrystals[i] = new TH1F(Name1, Name2, 3000, 0, 3000);
      fh1_EnergyCalCalifaCrystals[i]->GetXaxis()->SetTitle("Energy [keV]");
      fh1_EnergyCalCalifaCrystals[i]->GetYaxis()->SetTitle("counts");
      fh1_EnergyCalCalifaCrystals[i]->GetYaxis()->SetTitleOffset(1.4);
      fh1_EnergyCalCalifaCrystals[i]->GetXaxis()->CenterTitle(true);
      fh1_EnergyCalCalifaCrystals[i]->GetYaxis()->CenterTitle(true);
      fh1_EnergyCalCalifaCrystals[i]->Draw("col");
      calfolCalifa->Add(fh1_EnergyCalCalifaCrystals[i]);
    }
    mainfolCalifa->Add(calfolCalifa);

    // energy correlations CALIFA
    fh2_EnergyCorrelationsCrystals.resize(12);
    for (Int_t i=0; i<12; i++)
    {
      sprintf(Name1, "fh2_EnergyCorrelationsCrystals_%d", i);
      sprintf(Name2, "Energy Hit in Califa crystal");
      fh2_EnergyCorrelationsCrystals[i] = new TH2F(Name1, Name2, 3000, 0, 3000, 3000,0,3000);
      fh2_EnergyCorrelationsCrystals[i]->GetXaxis()->SetTitle("Energy [keV]");
      fh2_EnergyCorrelationsCrystals[i]->GetYaxis()->SetTitle("counts");
      fh2_EnergyCorrelationsCrystals[i]->GetYaxis()->SetTitleOffset(1.4);
      fh2_EnergyCorrelationsCrystals[i]->GetXaxis()->CenterTitle(true);
      fh2_EnergyCorrelationsCrystals[i]->GetYaxis()->CenterTitle(true);
      fh2_EnergyCorrelationsCrystals[i]->Draw("col");
      hitfolCalifa->Add(fh2_EnergyCorrelationsCrystals[i]);
    }
    mainfolCalifa->Add(hitfolCalifa);

    // Hit data Si
    fh2_EnergyHitVsStrip.resize(fNbDet);
    for (Int_t i = 0; i < fNbDet; i++)
    { // one histo per detector
        sprintf(Name1, "fh2_energyHit_vs_strip_det_%d", i + 1);
        sprintf(Name2, "Energy hit vs strip number for Si Det: %d", i + 1);
        fh2_EnergyHitVsStrip[i] = new TH2F(Name1, Name2, 32, -30, 30, 32, -30, 30);
        fh2_EnergyHitVsStrip[i]->GetXaxis()->SetTitle("Position X [mm]");
        fh2_EnergyHitVsStrip[i]->GetYaxis()->SetTitle("Position Y [mm]");
        fh2_EnergyHitVsStrip[i]->GetYaxis()->SetTitleOffset(1.4);
        fh2_EnergyHitVsStrip[i]->GetXaxis()->CenterTitle(true);
        fh2_EnergyHitVsStrip[i]->GetYaxis()->CenterTitle(true);
        fh2_EnergyHitVsStrip[i]->Draw("col");
        hitfolSi->Add(fh2_EnergyHitVsStrip[i]);
    }
    mainfolSi->Add(hitfolSi);

    // Looking for FairRunOnline
    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);
    run->AddObject(mainfolSi);
    run->AddObject(mainfolCalifa);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_CalifaJulich", Form("/Objects/%s/->Reset_CalifaJulich_Histo()", GetName()));

    return kSUCCESS;
}

void R3BCalifaJulichOnlineSpectra::Reset_CalifaJulich_Histo()
{
    LOG(INFO) << "R3BCalifaJulichOnlineSpectra::Reset_CalifaJulich_Histo";

    // Mapped data
    for (Int_t i = 0; i < fNbDet; i++)
    {
        fh2_EnergyVsStrip[i]->Reset();
        fh2_EnergyCalVsStrip[i]->Reset();
        fh2_EnergyHitVsStrip[i]->Reset();
    }

    for (Int_t i = 0; i < fNbCrystals; i++)
    {
        fh1_EnergyCalifaCrystals[i]->Reset();
        fh1_EnergyCalCalifaCrystals[i]->Reset();
        fh1_EnergyHitCalifaCrystals[i]->Reset();
    }

    for (Int_t i=0;i<12;i++) { fh2_EnergyCorrelationsCrystals[i] -> Reset(); }
}

void R3BCalifaJulichOnlineSpectra::Exec(Option_t* option)
{
    // Fill mapped data
    if (fMappedItemsSi && fMappedItemsSi->GetEntriesFast() > 0)
    {
        auto nHits = fMappedItemsSi->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BAmsMappedData* hit = (R3BAmsMappedData*)fMappedItemsSi->At(ihit);
            if (!hit)
                continue;
            fh2_EnergyVsStrip[hit->GetDetectorId()]->Fill(hit->GetStripId(), hit->GetEnergy());
        }
    }

    if (fMappedItemsCalifa && fMappedItemsCalifa->GetEntriesFast() > 0)
    {
        auto nHits = fMappedItemsCalifa->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BCalifaMappedData* hit = (R3BCalifaMappedData*)fMappedItemsCalifa->At(ihit);
            if (!hit)
                continue;
            fh1_EnergyCalifaCrystals[hit->GetCrystalId()-1]->Fill(hit->GetEnergy());
        }
    }

    // Fill Cal data
    if (fCalItemsSi && fCalItemsSi->GetEntriesFast() > 0)
    {
        auto nHits = fCalItemsSi->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BAmsStripCalData* hit = (R3BAmsStripCalData*)fCalItemsSi->At(ihit);
            if (!hit)
                continue;
            fh2_EnergyCalVsStrip[hit->GetDetId()]->Fill(hit->GetStripId()+32*hit->GetSideId(), hit->GetEnergy());
        }
    }

    if (fCalItemsCalifa && fCalItemsCalifa->GetEntriesFast() > 0)
    {
        auto nHits = fCalItemsCalifa->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BCalifaCrystalCalData* hit = (R3BCalifaCrystalCalData*)fCalItemsCalifa->At(ihit);
            if (!hit)
                continue;
            fh1_EnergyCalCalifaCrystals[hit->GetCrystalId()-1]->Fill(hit->GetEnergy());
        }
    }

    // Fill Hit data Si
    if (fHitItemsSi && fHitItemsSi->GetEntriesFast() > 0)
    {
        auto nHits = fHitItemsSi->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BAmsHitData* hit = (R3BAmsHitData*)fHitItemsSi->At(ihit);
            if (!hit)
                continue;
            fh2_EnergyHitVsStrip[hit->GetDetId()]->Fill(hit->GetPosLab().X(), hit->GetPosLab().Y());
        }
    }
    

    fNEvents += 1;
}

void R3BCalifaJulichOnlineSpectra::FinishEvent()
{
    if (fMappedItemsSi)
      {fMappedItemsSi->Clear();}

    if (fMappedItemsCalifa)
      {fMappedItemsCalifa->Clear();}

    if (fCalItemsCalifa)
      {fCalItemsCalifa->Clear();}

    if (fCalItemsSi)
      {fCalItemsSi->Clear();}

    if (fHitItemsCalifa)
      {fHitItemsCalifa->Clear();}

    if (fHitItemsSi)
      {fHitItemsSi->Clear();}
}

void R3BCalifaJulichOnlineSpectra::FinishTask()
{

   for (Int_t i = 0; i < fNbDet; i++)
      {
          if (fMappedItemsSi) { fh2_EnergyVsStrip[i]->Write(); }
          if (fCalItemsSi) { fh2_EnergyCalVsStrip[i]->Write(); }
          if (fHitItemsSi) { fh2_EnergyHitVsStrip[i]->Write(); }
      }

      for (Int_t i = 0; i < fNbCrystals; i++)
      {
          if (fMappedItemsCalifa){ fh1_EnergyCalifaCrystals[i]->Write(); }
          if (fCalItemsCalifa){ fh1_EnergyCalCalifaCrystals[i]->Write(); }
          if (fHitItemsCalifa){ fh1_EnergyHitCalifaCrystals[i]->Write(); }
      }

}

ClassImp(R3BCalifaJulichOnlineSpectra);

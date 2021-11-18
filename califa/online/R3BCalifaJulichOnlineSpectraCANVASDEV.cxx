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
// -----             R3BCalifaJulichOnlineSpectraCANVASDEV                 -----
// -----    Created 16/07/21  by J.L. Rodriguez-Sanchez   -----
// -----          Fill CalifaJulich online histograms             -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with CalifaJulich online data
 */

#include "R3BCalifaJulichOnlineSpectraCANVASDEV.h"
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

R3BCalifaJulichOnlineSpectraCANVASDEV::R3BCalifaJulichOnlineSpectraCANVASDEV()
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
    , fNbCrystals(16)
{
}

R3BCalifaJulichOnlineSpectraCANVASDEV::R3BCalifaJulichOnlineSpectraCANVASDEV(const TString& name, Int_t iVerbose)
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
    , fNbCrystals(16)
{
}

R3BCalifaJulichOnlineSpectraCANVASDEV::~R3BCalifaJulichOnlineSpectraCANVASDEV()
{
    LOG(DEBUG) << "R3BCalifaJulichOnlineSpectraCANVASDEV::Delete instance";
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

InitStatus R3BCalifaJulichOnlineSpectraCANVASDEV::Init()
{
    LOG(INFO) << "R3BCalifaJulichOnlineSpectraCANVASDEV::Init()";

    // Looking for FairRootManager
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(FATAL) << "R3BCalifaJulichOnlineSpectraCANVASDEV::FairRootManager not found";

    // Create histograms for all the detectors

    // Get access to Mapped data
    fMappedItemsSi = (TClonesArray*)mgr->GetObject("AmsMappedData");
    if (!fMappedItemsSi)
    {
        LOG(FATAL) << "R3BCalifaJulichOnlineSpectraCANVASDEV::AmsMappedData not found";
        return kFATAL;
    }

    // Get access to Mapped data
    fMappedItemsCalifa = (TClonesArray*)mgr->GetObject("CalifaMappedData");
    if (!fMappedItemsCalifa)
    {
        LOG(FATAL) << "R3BCalifaJulichOnlineSpectraCANVASDEV::CalifaMappedData not found";
        return kFATAL;
    }

    // Get access to Cal data
    fCalItemsCalifa = (TClonesArray*)mgr->GetObject("CalifaCrystalCalData");
    if (!fCalItemsCalifa)
        LOG(WARNING) << "R3BCalifaJulichOnlineSpectraCANVASDEV::CalifaCrystalCalData not found";

    fCalItemsSi = (TClonesArray*)mgr->GetObject("AmsStripCalData");
    if (!fCalItemsSi)
        LOG(WARNING) << "R3BCalifaJulichOnlineSpectraCANVASDEV::AmsStripCalData not found";

      fHitItemsSi = (TClonesArray*)mgr->GetObject("AmsHitData");
      if (!fHitItemsSi)
          LOG(WARNING) << "R3BCalifaJulichOnlineSpectraCANVASDEV::AmsHitData not found";


    // Energy range for strips
    Double_t binsE = 200;
    Double_t minE = 0;
    Double_t maxE = 3500;
    char Name1[255];
    char Name2[255];


    cSilicon = new TCanvas("cSilicon", "cSilicon", 10, 10, 500, 500);
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
    }


    // Mapped data Califa

    fh1_MultiplicityGamma = new TH1F("MultGamma", "MultGamma", 40, 0, 20);
    fh1_MultiplicityGamma->GetXaxis()->SetTitle("Multiplicity Gamma");
    fh1_MultiplicityGamma->GetYaxis()->SetTitleOffset(1.4);
    fh1_MultiplicityGamma->GetXaxis()->CenterTitle(true);
    fh1_MultiplicityGamma->GetYaxis()->CenterTitle(true);
    fh1_MultiplicityGamma->Draw("col");


    fh1_MultiplicityProton = new TH1F("MultProton", "MultProton", 40, 0, 20);
    fh1_MultiplicityProton->GetXaxis()->SetTitle("Multiplicity Proton");
    fh1_MultiplicityProton->GetYaxis()->SetTitleOffset(1.4);
    fh1_MultiplicityProton->GetXaxis()->CenterTitle(true);
    fh1_MultiplicityProton->GetYaxis()->CenterTitle(true);
    fh1_MultiplicityProton->Draw("col");


    fh1_EnergyCalifaCrystals.resize(fNbCrystals);
    for (Int_t i=0; i<fNbCrystals; i++)
    {
      if (i==0) {sprintf(Name1,  "fh1_energyBoxA_1_p"); }
      if (i==1) {sprintf(Name1,  "fh1_energyBoxA_2_p"); }
      if (i==2) {sprintf(Name1,  "fh1_energyBoxA_3_p"); }
      if (i==3) {sprintf(Name1,  "fh1_energyBoxA_4_p"); }
      if (i==4) {sprintf(Name1,  "fh1_energyBoxA_1_g"); }
      if (i==5) {sprintf(Name1,  "fh1_energyBoxA_2_g"); }
      if (i==6) {sprintf(Name1,  "fh1_energyBoxA_3_g"); }
      if (i==7) {sprintf(Name1,  "fh1_energyBoxA_4_g"); }
      if (i==8) {sprintf(Name1,  "fh1_energyBoxB_1_p"); }
      if (i==9) {sprintf(Name1,  "fh1_energyBoxB_2_p"); }
      if (i==10) {sprintf(Name1, "fh1_energyBoxB_3_p"); }
      if (i==11) {sprintf(Name1, "fh1_energyBoxB_4_p"); }
      if (i==12) {sprintf(Name1, "fh1_energyBoxB_1_g"); }
      if (i==13) {sprintf(Name1, "fh1_energyBoxB_2_g"); }
      if (i==14) {sprintf(Name1, "fh1_energyBoxB_3_g"); }
      if (i==15) {sprintf(Name1, "fh1_energyBoxB_4_g"); }
      //sprintf(Name1, "fh1_energy_CalifaCrystal_%d", i + 1);
      sprintf(Name2, "Energy in Califa crystal: %d", i + 1);
      fh1_EnergyCalifaCrystals[i] = new TH1F(Name1, Name2, 3000, 0, 30000);
      fh1_EnergyCalifaCrystals[i]->GetXaxis()->SetTitle("Energy [channels]");
      fh1_EnergyCalifaCrystals[i]->GetYaxis()->SetTitle("counts");
      fh1_EnergyCalifaCrystals[i]->GetYaxis()->SetTitleOffset(1.4);
      fh1_EnergyCalifaCrystals[i]->GetXaxis()->CenterTitle(true);
      fh1_EnergyCalifaCrystals[i]->GetYaxis()->CenterTitle(true);
      fh1_EnergyCalifaCrystals[i]->Draw("col");

    }


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
    }


    // Cal data CALIFA
    fh1_EnergyCalCalifaCrystals.resize(fNbCrystals);
    cEnergyCalCalifa_BoxA_p = new TCanvas("EnergyCalCalifa_BoxA_p", "EnergyCalCalifa_BoxA_p", 10, 10, 500, 500);
    cEnergyCalCalifa_BoxB_g = new TCanvas("EnergyCalCalifa_BoxB_g", "EnergyCalCalifa_BoxA_g", 10, 10, 500, 500);
    cEnergyCalCalifa_BoxA_g = new TCanvas("EnergyCalCalifa_BoxA_g", "EnergyCalCalifa_BoxA_g", 10, 10, 500, 500);
    cEnergyCalCalifa_BoxB_p = new TCanvas("EnergyCalCalifa_BoxB_p", "EnergyCalCalifa_BoxB_p", 10, 10, 500, 500);
    for (Int_t i=0; i<fNbCrystals; i++)
    {
      sprintf(Name2, "Energy calibrated in Califa crystal: %d", i + 1);
      Int_t pad = 0;
      if (i<4)
      {
        if (i==0) {sprintf(Name1,  "fh1_energyCalBoxA_1_p"); pad=2;}
        if (i==1) {sprintf(Name1,  "fh1_energyCalBoxA_2_p"); pad=4;}
        if (i==2) {sprintf(Name1,  "fh1_energyCalBoxA_3_p"); pad=3;}
        if (i==3) {sprintf(Name1,  "fh1_energyCalBoxA_4_p"); pad=1;}

        cEnergyCalCalifa_BoxA_p->Divide(2,2);
        cEnergyCalCalifa_BoxA_p->cd(pad);
        fh1_EnergyCalCalifaCrystals[i] = new TH1F(Name1, Name2, 30000, 0, 300000);
        fh1_EnergyCalCalifaCrystals[i]->Draw("col");
      }

      if (i<7 && i>=4)
      {
        if (i==4) {sprintf(Name1,  "fh1_energyCalBoxA_1_g"); pad=2; }
        if (i==5) {sprintf(Name1,  "fh1_energyCalBoxA_2_g"); pad=4; }
        if (i==6) {sprintf(Name1,  "fh1_energyCalBoxA_3_g"); pad=3; }
        if (i==7) {sprintf(Name1,  "fh1_energyCalBoxA_4_g"); pad=1; }

        cEnergyCalCalifa_BoxA_g->Divide(2,2);
        cEnergyCalCalifa_BoxA_g->cd(pad);
        fh1_EnergyCalCalifaCrystals[i] = new TH1F(Name1, Name2, 3000, 0, 3000);
        fh1_EnergyCalCalifaCrystals[i]->Draw("col");
      }

      if (i<11 && i>=7)
      {
        if (i==8) {sprintf(Name1,  "fh1_energyCalBoxB_1_p"); pad=2; }
        if (i==9) {sprintf(Name1,  "fh1_energyCalBoxB_2_p"); pad=4; }
        if (i==10) {sprintf(Name1, "fh1_energyCalBoxB_3_p"); pad=3; }
        if (i==11) {sprintf(Name1, "fh1_energyCalBoxB_4_p"); pad=1; }

        cEnergyCalCalifa_BoxB_p->Divide(2,2);
        cEnergyCalCalifa_BoxB_p->cd(pad);
        fh1_EnergyCalCalifaCrystals[i] = new TH1F(Name1, Name2, 30000, 0, 300000);
        fh1_EnergyCalCalifaCrystals[i]->Draw("col");
      }

      if (i>=11)
      {
        if (i==12) {sprintf(Name1, "fh1_energyCalBoxB_1_g"); pad=2;}
        if (i==13) {sprintf(Name1, "fh1_energyCalBoxB_2_g"); pad=4;}
        if (i==14) {sprintf(Name1, "fh1_energyCalBoxB_3_g"); pad=3;}
        if (i==15) {sprintf(Name1, "fh1_energyCalBoxB_4_g"); pad=1;}

       cEnergyCalCalifa_BoxB_g->Divide(2,2);
       cEnergyCalCalifa_BoxB_g->cd(pad);
       fh1_EnergyCalCalifaCrystals[i] = new TH1F(Name1, Name2, 3000, 0, 3000);
       }

      //sprintf(Name1, "fh1_energyCal_CalifaCrystal_%d", i + 1);
      fh1_EnergyCalCalifaCrystals[i]->GetXaxis()->SetTitle("Energy [keV]");
      fh1_EnergyCalCalifaCrystals[i]->GetYaxis()->SetTitle("counts");
      fh1_EnergyCalCalifaCrystals[i]->GetYaxis()->SetTitleOffset(1.4);
      fh1_EnergyCalCalifaCrystals[i]->GetXaxis()->CenterTitle(true);
      fh1_EnergyCalCalifaCrystals[i]->GetYaxis()->CenterTitle(true);
      //fh1_EnergyCalCalifaCrystals[i]->Draw("col");

    }


    // energy correlations CALIFA
    fh2_EnergyCorrelationsCrystals.resize(24);
    for (Int_t i=0; i<24; i++)
    {
      if (i==0) {sprintf(Name1, "fh2_EnergyCorr_BoxA_1_2_p"); }
      if (i==1) {sprintf(Name1, "fh2_EnergyCorr_BoxA_1_3_p"); }
      if (i==2) {sprintf(Name1, "fh2_EnergyCorr_BoxA_1_4_p"); }
      if (i==3) {sprintf(Name1, "fh2_EnergyCorr_BoxA_2_3_p"); }
      if (i==4) {sprintf(Name1, "fh2_EnergyCorr_BoxA_2_4_p"); }
      if (i==5) {sprintf(Name1, "fh2_EnergyCorr_BoxA_3_4_p"); }
      if (i==6) {sprintf(Name1, "fh2_EnergyCorr_BoxA_1_2_g"); }
      if (i==7) {sprintf(Name1, "fh2_EnergyCorr_BoxA_1_3_g"); }
      if (i==8) {sprintf(Name1, "fh2_EnergyCorr_BoxA_1_4_g"); }
      if (i==9) {sprintf(Name1, "fh2_EnergyCorr_BoxA_2_4_g"); }
      if (i==10) {sprintf(Name1, "fh2_EnergyCorr_BoxA_3_4_g"); }
      if (i==11) {sprintf(Name1, "fh2_EnergyCorr_BoxA_2_3_g"); }

      if (i==12) {sprintf(Name1, "fh2_EnergyCorr_BoxB_1_2_p"); }
      if (i==13) {sprintf(Name1, "fh2_EnergyCorr_BoxB_1_3_p"); }
      if (i==14) {sprintf(Name1, "fh2_EnergyCorr_BoxB_1_4_p"); }
      if (i==15) {sprintf(Name1, "fh2_EnergyCorr_BoxB_2_3_p"); }
      if (i==16) {sprintf(Name1, "fh2_EnergyCorr_BoxB_2_4_p"); }
      if (i==17) {sprintf(Name1, "fh2_EnergyCorr_BoxB_3_4_p"); }
      if (i==18) {sprintf(Name1, "fh2_EnergyCorr_BoxB_1_2_g"); }
      if (i==19) {sprintf(Name1, "fh2_EnergyCorr_BoxB_1_3_g"); }
      if (i==20) {sprintf(Name1, "fh2_EnergyCorr_BoxB_1_4_g"); }
      if (i==21) {sprintf(Name1, "fh2_EnergyCorr_BoxB_2_4_g"); }
      if (i==22) {sprintf(Name1, "fh2_EnergyCorr_BoxB_3_4_g"); }
      if (i==23) {sprintf(Name1, "fh2_EnergyCorr_BoxB_2_3_g"); }

      sprintf(Name2, "fh2_EnergyCorrelationsCrystals_%d", i);
      if (i<5 || (i>11 && i<18)) {fh2_EnergyCorrelationsCrystals[i] = new TH2F(Name1, Name2, 3000, 0, 300000, 3000,0,300000);}
      else {fh2_EnergyCorrelationsCrystals[i] = new TH2F(Name1, Name2, 3000, 0, 30000, 3000,0,30000);}
      fh2_EnergyCorrelationsCrystals[i]->GetYaxis()->SetTitleOffset(1.4);
      fh2_EnergyCorrelationsCrystals[i]->GetXaxis()->CenterTitle(true);
      fh2_EnergyCorrelationsCrystals[i]->GetYaxis()->CenterTitle(true);
      fh2_EnergyCorrelationsCrystals[i]->Draw("colz");

    }

    fh2_EnergyCorrelationsAlvProton = new TH2F("fh2_EnergyCorrelationsAlvProton", "Energies Alv ProtonRange", 300, 0, 30000, 300,0,30000);
    fh2_EnergyCorrelationsAlvProton->GetXaxis()->SetTitle("Energy Alv1");
    fh2_EnergyCorrelationsAlvProton->GetYaxis()->SetTitle("Energy Alv2");
    fh2_EnergyCorrelationsAlvProton->GetYaxis()->SetTitleOffset(1.4);
    fh2_EnergyCorrelationsAlvProton->GetXaxis()->CenterTitle(true);
    fh2_EnergyCorrelationsAlvProton->GetYaxis()->CenterTitle(true);
    fh2_EnergyCorrelationsAlvProton->Draw("col");


    fh2_EnergyCorrelationsAlvGamma = new TH2F("fh2_EnergyCorrelationsAlvGamma", "Energies Alv GammaRange", 300, 0, 30000, 300,0,30000);
    fh2_EnergyCorrelationsAlvGamma->GetXaxis()->SetTitle("Energy Alv1");
    fh2_EnergyCorrelationsAlvGamma->GetYaxis()->SetTitle("Energy Alv2");
    fh2_EnergyCorrelationsAlvGamma->GetYaxis()->SetTitleOffset(1.4);
    fh2_EnergyCorrelationsAlvGamma->GetXaxis()->CenterTitle(true);
    fh2_EnergyCorrelationsAlvGamma->GetYaxis()->CenterTitle(true);
    fh2_EnergyCorrelationsAlvGamma->Draw("COLZ");

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
    }

    cBoxEnergy = new TCanvas("BoxEnergy", "BoxEnergy", 10, 10, 500, 500);
    cBoxEnergy->Divide(2,2);
    fh1_EnergyBoxA_g = new TH1F("fh1_EnergyBoxA_g", "fh1_EnergyBoxA_g",3000,0,30000);
    cBoxEnergy->cd(1);
    fh1_EnergyBoxA_g->Draw("col");
    fh1_EnergyBoxA_p = new TH1F("fh1_EnergyBoxA_p", "fh1_EnergyBoxA_p",3000,0,300000);
    cBoxEnergy->cd(2);
    fh1_EnergyBoxA_p->Draw("col");
    fh1_EnergyBoxB_g = new TH1F("fh1_EnergyBoxB_g", "fh1_EnergyBoxB_g",3000,0,30000);
    cBoxEnergy->cd(3);
    fh1_EnergyBoxB_g->Draw("col");
    fh1_EnergyBoxB_p = new TH1F("fh1_EnergyBoxB_p", "fh1_EnergyBoxB_p",3000,0,300000);
    cBoxEnergy->cd(4);
    fh1_EnergyBoxB_p->Draw("col");

    // MAIN FOLDERs
    TFolder* mainfolCalifa = new TFolder("Califa", "Califa info");

    mainfolCalifa->Add(cEnergyCalCalifa_BoxA_p);
    mainfolCalifa->Add(cEnergyCalCalifa_BoxA_g);
    mainfolCalifa->Add(cEnergyCalCalifa_BoxB_p);
    mainfolCalifa->Add(cEnergyCalCalifa_BoxB_g);
    mainfolCalifa->Add(cBoxEnergy);

    TFolder* mainfolSi = new TFolder("Si", "Si info");


    // Folders for mapped data
    TFolder* mapfolCalifa = new TFolder("Mapped", "Califa Mapped");
    mapfolCalifa->Add(fh1_MultiplicityProton);
    mapfolCalifa->Add(fh1_MultiplicityGamma);
    for (Int_t i=0;i<fNbCrystals;i++)
    {
      mapfolCalifa->Add(fh1_EnergyCalifaCrystals[i]);
    }
    mainfolCalifa->Add(mapfolCalifa);

    TFolder* mapfolSi = new TFolder("MappedSi", "Si Mapped");
    for (Int_t i = 0; i < fNbDet; i++)
    {
      mapfolSi->Add(fh2_EnergyVsStrip[i]);
    }
    mainfolSi->Add(mapfolSi);

    // Folders for cal data
    TFolder* calfolCalifa = new TFolder("Cal", "Cal Califa info");
    for (Int_t i=0;i<fNbCrystals;i++)
    {
      calfolCalifa->Add(fh1_EnergyCalCalifaCrystals[i]);
    }
    mainfolCalifa->Add(calfolCalifa);

    TFolder* calfolSi = new TFolder("CalSi", "Cal Si info");
    for (Int_t i = 0; i < fNbDet; i++)
    {
      calfolSi->Add(fh2_EnergyCalVsStrip[i]);
    }
    mainfolSi->Add(calfolSi);

    // Folders for hit data and correlations
    TFolder* hitfolCalifa = new TFolder("Correlations", "Hit Califa info");
    hitfolCalifa->Add(fh2_EnergyCorrelationsAlvGamma);
    hitfolCalifa->Add(fh2_EnergyCorrelationsAlvProton);
    for (Int_t i=0;i<24;i++)
    {
      hitfolCalifa->Add(fh2_EnergyCorrelationsCrystals[i]);
    }
    mainfolCalifa->Add(hitfolCalifa);

    TFolder* hitfolSi = new TFolder("HitSi", "Hit Si info");
    for (Int_t i = 0; i < fNbDet; i++)
    {
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

void R3BCalifaJulichOnlineSpectraCANVASDEV::Reset_CalifaJulich_Histo()
{
    LOG(INFO) << "R3BCalifaJulichOnlineSpectraCANVASDEV::Reset_CalifaJulich_Histo";

    // Mapped data
    fh1_MultiplicityGamma->Reset();
    fh1_MultiplicityProton->Reset();

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

    for (Int_t i=0;i<24;i++) { fh2_EnergyCorrelationsCrystals[i] -> Reset(); }
}

void R3BCalifaJulichOnlineSpectraCANVASDEV::Exec(Option_t* option)
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
        Int_t multGamma=0; Int_t multProton=0;
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BCalifaMappedData* hit = (R3BCalifaMappedData*)fMappedItemsCalifa->At(ihit);
            if (!hit)
                continue;
            Int_t index=1000;
            if (hit->GetCrystalId()==324) {index=3;}
            if (hit->GetCrystalId()==325) {index=2;}
            if (hit->GetCrystalId()==326) {index=1;}
            if (hit->GetCrystalId()==327) {index=0;}

            if (hit->GetCrystalId()==340) {index=7;}
            if (hit->GetCrystalId()==341) {index=6;}
            if (hit->GetCrystalId()==342) {index=5;}
            if (hit->GetCrystalId()==343) {index=4;}

            if (hit->GetCrystalId()==356) {index=11;}
            if (hit->GetCrystalId()==357) {index=10;}
            if (hit->GetCrystalId()==358) {index=9;}
            if (hit->GetCrystalId()==359) {index=8;}

            if (hit->GetCrystalId()==372) {index=15;}
            if (hit->GetCrystalId()==373) {index=14;}
            if (hit->GetCrystalId()==374) {index=13;}
            if (hit->GetCrystalId()==375) {index=12;}

            if(index==1000) {return;}
            if (index<4 || (index>7 && index<12)) {multProton++;}
            else {multGamma++;}

            fh1_EnergyCalifaCrystals[index]->Fill(hit->GetEnergy());
        }

        fh1_MultiplicityGamma->Fill(multGamma);
        fh1_MultiplicityProton->Fill(multGamma);
    }

    //Fill Cal data
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
        Int_t crystals[16];
        Float_t energies[16];
        Float_t energy_alv[4]; for (Int_t i=0;i<4;i++) {energy_alv[i]=0.0;}
        for (Int_t i=0;i<16;i++) {crystals[i]=0; energies[i]=0.0;}

        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BCalifaCrystalCalData* hit = (R3BCalifaCrystalCalData*)fCalItemsCalifa->At(ihit);
            if (!hit)
                continue;
            Int_t index=1000;

            if (hit->GetCrystalId()==324) {index=3;}
            if (hit->GetCrystalId()==325) {index=2;}
            if (hit->GetCrystalId()==326) {index=1;}
            if (hit->GetCrystalId()==327) {index=0;}

            if (hit->GetCrystalId()==340) {index=7;}
            if (hit->GetCrystalId()==341) {index=6;}
            if (hit->GetCrystalId()==342) {index=5;}
            if (hit->GetCrystalId()==343) {index=4;}

            if (hit->GetCrystalId()==356) {index=11;}
            if (hit->GetCrystalId()==357) {index=10;}
            if (hit->GetCrystalId()==358) {index=9;}
            if (hit->GetCrystalId()==359) {index=8;}

            if (hit->GetCrystalId()==372) {index=15;}
            if (hit->GetCrystalId()==373) {index=14;}
            if (hit->GetCrystalId()==374) {index=13;}
            if (hit->GetCrystalId()==375) {index=12;}

            if (index==1000) {return;}

            fh1_EnergyCalCalifaCrystals[index]->Fill(hit->GetEnergy());
            crystals[index]=1;

            energies[index]=hit->GetEnergy();
        }


        for (Int_t i=0;i<4;i++)
        {
          energy_alv[0] =+ energies[i]; energy_alv[1] =+ energies[i+4];
          energy_alv[2] =+ energies[i+8]; energy_alv[3] =+ energies[i+12];
        }

        // BoxA ProtonRange
        if (crystals[0] && crystals[1])
        {
          fh2_EnergyCorrelationsCrystals[0]->Fill(energies[0],energies[1]);
          fh2_EnergyCorrelationsCrystals[0]->GetXaxis()->SetTitle("Energy [keV] crystal 1");
          fh2_EnergyCorrelationsCrystals[0]->GetYaxis()->SetTitle("Energy [keV] crystal 2");
        }

        if (crystals[0] && crystals[2])
        {
          fh2_EnergyCorrelationsCrystals[1]->Fill(energies[0],energies[2]);
          fh2_EnergyCorrelationsCrystals[1]->GetXaxis()->SetTitle("Energy [keV] crystal 1");
          fh2_EnergyCorrelationsCrystals[1]->GetYaxis()->SetTitle("Energy [keV] crystal 3");
        }

        if (crystals[0] && crystals[3])
        {
          fh2_EnergyCorrelationsCrystals[2]->Fill(energies[0],energies[3]);
          fh2_EnergyCorrelationsCrystals[2]->GetXaxis()->SetTitle("Energy [keV] crystal 1");
          fh2_EnergyCorrelationsCrystals[2]->GetYaxis()->SetTitle("Energy [keV] crystal 4");
        }

        if (crystals[1] && crystals[2])
        {
          fh2_EnergyCorrelationsCrystals[3]->Fill(energies[1],energies[2]);
          fh2_EnergyCorrelationsCrystals[3]->GetXaxis()->SetTitle("Energy [keV] crystal 2");
          fh2_EnergyCorrelationsCrystals[3]->GetYaxis()->SetTitle("Energy [keV] crystal 3");
        }

        if (crystals[1] && crystals[3])
        {
          fh2_EnergyCorrelationsCrystals[4]->Fill(energies[1],energies[3]);
          fh2_EnergyCorrelationsCrystals[4]->GetXaxis()->SetTitle("Energy [keV] crystal 2");
          fh2_EnergyCorrelationsCrystals[4]->GetYaxis()->SetTitle("Energy [keV] crystal 4");
        }

        if (crystals[2] && crystals[3])
        {
          fh2_EnergyCorrelationsCrystals[5]->Fill(energies[2],energies[3]);
          fh2_EnergyCorrelationsCrystals[5]->GetXaxis()->SetTitle("Energy [keV] crystal 3");
          fh2_EnergyCorrelationsCrystals[5]->GetYaxis()->SetTitle("Energy [keV] crystal 4");
        }

//BoxA GammaRange

      if (crystals[4] && crystals[5])
        {
          fh2_EnergyCorrelationsCrystals[6]->Fill(energies[4],energies[5]);
          fh2_EnergyCorrelationsCrystals[6]->GetXaxis()->SetTitle("Energy [keV] crystal 5");
          fh2_EnergyCorrelationsCrystals[6]->GetYaxis()->SetTitle("Energy [keV] crystal 6");
        }

        if (crystals[4] && crystals[6])
        {
          fh2_EnergyCorrelationsCrystals[7]->Fill(energies[4],energies[6]);
          fh2_EnergyCorrelationsCrystals[7]->GetXaxis()->SetTitle("Energy [keV] crystal 5");
          fh2_EnergyCorrelationsCrystals[7]->GetYaxis()->SetTitle("Energy [keV] crystal 7");
        }

        if (crystals[4] && crystals[7])
        {
          fh2_EnergyCorrelationsCrystals[8]->Fill(energies[4],energies[7]);
          fh2_EnergyCorrelationsCrystals[8]->GetXaxis()->SetTitle("Energy [keV] crystal 5");
          fh2_EnergyCorrelationsCrystals[8]->GetYaxis()->SetTitle("Energy [keV] crystal 8");
        }

        if (crystals[5] && crystals[7])
        {
          fh2_EnergyCorrelationsCrystals[9]->Fill(energies[5],energies[7]);
          fh2_EnergyCorrelationsCrystals[9]->GetXaxis()->SetTitle("Energy [keV] crystal 6");
          fh2_EnergyCorrelationsCrystals[9]->GetYaxis()->SetTitle("Energy [keV] crystal 8");
        }

        if (crystals[6] && crystals[7])
        {
          fh2_EnergyCorrelationsCrystals[10]->Fill(energies[6],energies[7]);
          fh2_EnergyCorrelationsCrystals[10]->GetXaxis()->SetTitle("Energy [keV] crystal 7");
          fh2_EnergyCorrelationsCrystals[10]->GetYaxis()->SetTitle("Energy [keV] crystal 8");
        }
        if (crystals[5] && crystals[6])
        {
          fh2_EnergyCorrelationsCrystals[11]->Fill(energies[5],energies[6]);
          fh2_EnergyCorrelationsCrystals[11]->GetXaxis()->SetTitle("Energy [keV] crystal 6");
          fh2_EnergyCorrelationsCrystals[11]->GetYaxis()->SetTitle("Energy [keV] crystal 7");
        }

//BoxB ProtonRange
        if (crystals[8] && crystals[9])
        {
          fh2_EnergyCorrelationsCrystals[12]->Fill(energies[8],energies[9]);
          fh2_EnergyCorrelationsCrystals[12]->GetXaxis()->SetTitle("Energy [keV] crystal 1");
          fh2_EnergyCorrelationsCrystals[12]->GetYaxis()->SetTitle("Energy [keV] crystal 2");
        }

        if (crystals[8] && crystals[10])
        {
          fh2_EnergyCorrelationsCrystals[13]->Fill(energies[8],energies[10]);
          fh2_EnergyCorrelationsCrystals[13]->GetXaxis()->SetTitle("Energy [keV] crystal 1");
          fh2_EnergyCorrelationsCrystals[13]->GetYaxis()->SetTitle("Energy [keV] crystal 3");
        }

        if (crystals[8] && crystals[11])
        {
          fh2_EnergyCorrelationsCrystals[14]->Fill(energies[8],energies[11]);
          fh2_EnergyCorrelationsCrystals[14]->GetXaxis()->SetTitle("Energy [keV] crystal 1");
          fh2_EnergyCorrelationsCrystals[14]->GetYaxis()->SetTitle("Energy [keV] crystal 4");
        }

        if (crystals[9] && crystals[10])
        {
          fh2_EnergyCorrelationsCrystals[15]->Fill(energies[9],energies[10]);
          fh2_EnergyCorrelationsCrystals[15]->GetXaxis()->SetTitle("Energy [keV] crystal 2");
          fh2_EnergyCorrelationsCrystals[15]->GetYaxis()->SetTitle("Energy [keV] crystal 3");
        }

        if (crystals[9] && crystals[11])
        {
          fh2_EnergyCorrelationsCrystals[16]->Fill(energies[9],energies[11]);
          fh2_EnergyCorrelationsCrystals[16]->GetXaxis()->SetTitle("Energy [keV] crystal 2");
          fh2_EnergyCorrelationsCrystals[16]->GetYaxis()->SetTitle("Energy [keV] crystal 4");
        }

        if (crystals[10] && crystals[11])
        {
          fh2_EnergyCorrelationsCrystals[17]->Fill(energies[10],energies[11]);
          fh2_EnergyCorrelationsCrystals[17]->GetXaxis()->SetTitle("Energy [keV] crystal 3");
          fh2_EnergyCorrelationsCrystals[17]->GetYaxis()->SetTitle("Energy [keV] crystal 4");
        }

//BoxB GammaRange

        if (crystals[12] && crystals[13])
        {
          fh2_EnergyCorrelationsCrystals[18]->Fill(energies[12],energies[13]);
          fh2_EnergyCorrelationsCrystals[18]->GetXaxis()->SetTitle("Energy [keV] crystal 5");
          fh2_EnergyCorrelationsCrystals[18]->GetYaxis()->SetTitle("Energy [keV] crystal 6");
        }

        if (crystals[12] && crystals[14])
        {
          fh2_EnergyCorrelationsCrystals[19]->Fill(energies[12],energies[14]);
          fh2_EnergyCorrelationsCrystals[19]->GetXaxis()->SetTitle("Energy [keV] crystal 5");
          fh2_EnergyCorrelationsCrystals[19]->GetYaxis()->SetTitle("Energy [keV] crystal 7");
        }

        if (crystals[12] && crystals[15])
        {
          fh2_EnergyCorrelationsCrystals[20]->Fill(energies[12],energies[15]);
          fh2_EnergyCorrelationsCrystals[20]->GetXaxis()->SetTitle("Energy [keV] crystal 5");
          fh2_EnergyCorrelationsCrystals[20]->GetYaxis()->SetTitle("Energy [keV] crystal 8");
        }

        if (crystals[13] && crystals[15])
        {
          fh2_EnergyCorrelationsCrystals[21]->Fill(energies[13],energies[15]);
          fh2_EnergyCorrelationsCrystals[21]->GetXaxis()->SetTitle("Energy [keV] crystal 6");
          fh2_EnergyCorrelationsCrystals[21]->GetYaxis()->SetTitle("Energy [keV] crystal 8");
        }

        if (crystals[14] && crystals[15])
        {
          fh2_EnergyCorrelationsCrystals[22]->Fill(energies[14],energies[15]);
          fh2_EnergyCorrelationsCrystals[22]->GetXaxis()->SetTitle("Energy [keV] crystal 7");
          fh2_EnergyCorrelationsCrystals[22]->GetYaxis()->SetTitle("Energy [keV] crystal 8");
        }
        if (crystals[13] && crystals[14])
        {
          fh2_EnergyCorrelationsCrystals[23]->Fill(energies[13],energies[14]);
          fh2_EnergyCorrelationsCrystals[23]->GetXaxis()->SetTitle("Energy [keV] crystal 6");
          fh2_EnergyCorrelationsCrystals[23]->GetYaxis()->SetTitle("Energy [keV] crystal 7");
        }

      fh2_EnergyCorrelationsAlvProton->Fill(energy_alv[0],energy_alv[2]);
      fh2_EnergyCorrelationsAlvGamma->Fill(energy_alv[1],energy_alv[3]);


      cBoxEnergy->cd(1); fh1_EnergyBoxA_p->Fill(energy_alv[0]);
      cBoxEnergy->cd(2); fh1_EnergyBoxA_g->Fill(energy_alv[1]);
      cBoxEnergy->cd(3); fh1_EnergyBoxB_p->Fill(energy_alv[2]);
      cBoxEnergy->cd(4); fh1_EnergyBoxB_g->Fill(energy_alv[3]);

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

void R3BCalifaJulichOnlineSpectraCANVASDEV::FinishEvent()
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

void R3BCalifaJulichOnlineSpectraCANVASDEV::FinishTask()
{

   for (Int_t i = 0; i < fNbDet; i++)
      {
          if (fMappedItemsSi) { fh2_EnergyVsStrip[i]->Write(); }
          if (fCalItemsSi) { fh2_EnergyCalVsStrip[i]->Write(); }
          if (fHitItemsSi) { fh2_EnergyHitVsStrip[i]->Write(); }
      }

      if(fMappedItemsCalifa) fh1_MultiplicityGamma->Write();

      for (Int_t i = 0; i < fNbCrystals; i++)
      {
          if (fMappedItemsCalifa){ fh1_EnergyCalifaCrystals[i]->Write(); }
          if (fCalItemsCalifa){ fh1_EnergyCalCalifaCrystals[i]->Write(); }
          if (fHitItemsCalifa){ fh1_EnergyHitCalifaCrystals[i]->Write(); }
      }

      if (fCalItemsCalifa)
      {
        for (Int_t i=0;i<24;i++)
        {
          fh2_EnergyCorrelationsCrystals[i]->Write();
        }
      }

      fh2_EnergyCorrelationsAlvProton->Write();
      fh2_EnergyCorrelationsAlvGamma->Write();

      cBoxEnergy->Write();
      fh1_EnergyBoxA_p->Write();
      fh1_EnergyBoxA_g->Write();
      fh1_EnergyBoxB_p->Write();
      fh1_EnergyBoxB_g->Write();

      cEnergyCalCalifa_BoxA_p->Write();
      cEnergyCalCalifa_BoxA_g->Write();
      cEnergyCalCalifa_BoxB_p->Write();
      cEnergyCalCalifa_BoxB_g->Write();
}

ClassImp(R3BCalifaJulichOnlineSpectraCANVASDEV);

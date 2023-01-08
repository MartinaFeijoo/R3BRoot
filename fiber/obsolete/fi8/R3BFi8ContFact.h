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

#ifndef R3BFI8CONTFACT_H
#define R3BFI8CONTFACT_H

#include "FairContFact.h"

class FairContainer;

class R3BFi8ContFact : public FairContFact
{
  private:
    void setAllContainers();

  public:
    R3BFi8ContFact();
    ~R3BFi8ContFact() {}
    FairParSet* createContainer(FairContainer*);
    void activateParIo(FairParIo* io);
    ClassDef(R3BFi8ContFact, 0) // Factory for all Fi4 parameter containers
};

#endif /* !R3BFI8CONTFACT_H */

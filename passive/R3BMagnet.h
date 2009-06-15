#ifndef R3BMAGNET_H
#define R3BMAGNET_H

#include "TNamed.h"
#include "TArrayI.h"
#include "TClonesArray.h"
#include "FairDetector.h"
#include "FairModule.h"
#include "FairModule.h"
#include "TGeoMatrix.h"

class R3BMagnet : public FairModule {

private:
   TGeoCombiTrans *gLobalPos;

public:
    R3BMagnet(const char * name, const char *Title="R3B Magnet");
    R3BMagnet();
    virtual ~R3BMagnet();
    void ConstructGeometry();
    void ConstructASCIIGeometry();
    Bool_t CheckIfSensitive(std::string name);
    ClassDef(R3BMagnet,1) //R3BMagnet




};

#endif //R3BMAGNET_H


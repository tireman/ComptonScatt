//******************************************************************************
// LICENSE DISCLAIMER:
// ===================

// Permission is granted to the public to copy and use this software
// without charge, provided that this Notice and any statement of
// authorship are reproduced on all copies.  Neither the Authors nor 
// Northern Michigan University makes any warranty, express or implied, 
// or assumes any liability or responsibility for the use of this software.
//
// SteppingAction.hh
//
//******************************************************************************
//

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class DetectorConstruction;
class EventAction;

class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction(DetectorConstruction* myDC,EventAction* myEA);
    virtual ~SteppingAction();

    virtual void UserSteppingAction(const G4Step*);
    
  private:
    DetectorConstruction* myDetector;
    EventAction* eventAction;
};

#endif

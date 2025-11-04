//******************************************************************************
// LICENSE DISCLAIMER:
// ===================

// Permission is granted to the public to copy and use this software
// without charge, provided that this Notice and any statement of
// authorship are reproduced on all copies.  Neither the Authors nor 
// Northern Michigan University makes any warranty, express or implied, 
// or assumes any liability or responsibility for the use of this software.
//
// SteppingAction.cc
//
//******************************************************************************
//

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "G4SteppingManager.hh"

SteppingAction::SteppingAction(DetectorConstruction* myDC, EventAction* myEA)
:myDetector(myDC), eventAction(myEA)
{ }

SteppingAction::~SteppingAction()
{ }

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // Accumulate energy deposits in detector volumes
  //
  G4double deStep = aStep->GetTotalEnergyDeposit();
  if (deStep > 0.) {
    G4String curVolName = 
      aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();

    if (curVolName == "D1") {
      eventAction->addD1dep(deStep);
    } else if (curVolName == "Target") {
      eventAction->addD2dep(deStep);
    } else if (curVolName == "detShield")
      eventAction->addD3dep(deStep);
  }
}

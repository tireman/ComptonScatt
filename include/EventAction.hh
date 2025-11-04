//******************************************************************************
// LICENSE DISCLAIMER:
// ===================

// Permission is granted to the public to copy and use this software
// without charge, provided that this Notice and any statement of
// authorship are reproduced on all copies.  Neither the Authors nor 
// Northern Michigan University makes any warranty, express or implied, 
// or assumes any liability or responsibility for the use of this software.
//
// EventAction.hh
//
//******************************************************************************
//

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;
class RunAction;

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
    ~EventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

    void addD1dep(G4double Edep) {D1dep += Edep;};
    void addD2dep(G4double Edep) {D2dep += Edep;}; 
    void addD3dep(G4double Edep) {D3dep += Edep;};
  

  private:
  G4double D1dep, D2dep, D3dep;
  
};

#endif

    

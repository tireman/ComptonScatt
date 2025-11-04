//******************************************************************************
// LICENSE DISCLAIMER:
// ===================

// Permission is granted to the public to copy and use this software
// without charge, provided that this Notice and any statement of
// authorship are reproduced on all copies.  Neither the Authors nor 
// Northern Michigan University makes any warranty, express or implied, 
// or assumes any liability or responsibility for the use of this software.
//
// RunManager.hh
//
//******************************************************************************
//

#ifndef RunManager_h
#define RunManager_h 1

#include "G4RunManager.hh"
#include "globals.hh"

class RunManager: public G4RunManager
{
public:
  RunManager();
  ~RunManager();

public:
  static RunManager* GetRunManager(); // Returns pointer to single instance.
  
  G4int GetnIncident(){ return nIncidentEvents; }
  
  G4bool GetEventTrigger(){ return eventTrigger; }
  
  void SetEventTrigger(){ eventTrigger = true; }

protected:
  static RunManager* fRunManager; // Pointer to single instance.
  
  virtual void DoEventLoop(G4int n_event, const char* macroFile=0,
			   G4int n_select=-1);

private:
  G4int nIncidentEvents;
  G4bool eventTrigger;
};

#endif


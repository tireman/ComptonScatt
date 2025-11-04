//******************************************************************************
// LICENSE DISCLAIMER:
// ===================

// Permission is granted to the public to copy and use this software
// without charge, provided that this Notice and any statement of
// authorship are reproduced on all copies.  Neither the Authors nor 
// Northern Michigan University makes any warranty, express or implied, 
// or assumes any liability or responsibility for the use of this software.
//
// RunManager.cc
//
//******************************************************************************
//


// On Sun, to prevent conflict with ObjectSpace, G4Timer.hh has to be
// loaded *before* globals.hh...
#include "G4Timer.hh"

#include "RunManager.hh"

#include "Randomize.hh"
#include "G4Run.hh"
#include "G4RunMessenger.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPhysicsList.hh"
#include "G4UserRunAction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeometryManager.hh"
#include "G4SDManager.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ApplicationState.hh"
#include "G4StateManager.hh"
#include "G4VPersistencyManager.hh"
#include "G4UImanager.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessTable.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"

#include "G4ios.hh"
#include <sstream>

//----------------------------------------------------------------------------//
RunManager* RunManager::fRunManager = 0;

//----------------------------------------------------------------------------//
RunManager::RunManager()
{
  if (fRunManager) { 
    G4Exception("G4RunManager","RunManager already Started",RunMustBeAborted,"Run Failed\n Another Run Manager is open"); 
  } else {
    fRunManager = this;
  }
}

//----------------------------------------------------------------------------//
RunManager::~RunManager()
{}

//----------------------------------------------------------------------------//
RunManager* RunManager::GetRunManager()
{ return fRunManager; }

//----------------------------------------------------------------------------//
// Override this function to loop over TRIGGERED events, not incident
// particles.
//----------------------------------------------------------------------------//
void RunManager::DoEventLoop(G4int n_event,const char* macroFile,G4int n_select)
{
  if(verboseLevel>0) 
  { timer->Start(); }

  G4String msg;
  if(macroFile!=0)
  { 
    if(n_select<0) n_select = n_event;
    msg = "/control/execute ";
    msg += macroFile;
  }
  else
  { n_select = -1; }

  nIncidentEvents = 0;      // Initilaize incident event counter for this run

  G4int i_event;
  for( i_event=0; i_event<n_event; i_event++ ) {

        eventTrigger = false;       // Keep generating until a trigger.
    do {                        // Event triggers are set in EventAction.
     nIncidentEvents++; 
     currentEvent = GenerateEvent(i_event);
     eventManager->ProcessOneEvent(currentEvent);
    if (eventTrigger) {
    AnalyzeEvent(currentEvent);
    if(i_event<n_select) G4UImanager::GetUIpointer()->ApplyCommand(msg);
    }
    StackPreviousEvent(currentEvent);
    currentEvent = 0;
    if(runAborted) break;
    } while (!eventTrigger);

      }

  if(verboseLevel>0)
  {
    timer->Stop();
    G4cout << "Run terminated." << G4endl;
    G4cout << "Run Summary" << G4endl;
    if(runAborted)
    { G4cout << "  Run Aborted after " << i_event << " events processed." << G4endl; }
    else
    { G4cout << "  Number of events processed : " << n_event << G4endl; }
    G4cout << "  "  << *timer << G4endl;
  }
}

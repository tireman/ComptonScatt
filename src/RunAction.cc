//******************************************************************************
// LICENSE DISCLAIMER:
// ===================

// Permission is granted to the public to copy and use this software
// without charge, provided that this Notice and any statement of
// authorship are reproduced on all copies.  Neither the Authors nor 
// Northern Michigan University makes any warranty, express or implied, 
// or assumes any liability or responsibility for the use of this software.
//
// RunAction.cc
//
//******************************************************************************

#include "RunAction.hh"

#include "G4Run.hh"
#include "RunManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4Timer.hh"
#include "G4ios.hh"
#include "Analysis.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

RunAction::RunAction()
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->SetVerboseLevel(1);

  man->CreateH1("TotalEnergy", "Total Energy in System", 100, 0., 2.0);

  // Create ntuple
  man->CreateNtuple("CSSData", "Energy Deposited");
  man->CreateNtupleIColumn("EventID");
  man->CreateNtupleDColumn("DetectorEnergy");
  man->CreateNtupleDColumn("TargetEnergy");
  man->CreateNtupleDColumn("ShieldEnergy");
  man->CreateNtupleDColumn("TotalEnergy");
  man->FinishNtuple();

  runIDcounter = 0;
  runTimer = new G4Timer(); // Setup the timer
}

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance(); 
}

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
 // Get/create analysis manager
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  
  // Open an output file
  man->OpenFile("CSSData");

  runTimer->Start(); // Start the timer so we can see how long it took

  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  G4cout << "      Event#   NaI(keV)   Pb(keV)    Cu(keV)   Etot(keV)" << G4endl;
  G4cout << "      -------  --------  ---------   --------  ---------" << G4endl;

  if (G4VVisManager::GetConcreteInstance()) {
    G4UImanager* UI = G4UImanager::GetUIpointer();
    UI->ApplyCommand("/vis/scene/notifyHandlers");
  } 
   
}

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  // Save histograms and ntuples
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();


  runTimer->Stop(); // stop the timer and print out time stuff

  G4cout << "### Run " << aRun->GetRunID() << " end." << G4endl;
  G4cout << "### Statistics: " << *(runTimer) << G4endl;
  G4cout << "                "
	 << RunManager::GetRunManager()->GetnIncident() 
	 << " incident photons" << G4endl;

  if (G4VVisManager::GetConcreteInstance()) {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }
}


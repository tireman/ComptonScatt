//******************************************************************************
// LICENSE DISCLAIMER:
// ===================

// Permission is granted to the public to copy and use this software
// without charge, provided that this Notice and any statement of
// authorship are reproduced on all copies.  Neither the Authors nor 
// Northern Michigan University makes any warranty, express or implied, 
// or assumes any liability or responsibility for the use of this software.
//
// EventAction.cc
//
//******************************************************************************
//

#include "EventAction.hh"
#include <fstream>
#include <iostream>
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4Run.hh"
#include "RunManager.hh"
#include "G4RunManager.hh"
#include "G4VVisManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "RunAction.hh"
#include "Analysis.hh"

EventAction::EventAction()
{
  // limit screen time. Can be altered in the macro files as well
  G4RunManager::GetRunManager()->SetPrintProgress(100); // limit screen time
}

EventAction::~EventAction()
{}

void EventAction::BeginOfEventAction(const G4Event*)
{
  D1dep = 0.;
  D2dep = 0.;
  D3dep = 0.;
}

void EventAction::EndOfEventAction(const G4Event* evt)
{
  // Here, decisions are made and output is generated to the screen
  // and to histograms and ntuple(s).
  
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  RunManager::GetRunManager()->SetEventTrigger();
  
  // let us be reasonable about the amount of "stuff" being posted to the 
  // screen.

  G4int event_id = evt->GetEventID();
  G4int output_num = G4RunManager::GetRunManager()->GetPrintProgress();
  if(event_id %output_num == 0){
    G4cout.width(12);
      G4cout << event_id;
      G4cout.width(10); G4cout.precision(4); G4cout.setf(std::ios::fixed);
      G4cout << D1dep/keV;
      G4cout.width(10); G4cout.precision(4); G4cout.setf(std::ios::fixed);
      G4cout << D2dep/keV;
      G4cout.width(10); G4cout.precision(4); G4cout.setf(std::ios::fixed);
      G4cout << D3dep/keV;
      G4cout.width(10); G4cout.precision(4); G4cout.setf(std::ios::fixed);
      G4cout << (D1dep+D2dep+D3dep)/keV << G4endl;   
    }

    // fill some histograms and ntuple columns
  
    man->FillH1(0,(D1dep+D2dep+D3dep)/MeV);

    man->FillNtupleIColumn(0, event_id);
    man->FillNtupleDColumn(1, D1dep);
    man->FillNtupleDColumn(2, D2dep);
    man->FillNtupleDColumn(3, D3dep);
    man->FillNtupleDColumn(3, (D1dep+D2dep+D3dep));
    man->AddNtupleRow(); 

    // Visualization of the event tracks; again, let's be reasonable
    
    if (G4VVisManager::GetConcreteInstance()) {
      if(event_id %output_num == 0){
      evt->Draw();
      }
    }    
}

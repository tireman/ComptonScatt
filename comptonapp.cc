//******************************************************************************
// LICENSE DISCLAIMER:
// ===================

// Permission is granted to the public to copy and use this software
// without charge, provided that this Notice and any statement of
// authorship are reproduced on all copies.  Neither the Authors nor 
// Northern Michigan University makes any warranty, express or implied, 
// or assumes any liability or responsibility for the use of this software.
//
// compton.cc
//
//******************************************************************************
//


// misc includes
//
#include <fstream>
#include <iostream>
#include <math.h>
#include "G4ios.hh"

// package includes
//
#include "RunManager.hh"
#include "DetectorConstruction.hh"
//#include "PhysicsList.hh"
#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

// geant4 includes
//
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "QGSP_BERT.hh"
#include "QGSP_BIC_HP.hh"
#include "QGSP_BIC.hh"
#include "FTFP_BERT.hh"

// visualization
//
//#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
//#endif

//------------------------------------------------------------------------------
int main(int argc,char** argv) {

  // Detect interactive mode (if no arguments) and define UI session
  //
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }
  
  // Cmd line args
  //--------------
  G4String macroFile = "";
  
  if (argc > 2) {
    G4cerr << "Usage: Compt [ <macro_file> ]"
	   << G4endl;
    return EXIT_FAILURE;
  }
  if (argc == 2) {
    macroFile = argv[1];
    std::ifstream mac(macroFile);
    if (!mac) {
      G4cout << "Cannot open macro file """ << macroFile << """" << G4endl;
      return 2;
    }
  }
  G4bool batch_mode = macroFile != ""; 

  // Run manager
  //------------
  RunManager * theRunManager = new RunManager;

  // UserInitialization classes
  //---------------------------
  DetectorConstruction* theDetector = new DetectorConstruction;

  theRunManager->SetUserInitialization(theDetector);
  theRunManager->SetUserInitialization(new FTFP_BERT);

  //#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  //-----------------------------------------
  // G4VisManager* theVisManager = new G4VisExecutive;
  //theVisManager->Initialize();
  //#endif
   
  // UserAction classes
  //-------------------
  theRunManager->SetUserAction(new RunAction);
  theRunManager->SetUserAction(new PrimaryGeneratorAction);

  EventAction* theEventAction = new EventAction;
  theRunManager->SetUserAction(theEventAction);
  theRunManager->SetUserAction(new SteppingAction(theDetector, theEventAction));

  // initialize Geant4 kernel
  theRunManager->Initialize();

  /*
  // User interactions
  //------------------
  G4UImanager* UI = G4UImanager::GetUIpointer();  
  
  // Define (G)UI terminal for interactive mode  
  //-------------------------------------------
  if(!batch_mode) { 
    // G4UIterminal is a (dumb) terminal.
    //-----------------------------------
    G4UIsession * theUIsession = new G4UIterminal(new G4UItcsh);
    UI->ApplyCommand("/control/execute init_vis.mac");    
    theUIsession->SessionStart();
    delete theUIsession;
  }
  // Batch mode
  //-----------
  else { 
    G4String command = "/control/execute ";
    UI->ApplyCommand(command+macroFile);
  }
  */

// Initialize visualization
  //
  G4VisManager* theVisManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  theVisManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( ! ui ) { 
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else { 
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }
  
  // Pause the program so I can look at the visual before it closes
  G4cout << "Press Return to continue" << G4endl;
  G4cin.get();
  
  //#ifdef G4VIS_USE
  delete theVisManager;
  //#endif
  delete theRunManager;
  return EXIT_SUCCESS;
  
}

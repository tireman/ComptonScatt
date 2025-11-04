//******************************************************************************
// LICENSE DISCLAIMER:
// ===================

// Permission is granted to the public to copy and use this software
// without charge, provided that this Notice and any statement of
// authorship are reproduced on all copies.  Neither the Authors nor 
// Northern Michigan University makes any warranty, express or implied, 
// or assumes any liability or responsibility for the use of this software.
//
// PrimaryGeneratorAction.cc
//
//******************************************************************************
//

#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
//#include "G4Particlegun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GeneralParticleSource.hh"

//----------------------------------------------------------------------------//
PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  particleGun = new G4GeneralParticleSource();
}

//----------------------------------------------------------------------------//
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
}

//----------------------------------------------------------------------------//
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
  particleGun->GeneratePrimaryVertex(anEvent);
}

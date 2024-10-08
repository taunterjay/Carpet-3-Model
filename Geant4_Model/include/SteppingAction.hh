#ifndef STEPPING_HH
#define STEPPING_HH

#include "G4UserSteppingAction.hh"
#include "G4Step.hh"
#include "G4VPhysicalVolume.hh"

#include "G4UserEventAction.hh"
#include "PhotosensorConstruction.hh"

class MySteppingAction : public G4UserSteppingAction{

	public:
		MySteppingAction(G4UserEventAction *EventAction);
		~MySteppingAction();
	
		virtual void UserSteppingAction(const G4Step*);
};

#endif

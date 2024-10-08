#include "SteppingAction.hh"

MySteppingAction::MySteppingAction(G4UserEventAction *EventAction){

}

MySteppingAction::~MySteppingAction(){
}

void MySteppingAction::UserSteppingAction(const G4Step *step){
	
	G4VPhysicalVolume *PreStepPointVolume = step->GetPreStepPoint()->GetPhysicalVolume();
	G4VPhysicalVolume *PostStepPointVolume = step->GetPostStepPoint()->GetPhysicalVolume(); 
	
	if (PostStepPointVolume != NULL){
		G4String PreVolName = PreStepPointVolume->GetName();
		G4String PostVolName = PostStepPointVolume->GetName();
	
		if ((PreVolName == "Floor") && (PostVolName == "Floor")){
		
			G4Track *track = step->GetTrack();
			track->SetTrackStatus(fStopAndKill);
	
		};
	};
	
}

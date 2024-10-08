#ifndef ACTION_HH
#define ACTION_HH

#include "G4VUserActionInitialization.hh"
//#include "G4Threading.hh"
//#include "G4MTRunManager.hh"

#include "Generator.hh"
#include "RunAction.hh"
//#include "EventAction.hh"
#include "G4UserEventAction.hh"
#include "G4UserRunAction.hh"
#include "SteppingAction.hh"

class Action : public G4VUserActionInitialization{

	public:
		Action();
		~Action();
    
		virtual void Build() const;
		virtual void BuildForMaster() const;
    
};

#endif

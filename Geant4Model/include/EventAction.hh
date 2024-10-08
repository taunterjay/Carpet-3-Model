#ifndef EVENT_HH
#define EVENT_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"

#include "RunAction.hh"

class MyEventAction : public G4UserEventAction{

	public:
		MyEventAction(MyRunAction* runAction);
		~MyEventAction();
    
		virtual void BeginOfEventAction(const G4Event*);
		virtual void EndOfEventAction(const G4Event*);
    
};

#endif

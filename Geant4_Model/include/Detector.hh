#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4SystemOfUnits.hh"
#include "G4VSensitiveDetector.hh"
#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
//#include "G4WorkerRunManager.hh"
#include "G4Threading.hh"
#include "G4OpticalPhoton.hh"

class SensitiveDetector : public G4VSensitiveDetector{

	public:
		SensitiveDetector(G4String);
		~SensitiveDetector();
		
		void Initialize(G4HCofThisEvent*);
		virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);
		void EndOfEvent(G4HCofThisEvent*);
    
	private:
		G4double CoverEnergy, RemoteOneEnergy, RemoteTwoEnergy, RemoteThreeEnergy, RemoteFourEnergy, CoverTime, RemoteOneTime, RemoteTwoTime, RemoteThreeTime, RemoteFourTime;
		G4bool CoverTrigger, RemoteOneTrigger, RemoteTwoTrigger, RemoteThreeTrigger, RemoteFourTrigger;
};

#endif

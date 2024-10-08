#ifndef MUONDETECTOR_HH
#define MUONDETECTOR_HH

#include "G4SystemOfUnits.hh"
#include "G4VSensitiveDetector.hh"
#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
//#include "G4WorkerRunManager.hh"
#include "G4Threading.hh"
#include "G4OpticalPhoton.hh"

class MuonSensitiveDetector : public G4VSensitiveDetector{

	public:
		MuonSensitiveDetector(G4String);
		~MuonSensitiveDetector();
		
		void Initialize(G4HCofThisEvent*);
		virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);
		void EndOfEvent(G4HCofThisEvent*);
    
	private:
		G4double SensorEnergy[410], SensorTime[410];
		G4bool SensorTrigger[410];

};

#endif

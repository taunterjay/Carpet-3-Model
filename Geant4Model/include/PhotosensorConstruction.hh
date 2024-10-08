#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

#include "G4SystemOfUnits.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4SubtractionSolid.hh"
#include "G4GDMLParser.hh"
#include "G4SDManager.hh"

#include "Detector.hh"
#include "MuonDetector.hh"

class CarpetConstruction : public G4VUserDetectorConstruction{

	public:
		CarpetConstruction();
		~CarpetConstruction();
    
		virtual G4VPhysicalVolume *Construct();
    
	private:
		G4LogicalVolume *PhotoLogical;
		G4LogicalVolume *MuonPhotoLogical;
    
		virtual void ConstructSDandField();
};

#endif

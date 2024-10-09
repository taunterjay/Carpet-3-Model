#include "RunAction.hh"

MyRunAction::MyRunAction(){
  
}

MyRunAction::~MyRunAction(){
  
}

void MyRunAction::BeginOfRunAction(const G4Run* run){

	std::time(&start);

	G4AnalysisManager *man = G4AnalysisManager::Instance();
	man->SetNtupleMerging(true); 
  
	G4int runNumber = run->GetRunID();
  
	std::stringstream strRunID;
	strRunID << runNumber;
 	
	G4int N_events = run->GetNumberOfEventToBeProcessed();
  	
	for (G4int i = 0; i < N_events; i++){
	
		std::fstream file;
		std::ostringstream oFileName;
  
		//oFileName << "/media/user/Elements\ SE/Read_showers/" << i << ".txt";
		oFileName << "../build/Showers/" << i << ".txt";
      
		file.open(oFileName.str ().c_str ());
  
		G4double E, theta, phi, x0, y0, Z_first;
  
		file >> E >> theta >> phi >> x0 >> y0 >> Z_first;
		
		//G4cout << E << G4endl;
  	
		file.close();
        
		std::ostringstream header, header_muon, header_time;
		header << E/1000 << "_" << theta << "_" << phi << "_" << x0 << "_" << y0 << "_" << Z_first;
		header_muon << "Muon_" << E/1000 << "_" << theta << "_" << phi << "_" << x0 << "_" << y0 << "_" << Z_first;
		header_time << "Time_" << E/1000 << "_" << theta << "_" << phi << "_" << x0 << "_" << y0 << "_" << Z_first;
	
		man->CreateH2(header.str ().c_str (), "", 70, -24.5*m, 24.5*m, 70, -24.5*m, 24.5*m, "m", "m");
		//man->CreateH2(header.str ().c_str (), "", 10, -7*m, 7*m, 10, -7*m, 7*m, "m", "m");
		man->CreateH2(header_muon.str ().c_str (), "", 11, -5.5*m, 5.5*m, 41, -20.5*m, 20.5*m, "m", "m");
		man->CreateNtuple(header_time.str ().c_str (), "RegistrationTimes");
		man->CreateNtupleFColumn("t_ns");
		man->FinishNtuple(i);
  	
	};
  	
	man->OpenFile("Output/output" + strRunID.str() + ".root");
}

void MyRunAction::EndOfRunAction(const G4Run*){
	
	G4AnalysisManager *man = G4AnalysisManager::Instance();
  
	man->Write();
	man->CloseFile();
  
	man->Clear();
	
	std::time(&end);
	if (IsMaster()) G4cout << G4endl << G4endl << G4endl << "The execution time equals " << (end - start)/60. << G4endl << G4endl << G4endl;
	
  
}

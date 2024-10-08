#include "Detector.hh"

//Detector отвечает за регистрацию фотонов, образующихся в ячейках жидкого сцинтиллятора и попадающих на окно ФЭУ

//Конструктор
SensitiveDetector::SensitiveDetector(G4String name) : G4VSensitiveDetector(name){
  
}

//Деструктор
SensitiveDetector::~SensitiveDetector(){
  
}

//Initialize выполняется в начале каждого события (1 событие - 1 ливень)

void SensitiveDetector::Initialize(G4HCofThisEvent*){
	
	//В начале каждого события обнуляем счётчик релятивистских частиц в Ковре и выносных пунктах. Триггер срабатывания порога переводим в ложное состояние
	
	CoverEnergy = 0., RemoteOneEnergy = 0., RemoteTwoEnergy = 0., RemoteThreeEnergy = 0., RemoteFourEnergy = 0., CoverTime = 0., RemoteOneTime = 0., RemoteTwoTime = 0., RemoteThreeTime = 0., RemoteFourTime = 0.;
	CoverTrigger = false, RemoteOneTrigger = false, RemoteTwoTrigger = false, RemoteThreeTrigger = false, RemoteFourTrigger = false;
}

//ProcessHits выполняется при каждом попадании частицы в детектор

G4bool SensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *R0hist){

	G4Track *track = aStep->GetTrack(); //достаём трэк частицы

	//Проверяем, что мы зарегистрировали именно сцинтилляционный фотон
	
	if ((track->GetDefinition() == G4OpticalPhoton::Definition()) && (track->GetCreatorProcess()->GetProcessName() == "Scintillation")){
  	
		const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable(); //Смотрим, в какой детектор попала частица
  	
		G4VPhysicalVolume *physVol = touchable->GetVolume(); //Получаем информацию об объёме детектора
		G4ThreeVector posDetector = physVol->GetTranslation(); //Получаем информацию о положении детектора
  	
		G4int EvtID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID(); //Достаём номер события
  	
		G4AnalysisManager *man = G4AnalysisManager::Instance();  //Вызываем AnalysisManager, чтобы записать данные в гистограмму
  	
		man->FillH2(EvtID*2, posDetector[0], posDetector[1], 1./666.); //Записываем данные в гистограмму с номером EvtID*2 в бин соответствующий положению детектора. 1/666 - нормировка
	
		track->SetTrackStatus(fStopAndKill); //Убиваем зарегистрированную частицу
	
		G4String DetName = physVol->GetName(); //Достаём имя детектора
		
		
		//Здесь мы проверяем, в каком детекторе зарегистрировался фотон: в Ковре или в выносном пункте. Также проверяется срабатывание триггера (15 р.ч. - для Ковра; 0.5 р.ч. - для выносных пунктов). Время срабатывания триггера фиксируется. Оно отсчитывается от начала запуска ливня
		if (DetName == "PhotoSensor"){
			if (not CoverTrigger){
				CoverEnergy += 1./666.;
				if (CoverEnergy > 15){
					CoverTrigger = true;
					CoverTime = aStep->GetPreStepPoint()->GetGlobalTime();
				}
			}
		
		}
		else if (DetName == "Remote1Sensor"){
			if (not RemoteOneTrigger){
				RemoteOneEnergy += 1./666.;
				if (RemoteOneEnergy > 0.5){
					RemoteOneTrigger = true;
					RemoteOneTime = aStep->GetPreStepPoint()->GetGlobalTime();
				}
			}
		
		}
		else if (DetName == "Remote2Sensor"){
			if (not RemoteTwoTrigger){
				RemoteTwoEnergy += 1./666.;
				if (RemoteTwoEnergy > 0.5){
					RemoteTwoTrigger = true;
					RemoteTwoTime = aStep->GetPreStepPoint()->GetGlobalTime();
				}
			}
		
		}
		else if (DetName == "Remote3Sensor"){
			if (not RemoteThreeTrigger){
				RemoteThreeEnergy += 1./666.;
				if (RemoteThreeEnergy > 0.5){
					RemoteThreeTrigger = true;
					RemoteThreeTime = aStep->GetPreStepPoint()->GetGlobalTime();
				}
			}
		
		}
		else if (DetName == "Remote4Sensor"){
			if (not RemoteFourTrigger){
				RemoteFourEnergy += 1./666.;
				if (RemoteFourEnergy > 0.5){
					RemoteFourTrigger = true;
					RemoteFourTime = aStep->GetPreStepPoint()->GetGlobalTime();
				}
			}
		
		}
  	
	}
  
	return true;
  
}

//EndOfEvent выполняется в конце каждого события

void SensitiveDetector::EndOfEvent(G4HCofThisEvent*){

	G4AnalysisManager *man = G4AnalysisManager::Instance(); //Вызываем AnalysisManager, чтобы записать данные в гистограмму
	G4int EvtID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID(); //Достаём номер события
	
	// Заполняем дерево времён срабатывания триггера
	
	man->FillNtupleFColumn(EvtID, 0, CoverTime);
	man->AddNtupleRow(EvtID);
	man->FillNtupleFColumn(EvtID, 0, RemoteOneTime);
	man->AddNtupleRow(EvtID);
	man->FillNtupleFColumn(EvtID, 0, RemoteTwoTime);
	man->AddNtupleRow(EvtID);
	man->FillNtupleFColumn(EvtID, 0, RemoteThreeTime);
	man->AddNtupleRow(EvtID);
	man->FillNtupleFColumn(EvtID, 0, RemoteFourTime);
	man->AddNtupleRow(EvtID);

}

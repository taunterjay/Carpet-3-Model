#include "MuonDetector.hh"

//MuonDetector отвечает за регистрацию фотонов, образующихся в ячейках пластического сцинтиллятора и попадающих на окно ФЭУ


//Конструктор
MuonSensitiveDetector::MuonSensitiveDetector(G4String name) : G4VSensitiveDetector(name){
  
}

//Деструктор
MuonSensitiveDetector::~MuonSensitiveDetector(){
  
}

//Initialize выполняется в начале каждого события (1 событие - 1 ливень)

void MuonSensitiveDetector::Initialize(G4HCofThisEvent*){

	//В начале каждого события обнуляем счётчик релятивистских частиц в каждой ячейке мюонного детектора. Триггер срабатывания порога переводим в ложное состояние

	for(G4int i = 0; i < 410; i++) SensorEnergy[i] = 0, SensorTime[i] = 0, SensorTrigger[i] = false;
	
}

//ProcessHits выполняется при каждом попадании частицы в детектор

G4bool MuonSensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *R0hist){

	G4Track *track = aStep->GetTrack(); //достаём трэк частицы
	
	//Проверяем, что мы зарегистрировали именно сцинтилляционный фотон
  
	if ((track->GetDefinition() == G4OpticalPhoton::Definition()) && (track->GetCreatorProcess()->GetProcessName() == "Scintillation")){
  	
		const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable(); //Смотрим, в какой детектор попала частица
  	
		G4VPhysicalVolume *physVol = touchable->GetVolume(); //Получаем информацию об объёме детектора
		G4ThreeVector posDetector = physVol->GetTranslation(); //Получаем информацию о положении детектора
  	
		G4int EvtID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID(); //Достаём номер события
  	
		G4AnalysisManager *man = G4AnalysisManager::Instance(); //Вызываем AnalysisManager, чтобы записать данные в гистограмму
  	
		man->FillH2(2*EvtID + 1, posDetector[0] + 45000.*mm, posDetector[1] + 3950.*mm, 1./524.4); //Записываем данные в гистограмму с номером EvtID*2 + 1 в бин соответствующий положению детектора. 1/524.4 - нормировка
	
		track->SetTrackStatus(fStopAndKill); //Убиваем зарегистрированную частицу
	
		G4int copyN = touchable->GetCopyNumber(); //Достаём номер детектора
		
		//Здесь мы проверяем срабатывание триггера 0.5 р.ч. для соответствующего номера детектора. Время срабатывания триггера фиксируется. Оно отсчитывается от начала запуска ливня
	
		if (not SensorTrigger[copyN]){
			SensorEnergy[copyN] += 1./524.4;
			if (SensorEnergy[copyN] > 0.5){
	   
				SensorTrigger[copyN] = true;
				SensorTime[copyN] = aStep->GetPreStepPoint()->GetGlobalTime();
			}
		}
  
  }
  
  return true;

}

//EndOfEvent выполняется в конце каждого события

void MuonSensitiveDetector::EndOfEvent(G4HCofThisEvent*){

	G4AnalysisManager *man = G4AnalysisManager::Instance(); //Вызываем AnalysisManager
	G4int EvtID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID(); //Достаём номер события
		
	// Заполняем дерево времён срабатывания триггера
	
	for(G4int i = 0; i < 410; i++){
		man->FillNtupleFColumn(EvtID, 0, SensorTime[i]);
		man->AddNtupleRow(EvtID);
	}
}

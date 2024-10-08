//#include <iostream>

#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"

#include "PhotosensorConstruction.hh"
#include "Physics.hh"
#include "Action.hh"

#include "Randomize.hh"

int main(int argc, char** argv){
  
	G4long seed = time(NULL); //Генерируем рандомный сид
	G4Random::setTheSeed(seed); //Генерируем рандомный сид
  
	#ifdef G4MULTITHREADED
		G4MTRunManager *runManager = new G4MTRunManager(); //запускаем программу в режиме многопоточности
	#else
		G4RunManager *runManager = new G4RunManager();
	#endif
  
	G4VModularPhysicsList* physics= new QGSP_BERT(); //Используем стандартный physicslist в связке с оптическими процессами (сцинтилляцией)
	physics->RegisterPhysics(new G4OpticalPhysics()); //Используем стандартный physicslist в связке с оптическими процессами (сцинтилляцией)
  
	runManager->SetUserInitialization(new CarpetConstruction); //Объявляем геометрию детектора
	runManager->SetUserInitialization(physics); //Обявляем PhysicsList
	runManager->SetUserInitialization(new Action); //Объявляем пользовательские действия
  
	G4UIExecutive *ui = 0; //Создаём пользовательский интерфейс, если не подаём на вход mac файлы
	if (argc == 1) ui = new G4UIExecutive(argc, argv); //Создаём пользовательский интерфейс, если не подаём на вход mac файлы
  
  
	G4VisManager *visManager = new G4VisExecutive(); //Создаём VisManager и UImanager
	visManager->Initialize(); //Создаём VisManager и UImanager
  
	G4UImanager *UImanager = G4UImanager::GetUIpointer(); //Создаём VisManager и UImanager
  
  	//Запускаем визуализацию из vis.mac или же исполняем run.mac файл
  
	if (ui){
		UImanager->ApplyCommand("/control/execute vis.mac");
		ui->SessionStart();
	}
	else{
 	G4String command = "/control/execute ";
 	G4String fileName = argv[1];
 	UImanager->ApplyCommand(command + fileName);
	}
  
	delete visManager; //отчистка памяти
	delete runManager; //отчистка памяти
	delete ui; //отчистка памяти
  
	return 0;
}

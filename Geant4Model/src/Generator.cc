#include "Generator.hh"

//В Generator задаются параметры первичных частиц: их положение, импульс и энергия 


//PrimaryGenerator - здесь мы задаём некторые значения, которые будут использоваться в случае, если далее программе ничего дополнительного не указать. В нашем случае по умолчанию генерируется вертикальный ГэВный мюон, расположенный над мюонным детектором
PrimaryGenerator::PrimaryGenerator(){

	fParticleGun = new G4ParticleGun(1);
  
	G4ParticleTable *table = G4ParticleTable::GetParticleTable();
	G4String particleName = "mu-";
	G4ParticleDefinition *particle = table->FindParticle(particleName);

	G4ThreeVector position(-45000.*mm, 16177.*mm,10.*m);
	G4ThreeVector m(0.,0.,-1.);
	G4double energy = 1*GeV;

	fParticleGun->SetParticlePosition(position);
	fParticleGun->SetParticleMomentumDirection(m);
	fParticleGun->SetParticleEnergy(energy);
	fParticleGun->SetParticleDefinition(particle);
}

//Деструктор
PrimaryGenerator::~PrimaryGenerator(){
	delete fParticleGun;
}

//В GeneratePrimaries происходит считывание корсиканского ливня и создание события для дальнейшей обработки  
void PrimaryGenerator::GeneratePrimaries(G4Event *Event){
  
	G4int EventID = Event->GetEventID(); //Достаём номер события
  
	G4ParticleTable *table = G4ParticleTable::GetParticleTable(); //Достаём список частиц и их свойств
  
	std::fstream file; //переменная для работы с корсиканским файлом
	std::ostringstream oFileName; //переменная имя файла
  
	//oFileName << "/media/user/Elements\ SE/Read_showers/" << EventID << ".txt"; //Считываем название файла в переменную
	oFileName << "/home/user/Showers/Result/" << EventID << ".txt";
	
      
	file.open(oFileName.str ().c_str ()); //Открываем соответствующий файл
  
	G4double E, theta, phi, x0, y0, Z_first; //Переменные для энергии, углов, положения оси и высоты первого взаимодействия ливня
  
	file >> E >> theta >> phi >> x0 >> y0 >> Z_first; //Считываем параметры ливня из файла

	G4String id; //переменная имя частицы
	G4double x, y, e, px_dir, py_dir, pz_dir, t; //переменные координаты, энергия, компоненты импульса частицы
  
	while(1){
		file >> id >> x >> y >> e >> px_dir >> py_dir >> pz_dir >> t; //считываем ливень почастично
		if(file.eof()) break;
    	
		G4ParticleDefinition *particle = table->FindParticle(id); //Находим частицу по её имени
    	
		G4ThreeVector position(x*cm + x0*m, y*cm + y0*m, 10*m); //вектор положение частицы
	
		G4ThreeVector m(px_dir, py_dir, pz_dir); //вектор направление импульса частицы

		fParticleGun->SetParticlePosition(position); //Задаём положение частицы
		fParticleGun->SetParticleMomentumDirection(m); //Задаём направление импульса частицы
		fParticleGun->SetParticleEnergy(e*GeV); //Задаём энергию частицы
		fParticleGun->SetParticleTime(t*ns);
		fParticleGun->SetParticleDefinition(particle); //Задаём частицу
    	
		fParticleGun->GeneratePrimaryVertex(Event); //Добавляем её в событие
	}
  
	file.close(); //Закрываем файл

	G4cout << "Generated event " << EventID << G4endl; //Для удобства выводим номер сгенерированного ливня

}

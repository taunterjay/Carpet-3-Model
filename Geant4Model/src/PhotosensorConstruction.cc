#include "PhotosensorConstruction.hh"

//PhotosensorConstruction содержит описание геометрии установки, свойства материалов, полей и т.д.


//Конструктор
CarpetConstruction::CarpetConstruction(){

}

//Деструктор
CarpetConstruction::~CarpetConstruction(){

}

//Construct - описание геометрии и материалов
G4VPhysicalVolume *CarpetConstruction::Construct(){

	//Задаём параметры: Длина стороны ковра, длина/ширина выносных пунктов, длина массива энергий для описания сцинтилляции, длина/ширина мюонного детектора, корень количества пластиковых сцинтилляторов в одном кожухе, количество тоннелей мюонного детектора, переменная для работы с составными материалами, массив энергий для описания сцинтилляции, координаты выносных пунктов, переменные для работы с плотностью и химическими соединениями, переменная имя
  
	G4int CarpetCellNumber = 20, Remote_n = 3, Remote_m = 6, nE = 3, N_rows = 41, N_columns = 5, N_scint = 2, N_rooms = 3, ncomponents;
	G4double Energy[nE] = {7.00*eV, 7.07*eV, 7.14*eV};
	G4double x_move1, y_move1, x_move2, y_move2, x_move3, y_move3, x_move4, y_move4, density, fractionmass;
	G4String name;
  
	auto nist = G4NistManager::Instance(); //NistManager используется для поиска маетриалов в база Geant4

	G4GDMLParser parser; //GDMLParser позволяет читать GDML-файлы и использовать их в моделях
    
	parser.SetOverlapCheck(true); //проверка пересечений с другими объёмами
	parser.Read("world.gdml", false); //считываем файл с крышой
  
	G4Material *wmt = nist->FindOrBuildMaterial("G4_AIR"); //Создаём материал "Воздух" для заполнения мирового объёма
  
	G4double RindexWorld[nE] = {1., 1., 1.}; //показатель преломления воздуха
  
	G4MaterialPropertiesTable *MPTworld = new G4MaterialPropertiesTable(); //Создаём таблицу свойств воздуха
	MPTworld->AddProperty("RINDEX", Energy, RindexWorld, nE); //Задаём показатель преломления воздуха
  
	wmt->SetMaterialPropertiesTable(MPTworld); //Присваиваем воздуху его показатель преломления
	MPTworld->DumpTable(); //Отчищаем память
    
	G4double world_x = 65*m; //размеры мирового объёма
	G4double world_y = 50*m; //размеры мирового объёма
	G4double world_z = 20*m; //размеры мирового объёма
  
	G4Box *WorldBox = new G4Box(name="World", world_x, world_y, world_z); //Создаём форму мирового объёма
  
	G4LogicalVolume *WorldLogical = new G4LogicalVolume(WorldBox, wmt, name="World"); //Создаём логический мировой объём, состоящий из воздуха
	G4VPhysicalVolume *WorldPhysical = new G4PVPlacement(0, G4ThreeVector(), WorldLogical, name="World", nullptr, false, 0, true); //Создаём мировой объём физически 
  
	G4LogicalVolume *RoofLogical = parser.GetVolume("CarpetRoof_Concrete"); //считываем из GDML-файла логический объём крыши
  
	G4RotationMatrix* zRot = new G4RotationMatrix; //Задаём матрицу попворота, чтобы правильно поставить крышу
	zRot->rotateZ(M_PI/2.*rad); //Задаём матрицу попворота, чтобы правильно поставить крышу
  
	G4VPhysicalVolume *RoofPhysical = new G4PVPlacement(zRot, G4ThreeVector(-990*cm, 920*cm, 0), RoofLogical, name="Roof", WorldLogical, false, 0, true); //Создаём крышу физически
    
	G4Material *Al = nist->FindOrBuildMaterial("G4_Al"); //Создаём материал "Алюминий" для кожухов сцинтилляторов
  
	G4double RindexAl[nE] = {1.44, 1.44, 1.44}; //показатель преломления алюминия
  
	G4MaterialPropertiesTable *MPTal = new G4MaterialPropertiesTable(); //Создаём таблицу свойств алюминия
	MPTal->AddProperty("RINDEX", Energy, RindexAl, nE); //Задаём показатель преломления алюминия
  
	Al->SetMaterialPropertiesTable(MPTal); //Присваиваем алюминию показатель преломления
	MPTal->DumpTable(); //Отчищаем память
  
	G4double al_x = 0.35*m; //размер кожухов для жидкого сцинтиллятора по х
	G4double al_y = 0.35*m; //размеры кожухов для жидкого сцинтиллятора по у
	G4double al_z = 0.15*m; //размерыкожухов для жидкого сцинтиллятора по z
	G4double al_width = 3*mm; //толщина кожухов для жидкого сцинтиллятора
  
	G4double hole_Rmin = 0.*m; //внутренний радиус отверстия в кожухе под окно ФЭУ
	G4double hole_Rmax = 7.5*cm; //внешний радиус отверстия в кожухе под окно ФЭУ
	G4double hole_PhiMin = 0.; //стартовый угол отверстия в кожухе под окно ФЭУ
	G4double hole_PhiMax = CLHEP::twopi; //конечный угол отверстия в кожухе под окно ФЭУ
  
	G4RotationMatrix rotm = G4RotationMatrix(); //Задаём матрицу поворота для отверстия
  
	G4Box *AlBig = new G4Box(name="AlBig", al_x, al_y, al_z); //Создаём форму кожуха
	G4Box *AlLil = new G4Box(name="AlLil", al_x - al_width, al_y - al_width, al_z - al_width); //Создаём форму полости, в которой будет находится жидкий сцинтиллятор
	G4Tubs *PhotoHole = new G4Tubs(name="PhotoHole", hole_Rmin, hole_Rmax, al_width/2, hole_PhiMin, hole_PhiMax); //Создаём форму отверстия в кожухе под окно ФЭУ
  
	G4ThreeVector PhotoPosition = G4ThreeVector(0., 0., al_z - al_width/2); //Задаём координаты отверстия в кожухе под окно ФЭУ
	G4Transform3D PhotoTr = G4Transform3D(rotm, PhotoPosition); //Задаём комбинацию поворота и координат отверстия в кожухе под окно ФЭУ
  
	G4SubtractionSolid *CoverSolid =  new G4SubtractionSolid(name="Cover", AlBig, AlLil); //Вырезаем в форме кожуха отверстие под сцинтиллятор
	G4SubtractionSolid *FinalCoverSolid = new G4SubtractionSolid(name="FinalCover", CoverSolid, PhotoHole, PhotoTr); //Вырезаем в форме кожуха отверстие под окно ФЭУ
  
	G4LogicalVolume *CoverLogical = new G4LogicalVolume(FinalCoverSolid, Al, name="CoverLogical"); //Задаём логический объём кожуха из алюминия
  
	G4OpticalSurface *AlSurface = new G4OpticalSurface(name="AlSurface"); //OpticalSurface используется для симуляции отражающей эмали
  
	AlSurface->SetType(dielectric_dielectric); //Выбираем тип взимодейтсвия диэлектрик-диэлектрик
	AlSurface->SetModel(unified); //Задаём модель
	AlSurface->SetFinish(groundfrontpainted); //Указываем, что кожух смазан эмалью изнутри
  
	G4LogicalSkinSurface *Surface = new G4LogicalSkinSurface(name="Surface", CoverLogical, AlSurface); //Присваиваем свойства отражающей поверхности кожуху
  
	std::vector<G4double> pp = {0.5*eV, 4.144*eV}; //массив энергий
	std::vector<G4double> reflectivity = {.95, .95}; //коэффициент отражения
	std::vector<G4double> rindex = {1.5, 1.5}; //показатель преломления
  
	G4MaterialPropertiesTable* AlSurfaceProperty = new G4MaterialPropertiesTable(); //Создаём маблицу свойств поверхности
	AlSurfaceProperty->AddProperty("REFLECTIVITY", pp, reflectivity); //Задаём вероятность отражения
	AlSurfaceProperty->AddProperty("RINDEX", pp, rindex); //Задаём показатель преломления
  
	AlSurface->SetMaterialPropertiesTable(AlSurfaceProperty); //Присваиваем свойства повверхности
  
	//Помещаем 400 кожухов в мировой объём
	for(G4int i = 0; i < CarpetCellNumber; i++)
		for(G4int j = 0; j < CarpetCellNumber; j++){
			G4double x = ((-7.+0.35) + i*0.7)*m;
			G4double y = ((7.-0.35) - j*0.7)*m;
			G4VPhysicalVolume *CoverPhysical = new G4PVPlacement(0, G4ThreeVector(x, y, al_z), CoverLogical, name="CoverPhysical", WorldLogical, false, i + j*CarpetCellNumber, true);
		}
  
	//Carpet Cell
   
   
	//Задаём материалы жидкого сцинтиллятора
	G4Material *elH = nist->FindOrBuildMaterial("G4_H");
	G4Material *elC = nist->FindOrBuildMaterial("G4_C");
  
	density = 0.78*g/cm3; //плотность жидкого сцинтиллятора
	G4Material *cmt = new G4Material(name="Liquid Scintillator", density, ncomponents=2); //Задаём материал жидкого сцинтиллятора
	cmt->AddMaterial(elH, fractionmass=0.84375); //Задаём материал жидкого сцинтиллятора
	cmt->AddMaterial(elC, fractionmass=0.15625); //Задаём материал жидкого сцинтиллятора
  
	G4double RindexCell[nE] = {1.58, 1.58, 1.58}; //показатель преломления жидкого сцинтиллятора
	G4double AbsLength[nE] = {2000.*cm, 2000.*cm, 2000.*cm}; //absorption length
	G4double ScintCell[nE] = {0.1, 1., 0.1}; //энергетический спектр сцинтилляционных фотонов
  
	G4MaterialPropertiesTable *MPTcell = new G4MaterialPropertiesTable(); //Создаём таблицу свойств жидкого сцинтиллятора
  
	MPTcell->AddProperty("SCINTILLATIONCOMPONENT1", Energy, ScintCell, nE); //Задаём энергетический спектр сцинтилляционных фотонов
	MPTcell->AddConstProperty("SCINTILLATIONYIELD", 100./MeV); //Задаём количество фотонов на МэВ ионизации
	MPTcell->AddProperty("RINDEX", Energy, RindexCell, nE); //Задаём показатель преломления
	MPTcell->AddProperty("ABSLENGTH", Energy, AbsLength, nE); //Задаём absorption length
	MPTcell->AddConstProperty("RESOLUTIONSCALE", 0.); //Задаём нулевую флуктуацию числа фотонов, образующихся на каждом шаге
	MPTcell->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 20.*ns); //Задаём время высвечивания?
  
	cmt->SetMaterialPropertiesTable(MPTcell); //Передаём установленные свойства материалу жидкого сцинтиллятора
	MPTcell->DumpTable(); //Отчистка памяти
  
	G4double cell_x = al_x - al_width; //Размер жидкого сцинтиллятора х
	G4double cell_y = al_y - al_width; //Размер жидкого сцинтиллятора y
	G4double cell_z = al_z - al_width; //Размер жидкого сцинтиллятора z
  
	G4Box *CellSolid = new G4Box(name="Carpet Cell", cell_x, cell_y, cell_z); //Создаём форму ячейки жидкого сцинтиллятора
  
	G4LogicalVolume *CellLogical = new G4LogicalVolume(CellSolid, cmt, name="Carpet Cell"); //Задаём логический объём ячейки жидкого сцинтиллятора
	
	//Создаём 400 ячеек жидкого сцинтиллятора Ковра

	for(G4int i = 0; i < CarpetCellNumber; i++)
		for(G4int j = 0; j < CarpetCellNumber; j++){
			G4double x = ((-7.+0.35) + i*0.7)*m;
			G4double y = ((7.-0.35) - j*0.7)*m;
			G4VPhysicalVolume *CellPhysical = new G4PVPlacement(0, G4ThreeVector(x, y, al_z), CellLogical, name="Carpet Cell", WorldLogical, false, i + j*CarpetCellNumber, true);
		}
    
	//Photocatode
  
	G4double photo_Rmin = hole_Rmin; //стартовый радиус окна ФЭУ
	G4double photo_Rmax = hole_Rmax; //конечный радиус окна ФЭУ
	G4double photo_PhiMin = 0.; //стартовый угол окна ФЭУ
	G4double photo_PhiMax = CLHEP::twopi; //конечный угол окна ФЭУ
	G4double photo_z = 0.006*m; //полувысота окна ФЭУ
  
	G4Tubs *PhotoSolid = new G4Tubs(name="PhotoSolid", photo_Rmin, photo_Rmax, photo_z, photo_PhiMin, photo_PhiMax); //Создаём форму окна ФЭУ
  
	PhotoLogical = new G4LogicalVolume(PhotoSolid, wmt, name="PhotoLogical"); //Создаём логический объём окна
	
	//Выставляем 400 окон ФЭУ Ковра
  
	for(G4int i = 0; i < CarpetCellNumber; i++)
		for(G4int j = 0; j < CarpetCellNumber; j++){
		G4double x = ((-7.+0.35) + i*0.7)*m;
		G4double y = ((7.-0.35) - j*0.7)*m;
		G4VPhysicalVolume *PhotoPhysical = new G4PVPlacement(0, G4ThreeVector(x, y, al_z + cell_z + photo_z), PhotoLogical, name="PhotoSensor", WorldLogical, false, i + j*CarpetCellNumber, true);
		}
    
	//Remote points
    
	//1
    
	x_move1 = -21.29 - Remote_n*0.35; //координата х первого выносного пункта
	y_move1 = 21.68 + Remote_m*0.35; //координата у первого выносного пункта
	
	//Выставляем 18 кожухов, сцинтилляторов и окон ФЭУ
    
	for(G4int i = 0; i < Remote_n; i++)
		for(G4int j = 0; j < Remote_m; j++){
			G4double x = ((x_move1 + 0.35) + i*0.7)*m;
			G4double y = ((y_move1 - 0.35) - j*0.7)*m;
			G4VPhysicalVolume *RemoteOnePhysical = new G4PVPlacement(0, G4ThreeVector(x, y, al_z), CoverLogical, name="Remote1Physical", WorldLogical, false, i*Remote_m + j, true);
			G4VPhysicalVolume *CellPhysical = new G4PVPlacement(0, G4ThreeVector(x, y, al_z), CellLogical, name="Remote1Cell", WorldLogical, false, i*Remote_m + j, true);
			G4VPhysicalVolume *PhotoPhysical = new G4PVPlacement(0, G4ThreeVector(x, y, al_z + cell_z + photo_z), PhotoLogical, name="Remote1Sensor", WorldLogical, false, i*Remote_m + j, true);
		}
      
	//2
    
	x_move2 = -21.63 - Remote_n*0.35; //координата х второго выносного пункта
	y_move2 = -21.13 + Remote_m*0.35; //координата у второго выносного пункта
	
	//Выставляем 18 кожухов, сцинтилляторов и окон ФЭУ
    
	for(G4int i = 0; i < Remote_n; i++)
		for(G4int j = 0; j < Remote_m; j++){
			G4double x = ((x_move2 + 0.35) + i*0.7)*m;
			G4double y = ((y_move2 - 0.35) - j*0.7)*m;
			G4VPhysicalVolume *RemoteTwoPhysical = new G4PVPlacement(0, G4ThreeVector(x, y, al_z), CoverLogical, name="Remote2Physical", WorldLogical, false, i*Remote_m + j, true);
			G4VPhysicalVolume *CellPhysical = new G4PVPlacement(0, G4ThreeVector(x, y, al_z), CellLogical, name="Remote2Cell", WorldLogical, false, i*Remote_m + j, true);
			G4VPhysicalVolume *PhotoPhysical = new G4PVPlacement(0, G4ThreeVector(x, y, al_z + cell_z + photo_z), PhotoLogical, name="Remote2Sensor", WorldLogical, false, i*Remote_m + j, true);
		}

	//3
    
	x_move3 = 17.5 - Remote_n*0.35; //координата х третьего выносного пункта
	y_move3 = 9.34 + Remote_m*0.35; //координата у третьего выносного пункта
	
	//Выставляем 18 кожухов, сцинтилляторов и окон ФЭУ
    
	for(G4int i = 0; i < Remote_n; i++)
		for(G4int j = 0; j < Remote_m; j++){
			G4double x = ((x_move3 + 0.35) + i*0.7)*m;
			G4double y = ((y_move3 - 0.35) - j*0.7)*m;
			G4VPhysicalVolume *RemoteThreePhysical = new G4PVPlacement(0, G4ThreeVector(x, y, al_z), CoverLogical, name="Remote3Physical", WorldLogical, false, i*Remote_m + j, true);
			G4VPhysicalVolume *CellPhysical = new G4PVPlacement(0, G4ThreeVector(x, y, al_z), CellLogical, name="Remote3Cell", WorldLogical, false, i*Remote_m + j, true);
			G4VPhysicalVolume *PhotoPhysical = new G4PVPlacement(0, G4ThreeVector(x, y, al_z + cell_z + photo_z), PhotoLogical, name="Remote3Sensor", WorldLogical, false, i*Remote_m + j, true);
		}
      
	//4
    
	x_move4 = 21.16 - Remote_n*0.35; //координата х четвёртого выносного пункта
	y_move4 = -21.05 + Remote_m*0.35; ////координата у четвёртого выносного пункта
	
	//Выставляем 18 кожухов, сцинтилляторов и окон ФЭУ
    
	for(G4int i = 0; i < Remote_n; i++)
		for(G4int j = 0; j < Remote_m; j++){
			G4double x = ((x_move4 + 0.35) + i*0.7)*m;
			G4double y = ((y_move4 - 0.35) - j*0.7)*m;
			G4VPhysicalVolume *RemoteFourPhysical = new G4PVPlacement(0, G4ThreeVector(x, y, al_z), CoverLogical, name="Remote4Physical", WorldLogical, false, i*Remote_m + j, true);
			G4VPhysicalVolume *CellPhysical = new G4PVPlacement(0, G4ThreeVector(x, y, al_z), CellLogical, name="Remote4Cell", WorldLogical, false, i*Remote_m + j, true);
			G4VPhysicalVolume *PhotoPhysical = new G4PVPlacement(0, G4ThreeVector(x, y, al_z + cell_z + photo_z), PhotoLogical, name="Remote4Sensor", WorldLogical, false, i*Remote_m + j, true);
		}
    
	//Muon Detector
    
	G4Material *Concrete = nist->FindOrBuildMaterial("G4_CONCRETE"); //материал стен, пола и крыши мюонного детектора
    
	G4double MuonXmove = 48000.*mm; //сдвиг мюонного детектора по х
	G4double MuonYmove = 4000*mm; //сдвиг мюонного детектора по у
    
	G4double Room_x = 6000./2.*mm; //размер комнаты мюонного детектора по х
	G4double Room_y = 43000./2.*mm; //размер комнаты мюонного детектора по у
	G4double Room_z = 2400./2.*mm; //размер комнаты мюонного детектора по z
	G4double Room_width = (6000.-5370.)/2.*mm; //толщина стен команты мюонного детектора
	G4double Room_width_z = 300.*mm; //толщина пола и потолка мюонного детектора
  
	G4double Room_depth = 1800.*mm; //глубина залегания команты мюонного детектора 
  
	G4Box *RoomBig = new G4Box(name="RoomBig", Room_x, Room_y, Room_z + Room_width_z); //Форма комнаты мюонного детектора
	G4Box *RoomLil = new G4Box(name="RoomLil", Room_x - Room_width, Room_y - Room_width, Room_z); //Форма полости в комнате мюонного детектора
  
	G4SubtractionSolid *RoomSolid =  new G4SubtractionSolid(name="MuonRoom", RoomBig, RoomLil); //Создаём комнату с полостью внутри
	G4LogicalVolume *RoomLogical = new G4LogicalVolume(RoomSolid, Concrete, name="MuonRoomLogical"); //Логический объём комнаты мюонного детектора
	
	//Расставляем все три комнаты
	
	for(G4int i = 0; i < N_rooms; i++){
		G4double x = (i - 1)*2*Room_x;
		G4VPhysicalVolume *RoomPhysical = new G4PVPlacement(0, G4ThreeVector(x - MuonXmove, -MuonYmove, -Room_depth + Room_z), RoomLogical, name="MuonRoomPhysical", WorldLogical, false, i, true);
	}
    
	G4Material *Fe = nist->FindOrBuildMaterial("G4_Fe"); //материал обивки стен мюонного детектора
    
	G4double Fe_RoomWidth = 5.*mm; //толщина слоя железа
  
	G4Box *RoomFe = new G4Box(name="RoomFe", Room_x - Room_width - Fe_RoomWidth, Room_y - Room_width - Fe_RoomWidth, Room_z - Fe_RoomWidth); //Форма комнаты с учётом слоя железа
  
	G4SubtractionSolid *RoomFeSolid =  new G4SubtractionSolid(name="MuonRoomFe", RoomLil, RoomFe); //Форма железной обивки
	G4LogicalVolume *RoomFeLogical = new G4LogicalVolume(RoomFeSolid, Fe, name="MuonRoomFeLogical"); //Логический объём железной обивки
	
	//Расставляем обивку в каждую комнату
	
	for(G4int i = 0; i < N_rooms; i++){
		G4double x = (i - 1)*2*Room_x;
		G4VPhysicalVolume *RoomFePhysical = new G4PVPlacement(0, G4ThreeVector(x - MuonXmove, -MuonYmove, -Room_depth + Room_z), RoomFeLogical, name="MuonRoomFePhysical", WorldLogical, false, i, true);
	}
  
  	//Самая левая и самая правая стена мюонного детектора толще остальных. Поэтому слева и справа добавлены дополнительные "края"
  
	G4double Edge_x = (6000.-5370.)/4.*mm; //размер "края" по х
	G4double Edge_y = 43000./2.*mm; //размер "края" по у
	G4double Edge_z = 2400./2.*mm; //размер "края" по z
  
	G4Box *RoomEdge = new G4Box(name="RoomEdge", Edge_x, Edge_y, Edge_z + Room_width_z); //Форма "краёв"
	G4LogicalVolume *RoomEdgeLogical = new G4LogicalVolume(RoomEdge, Concrete, name="MuonRoomEdgeLogical"); //Логический объём "краёв"
	
	//Расставляем "края"
	for(G4int i = 0; i < 2; i++){
		G4double x = (2*i - 1)*(3*Room_x + Edge_x);
		G4VPhysicalVolume *RoomEdgePhysical = new G4PVPlacement(0, G4ThreeVector(x - MuonXmove, -MuonYmove, -Room_depth + Room_z), RoomEdgeLogical, name="MuonRoomEdgePhysical", WorldLogical, false, i, true);
	}
    
	G4double ConcreteRoof_x = 3*Room_x + 2*Edge_x; //Размер бетонной стяжки по х
	G4double ConcreteRoof_y = 43000./2.*mm; //Размер бетонной стяжки по у
	G4double ConcreteRoof_z = (800*mm - Room_width_z)/2.; //Размер бетонной стяжки по z
    
	G4Box *ConcreteMuonRoof = new G4Box(name="ConcreteMuonRoof", ConcreteRoof_x, ConcreteRoof_y, ConcreteRoof_z); //Форма бетонной стяжки
	G4LogicalVolume *ConcreteMuonRoofLogical = new G4LogicalVolume(ConcreteMuonRoof, Concrete, name="ConcreteMuonRoofLogical"); //Логический объём бетонной стяжки 
	G4VPhysicalVolume *ConcreteMuonRoofPhysical = new G4PVPlacement(0, G4ThreeVector(-MuonXmove, -MuonYmove, -Room_depth + 2*Room_z+Room_width_z + ConcreteRoof_z), ConcreteMuonRoofLogical, name="ConcreteMuonRoofPhysical", WorldLogical, false, 0, true); //Физический объём бетонной стяжки
    
    
	G4Material* SiO2 = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE"); //Задаём материалы, из которых состоит гранит
	G4Material* Al2O3 = nist->FindOrBuildMaterial("G4_ALUMINUM_OXIDE"); //Задаём материалы, из которых состоит гранит
	G4Material* K2O = nist->FindOrBuildMaterial("G4_POTASSIUM_OXIDE"); //Задаём материалы, из которых состоит гранит
	G4Material* Na2O = nist->FindOrBuildMaterial("G4_SODIUM_MONOXIDE"); //Задаём материалы, из которых состоит гранит
	G4Material* CaO = nist->FindOrBuildMaterial("G4_CALCIUM_OXIDE"); //Задаём материалы, из которых состоит гранит
	G4Material* FeO = nist->FindOrBuildMaterial("G4_FERROUS_OXIDE"); //Задаём материалы, из которых состоит гранит
	G4Material* Fe2O3 = nist->FindOrBuildMaterial("G4_FERRIC_OXIDE"); //Задаём материалы, из которых состоит гранит
	G4Material* MgO = nist->FindOrBuildMaterial("G4_MAGNESIUM_OXIDE"); //Задаём материалы, из которых состоит гранит
	G4Material* TiO2 = nist->FindOrBuildMaterial("G4_TITANIUM_DIOXIDE"); //Задаём материалы, из которых состоит гранит
  
	G4Material* P = nist->FindOrBuildMaterial("G4_P"); //Задаём материалы, из которых состоит гранит
	G4Material* O = nist->FindOrBuildMaterial("G4_O"); //Задаём материалы, из которых состоит гранит
  
	density = 2.39*g/cm3; //Задаём материалы, из которых состоит гранит
	G4Material* P2O5 = new G4Material(name="P2O5", density, ncomponents=2); //Задаём материалы, из которых состоит гранит
	P2O5->AddMaterial(P, fractionmass=0.42857); //Задаём материалы, из которых состоит гранит
	P2O5->AddMaterial(O , fractionmass=0.57143); //Задаём материалы, из которых состоит гранит
  
	G4Material* Mn = nist->FindOrBuildMaterial("G4_Mn"); //Задаём материалы, из которых состоит гранит
	density = 5.18*g/cm3;
	G4Material* MnO = new G4Material(name="MnO", density, ncomponents=2); //Задаём материалы, из которых состоит гранит
	MnO->AddMaterial(Mn, fractionmass=0.757575); //Задаём материалы, из которых состоит гранит
	MnO->AddMaterial(O , fractionmass=0.242425); //Задаём материалы, из которых состоит гранит
  
	density = 0.8*2.6*g/cm3; //плотность гранита с поправкой на полости
	G4Material* Granite = new G4Material(name="Granite", density, ncomponents=11); //Создаём гранит
	Granite->AddMaterial(SiO2, fractionmass=0.7189);  //Добавляем материалы в соответствующей пропорции
	Granite->AddMaterial(Al2O3, fractionmass=0.144); //Добавляем материалы в соответствующей пропорции
	Granite->AddMaterial(K2O, fractionmass=0.0412); //Добавляем материалы в соответствующей пропорции
	Granite->AddMaterial(Na2O, fractionmass=0.0369); //Добавляем материалы в соответствующей пропорции
	Granite->AddMaterial(CaO, fractionmass=0.0182); //Добавляем материалы в соответствующей пропорции
	Granite->AddMaterial(FeO, fractionmass=0.0168); //Добавляем материалы в соответствующей пропорции
	Granite->AddMaterial(Fe2O3, fractionmass=0.0122); //Добавляем материалы в соответствующей пропорции
	Granite->AddMaterial(MgO, fractionmass=0.0071); //Добавляем материалы в соответствующей пропорции
	Granite->AddMaterial(TiO2, fractionmass=0.0030); //Добавляем материалы в соответствующей пропорции
	Granite->AddMaterial(P2O5, fractionmass=0.0012); //Добавляем материалы в соответствующей пропорции
	Granite->AddMaterial(MnO, fractionmass=0.0005); //Добавляем материалы в соответствующей пропорции
  
	G4double Earth_x1 = 3*Room_x + 2*Edge_x + 10000./2.*mm; //Размер нижней плоскости гранитной насыпи по х
	G4double Earth_x2 = 3*Room_x + 2*Edge_x; //Размер верхней плоскости гранитной насыпи по х
	G4double Earth_y1 = (43000. + 10000.)/2.*mm; //Размер нижней плоскости гранитной насыпи по у
	G4double Earth_y2 = 43000./2.*mm; //Размер верхней плоскости гранитной насыпи по у
	G4double Earth_z = 3400./2.*mm; //Размер гранитной насыпи по z
  
	G4double Hole_z = Room_z - Room_depth/2 + Room_width_z/2 + ConcreteRoof_z; //Высота полости в насыпи. Эту часть занимает мюонный детектор
  
	G4Trd *MuonEarth = new G4Trd(name="MuonEarth", Earth_x1, Earth_x2, Earth_y1, Earth_y2, Earth_z); //Форма насыпи без полости
	G4Box *EarthHole = new G4Box(name="EarthHole", Earth_x2, Earth_y2, Hole_z); //Форма полости
  
	G4ThreeVector HolePosition = G4ThreeVector(0., 0., -Earth_z + Hole_z); //Координаты полости
	G4Transform3D HoleTr = G4Transform3D(rotm, HolePosition); //Задаём комбинацию поворота и координат полости
	G4SubtractionSolid *EarthCombined = new G4SubtractionSolid(name="EarthCombined", MuonEarth, EarthHole, HoleTr); //Вырезаем в насыпи полость
  
	G4LogicalVolume *EarthLogical = new G4LogicalVolume(EarthCombined, Granite, name="EarthLogical"); //Логический объём насыпи
  
	G4VPhysicalVolume *EarthPhysical = new G4PVPlacement(0, G4ThreeVector(-MuonXmove, -MuonYmove, Earth_z), EarthLogical, name="EarthPhysical", WorldLogical, false, 0, true); //Помещаем насыпь
    
	//Aluminum Muon Cover
    
	G4double R  = 1005./2. *mm;           // Ширина короба сцинтиллятора
	G4double H  = 55   *mm;           // Высота короба сцинтиллятора
	G4double r1 = 960./2.  *mm;           // Ширина основания диффузора
	G4double r2 = 260./2.  *mm;           // Ширина вершины диффузора
	G4double h  = 510  *mm;           // Высота диффузора+
	G4double dr = 0.55 *mm;           // Толщина стенок детектора+
	G4double upper_layer_height = 5 *mm;
	//
	G4double phiStart = 45 *deg;              // phiStart
	G4double phiTotal = 360*deg;              // phiTotal
	G4int    sides    = 4;                    // numbers of sides of each polygon in the x-y plane
	G4int    nZplanes = 6;                    // numbers of planes perpendicular to the z axis 
	G4double zval[nZplanes]  = {0 , 0+dr , H   , H+dr , h+H+dr , h+H+dr+upper_layer_height}; // zcoordinates of each plane
	G4double rmin[nZplanes]  = {0 , R    , R   , r1   , r2    , 0}; // Radii of inner and    
	G4double rmax[nZplanes]  = {R , R+dr , R+dr, r1+dr, r2+dr, r2+dr}; // outer polygon at each plane
	//
	G4Polyhedra* AluminumMuonCover = new G4Polyhedra(name="Polyhedra", phiStart, phiTotal, sides, nZplanes, zval, rmin, rmax); //Форма кожуха пластического сцинтиллятора
	//

	G4Tubs *PlasticHole = new G4Tubs(name="PlasticHole", hole_Rmin, hole_Rmax, upper_layer_height/2, hole_PhiMin, hole_PhiMax); //Форма отверстия под окно ФЭУ

	G4ThreeVector PlasticHolePosition = G4ThreeVector(0., 0., h+H+dr+upper_layer_height - upper_layer_height/2); //Координаты отверстия под окно ФЭУ
	G4Transform3D PlasticHoleTr = G4Transform3D(rotm, PlasticHolePosition); //Задаём комбинацию поворота и координат отверстия в кожухе под окно ФЭУ
  
	G4SubtractionSolid *AluminumCoverSolid = new G4SubtractionSolid(name="PlasticCover", AluminumMuonCover, PlasticHole, PlasticHoleTr); //Вырезаем отверстие под окно ФЭУ
    
	G4LogicalVolume* AluminumCoverLogical = new G4LogicalVolume(AluminumCoverSolid, Al, name="AluminumCoverLogical"); //Логический объём кожуха пластического сцинтиллятора
	G4LogicalSkinSurface *MuonCoverSurfaceSurface = new G4LogicalSkinSurface(name="MuonCoverSurface", AluminumCoverLogical, AlSurface); //Добавляем внутренней сотороне кожуха пластического сцинтиллятора такие же отражающие свойства, как у внутренней стороны кожуха жидкого сцинтиллятора
    
	//Plastic Scintillators
    
	density = 1.05*g/cm3; //плотность пластического сцинтиллятора
	G4Material *PlasticMaterial = new G4Material(name="Plastic Scintillator", density, ncomponents=2); //Создаём материал пластического сцинтиллятора
	PlasticMaterial->AddMaterial(elH, fractionmass=0.083333); //Добавляем водород
	PlasticMaterial->AddMaterial(elC, fractionmass=0.916667); //Добавляем кремний
  
	G4double RindexPlastic[nE] = {1.59, 1.59, 1.59}; //показатель преломления пластика
	G4double AbsLengthPlastic[nE] = {2000.*cm, 2000.*cm, 2000.*cm}; //absorption length пластика
	G4double ScintPlastic[nE] = {0.1, 1., 0.1}; //энергетический спектр сцинтилляционных фотонов
  
	G4MaterialPropertiesTable *MPTplastic = new G4MaterialPropertiesTable(); //таблица свойств пластика
  
	MPTplastic->AddProperty("SCINTILLATIONCOMPONENT1", Energy, ScintPlastic, nE); //энергетический спектр сцинтилляционных фотонов
	MPTplastic->AddConstProperty("SCINTILLATIONYIELD", 500./MeV); //количество фотонов на МэВ ионизации
	MPTplastic->AddProperty("RINDEX", Energy, RindexPlastic, nE); //показатель преломления
	MPTplastic->AddProperty("ABSLENGTH", Energy, AbsLengthPlastic, nE); //absorption length
	MPTplastic->AddConstProperty("RESOLUTIONSCALE", 0.); //Задаём нулевую флуктуацию числа фотонов, образующихся на каждом шаге
	MPTplastic->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 20.*ns); //время высвечивания
  
	PlasticMaterial->SetMaterialPropertiesTable(MPTplastic); //Присваиваем пластику указанные свойства
	MPTplastic->DumpTable(); //Отчистка памяти
    
	G4double plastic_x = R/2; //размер плстического сцинтиллятора по х
	G4double plastic_y = R/2; //размер плстического сцинтиллятора по у
	G4double plastic_z = (H - dr)/2; //размер плстического сцинтиллятора по z
  
	G4Box *PlasticSolid = new G4Box(name="Plastic Cell", plastic_x, plastic_y, plastic_z); //форма пластического сцинтиллятора
  
	G4LogicalVolume *PlasticLogical = new G4LogicalVolume(PlasticSolid, PlasticMaterial, name="Plastic Cell"); //логический объём пластического сцинтиллятора
    
	G4RotationMatrix* yRot = new G4RotationMatrix; //матрица поворота
	yRot->rotateY(M_PI*rad); //поворот кожуха на 180*
    
	MuonPhotoLogical = new G4LogicalVolume(PlasticHole, wmt, name="MuonPhotoLogical");  //Логический объём окна ФЭУ
	
	//Зполняем правый и средний тоннели мюонного детектора
    
	for(G4int w = 0; w < 2; w++)
		for(G4int i = 0; i < N_rows; i++)
			for(G4int j = 0; j < N_columns; j++){
        	
				G4double x = -(Room_x - Room_width - Fe_RoomWidth - (R+dr)) + 335./2.*mm + ((N_columns-1)*2*(R+dr) - j*2*(R+dr)) + (1-w)*2*Room_x - MuonXmove;
				//G4cout << G4endl << G4endl << x << G4endl << G4endl;
				G4double y = (Room_y - Room_width - Fe_RoomWidth - (R+dr)) - 1000./2.*mm - i*2*(R+dr) - MuonYmove;
				//G4cout << G4endl << G4endl << y << G4endl << G4endl;           
				G4double z = 2*Room_z - Room_depth - Fe_RoomWidth;           
				G4ThreeVector pos = G4ThreeVector(x, y ,z);
				G4VPhysicalVolume *AluminumCoverPhysical = new G4PVPlacement(yRot, pos, AluminumCoverLogical , name="AluminumCoverPhysical", WorldLogical, false, w*N_rows*N_columns + i*N_columns + j, true);
				G4VPhysicalVolume *PhotoMuonPhysical = new G4PVPlacement(0, G4ThreeVector(x, y, z - (h+H+dr+upper_layer_height)/2.), MuonPhotoLogical, name="MuonPhotoSensor", WorldLogical, false, w*N_rows*N_columns + i*N_columns + j, true);
      
				for(G4int k = 0; k < N_scint; k++)
					for(G4int l = 0; l < N_scint; l++){
  
						G4double x_pl = plastic_x * (2*k - 1);
						G4double y_pl = plastic_y * (2*l - 1);
						G4VPhysicalVolume *PlasticPhysical = new G4PVPlacement(0, G4ThreeVector(x + x_pl, y + y_pl, z-(dr + plastic_z)), PlasticLogical, name="PlasticPhysical", WorldLogical, false, (w*N_rows*N_columns + i*N_columns + j)*N_scint*N_scint + k+l*2, true);
  
					}
  
			} 
    
	//Floor

	ncomponents= 7;
	density= 2.5*g/cm3;
	G4Material *Basalt= new G4Material(name="Basalt", density, ncomponents); //материал земли
    
	G4Element *elO = nist->FindOrBuildElement("O"); //задаём материалы базальта
	G4Element *elNa = nist->FindOrBuildElement("Na"); //задаём материалы базальта
	G4Element *elMg = nist->FindOrBuildElement("Mg"); //задаём материалы базальта 
	G4Element *elAl = nist->FindOrBuildElement("Al"); //задаём материалы базальта
	G4Element *elSi = nist->FindOrBuildElement("Si"); //задаём материалы базальта
	G4Element *elCa = nist->FindOrBuildElement("Ca"); //задаём материалы базальта
	G4Element *elFe = nist->FindOrBuildElement("Fe"); //задаём материалы базальта
    
	Basalt->AddElement(elO ,fractionmass=0.600); //Добавляем материалы в соответствующей пропорции
	Basalt->AddElement(elNa,fractionmass=0.020); //Добавляем материалы в соответствующей пропорции
	Basalt->AddElement(elMg,fractionmass=0.025); //Добавляем материалы в соответствующей пропорции
	Basalt->AddElement(elAl,fractionmass=0.060); //Добавляем материалы в соответствующей пропорции
	Basalt->AddElement(elSi,fractionmass=0.175); //Добавляем материалы в соответствующей пропорции
	Basalt->AddElement(elCa,fractionmass=0.050); //Добавляем материалы в соответствующей пропорции
	Basalt->AddElement(elFe,fractionmass=0.070); //Добавляем материалы в соответствующей пропорции
    
	G4double RindexFloor[nE] = {1.66, 1.66, 1.66}; //показатель преломления базальта
    
	G4MaterialPropertiesTable *MPTfloor = new G4MaterialPropertiesTable(); //таблица свойств базальта
	MPTfloor->AddProperty("RINDEX", Energy, RindexFloor, nE); //показатель преломления
  
	Basalt->SetMaterialPropertiesTable(MPTfloor); //присваиваем базальту свойства, описанные выше
	MPTfloor->DumpTable(); //Отчитска памяти
    
	G4double floor_x = world_x; //размер земли по х
	G4double floor_y = world_y; //размер земли по у
	G4double floor_z = (world_z)/2; //размер земли по z
    
	G4double floor_hole_x = 3*Room_x + 2*Edge_x; //размер полости под мюонный детектор по х
	G4double floor_hole_y = 43000./2.*mm; //размер полости под мюонный детектор по у
	G4double floor_hole_z = (Room_depth + Room_width_z)/2.; //размер полости под мюонный детектор по z
    
	G4Box *FloorSolid = new G4Box(name="FloorSoild", floor_x, floor_y, floor_z); //форма земли
	G4Box *FloorMuonHole = new G4Box(name="FloorMuonHole", floor_hole_x, floor_hole_y, floor_hole_z); //форма полости под мюонный детектор
    
	G4ThreeVector FloorHolePosition = G4ThreeVector(-MuonXmove, -MuonYmove, floor_z - floor_hole_z); //координаты полости под мюонный детектор
	G4Transform3D FloorHoleTr = G4Transform3D(rotm, FloorHolePosition); //Задаём комбинацию поворота и координат полости под мюонный детектор
	G4SubtractionSolid *FloorCombined = new G4SubtractionSolid(name="FloorCombined", FloorSolid, FloorMuonHole, FloorHoleTr); //Вырезаем полость в земле
    
	G4LogicalVolume *FloorLogical = new G4LogicalVolume(FloorCombined, Basalt, name="FloorLogical"); //логический объём земли
    
	G4VPhysicalVolume *FloorPhysical = new G4PVPlacement(0, G4ThreeVector(0, 0, -floor_z), FloorLogical, name="Floor", WorldLogical, false, 0, true); //физический объём земли
    
	return WorldPhysical;
}

//В ConstructSDandField задаются чувствительные детекторы и параметры электромагнитного поля 

void CarpetConstruction::ConstructSDandField(){
  
  SensitiveDetector *PhotoDetector = new SensitiveDetector("PhotoDetector"); //Задаём параметры окон ФЭУ жидкого сцинтиллятора
  G4SDManager::GetSDMpointer()->AddNewDetector(PhotoDetector); //Задаём параметры окон ФЭУ жидкого сцинтиллятора
  PhotoLogical->SetSensitiveDetector(PhotoDetector); //Задаём параметры окон ФЭУ жидкого сцинтиллятора
  
  MuonSensitiveDetector *MuonPhotoDetector = new MuonSensitiveDetector("MuonPhotoDetector"); //Задаём параметры окон ФЭУ пластического сцинтиллятора
  G4SDManager::GetSDMpointer()->AddNewDetector(MuonPhotoDetector); //Задаём параметры окон ФЭУ пластического сцинтиллятора
  MuonPhotoLogical->SetSensitiveDetector(MuonPhotoDetector); //Задаём параметры окон ФЭУ пластического сцинтиллятора
}

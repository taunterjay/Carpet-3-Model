

//В блоке Action задаются действия, которые совершаются в каждом прогоне, событии и шаге


#include "Action.hh"

//Конструктор
Action::Action(){

}

//Деструктор
Action::~Action(){

}

//BuildForMaster позволяет работать в режиме многопоточности. Здесь указывается только действие, выполняемое в каждом прогоне.
   
void Action::BuildForMaster() const {
  

	MyRunAction *runAction = new MyRunAction(); //Создаём действие, выполняемое в каждом прогоне
	SetUserAction(runAction); //Задаём действие, выполняемое в каждом прогоне	
 
}

//Build запускается на каждом ядре. Здесь задаются все нужные нам действия
   
void Action::Build() const {
  
	PrimaryGenerator *generator = new PrimaryGenerator(); //Создаём генератор первичынх частиц
	SetUserAction(generator); //Задаём генератор первичынх частиц
  
	MyRunAction *runAction = new MyRunAction(); //Создаём действие, выполняемое в каждом прогоне
	SetUserAction(runAction); //Задаём действие, выполняемое в каждом прогоне
  
	MySteppingAction *steppingAction = new MySteppingAction(new G4UserEventAction()); //Создаём действие, выполняемое на каждом шаге
	SetUserAction(steppingAction); //Задаём действие, выполняемое на каждом шаге
    
}

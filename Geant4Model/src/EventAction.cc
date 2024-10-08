#include "EventAction.hh"

//В Event можно указать действие, выполняемое в каждом событии. Мне не пригодилось.


//Конструктор
MyEventAction::MyEventAction(MyRunAction* runAction){

}


//Деструктор
MyEventAction::~MyEventAction(){

}

//BeginOfEventAction выполняется в начале события

void MyEventAction::BeginOfEventAction(const G4Event* evt){
  
}

//BeginOfEventAction выполняется в конце события

void MyEventAction::EndOfEventAction(const G4Event*){

}

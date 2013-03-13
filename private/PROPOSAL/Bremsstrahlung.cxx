#include "PROPOSAL/Bremsstrahlung.h"

Bremsstrahlung::Bremsstrahlung(){

}
//----------------------------------------------------------------------------//

Bremsstrahlung::Bremsstrahlung(const Bremsstrahlung &brems)
{
    *this = brems;
}
//----------------------------------------------------------------------------//

Bremsstrahlung& Bremsstrahlung::operator=(const Bremsstrahlung &brems){
    return *this;
}

//----------------------------------------------------------------------------//

void Bremsstrahlung::SetIntegralLimits(){
}

//----------------------------------------------------------------------------//

double Bremsstrahlung::CalculatedEdx(){
    return 0;
}
//----------------------------------------------------------------------------//

double Bremsstrahlung::CalculatedNdx(){
    return 0;
}
//----------------------------------------------------------------------------//

double Bremsstrahlung::CalculatedNdx(double rnd){
    return 0;
}
//----------------------------------------------------------------------------//

double Bremsstrahlung::CalculateStochasticLoss(){
    return 0;
}
//----------------------------------------------------------------------------//

void Bremsstrahlung::EnableStochasticInerpolation(){
    doStochasticInterpolation_=true;
}
//----------------------------------------------------------------------------//

void Bremsstrahlung::EnableContinuousInerpolation(){
    doContinuousInterpolation_=true;
}

//----------------------------------------------------------------------------//

double Bremsstrahlung::FunctionToContinuousIntegral(double variable){
    return 0;
}

//----------------------------------------------------------------------------//

double Bremsstrahlung::FunctionToStochasticalIntegral(double variable){
    return 0;
}
//----------------------------------------------------------------------------//

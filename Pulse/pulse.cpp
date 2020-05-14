#include "pulse.hpp"

Pulse::Pulse(double w_, double width = 10.0, double dt_ = 0.01, int numSteps_ = 15000):{
    tu = 1.05457/1.60218; // in unit of fs
    Eu = 0.160218/3/1.055; // in unit of 10^8 V/m
    Au = Eu*tu*1e-7; // in unit of V/m*s
    Lu = 300.0*tu; // nm
    a = 0.2429/Lu;

    params.resize(4);
    params.at(0) = 10.0/Eu; // default E0 10*10^8 V/m = 0.1V/am
    params.at(1) = width * std::sqrt(2)/tu;
    params.at(2) = w_;
    params.at(3) = PI/2.0;
    dt = dt_;
    numSteps = numSteps_;

    tc = dt*numSteps/2.0;
    ti = 0.0-tc;
    tf = dt*numSteps-tc;

    FuncE.function = &GaussPulse;
    FuncE.params = params.data();

    FuncE2.function = &GaussPulse2;
    FuncE2.params = params.data();
    epsabs = 1e-10;
    epsrel = 1e-6;

    Fluence = computeFluence();

}

double Pulse::getE(int stepIdx){
    double t = stepIdx*dt-tc;
    return GaussPulse(t,params.data());
}

double Pulse::getA(int stepIdx){
    double t = stepIdx*dt-tc;
    double result, abserr;
    size_t neval;
    gsl_integration_qng(&FuncE, ti, t, epsabs, epsrel, &result, &abserr, &neval);
    return result;
}

double Pulse::computeFluence(){
    double result, abserr;
    size_t neval;
    gsl_integration_qng(&FuncE2, ti, tf, epsabs, epsrel, &result, &abserr, &neval);
    return 3.0*8.85*tu*Eu*Eu*result*1e-7;
}

void Pulse::setFluence(double Flu){
    params.at(0) *= std::sqrt(Flu/Fluence); 
    double Flu_ckeck = computeFluence();
    assert_msg(std::abs(Flu-Flu_check)<1e-6,"Err bigger than 1e-6 in Pulse::setFluence!");
    Fluence = Flu;
}

void Pulse::print() const {
    std::cout<<"Pulse Info:\n";
    std::cout<<"E0:"<<params.at(0)*Eu<<"*10^8 V/m\n"\
             <<"sigma:"<<params.at(1)*tu<<"fs\n"\
             <<"w:"<<params.at(2)<<"eV\n"\
             <<"phase:"<<params.at(3)<<"\n"\
             <<"Fluence:"<<Fluence<<"mJ/cm^2\n";
}
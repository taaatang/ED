#include "pulse.hpp"

double GaussPulse(double t, void* params){
    double *Params = (double *) params;
    double amp = Params[0], sigma = Params[1], freq = Params[2], phase = Params[3];
    return amp*std::exp(-t*t/sigma/sigma)*std::sin(freq*t+phase);
}

double GaussPulse2(double t, void* params){
    double E = GaussPulse(t, params);
    return E*E;
}

Pulse::Pulse(double w_, double width, double dt_ , int numSteps_){
    A = 0.0;

    tu = 1.05457/1.60218; // in unit of fs
    Eu = 0.160218/3/1.055; // in unit of 10^8 V/m
    Au = Eu*tu*1e-7; // in unit of V/m*s
    Lu = 300.0*tu; // nm
    a = 0.378/Lu; // Cu-Ox/y = 1.89A, Cu-Oz = 2.429A

    params = VecD(4, 0.0);
    params.at(0) = 10.0/Eu; // default E0 10*10^8 V/m = 0.1V/am
    params.at(1) = width * std::sqrt(2)/tu;
    params.at(2) = w_;
    params.at(3) = PI/2.0;
    dt = dt_;
    numSteps = numSteps_;

    tc = dt*numSteps/2.0;
    ti = 0.0-tc;
    tf = dt*numSteps-tc;
    t = ti;

    FuncE.function = &GaussPulse;
    FuncE.params = (void *)params.data();

    FuncE2.function = &GaussPulse2;
    FuncE2.params = (void *)params.data();
    epsabs = 1e-8;
    epsrel = 1e-8;

    Fluence = computeFluence();
}

void Pulse::setFuncPara( ) {
    FuncE.params = (void *)params.data();
    FuncE2.params = (void *)params.data(); 
}

double Pulse::getE(int stepIdx) const {
    double t = stepIdx*dt-tc;
    return GaussPulse(t,(void *)params.data());
}

double Pulse::getA() const {
    double result, abserr;
    size_t neval;
    gsl_integration_qng(&FuncE, t, t+dt, epsabs, epsrel, &result, &abserr, &neval);
    A += result;
    return A;
}

double Pulse::computeFluence() const {
    double result_sum=0.0;
    double intval = 2*PI/params.at(2);
    double result, abserr;
    size_t neval;
    for(double t_start = ti; t_start < tf; t_start += intval){
        double t_next = t_start + intval;
        double t_end = t_next<tf?t_next:tf;
        gsl_integration_qng(&FuncE2, t_start, t_end, epsabs, epsrel, &result, &abserr, &neval);
        result_sum += result;
    }
    return 3.0*8.85*tu*Eu*Eu*result_sum*1e-4;
}

void Pulse::setFluence(double Flu){
    params.at(0) *= std::sqrt(Flu/Fluence); 
    double Flu_check = computeFluence();
    // std::cout<<"E0:"<<params.at(0)<<".\n";
    // std::cout<<"Original Fluence:"<<Fluence<<", Traget:"<<Flu<<", result:"<<Flu_check<<"\n";
    assert_msg(std::abs(Flu-Flu_check)<1e-6,"Err bigger than 1e-6 in Pulse::setFluence!");
    Fluence = Flu;
}

void Pulse::print(std::ostream& os) const {
    os<<"Pulse Info:\n";
    os<<"E0:"<<params.at(0)*Eu<<"*10^8 V/m\n"\
      <<"sigma:"<<params.at(1)*tu<<" fs\n"\
      <<"w:"<<params.at(2)<<" eV\n"\
      <<"phase:"<<params.at(3)<<"\n"\
      <<"Fluence:"<<Fluence<<" mJ/cm^2\n";
}
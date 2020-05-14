//
//  pulse.hpp
//  ED
//
//  Created by tatang on 10/27/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#ifndef pulse_hpp
#define pulse_hpp

#include "../Global/globalPara.hpp"
#include "../Utils/utils.hpp"
#include <gsl/gsl_integration.h>
#include <cmath>

double GaussPulse(double t, void* params){
    double *Params = (double *) params;
    double amp = Params[0], sigma = Params[1], freq = Params[2], phase = Params[3];
    return amp*std::exp(-t*t/sigma/sigma)*std::sin(freq*t+phase);
}

double GaussPulse2(double t, void* params){
    double E = GaussPulse(t, params);
    return E*E;
}

class Pulse{
/*
    Gaussian shape elextric field with frequency w.
    return vector potential A(t) = \int_{ti}^{t} E(t).
    unit: hbar = 1, eV = 1, c = 1, e = 1
          tu = hbar/eV, Lu = c*tu, Eu = eV/e/Lu
*/
public:
    Pulse(double w_/*eV*/, double width = 10.0/*fs*/, double dt_ = 0.01/*tu*/, int numSteps_ = 15000);
    ~Pulse(){};

    double getE(int stepIdx) const;
    double getFluence() const {return Fluence;}
    double computeFluence() const;
    double getA(int stepIdx) const;
    double getAa(int stepIdx) const {return a*getA(stepIdx);} // return A*a

    setE0(double E0){params.at(0)=E0;}
    setWidth(double width){params.at(1)=width*std::sqrt(2);}
    setW(double w){params.at(2)=w;}
    setPhase(double phase){params.at(3)=phase;}
    void setFluence(double Flu /* mJ/cm^2 */);

    void seta(double a_/*nm*/){a = a_/Lu;}

    void print() const;

private:
    double E, A;

    double Eu, Au; // electric field unit
    std::vector<double> params; // params = {E0, sigma, w, phase}, sigma = sqrt(2)*width, where I(width) = I(0)/e
    double Fluence;

    double Lu, a;

    double dt;
    int numSteps;
    double tu; // time unit: hbar/eV, fs
    double ti,tc,tf; 

    gsl_function FuncE;
    gsl_function FuncE2;
    double epsabs, epsrel;
};

#endif // pulse_hpp
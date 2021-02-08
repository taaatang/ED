//
//  pulse.hpp
//  ED
//
//  Created by tatang on 10/27/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#ifndef pulse_hpp
#define pulse_hpp

#include <gsl/gsl_integration.h>
#include <cmath>
#include <iostream>

#include "Global/globalPara.hpp"
#include "Utils/utils.hpp"

double GaussPulse(double t, void* params);

double GaussPulse2(double t, void* params);

class Pulse{
/*
    Gaussian shape elextric field with frequency w.
    return vector potential A(t) = \int_{ti}^{t} E(t).
    unit: hbar = 1, eV = 1, c = 1, e = 1
          tu = hbar/eV, Lu = c*tu, Eu = eV/e/Lu
*/
public:
    Pulse(double w_ = 0.0/*eV*/, double width = 10.0/*fs*/, double dt_ = 0.01/*tu*/, int numSteps_ = 15000);
    ~Pulse(){};

    double getE(int stepIdx) const;
    double getFluence( ) const { return Fluence; }
    double computeFluence( ) const;
    double getA( ) const;
    double getAa( ) const { return a*getA(); } // return A*a
    VecD getPol( ) const { return pol; }

    void setE0(double E0) { params.at(0)=E0; }
    void setWidth(double width) { params.at(1)=width*std::sqrt(2); }
    void setW(double w) { params.at(2)=w; }
    void setPhase(double phase) { params.at(3)=phase; }
    void setFluence(double Flu /* mJ/cm^2 */);
    void setPol(VecD pol) { this->pol = pol;}
    void setFuncPara( );

    void seta(double a_/*nm*/) { a = a_/Lu; }

    bool next( ) { ++count; t += dt; return count<numSteps;}

    void print(std::ostream& os = std::cout) const;

private:
    mutable double E, A;

    double Eu, Au; // electric field unit
    VecD pol; // polarization in x-y-z coordinates
    VecD params; // params = {E0, sigma, w, phase}, sigma = sqrt(2)*width, where I(width) = I(0)/e
    double Fluence;

    double Lu, a;

    double dt;
    int count{0};
    int numSteps{0};
    double tu; // time unit: hbar/eV, fs
    double ti,tc,tf; 
    mutable double t;

    gsl_function FuncE;
    gsl_function FuncE2;
    double epsabs, epsrel;
};

#endif // pulse_hpp
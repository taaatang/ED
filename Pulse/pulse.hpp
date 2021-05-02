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
#include "Geometry/Vec3.hpp"
#include "Utils/io.hpp"

double GaussPulse(double t, void* params);

double GaussPulse2(double t, void* params);

double GaussPulseDerivative(double t, void* params);

double GaussPulseDerivative2(double t, void* params);

class Pulse{
/*
    Gaussian shape elextric field with frequency w.
    return vector potential A(t) = \int_{ti}^{t} E(t).
    unit: hbar = 1, eV = 1, c = 1, e = 1
          tu = hbar/eV, Lu = c*tu, Eu = eV/e/Lu
*/
public:
    Pulse(double w_ = 0.0/*eV*/, double sigma = 50.0/*fs*/, double dt_ = 0.01/*tu*/, double duration_ = 10.0/* in unit of sqrt(2)*sigma */, bool useA_ = true);
    ~Pulse(){};

    double getE(int stepIdx) const;
    double getw( ) const { return params[2]; }
    double getFluence( ) const { return Fluence; }
    double computeFluence( ) const;
    double getA( ) const;
    double getA(double t_) const;
    double getAa( ) const { return a * getA(); } // return A*a
    Vec3d getPol( ) const { return pol; }
    double getdt( ) const { return dt; }
    int getCount( ) const { return count; }
    int getStepNum( ) const { return numSteps; }

    void setAmp(double A0) { params[0] = A0; }
    void setWidth(double sigma_, double duration_);
    void setFreq(double w) { params[2] = w; }
    void setPhase(double phase) { params[3] = phase * 2.0 * PI; }
    void setFluence(double Flu /* mJ/cm^2 */);
    void setPol(Vec3d pol) { this->pol = pol;}
    void setFuncPara( );

    void seta(double a_/*nm*/) { a = a_/Lu; }

    bool next( ) { ++count; t += dt; return count <= numSteps;}

    void progressBar(int total_) const;

    void progress( ) const;

    void profile(int n = 21) const;

    void print(std::ostream& os = std::cout) const;

private:
    bool useA;
    mutable double E, A;

    double Eu, Au; // electric field unit
    Vec3d pol; // polarization in x-y-z coordinates
    double params[4]; // params = {E0, sigma, w, phase}, sigma = sqrt(2)*width, where I(width) = I(0)/e
    double Fluence;

    double Lu, a;

    double dt;
    int count{0};
    int numSteps{0};
    mutable int total{1};
    double tu; // time unit: hbar/eV, fs
    double ti, tc, tf; 
    mutable double t;

    gsl_function FuncE;
    gsl_function FuncE2;
    double epsabs{1e-8}, epsrel{1e-8};
};

#endif // pulse_hpp
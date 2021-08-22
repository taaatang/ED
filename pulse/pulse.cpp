#include <cmath>
#include "pulse.hpp"
#include "utils/runtimeCheck.hpp"

/**
 * @brief A0 * exp(-t^2 / sigma^2) * sin(fre * t + phase)
 * 
 * @param t 
 * @param params [A0, sigma, fre, phase]. defined as void* for gsl function object
 * @return double 
 */
double GaussPulse(double t, void* params){
    double *Params = (double *) params;
    double amp = Params[0], sigma = Params[1], freq = Params[2], phase = Params[3];
    return amp * std::exp(-t * t / sigma / sigma) * std::sin(freq * t + phase);
}

double GaussPulse2(double t, void* params){
    double E = GaussPulse(t, params);
    return E * E;
}

/**
 * @brief -d/dt (GaussPulse). E  = -dA/dt
 * 
 * @param t 
 * @param params GaussPulse params
 * @return double 
 */
double GaussPulseDerivative(double t, void* params) {
    double *Params = (double *) params;
    double amp = Params[0], sigma = Params[1], fre = Params[2], phase = Params[3];
    return amp * std::exp(-t * t / sigma / sigma) * ((2.0 * t / sigma / sigma) * std::sin(fre * t + phase) - fre * std::cos(fre * t + phase));
}

double GaussPulseDerivative2(double t, void* params) {
    double amp = GaussPulseDerivative(t, params);
    return amp * amp;
}

Pulse::Pulse(double w_, double sigma_, double dt_ , double duration_, bool useA_){
    useA = useA_;
    A = 0.0;

    tu = 1.05457/1.60218; // in unit of fs
    Eu = 0.160218/3/1.055; // in unit of 0.01 V/am
    Au = Eu*tu; // in unit of 0.01 V/am * fs
    Lu = 300.0*tu; // nm
    a = 0.378/Lu; // Cu-Ox/y = 1.89am, Cu-Oz = 2.429am

    if (useA) {
        params[0] = 1.0 / Au;
    } else {
        params[0] = 1.0 / Eu; // default E0 10^8 V/m = 0.01V/am
    }
    params[1] = sigma_ * std::sqrt(2) / tu;
    params[2] = w_;
    params[3] = PI / 2.0;
    dt = dt_;
    numSteps = std::ceil(duration_ * params[1] / dt_);

    tc = dt * numSteps / 2.0;
    ti = 0.0 - tc;
    tf = dt * numSteps - tc;
    t = ti;

    if (useA) {
        FuncE.function = &GaussPulseDerivative;
        FuncE2.function = &GaussPulseDerivative2;
    } else {
        FuncE.function = &GaussPulse;
        FuncE2.function = &GaussPulse2;
    }
    FuncE.params = params;
    FuncE2.params = params;
    Fluence = computeFluence();
}

void Pulse::setFuncPara( ) {
    FuncE.params = params;
    FuncE2.params = params; 
}

void Pulse::profile(int n) const {
    double step = (tf - ti) / n;
    for (int i = 0; i <= n; ++i) {
        std::cout<<"step "<<i<<", A = "<<getA(ti + i * step) * Au * 0.01<<" V/am*fs\n";
    }
}

double Pulse::getE(int stepIdx) const {
    double t = stepIdx*dt-tc;
    if (useA) {
        return GaussPulseDerivative(t, (void *)params);
    } else {
        return GaussPulse(t, (void *)params);
    }
}

double Pulse::getA() const {
    if (useA) {
        return GaussPulse(t, (void *)params);
    } else {
        double result, abserr;
        size_t neval;
        gsl_integration_qng(&FuncE, t, t+dt, epsabs, epsrel, &result, &abserr, &neval);
        A += result;
        return A;
    }
}

double Pulse::getA(double t_) const {
    if (useA) {
        return GaussPulse(t_, (void *)params);
    } else {
        double result, abserr;
        size_t neval;
        double A_ = 0.0;
        double t1 = ti;
        for (; t1 < t_-dt; t1 += dt) {
            gsl_integration_qng(&FuncE, t1, t1+dt, epsabs, epsrel, &result, &abserr, &neval);
            A_ += result;
        }
        gsl_integration_qng(&FuncE, t1, t_, epsabs, epsrel, &result, &abserr, &neval);
        A_ += result;
        return A_;
    } 
}

void Pulse::setWidth(double sigma_, double duration_) {
    params[1] = sigma_ * std::sqrt(2);
    numSteps = std::ceil(duration_ * params[1] / dt);

    tc = dt * numSteps / 2.0;
    ti = 0.0 - tc;
    tf = dt * numSteps - tc;
    t = ti;
}

double Pulse::computeFluence() const {
    double result_sum=0.0;
    double intval = 2*PI/params[2];
    double result, abserr;
    size_t neval;
    for(double t_start = ti; t_start < tf; t_start += intval){
        double t_next = t_start + intval;
        double t_end = t_next  < tf ?  t_next : tf;
        gsl_integration_qng(&FuncE2, t_start, t_end, epsabs, epsrel, &result, &abserr, &neval);
        result_sum += result;
    }
    return 3.0 * 8.85 * tu * Eu * Eu * result_sum * 1e-4;
}

void Pulse::setFluence(double Flu){
    params[0] *= std::sqrt(Flu/Fluence); 
    double Flu_check = computeFluence();
    // std::cout<<"E0:"<<params.at(0)<<".\n";
    // std::cout<<"Original Fluence:"<<Fluence<<", Traget:"<<Flu<<", result:"<<Flu_check<<"\n";
    assert_msg(std::abs(Flu-Flu_check)<1e-6,"Err bigger than 1e-6 in Pulse::setFluence!");
    Fluence = Flu;
}

void Pulse::print(std::ostream& os) const {
    if (useA) {
        os<<"Pulse A0";
    } else {
        os<<"Pulse E0";
    }
    os<<" * exp(-(t^2)/2/sigma^2) * sin(w * t + phase):\n";
    if (useA) {
        os<<"A0: "<<params[0] * Au * 0.01<<" V/am*fs\n";
    } else  {
        os<<"E0: "<<params[0] * Eu * 0.01<<" V/am\n";
    }
    os<<"sigma: "<<params[1] * tu / std::sqrt(2)<<" fs\n"\
      <<"w: "<<params[2]<<" eV\n"\
      <<"phase: "<<params[3]<<"\n"\
      <<"Fluence: "<<Fluence<<" mJ/cm^2\n";
}
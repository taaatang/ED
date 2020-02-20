//
// HelperClass.hpp
// ED
//
// Created by tatang on 9/2/20
//

#ifndef HelperClass_hpp
#define HelperClass_hpp

#include <chrono>

class Timer{
    std::chrono::system_clock::time_point tik_, tok_;
    std::chrono::duration<double> duration_;
    bool is_tok;
public:
    Timer():is_tok(false){tik();}
    ~Timer(){}

    void tik(){tik_ = std::chrono::system_clock::now(); is_tok = false;}
    void tok(){tok_ = std::chrono::system_clock::now(); is_tok = true;}
    double elapse(){
        if (!is_tok) tok();
        duration_ = tok_ - tik_;
        return duration_.count()*1000.0;
    }
};

#endif
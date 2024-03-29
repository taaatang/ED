/**
 * @file timer.hpp
 * @author Ta Tang (tatang.physics@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2020-09-02
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#pragma once

#include <chrono>
#include <iostream>
#include <string>

class Timer{
    std::chrono::system_clock::time_point tik_, tok_;
    std::chrono::duration<double> duration_;
    bool is_tok{false};
    bool talk{true};
public:
    Timer( ) { }
    Timer(bool isTalk):talk(isTalk){ }
    ~Timer( ) { }

    void tik( ) { tik_ = std::chrono::system_clock::now(); is_tok = false; }
    void tok( ) { tok_ = std::chrono::system_clock::now(); is_tok = true; }
    // return duration in unit of second
    double elapse( ) {
        if (!is_tok) tok();
        duration_ = tok_ - tik_;
        is_tok = false;
        return duration_.count();
    }
    void print(std::string event, std::ostream& os = std::cout) {
        if (talk) os << event << " time:" << elapse() << "s." << std::endl;
    }
};
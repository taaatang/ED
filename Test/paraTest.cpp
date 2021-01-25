#include "../Utils/paras.hpp"

int main(){
    Parameters para;
    para.read("../Input/lattice.txt");
    para.read("../Input/HubbardMultiband.txt");
    para.print(std::cout);
    std::cout<<3.0<<"\n";
    return 0;
}
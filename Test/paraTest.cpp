#include "../Utils/paras.hpp"

int main(){
    Parameters para;
    para.read("../Input/lattice.txt");
    para.read("../Input/HubbardMultiband.txt");
    std::ofstream outfile("para.txt");
    para.print(outfile);
    return 0;
}
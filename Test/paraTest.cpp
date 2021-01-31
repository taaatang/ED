#include "../Utils/paras.hpp"

int main(){
    Parameters para("../Input/config.txt");
    std::ofstream outfile("para.txt");
    para.print(outfile);
    return 0;
}
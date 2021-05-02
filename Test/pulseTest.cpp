#include "Pulse/pulse.hpp"
using namespace std;
int main( ) {
    /**
     * @brief 
     * 
     */
    Pulse pulse(1.5, 10.0, 0.01, 10.0, true);
    pulse.setPhase(PI / 2.0);
    pulse.setFluence(1.0);
    pulse.profile(31);
    pulse.print();
    cout<<5/2<<"\n";
    cout<<5%2<<"\n";
    return 0;
}
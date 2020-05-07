//
//  geometryTest.cpp
//  ED
//
//  Created by tatang on 8/4/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#include "../Global/globalPara.hpp"
#include "../Utils/utils.hpp"
#include "../Geometry/Geometry.hpp"

#include <iostream>

int main(int argc, const char * argv[]) {
    // int N = 12;
    // TriAngLattice Lattice(N);
    // Lattice.addOrb({});
    // int Nx=2, Ny=2;
    SquareLattice Lattice(8);
    Lattice.addOrb({ORBITAL::Dx2y2,0,{0.0,0.0,0.0}}).addOrb({ORBITAL::Px,1,{0.5,0.0,0.0}}).addOrb({ORBITAL::Py,2,{0.0,0.5,0.0}});
    Lattice.addOrb({ORBITAL::Pzu,3,{0.0,0.0,0.5}}).addOrb({ORBITAL::Pzd,4,{0.0,0.0,-0.5}});
    Lattice.construct();
    Lattice.print();
    return 0;
}
//
//  geometryTest.cpp
//  ED
//
//  Created by tatang on 8/4/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#include "../global/globalPara.hpp"
#include "../geometry/geometry.hpp"
#include "../operator/operatorsBase.hpp"

#include <iostream>

int main(int argc, const char * argv[]) {
    int N = 12;
    TriAngLattice Lattice(N);
    // Lattice.addOrb({});
    // int Nx=2, Ny=2;
    // SquareLattice Lattice(8);
    // Lattice.addOrb({ORBITAL::Dx2y2,0,{0.0,0.0,0.0}}).addOrb({ORBITAL::Px,1,{0.5,0.0,0.0}}).addOrb({ORBITAL::Py,2,{0.0,0.5,0.0}});
    // Lattice.addOrb({ORBITAL::Pzu,3,{0.0,0.0,0.5}}).addOrb({ORBITAL::Pzd,4,{0.0,0.0,-0.5}});
    Lattice.addOrb({});
    Lattice.construct();
    Lattice.print();

    dataType Jk = 1.0;
    Link<dataType> JkLink(LINK_TYPE::CHIRAL_K, {ORBITAL::SINGLE, ORBITAL::SINGLE, ORBITAL::SINGLE},Jk);
    JkLink.addLinkVec({0.0,1.0,0.0}).addLinkVec({1.0,0.0,0.0}).\
        addLinkVec({1.0,0.0,0.0}).addLinkVec({1.0,-1.0,0.0});
    JkLink.genLinkMaps(&Lattice);
    JkLink.print();
    return 0;
}
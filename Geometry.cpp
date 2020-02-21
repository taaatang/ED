//
//  Geometry.cpp
//  ED
//
//  Created by tatang on 10/26/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#include "Geometry.hpp"


/*
   ***************************
   *   Geometry Base Class   *
   ***************************
*/

bool Geometry::coordToOrbid(double* coord, int &orbid, double tol) const {
    bool cond;
    for (auto it = enlgOrbs.begin(); it != enlgOrbs.end(); it++){
        cond = true;
        for (int i = 0; i < dim; i++){
            if(std::abs(coord[i]-(*it).coord[i])>tol) {cond = false; break;}
        }
        if (cond) {orbid = (*it).id; return cond;}
        else{continue;}
    }
    return false;
}

void Geometry::genTransList(){
    VecD coordi(3), coordr(3), coordf(3);
    TransList.clear();
    std::vector<int> tmp;
    for (int r = 0; r < getSiteNum(); r++){
        tmp.clear();
        getSiteR(r,coordr.data());
        for(int orbidi = 0; orbidi < getOrbNum(); orbidi++){
            getOrbR(orbidi,coordi.data());
            vecXAdd(1.0, coordi.data(), 1.0, coordr.data(), coordf.data(), dim);
            int orbidf;
            if (coordToOrbid(coordf.data(), orbidf)){
                tmp.push_back(orbidf);
            }else{
                std::cout<<"translation position not found for orbid = "<<orbidi<<", transVecid = "<<r<<std::endl;
                exit(1);
            }
        }
        TransList.push_back(tmp);
    }
}

void Geometry::construct(){
    Norb = Nsite * unitSite.size();
    Norb_enlg = is_PBC?Norb*TranVecs.size():Norb;
    assert(Nsite>0 and Norb>0 and Norb_enlg>=Norb);
    assert(xlist.size()==getSiteNum() and ylist.size()==getSiteNum() and zlist.size()==getSiteNum());
    VecD vsite(3,0.0);
    VecD vorb(3,0.0);
    int id = 0;
    // construc Lattice and orbs
    for (int siteid = 0; siteid < getSiteNum(); siteid++){
        vecXAdd(xlist.at(siteid), a1.data(), ylist.at(siteid), a2.data(), zlist.at(siteid), a3.data(), vsite.data(), getDim());
        Lattice.push_back(Site{siteid, vsite});
        for (int i = 0; i < unitSite.size(); i++){
            vecXAdd(1.0, vsite.data(), 1.0, unitSite.at(i).coord.data(), vorb.data(), getDim());
            orbs.push_back(Orbital{unitSite.at(i).orb,unitSite.at(i).orbid,vorb});
            orbs.at(id).siteid = siteid;
            orbs.at(id).id = id;
            id++;
        }
    }
    // construc enlarged orbs (periodic boundary condition)
    int count = 0;
    for (int R = 0; R < TranVecs.size(); R++){
        for (int r = 0; r < orbs.size(); r++){
            vecXAdd(1.0, TranVecs.at(R).data(), 1.0, orbs.at(r).coord.data(), vorb.data(), getDim());
            enlgOrbs.push_back(Orbital{orbs.at(r).orb, orbs.at(r).orbid, vorb});
            enlgOrbs.at(count).id = orbs.at(r).id;
            enlgOrbs.at(count).siteid = orbs.at(r).siteid;
            count++;
        }
    }
    assert(enlgOrbs.size()==Norb_enlg);
    // construct corresponding k-space lattice
    if(is_PBC){
        assert(kxlist.size()==getSiteNum() and kylist.size()==getSiteNum() and kzlist.size()==getSiteNum());
        for (int siteid = 0; siteid < getSiteNum(); siteid++){
            vecXAdd(kxlist.at(siteid), b10.data(), kylist.at(siteid), b20.data(), kzlist.at(siteid), b30.data(), vsite.data(), getDim());
            KLattice.push_back(Site{siteid, vsite});
        }
        // generate transformation matrix for translation operation
        genTransList();
    } 

    // return resources no longer needed
    xlist.clear(); ylist.clear(); zlist.clear();
    kxlist.clear(); kylist.clear(); kzlist.clear();
}

void Geometry::printLattice() const {
    // print Lattice
    std::cout<<"Lattice Unit Cell:"<<std::endl;
    for (int i = 0; i < Nsite; i ++){
        std::cout<<"siteID:"<<Lattice[i].id<<", coord:["<<Lattice[i].coord[0]<<", "<<Lattice[i].coord[1]<<", "<<Lattice[i].coord[2]<<"]"<<std::endl;
    }
    std::cout<<std::endl;
}
void Geometry::printOrbs() const {
    std::cout<<"Lattice Orbitals:"<<std::endl;
    for (int i = 0; i < Norb; i ++){
        std::cout<<"orbID:"<<orbs[i].id<<", from site"<<orbs[i].siteid<<", coord:["<<orbs[i].coord[0]<<", "<<orbs[i].coord[1]<<", "<<orbs[i].coord[2]<<"]"<<std::endl;
    }
    std::cout<<std::endl;
}
void Geometry::printKLattice() const {
    // print kLattice
    std::cout<<"k-Lattice:"<<std::endl;
    for (int i = 0; i < Nsite; i ++){
        std::cout<<"siteID:"<<KLattice[i].id<<", coord:["<<KLattice[i].coord[0]<<", "<<KLattice[i].coord[1]<<", "<<KLattice[i].coord[2]<<"]"<<std::endl;
    }
    std::cout<<std::endl;
}
void Geometry::printTrans() const {
    std::cout<<"TransList:"<<std::endl;
    for (int r = 0; r < Nsite; r++){
        std::cout<<"r("<<r<<"):"<<std::endl;
        for (int i = 0; i < Norb; i++){
            std::cout<<getOrbTran(r,i)<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}
void Geometry::print() const {
    printLattice();
    printOrbs();
    printTrans();
    printKLattice();
}

/*
    **********************
    * Triangular Lattice *
    **********************
*/

//constructor with T x D6 Symm
TriAngLattice::TriAngLattice(int numSites, bool PBC){
    is_PBC = PBC;
    Nsite = numSites;
    name = "TriAng_D6_N"+std::to_string(Nsite);
    ax = VecD {1.0, 0.0, 0.0};
    ay = VecD {0.5, std::sqrt(3.0)/2.0, 0.0};
    az = VecD {0.0, 0.0, 1.0};
    bx = VecD {2*PI, -2*PI/std::sqrt(3.0), 0.0};
    by = VecD {0.0, 4*PI/std::sqrt(3.0), 0.0};
    bz = VecD {0.0, 0.0, 2*PI};
    // unitSite.push_back(Orbital{0,0,VecD{0.0,0.0,0.0},ORBITAL::SINGLE});
    switch(numSites){
        case 9:{
            vecXAdd(3.0, a1.data(), -3.0, a2.data(), R1.data(), getDim());
            vecXAdd(3.0, a1.data(), 0.0, a2.data(), R2.data(), getDim());
            vecXAdd(1.0/3.0, b1.data(), 1.0/3.0, b2.data(), b10.data(), getDim());
            vecXAdd(0.0, b1.data(), 1.0/3.0, b2.data(), b20.data(), getDim());
            xlist = VecD {+0, +1, -1, 0, 1, -1, 0, 1, -1};
            ylist = VecD {-1, -1, +0, 0, 0, +1, 1, 1, +2};
            kxlist = VecD {0.0, 1.0, -1.0, -1.0, 0.0, 1.0, +0.0, +1.0, +2.0};
            kylist = VecD {0.0, 0.0, +0.0, +1.0, 1.0, 1.0, -1.0, -1.0, -1.0};
            break;
        }
        case 12:{
            vecXAdd(4.0, a1.data(), -2.0, a2.data(), R1.data(), getDim());
            vecXAdd(2.0, a1.data(), 2.0, a2.data(), R2.data(), getDim());
            vecXAdd(1.0/6.0, b1.data(), -1.0/6.0, b2.data(), b10.data(), getDim());
            vecXAdd(1.0/6.0, b1.data(), 2.0/6.0, b2.data(), b20.data(), getDim());
            xlist = VecD {+0, +1, +2, -1, 0, 1, 2, -1, 0, 1, -1, 0};
            ylist = VecD {-1, -1, -1, +0, 0, 0, 0, +1, 1, 1, +2, 2};
            kxlist = VecD {0.0, 1.0, -1.0, -1.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, -1.0, +0.0};
            kylist = VecD {0.0, 0.0, +0.0, +1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, -1.0, -1.0};
            break;
        }
        case 21:{
            vecXAdd(5.0, a1.data(), -4.0, a2.data(), R1.data(), getDim());
            vecXAdd(4.0, a1.data(), 1.0, a2.data(), R2.data(), getDim());
            vecXAdd(1.0/21.0, b1.data(), -4.0/21.0, b2.data(), b10.data(), getDim());
            vecXAdd(4.0/21.0, b1.data(), 5.0/21.0, b2.data(), b20.data(), getDim());
            xlist = VecD {+0, +1, +2, -1, +0, +1, +2, -2, -1, 0, 1, 2, -2, -1, 0, 1, -2, -1, 0, 1, -2};
            ylist = VecD {-2, -2, -2, -1, -1, -1, -1, +0, +0, 0, 0, 0, +1, +1, 1, 1, +2, +2, 2, 2, +3};
            kxlist = VecD {0.0, 1.0, 2.0, -2.0, -1.0, -1.0, 0.0, 1.0, 2.0, -1.0, 0.0, 1.0, 2.0, 2.0, -2.0, -1.0, +0.0, +1.0, -2.0, -1.0, +0.0};
            kylist = VecD {0.0, 0.0, 0.0, +0.0, +0.0, +1.0, 1.0, 1.0, 1.0, +2.0, 2.0, 2.0, 2.0, 3.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -2.0};
            break;
        }
        case 27:{
            vecXAdd(6.0, a1.data(), -3.0, a2.data(), R1.data(), getDim());
            vecXAdd(3.0, a1.data(), 3.0, a2.data(), R2.data(), getDim());
            vecXAdd(1.0/9.0, b1.data(), -1.0/9.0, b2.data(), b10.data(), getDim());
            vecXAdd(1.0/9.0, b1.data(), 2.0/9.0, b2.data(), b20.data(), getDim());
            xlist = VecD {+0, +1, +2, +3, -1, +0, +1, +2, +3, -2, -1, 0, 1, 2, 3, -2, -1, 0, 1, 2, -2, -1, 0, 1, -2, -1, 0};
            ylist = VecD {-2, -2, -2, -2, -1, -1, -1, -1, -1, +0, +0, 0, 0, 0, 0, +1, +1, 1, 1, 1, +2, +2, 2, 2, +3, +3, 3};
            kxlist = VecD {0.0, 1.0, 2.0, -2.0, -1.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, -1.0, 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, -2.0, -1.0, +0.0, +1.0, -2.0, -1.0, +0.0};
            kylist = VecD {0.0, 0.0, 0.0, +0.0, +0.0, +1.0, +1.0, 1.0, 1.0, 1.0, 1.0, +2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -2.0};
            break;
        }
        case 36:{
            vecXAdd(6.0, a1.data(), 0.0, a2.data(), R1.data(), getDim());
            vecXAdd(0.0, a1.data(), 6.0, a2.data(), R2.data(), getDim());
            vecXAdd(1.0/6.0, b1.data(), 0.0, b2.data(), b10.data(), getDim());
            vecXAdd(0.0, b1.data(), 1.0/6.0, b2.data(), b20.data(), getDim());
            xlist = VecD {+0, +1, +2, -2, -1, +0, +1, +2, +3, -2, -1, +0, +1, +2, +3, -3, -2, -1, +0, +1, +2, -3, -2, -1, +0, +1, +2, -4, -3, -2, -1, 0, 1, -3, -2, -1};
            ylist = VecD {-3, -3, -3, -2, -2, -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, +0, +0 ,+0 ,+0, +0, +0, +1, +1, +1, +1, +1, +1, +2, +2, +2, +2, 2, 2, +3, +3, +3};
            kxlist = VecD {0.0, 1.0, 2.0, 3.0, -2.0, -1.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 2.0, -3.0, -2.0, -1.0, +0.0, +1.0, +2.0, -3.0, -2.0, -1.0, +0.0, +1.0, -2.0, -1.0};
            kylist = VecD {0.0, 0.0, 0.0, 0.0, +0.0, +0.0, +1.0, +1.0, 1.0, 1.0, 1.0, 1.0, +2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 4.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -2.0, -2.0, -2.0, -3.0, -3.0};
            break;
        }
        default: 
            std::cout<<"N="<<Nsite<<" must be one of the values: 9, 12, 21, 27, 36!"<<std::endl;
            exit(1);
            break;
    }
    zlist.resize(Nsite,0.0);
    kzlist.resize(Nsite,0.0);
    VecD vtmp(3);
    vecXAdd(0.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
    vecXAdd(1.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
    vecXAdd(0.0, R1.data(), 1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
    vecXAdd(-1.0, R1.data(), 1.0 , R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
    vecXAdd(-1.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
    vecXAdd(0.0, R1.data(), -1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
    vecXAdd(1.0, R1.data(), -1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
    
    // construct();
};

//constructor with T without D6 (possible C6?)
TriAngLattice::TriAngLattice(int N1, int N2, bool PBC){
    // default PBC condition
    is_PBC = PBC;
    Nsite = N1 * N2;
    name = "TriAng"+std::to_string(N1)+"x"+std::to_string(N2);

    ax = VecD {1.0, 0.0, 0.0};
    ay = VecD {0.5, std::sqrt(3.0)/2.0, 0.0};
    az = VecD {0.0, 0.0, 1.0};
    bx = VecD {2*PI, -2*PI/std::sqrt(3.0), 0.0};
    by = VecD {0.0, 4*PI/std::sqrt(3.0), 0.0};
    bz = VecD {0.0, 0.0, 2*PI};
    
    // unitSite.push_back(Orbital{0,0,VecD{0.0,0.0,0.0},ORBITAL::SINGLE});

    vecXAdd((double)N1, a1.data(), 0.0, a2.data(), R1.data(), getDim());
    vecXAdd(0.0, a1.data(), (double)N2, a2.data(), R2.data(), getDim());
    vecXAdd(1.0/(double)N1, b1.data(), 0.0, b2.data(), b10.data(), getDim());
    vecXAdd(0.0, b1.data(), 1.0/(double)N2, b2.data(), b20.data(), getDim());
    for (int y = 0; y < N2; y++){
        for (int x = 0; x < N1; x++){
            xlist.push_back(x);
            ylist.push_back(y);
            kxlist.push_back(x);
            kylist.push_back(y);
        }
    }
    zlist.resize(Nsite,0.0);
    kzlist.resize(Nsite,0.0);
    VecD vtmp(3);
    vecXAdd(0.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
    if (N1==1 and N2>1){
        vecXAdd(0.0, R1.data(), 1.0 , R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(0.0, R1.data(), -1.0 , R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
    }
    else if (N1>1 and N2==1){
        vecXAdd(1.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(-1.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
    }
    else if (N1>1 and N2>1){
        vecXAdd(1.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(1.0, R1.data(), 1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(0.0, R1.data(), 1.0 , R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(-1.0, R1.data(), 1.0 , R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(-1.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(-1.0, R1.data(), -1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(0.0, R1.data(), -1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(1.0, R1.data(), -1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
    }
    
    // construct();
};


/*
    ******************
    * Square Lattice *
    ******************
*/
SquareLattice::SquareLattice(int N1, int N2, bool PBC){
    is_PBC = PBC;
    Nsite = N1 * N2;
    name = "Square"+std::to_string(N1)+"x"+std::to_string(N2);
    ax = VecD {1.0, 0.0, 0.0};
    ay = VecD {0.0, 1.0, 0.0};
    az = VecD {0.0, 0.0, 1.0};
    bx = VecD {2*PI, 0.0, 0.0};
    by = VecD {0.0, 2*PI, 0.0};
    bz = VecD {0.0, 0.0, 2*PI};

    vecXAdd((double)N1, a1.data(), 0.0, a2.data(), R1.data(), getDim());
    vecXAdd(0.0, a1.data(), (double)N2, a2.data(), R2.data(), getDim());
    vecXAdd(1.0/(double)N1, b1.data(), 0.0, b2.data(), b10.data(), getDim());
    vecXAdd(0.0, b1.data(), 1.0/(double)N2, b2.data(), b20.data(), getDim());
    for (int y = 0; y < N2; y++){
        for (int x = 0; x < N1; x++){
            xlist.push_back(x);
            ylist.push_back(y);
            kxlist.push_back(x);
            kylist.push_back(y);
        }
    }
    zlist.resize(Nsite,0.0);
    kzlist.resize(Nsite,0.0);
    VecD vtmp(3);
    vecXAdd(0.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
    if (N1==1 and N2>1){
        vecXAdd(0.0, R1.data(), 1.0 , R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(0.0, R1.data(), -1.0 , R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
    }
    else if (N1>1 and N2==1){
        vecXAdd(1.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(-1.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
    }
    else if (N1>1 and N2>1){
        vecXAdd(1.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(1.0, R1.data(), 1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(0.0, R1.data(), 1.0 , R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(-1.0, R1.data(), 1.0 , R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(-1.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(-1.0, R1.data(), -1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(0.0, R1.data(), -1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(1.0, R1.data(), -1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
    }    
}
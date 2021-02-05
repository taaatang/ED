#ifndef __LINKS_H__
#define __LINKS_H__

#include <vector>
#include "Operator/OperatorBasics.hpp"
#include "Global/globalPara.hpp"

/**********************
 *   SQUARE LATTICE   *
 * SINGLEBAND HUBBARD *
 **********************/
Link<dataType> tnn_px(LINK_TYPE::HOPPING_T, {ORBITAL::SINGLE, ORBITAL::SINGLE}, -1.0); tnn_px.addLinkVec({1.0, 0.0, 0.0});
Link<dataType> tnn_py(LINK_TYPE::HOPPING_T, {ORBITAL::SINGLE, ORBITAL::SINGLE}, -1.0); tnn_py.addLinkVec({0.0, 1.0, 0.0});
Link<dataType> tnnn_pxpy(LINK_TYPE::HOPPING_T, {ORBITAL::SINGLE, ORBITAL::SINGLE}, 1.0); tnnn_pxpy.addLinkVec({1.0, 1.0, 0.0});
Link<dataType> tnnn_pxny(LINK_TYPE::HOPPING_T, {ORBITAL::SINGLE, ORBITAL::SINGLE}, 1.0); tnnn_pxny.addLinkVec({1.0, -1.0, 0.0});

/*********************
 *  SQUARE LATTICE   *
 * MULTIBAND HUBBARD *
 *      HOPPING      *
 *********************/

Link<dataType> tdpx(LINK_TYPE::HOPPING_T, {ORBITAL::Dx2y2, ORBITAL::Px}, -1.0); tdpx.addLinkVec({0.5,0.0,0.0});
Link<dataType> tdpy(LINK_TYPE::HOPPING_T, {ORBITAL::Dx2y2, ORBITAL::Py}, 1.0); tdpy.addLinkVec({0.0,0.5,0.0});
Link<dataType> tpxd(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Dx2y2}, 1.0); tpxd.addLinkVec({0.5,0.0,0.0});
Link<dataType> tpxpy(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Py}, -1.0); tpxpy.addLinkVec({0.5,0.5,0.0});
Link<dataType> tpypx(LINK_TYPE::HOPPING_T, {ORBITAL::Py, ORBITAL::Px}, -1.0); tpypx.addLinkVec({0.5,0.5,0.5});
Link<dataType> tpxpyp(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Py}, 1.0); tpxpyp.addLinkVec({0.5,-0.5,0.0});
Link<dataType> tpypxp(LINK_TYPE::HOPPING_T, {ORBITAL::Py, ORBITAL::Px}, 1.0); tpypxp.addLinkVec({0.5,-0.5,0.0});
Link<dataType> tpyd(LINK_TYPE::HOPPING_T, {ORBITAL::Py, ORBITAL::Dx2y2}, -1.0); tpyd.addLinkVec({0.0,0.5,0.0});

Link<dataType> tdzpx(LINK_TYPE::HOPPING_T, {ORBITAL::Dz2, ORBITAL::Px}, 1.0); tdzpx.addLinkVec({0.5,0.0,0.0});
Link<dataType> tdzpy(LINK_TYPE::HOPPING_T, {ORBITAL::Dz2, ORBITAL::Py}, 1.0); tdzpy.addLinkVec({0.0,0.5,0.0});
Link<dataType> tdzpzu(LINK_TYPE::HOPPING_T, {ORBITAL::Dz2, ORBITAL::Pzu}, -1.0); tdzpzu.addLinkVec({0.0,0.0,0.5});
Link<dataType> tpxdz(LINK_TYPE::HOPPING_T, {ORBITAL::Dz2, ORBITAL::Px}, -1.0); tpxdz.addLinkVec({-0.5,0.0,0.0});
Link<dataType> tpydz(LINK_TYPE::HOPPING_T, {ORBITAL::Dz2, ORBITAL::Py}, -1.0); tpydz.addLinkVec({0.0,-0.5,0.0});
Link<dataType> tdzpzd(LINK_TYPE::HOPPING_T, {ORBITAL::Dz2, ORBITAL::Pzd}, 1.0); tdzpzd.addLinkVec({0.0,0.0,-0.5});

Link<dataType> tpxpz(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Pzu}, 1.0); tpxpz.addLinkVec({-0.5,0.0,0.5});
Link<dataType> tpzpx(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Px}, -1.0); tpzpx.addLinkVec({0.5,0.0,0.5});
Link<dataType> tpypz(LINK_TYPE::HOPPING_T, {ORBITAL::Py, ORBITAL::Pzu}, 1.0); tpypz.addLinkVec({0.0,-0.5,0.5});
Link<dataType> tpzpy(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Py}, -1.0); tpzpy.addLinkVec({0.0,0.5,0.5});

Link<dataType> tpxpzp(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Pzu}, -1.0); tpxpzp.addLinkVec({0.5,0.0,0.5});
Link<dataType> tpzpxp(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Px}, 1.0); tpzpxp.addLinkVec({-0.5,0.0,0.5});
Link<dataType> tpypzp(LINK_TYPE::HOPPING_T, {ORBITAL::Py, ORBITAL::Pzu}, -1.0); tpypzp.addLinkVec({0.0,0.5,0.5});
Link<dataType> tpzpyp(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Py}, 1.0); tpzpyp.addLinkVec({0.0,-0.5,0.5});

// interband Coulomb interaction
Link<dataType> ndnpx(LINK_TYPE::HUBBARD_U, {ORBITAL::Dx2y2, ORBITAL::Px}, 1.0); ndnpx.addLinkVec({0.5,0.0,0.0}).addLinkVec({-0.5,0.0,0.0});
Link<dataType> ndnpy(LINK_TYPE::HUBBARD_U, {ORBITAL::Dx2y2, ORBITAL::Py}, 1.0); ndnpy.addLinkVec({0.0,0.5,0.0}).addLinkVec({0.0,-0.5,0.0});

Link<dataType> ndnpzu(LINK_TYPE::HUBBARD_U, {ORBITAL::Dx2y2, ORBITAL::Pzu}, 1.0); ndnpzu.addLinkVec({0.0,0.0,0.5});
Link<dataType> ndnpzd(LINK_TYPE::HUBBARD_U, {ORBITAL::Dx2y2, ORBITAL::Pzd}, 1.0); ndnpzd.addLinkVec({0.0,0.0,-0.5});

Link<dataType> ndndz(LINK_TYPE::HUBBARD_U, {ORBITAL::Dx2y2, ORBITAL::Dz2}, 1.0); ndndz.addLinkVec({0.0,0.0,0.0});
// exchange
Link<dataType> ex_ddz(LINK_TYPE::EXCHANGE_J, {ORBITAL::Dx2y2, ORBITAL::Dz2}, 1.0); ex_ddz.addLinkVec({0.0,0.0,0.0});
// pair hoping
Link<dataType> ph_ddz(LINK_TYPE::PAIR_HOPPING_J, {ORBITAL::Dx2y2, ORBITAL::Dz2}, 1.0); ph_ddz.addLinkVec({0.0,0.0,0.0});

/**********************
 * TRIANGULAR LATTICE *
 *      HOPPING       *
 **********************/

Link<dataType> t1Link(LINK_TYPE::HOPPING_T, {ORBITAL::SINGLE, ORBITAL::SINGLE}, t1);
Link<dataType> t2Link(LINK_TYPE::HOPPING_T, {ORBITAL::SINGLE, ORBITAL::SINGLE}, t2);
t1Link.addLinkVec(VecD{1.0,0.0,0.0}).addLinkVec(VecD{0.0,1.0,0.0}).addLinkVec(VecD{1.0,-1.0,0.0});
t2Link.addLinkVec(VecD{1.0,1.0,0.0}).addLinkVec(VecD{-1.0,2.0,0.0}).addLinkVec(VecD{2.0,-1.0,0.0});

/**********************
 * TRIANGULAR LATTICE *
 *   SUPEREXCHANGE    *
 **********************/

Link<dataType> J1Link(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE}, J1);
Link<dataType> J2Link(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE}, J2, true);
Link<dataType> JkLink(LINK_TYPE::CHIRAL_K, {ORBITAL::SINGLE, ORBITAL::SINGLE, ORBITAL::SINGLE},Jk, true);

J1Link.addLinkVec(VecD{1.0,0.0,0.0}).addLinkVec(VecD{0.0,1.0,0.0}).addLinkVec(VecD{1.0,-1.0,0.0});
J2Link.addLinkVec(VecD{1.0,1.0,0.0}).addLinkVec(VecD{-1.0,2.0,0.0}).addLinkVec(VecD{2.0,-1.0,0.0});
JkLink.addLinkVec({0.0,1.0,0.0}).addLinkVec({1.0,0.0,0.0}).addLinkVec({1.0,0.0,0.0}).addLinkVec({1.0,-1.0,0.0});

std::vector<Link<dataType>> HeisenbergLink( ) {
    return std::vector<Link<dataType>{J1Link, J2Link, JkLink};
}

std::vector<Link<dataType>> HubbardSingleBandLink( ) {
    return std::vector<Link<dataType>>{tnn_px, tnn_py, tnnn_pxpy, tnnn_pxny};
}

std::vector<Link<dataType>> HubbardMultiBandLink( ) {
    std::vector<Link<dataType>> links_3band{tdpx, tdpy, tpxd, tpyd, tpxpy, tpypx, tpxpyp, tpypxp, ndnpx, ndnpy};
    std::vector<Link<dataType>> links_pz{tpxpz, tpzpx, tpypz, tpzpy, tpxpzp, tpzpxp, tpypzp, tpzpyp, ndnpzu, ndnpzd};
    std::vector<Link<dataType>> links_dz2{tdzpx, tdzpy, tpxdz, tpydz, tdzpzu, tdzpzd, ndndz, ex_ddz, ph_ddz};
    std::vector<Link<dataType>> links;
    links.insert(links.end(), links_3band.begin(), links_3band.end());
    links.insert(links.end(), links_pz.begin(), links_pz.end());
    links.insert(links.end(), links_dz2.begin(), links_dz2.end());
    return links;
}

std::vector<Link<dataType>> HubbardLink( ) {
    auto single = HubbardSingleBandLink();
    auto multi = HubbardMultiBandLink();
    std::vector<Link<dataType>> links;
    links.insert(links.end(), single.begin(), single.end());
    links.insert(links.end(), multi.begin(), multi.end());
    return links;
}
#endif // __LINKS_H__
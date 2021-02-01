#ifndef __LINKS_H__
#define __LINKS_H__

#include "../Operator/OperatorBasics.hpp"
#include "../Global/globalPara.hpp"

/*********************
 * MULTIBAND HUBBARD *
 *********************/

Link<dataType> tdpx(LINK_TYPE::HOPPING_T, {ORBITAL::Dx2y2, ORBITAL::Px}, -tdp_val); tdpx.addLinkVec({0.5,0.0,0.0});
Link<dataType> tdpy(LINK_TYPE::HOPPING_T, {ORBITAL::Dx2y2, ORBITAL::Py}, tdp_val); tdpy.addLinkVec({0.0,0.5,0.0});
Link<dataType> tpxd(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Dx2y2}, tdp_val); tpxd.addLinkVec({0.5,0.0,0.0});
Link<dataType> tpxpy(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Py}, -tpp_val); tpxpy.addLinkVec({0.5,0.5,0.0}).addLinkVec({-0.5,-0.5,0.0});
Link<dataType> tpxpyp(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Py}, tpp_val); tpxpyp.addLinkVec({0.5,-0.5,0.0}).addLinkVec({-0.5,0.5,0.0});
Link<dataType> tpyd(LINK_TYPE::HOPPING_T, {ORBITAL::Py, ORBITAL::Dx2y2}, -tdp_val); tpyd.addLinkVec({0.0,0.5,0.0});

Link<dataType> tdzpx(LINK_TYPE::HOPPING_T, {ORBITAL::Dz2, ORBITAL::Px}, tdzp_val); tdzpx.addLinkVec({0.4,0.0,0.0});
Link<dataType> tdzpy(LINK_TYPE::HOPPING_T, {ORBITAL::Dz2, ORBITAL::Py}, tdzp_val); tdzpy.addLinkVec({-0.1,0.5,0.0});
Link<dataType> tdzpzu(LINK_TYPE::HOPPING_T, {ORBITAL::Dz2, ORBITAL::Pzu}, -tdzpz_val); tdzpzu.addLinkVec({-0.1,0.0,0.5});
Link<dataType> tpxdz(LINK_TYPE::HOPPING_T, {ORBITAL::Dz2, ORBITAL::Px}, -tdzp_val); tpxdz.addLinkVec({-0.6,0.0,0.0});
Link<dataType> tpydz(LINK_TYPE::HOPPING_T, {ORBITAL::Dz2, ORBITAL::Py}, -tdzp_val); tpydz.addLinkVec({-0.1,-0.5,0.0});
Link<dataType> tdzpzd(LINK_TYPE::HOPPING_T, {ORBITAL::Dz2, ORBITAL::Pzd}, tdzpz_val); tdzpzd.addLinkVec({-0.1,0.0,-0.5});

// density interaction
Link<dataType> ndndz(LINK_TYPE::HUBBARD_U, {ORBITAL::Dx2y2, ORBITAL::Dz2}, Udz); ndndz.addLinkVec({0.1,0.0,0.0});
// exchange
Link<dataType> ex_ddz(LINK_TYPE::EXCHANGE_J, {ORBITAL::Dx2y2, ORBITAL::Dz2}, Jex); ex_ddz.addLinkVec({0.1,0.0,0.0});
// pair hoping
Link<dataType> ph_ddz(LINK_TYPE::PAIR_HOPPING_J, {ORBITAL::Dx2y2, ORBITAL::Dz2}, Jex); ph_ddz.addLinkVec({0.1,0.0,0.0});

bool NotConst=true, isOrdered=false;
Link<dataType> tpxpz(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Pzu}, tppz_val, NotConst, isOrdered); tpxpz.addLinkVec({-0.5,0.0,0.5});
Link<dataType> tpzpx(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Px}, -tppz_val, NotConst, isOrdered); tpzpx.addLinkVec({0.5,0.0,0.5});
Link<dataType> tpypz(LINK_TYPE::HOPPING_T, {ORBITAL::Py, ORBITAL::Pzu}, tppz_val, NotConst, isOrdered); tpypz.addLinkVec({0.0,-0.5,0.5});
Link<dataType> tpzpy(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Py}, -tppz_val, NotConst, isOrdered); tpzpy.addLinkVec({0.0,0.5,0.5});

Link<dataType> tpxpzp(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Pzu}, -tppz_val, NotConst, isOrdered); tpxpzp.addLinkVec({0.5,0.0,0.5});
Link<dataType> tpzpxp(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Px}, tppz_val, NotConst, isOrdered); tpzpxp.addLinkVec({-0.5,0.0,0.5});
Link<dataType> tpypzp(LINK_TYPE::HOPPING_T, {ORBITAL::Py, ORBITAL::Pzu}, -tppz_val, NotConst, isOrdered); tpypzp.addLinkVec({0.0,0.5,0.5});
Link<dataType> tpzpyp(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Py}, tppz_val, NotConst, isOrdered); tpzpyp.addLinkVec({0.0,-0.5,0.5});

Link<dataType> ndnpx(LINK_TYPE::HUBBARD_U, {ORBITAL::Dx2y2, ORBITAL::Px}, Upd); ndnpx.addLinkVec({0.5,0.0,0.0}).addLinkVec({-0.5,0.0,0.0});
Link<dataType> ndnpy(LINK_TYPE::HUBBARD_U, {ORBITAL::Dx2y2, ORBITAL::Py}, Upd); ndnpy.addLinkVec({0.0,0.5,0.0}).addLinkVec({0.0,-0.5,0.0});
Link<dataType> ndnpzu(LINK_TYPE::HUBBARD_U, {ORBITAL::Dx2y2, ORBITAL::Pzu}, Uda); ndnpzu.addLinkVec({0.0,0.0,0.5});
Link<dataType> ndnpzd(LINK_TYPE::HUBBARD_U, {ORBITAL::Dx2y2, ORBITAL::Pzd}, Uda); ndnpzd.addLinkVec({0.0,0.0,-0.5});


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

#endif // __LINKS_H__
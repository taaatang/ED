#include "Operator/links.hpp"

/**********************
 *   SQUARE LATTICE   *
 * SINGLEBAND HUBBARD *
 **********************/
Link<dataType> tnn_px(LINK_TYPE::HOPPING_T, {ORBITAL::SINGLE, ORBITAL::SINGLE}, -1.0, {{1.0, 0.0, 0.0}});
Link<dataType> tnn_py(LINK_TYPE::HOPPING_T, {ORBITAL::SINGLE, ORBITAL::SINGLE}, -1.0, {{0.0, 1.0, 0.0}});
Link<dataType> tnnn_pxpy(LINK_TYPE::HOPPING_T, {ORBITAL::SINGLE, ORBITAL::SINGLE}, 1.0, {{1.0, 1.0, 0.0}});
Link<dataType> tnnn_pxny(LINK_TYPE::HOPPING_T, {ORBITAL::SINGLE, ORBITAL::SINGLE}, 1.0, {{1.0, -1.0, 0.0}});

/*********************
 *  SQUARE LATTICE   *
 * MULTIBAND HUBBARD *
 *      HOPPING      *
 *********************/

Link<dataType> tdpx(LINK_TYPE::HOPPING_T, {ORBITAL::Dx2y2, ORBITAL::Px}, -1.0, {{0.5,0.0,0.0}});
Link<dataType> tdpy(LINK_TYPE::HOPPING_T, {ORBITAL::Dx2y2, ORBITAL::Py}, 1.0, {{0.0,0.5,0.0}});
Link<dataType> tpxd(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Dx2y2}, 1.0, {{0.5,0.0,0.0}});
Link<dataType> tpxpy(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Py}, -1.0, {{0.5,0.5,0.0}});
Link<dataType> tpypx(LINK_TYPE::HOPPING_T, {ORBITAL::Py, ORBITAL::Px}, -1.0, {{0.5,0.5,0.0}});
Link<dataType> tpxpyp(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Py}, 1.0, {{0.5,-0.5,0.0}});
Link<dataType> tpypxp(LINK_TYPE::HOPPING_T, {ORBITAL::Py, ORBITAL::Px}, 1.0, {{0.5,-0.5,0.0}});
Link<dataType> tpyd(LINK_TYPE::HOPPING_T, {ORBITAL::Py, ORBITAL::Dx2y2}, -1.0, {{0.0,0.5,0.0}});

Link<dataType> tdzpx(LINK_TYPE::HOPPING_T, {ORBITAL::Dz2, ORBITAL::Px}, 1.0, {{0.5,0.0,0.0}});
Link<dataType> tdzpy(LINK_TYPE::HOPPING_T, {ORBITAL::Dz2, ORBITAL::Py}, 1.0, {{0.0,0.5,0.0}});
Link<dataType> tdzpzu(LINK_TYPE::HOPPING_T, {ORBITAL::Dz2, ORBITAL::Pzu}, -1.0, {{0.0,0.0,0.5}});
Link<dataType> tpxdz(LINK_TYPE::HOPPING_T, {ORBITAL::Dz2, ORBITAL::Px}, -1.0, {{-0.5,0.0,0.0}});
Link<dataType> tpydz(LINK_TYPE::HOPPING_T, {ORBITAL::Dz2, ORBITAL::Py}, -1.0, {{0.0,-0.5,0.0}});
Link<dataType> tdzpzd(LINK_TYPE::HOPPING_T, {ORBITAL::Dz2, ORBITAL::Pzd}, 1.0, {{0.0,0.0,-0.5}});

Link<dataType> tpxpz(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Pzu}, 1.0, {{-0.5,0.0,0.5}});
Link<dataType> tpzpx(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Px}, -1.0, {{0.5,0.0,0.5}});
Link<dataType> tpypz(LINK_TYPE::HOPPING_T, {ORBITAL::Py, ORBITAL::Pzu}, 1.0, {{0.0,-0.5,0.5}});
Link<dataType> tpzpy(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Py}, -1.0, {{0.0,0.5,0.5}});

Link<dataType> tpxpzp(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Pzu}, -1.0, {{0.5,0.0,0.5}});
Link<dataType> tpzpxp(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Px}, 1.0, {{-0.5,0.0,0.5}});
Link<dataType> tpypzp(LINK_TYPE::HOPPING_T, {ORBITAL::Py, ORBITAL::Pzu}, -1.0, {{0.0,0.5,0.5}});
Link<dataType> tpzpyp(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Py}, 1.0, {{0.0,-0.5,0.5}});

// interband Coulomb interaction
Link<dataType> ndnpx(LINK_TYPE::HUBBARD_U, {ORBITAL::Dx2y2, ORBITAL::Px}, 1.0, { {0.5,0.0,0.0}, {-0.5,0.0,0.0} });
Link<dataType> ndnpy(LINK_TYPE::HUBBARD_U, {ORBITAL::Dx2y2, ORBITAL::Py}, 1.0, { {0.0,0.5,0.0}, {0.0,-0.5,0.0} });

Link<dataType> ndnpzu(LINK_TYPE::HUBBARD_U, {ORBITAL::Dx2y2, ORBITAL::Pzu}, 1.0, {{0.0,0.0,0.5}});
Link<dataType> ndnpzd(LINK_TYPE::HUBBARD_U, {ORBITAL::Dx2y2, ORBITAL::Pzd}, 1.0, {{0.0,0.0,-0.5}});

Link<dataType> ndndz(LINK_TYPE::HUBBARD_U, {ORBITAL::Dx2y2, ORBITAL::Dz2}, 1.0, {{0.0,0.0,0.0}});
// exchange
Link<dataType> ex_ddz(LINK_TYPE::EXCHANGE_J, {ORBITAL::Dx2y2, ORBITAL::Dz2}, 1.0, {{0.0,0.0,0.0}});
// pair hoping
Link<dataType> ph_ddz(LINK_TYPE::PAIR_HOPPING_J, {ORBITAL::Dx2y2, ORBITAL::Dz2}, 1.0, {{0.0,0.0,0.0}});

/**********************
 * TRIANGULAR LATTICE *
 *      HOPPING       *
 **********************/

Link<dataType> t1Link(LINK_TYPE::HOPPING_T, {ORBITAL::SINGLE, ORBITAL::SINGLE}, 1.0, {{1.0,0.0,0.0}, {0.0,1.0,0.0}, {1.0,-1.0,0.0}});
Link<dataType> t2Link(LINK_TYPE::HOPPING_T, {ORBITAL::SINGLE, ORBITAL::SINGLE}, -1.0, {{1.0,1.0,0.0}, {-1.0,2.0,0.0}, {2.0,-1.0,0.0}});

/**********************
 * TRIANGULAR LATTICE *
 *   SUPEREXCHANGE    *
 **********************/

Link<dataType> J1Link(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE}, 1.0, {{1.0,0.0,0.0}, {0.0,1.0,0.0}, {1.0,-1.0,0.0}});
Link<dataType> J2Link(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE}, 1.0, {{1.0,1.0,0.0}, {-1.0,2.0,0.0}, {2.0,-1.0,0.0}}, true);
Link<dataType> JkLink(LINK_TYPE::CHIRAL_K, {ORBITAL::SINGLE, ORBITAL::SINGLE, ORBITAL::SINGLE}, 1.0, {{0.0,1.0,0.0}, {1.0,0.0,0.0}, {1.0,0.0,0.0}, {1.0,-1.0,0.0}},  true);

std::vector<Link<dataType>> HeisenbergLink( ) {
    return std::vector<Link<dataType>>{J1Link, J2Link, JkLink};
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
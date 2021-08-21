#pragma once

#include <memory>

#include "operator/operators.hpp"
#include "utils/path.hpp"
#include "utils/paras.hpp"

void setpath(Parameters&);

Orbital stringToOrb(std::string name, int id);

LATTICE getLattType(const Parameters& para);

MODEL getModel(const Parameters& para);

void setlatt(const Parameters&, std::unique_ptr<Geometry>& latt);

template <IsBasisType B>
void setbasis(const Parameters&, std::unique_ptr<Basis<B>>&, Geometry*);

template <IsBasisType B>
void setbasis(std::unique_ptr<Basis<B>>&, MODEL, Geometry*, int nuf, int ndf, int kf, int pf);

template <typename T, IsBasisType B>
void setham(const Parameters&, std::unique_ptr<Hamiltonian<T, B>>&, Geometry*, Basis<B>*);

template <typename T, IsBasisType B>
void setBasics(const Parameters&, std::unique_ptr<Geometry>& latt, std::unique_ptr<Basis<B>>& Basis, std::unique_ptr<Hamiltonian<T, B>>& H);

void setmeasure(const Parameters&);

void setpulse(const Parameters&, Pulse&);

template <IsBasisType B>
void setbasis(const Parameters& para, std::unique_ptr<Basis<B>>& ba, Geometry* latt){
    auto kid = para.template get<int>("kidx").value();
    auto pid = para.template get<int>("pidx").value();
    auto nu = para.template get<int>("nu").value();
    auto nd = para.template get<int>("nd").value();
	auto nPho = para.template get<int>("maxPhPerSite");
    ba = std::unique_ptr<Basis<B>>(new Basis<B>(latt, nu, nd, nPho, kid, pid));
}


template <IsBasisType B>
void setbasis(std::unique_ptr<Basis<B>>& ba, MODEL model, Geometry* latt, int nuf, int ndf, int nPhof, int kf, int pf) {
    ba = std::unique_ptr<Basis<B>>(new Basis<B>(latt, nuf, ndf, nPhof, kf, pf));
}

template <typename T, IsBasisType B>
void setham(const Parameters& para, std::unique_ptr<Hamiltonian<T, B>>& H, Geometry* latt, Basis<B>* b) {
	H = std::unique_ptr<Hamiltonian<T, B>>(new Hamiltonian<T, B>(latt, b, b, true, true, 1, 1));
    if constexpr (ContainCharge<B>) {
        auto unitcell = latt->getUnitCell();
        if (unitcell.size() == 1) {/* single band */
            H->pushU({ORBITAL::SINGLE}, para.get<double>("Uss").value());
            auto links = HubbardSingleBandLink(*latt);
            auto tnn = para.get<double>("tnn").value();
            auto tnnn = para.get<double>("tnnn").value();
            if (std::abs(tnn) > INFINITESIMAL) {/* nearest neighbor hopping */
                for (int id = 0; id < 2; ++id) {
                    links.at(id).setVal(links.at(id).getVal() * tnn);
                    H->pushLink(links.at(id), 0);
                }
            }
            if (std::abs(tnnn) > INFINITESIMAL) {/* next nearest neighbor hopping */
                for (int id = 2; id < 4; ++id) {
                    links.at(id).setVal(links.at(id).getVal() * tnnn);
                    H->pushLink(links.at(id), 0);
                }
            }
        } else {/* multi band */
            for (auto orb : unitcell) {
                auto ids = latt->getOrbID(orb.orb);
                auto id = ids.at(0);
                auto v = (*para.get<ArrD>("t")).at(id).at(id);
                auto u = (*para.get<ArrD>("U")).at(id).at(id);
                H->pushV({orb.orb}, v);
                H->pushU({orb.orb}, u);
            }
            auto links = HubbardMultiBandLink(*latt);
            for (auto& link : links) {
                auto orbs = link.getOrbs();
                bool validLink = true;
                VecI orbids;
                for (auto orb : orbs) {
                    auto ids = latt->getOrbID(orb);
                    if (ids.empty()) {
                        validLink = false;
                        break;
                    }
                    orbids.push_back(ids.at(0));
                }
                if (!validLink) {
                    continue;
                }
                assert(orbids.size()==2 and orbids[0]!=orbids[1]);
                std::optional<ArrD> arrdOpt = std::nullopt;
                switch (link.getLinkType()) {
                    case LINK_TYPE::HOPPING_T:
                        arrdOpt = para.get<ArrD>("t");
                        break;
                    case LINK_TYPE::HUBBARD_U:
                        arrdOpt = para.get<ArrD>("U");
                        break;
                    case LINK_TYPE::EXCHANGE_J:
                        arrdOpt = para.get<ArrD>("J");
                        break;
                    case LINK_TYPE::PAIR_HOPPING_J:
                        arrdOpt = para.get<ArrD>("J");
                        break;
                    default:
                        break;
                }
                if (arrdOpt) {
                    auto val = arrdOpt->at(orbids[0]).at(orbids[1]);
                    if (std::abs(val) > INFINITESIMAL) {
                        link.setVal(link.getVal() * val);
                        H->pushLink(link, 0);
                    }
                }
            }
        }
    }
	
	if constexpr (ContainSpin<B>) {
        auto J1 = para.get<double>("J1");
        if (std::abs(J1.value_or(0.0)) > INFINITESIMAL) {
            auto link = HeisenbergLink("J1", *latt);
            link.setVal(link.getVal() * J1.value_or(0.0));
            H->pushLink(link, 0);
        }
        auto J2 = para.get<double>("J2");
        if (std::abs(J2.value_or(0.0)) > INFINITESIMAL) {
            auto link = HeisenbergLink("J2", *latt);
            link.setVal(link.getVal() * J2.value_or(0.0));
            H->pushLink(link, 0);
        }
        auto Jk = para.get<double>("Jk");
        if (std::abs(Jk.value_or(0.0)) > INFINITESIMAL) {
            auto link = HeisenbergLink("Jk", *latt);
            link.setVal(link.getVal() * Jk.value_or(0.0));
            H->pushLink(link, 0);
        }

    }
    H->transform();
}

template <typename T, IsBasisType B>
void setBasics(const Parameters& para, std::unique_ptr<Geometry>& latt, std::unique_ptr<Basis<B>>& b, std::unique_ptr<Hamiltonian<T, B>>& H) {
    setlatt(para, latt);
    setbasis(para, b, latt.get());
    setham(para, H, latt.get(), b.get());
}
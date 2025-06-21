#ifndef NEEBC_H
#define NEEBC_H
#include "../../src/palabos2D.h"
#include "../../src/palabos2D.hh"

namespace plb
{
    template <typename T, template <typename U> class Descriptor, int direction, int orientation>
    class NeeFlatAdiabaticBoundaryFunctional2D : public BoxProcessingFunctional2D_L<T, Descriptor>
    {
    public:
        virtual plint extent() const { return 2; }
        virtual plint extent(int whichDirection) const { return 2; }
        virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
        {
            plint velOffset = Descriptor<T>::ExternalField::velocityBeginsAt;
            for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            {
                for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                {
                    plint iX_prev = iX + ((direction == 0) ? (-orientation) : 0);
                    plint iY_prev = iY + ((direction == 1) ? (-orientation) : 0);

                    auto cell1 = lattice.get(iX_prev, iY_prev);
                    auto *cell0 = &(lattice.get(iX, iY));
                    T temperature_1 = lattice.get(iX_prev, iY_prev).computeDensity();

                    T C1 = cell1.computeDensity();
                    T C0 = C1;
                    T u0_[Descriptor<T>::d];
                    T u1_[Descriptor<T>::d];
                    cell0->getExternalField(velOffset, Descriptor<T>::d, u0_);
                    cell1.getExternalField(velOffset, Descriptor<T>::d, u1_);

                    Array<T, Descriptor<T>::d> u0;
                    Array<T, Descriptor<T>::d> u1;
                    for (plint i = 0; i < Descriptor<T>::d; ++i)
                    {
                        u0[i] = u0_[i];
                        u1[i] = u1_[i];
                    }

                    Array<T, Descriptor<T>::d> jEq0(C0 * u0);
                    Array<T, Descriptor<T>::d> jEq1(C1 * u1);

                    auto jSqr = 0;

                    for (pluint iPop = 0; iPop < Descriptor<T>::q; ++iPop)
                    {
                        (*cell0)[iPop] =
                            // cell0->[iPop] =
                            (*cell0).computeEquilibrium(iPop, C0, jEq0, jSqr) + cell1[iPop] - cell1.computeEquilibrium(iPop, C1, jEq1, jSqr);
                    }

                    lattice.get(iX, iY).defineDensity(temperature_1);
                }
            }
        }
        virtual NeeFlatAdiabaticBoundaryFunctional2D<T, Descriptor, direction, orientation> *clone() const
        {
            return new NeeFlatAdiabaticBoundaryFunctional2D<T, Descriptor, direction, orientation>(*this);
        }

        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
        {
            modified[0] = modif::dynamicVariables;
        }
        virtual BlockDomain::DomainT appliesTo() const
        {
            return BlockDomain::bulk;
        }
    };

    template <typename T, template <typename U> class Descriptor, int direction, int orientation>
    class NeeFlatDirichletBoundaryFunctional2D : public BoxProcessingFunctional2D_L<T, Descriptor>
    {
    public:
        virtual plint extent() const { return 2; }
        virtual plint extent(int whichDirection) const { return 2; }

        NeeFlatDirichletBoundaryFunctional2D<T, Descriptor, direction, orientation>(T bcValue_) : bcValue(bcValue_) {}

        virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
        {
            plint velOffset = Descriptor<T>::ExternalField::velocityBeginsAt;

            for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            {
                plint iX_prev = iX + ((direction == 0) ? (-orientation) : 0);
                for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                {
                    plint iY_prev = iY + ((direction == 1) ? (-orientation) : 0);

                    auto cell1 = lattice.get(iX_prev, iY_prev);
                    auto *cell0 = &(lattice.get(iX, iY));

                    T C1 = cell1.computeDensity();
                    T C0 = bcValue;

                    T u0_[Descriptor<T>::d];
                    T u1_[Descriptor<T>::d];
                    cell0->getExternalField(velOffset, Descriptor<T>::d, u0_);
                    cell1.getExternalField(velOffset, Descriptor<T>::d, u1_);

                    Array<T, Descriptor<T>::d> u0;
                    Array<T, Descriptor<T>::d> u1;
                    for (plint i = 0; i < Descriptor<T>::d; ++i)
                    {
                        u0[i] = u0_[i];
                        u1[i] = u1_[i];
                    }
                    Array<T, Descriptor<T>::d> jEq0(C0 * u0);
                    Array<T, Descriptor<T>::d> jEq1(C1 * u1);
                    auto jSqr = 0;
                    for (pluint iPop = 0; iPop < Descriptor<T>::q; ++iPop)
                    {
                        (*cell0)[iPop] =
                            (*cell0).computeEquilibrium(iPop, C0, jEq0, jSqr) //
                            + cell1[iPop] - cell1.computeEquilibrium(iPop, C1, jEq1, jSqr);
                    }
                }
            }
        }
        virtual NeeFlatDirichletBoundaryFunctional2D<T, Descriptor, direction, orientation> *clone() const
        {
            return new NeeFlatDirichletBoundaryFunctional2D<T, Descriptor, direction, orientation>(*this);
        }

        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
        {
            modified[0] = modif::dynamicVariables;
        }
        virtual BlockDomain::DomainT appliesTo() const
        {
            return BlockDomain::bulk;
        }

    private:
        T bcValue;
    };

}

#endif
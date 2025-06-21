#ifndef POISSONEXTRA_H
#define POISSONEXTRA_H

#include "../../src/palabos3D.h"
#include "../../src/palabos3D.hh"
#include "../../src/palabos2D.h"
#include "../../src/palabos2D.hh"

#include "poisson.h"
#include "poisson.hh"
namespace plb
{

    template <typename T, template <typename U> class Descriptor>
    void iniCellAtEquilibriumPoisson(Cell<T, Descriptor> &cell, T density)
    {
        Array<T, Descriptor<T>::d> jEq(0, 0);

        for (plint iPop = 0; iPop < Descriptor<T>::numPop; ++iPop)
        {
            cell[iPop] = cell.computeEquilibrium(iPop, density, jEq, T());
        }
    }

    template <typename T, template <typename U> class Descriptor>
    void iniCellAtEquilibriumPoisson3D(Cell<T, Descriptor> &cell, T density)
    {
        Array<T, Descriptor<T>::d> jEq(0, 0, 0);

        for (plint iPop = 0; iPop < Descriptor<T>::numPop; ++iPop)
        {
            cell[iPop] = cell.computeEquilibrium(iPop, density, jEq, T());
        }
    }

    template <typename T, template <typename U> class Descriptor>
    class IniConstEquilibriumFunctional2DPoisson : public BoxProcessingFunctional2D_L<T, Descriptor>
    {
    public:
        IniConstEquilibriumFunctional2DPoisson(T density) : rho(density) {}

        virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
        {
            Array<T, Descriptor<T>::d> jEq(0.0, 0.0);
            for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            {
                for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                {
                    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop)
                    {
                        lattice.get(iX, iY)[iPop] =
                            lattice.get(iX, iY).computeEquilibrium(iPop, rho, jEq, T());
                    }
                }
            }
        }

        virtual IniConstEquilibriumFunctional2DPoisson<T, Descriptor> *clone() const
        {
            return new IniConstEquilibriumFunctional2DPoisson<T, Descriptor>(*this);
        }

        virtual BlockDomain::DomainT appliesTo() const
        {
            return BlockDomain::bulkAndEnvelope;
        }
        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
        {
            modified[0] = modif::staticVariables;
        }

    private:
        T rho;
    };

    template <typename T, template <typename U> class Descriptor>
    class IniConstEquilibriumFunctional3DPoisson : public BoxProcessingFunctional3D_L<T, Descriptor>
    {
    public:
        IniConstEquilibriumFunctional3DPoisson(T density) : rho(density) {}

        virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
        {
            Array<T, Descriptor<T>::d> jEq(0.0, 0.0, 0.0);
            for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            {
                for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                {
                    for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ)
                    {
                        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop)
                        {
                            lattice.get(iX, iY, iZ)[iPop] =
                                lattice.get(iX, iY, iZ).computeEquilibrium(iPop, rho, jEq, T());
                        }
                    }
                }
            }
        }

        virtual IniConstEquilibriumFunctional3DPoisson<T, Descriptor> *clone() const
        {
            return new IniConstEquilibriumFunctional3DPoisson<T, Descriptor>(*this);
        }

        virtual BlockDomain::DomainT appliesTo() const
        {
            return BlockDomain::bulkAndEnvelope;
        }
        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
        {
            modified[0] = modif::staticVariables;
            modified[1] = modif::staticVariables;
        }

    private:
        T rho;
    };

    template <typename T, template <class U> class Descriptor>
    void initializeAtEquilibrium(MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain, T rho)
    {
        applyProcessingFunctional(new IniConstEquilibriumFunctional2DPoisson<T, Descriptor>(rho),
                                  domain, lattice);
    }

    template <typename T, template <class U> class Descriptor>
    void initializeAtEquilibrium(MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, T rho)
    {
        applyProcessingFunctional(new IniConstEquilibriumFunctional3DPoisson<T, Descriptor>(rho),
                                  domain, lattice);
    }
}
#endif
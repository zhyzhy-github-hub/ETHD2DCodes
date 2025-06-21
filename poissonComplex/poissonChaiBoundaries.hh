#ifndef POISSONCHAIBOUNDARIES_HH
#define POISSONCHAIBOUNDARIES_HH

#include "poissonChaiBoundaries.h"
#include "poissonChaiBoundaries.hh"

#include "../../src/core/util.h"
#include "../../src/complexDynamics/utilAdvectionDiffusion.h"
#include "../../src/latticeBoltzmann/advectionDiffusionLattices.h"
#include "../../src/latticeBoltzmann/advectionDiffusionDynamicsTemplates.h"
#include "../../src/latticeBoltzmann/indexTemplates.h"

namespace plb
{
    template <typename T, template <typename U> class Descriptor,
              int direction, int orientation>
    int PoissonLocalBoundaryDynamics<T, Descriptor, direction, orientation>::id =
        meta::registerGeneralDynamics<T, Descriptor, PoissonLocalBoundaryDynamics<T, Descriptor, direction, orientation>>(
            std::string("Boundary_Poisson") + util::val2str(direction) +
            std::string("_") + util::val2str(orientation));

    /// Constructor
    template <typename T, template <typename U> class Descriptor, int direction, int orientation>
    PoissonLocalBoundaryDynamics<T, Descriptor, direction, orientation>::
        PoissonLocalBoundaryDynamics(Dynamics<T, Descriptor> *baseDynamics,
                                     bool automaticPrepareCollision_)
        : StoreDensityDynamics<T, Descriptor>(baseDynamics, automaticPrepareCollision_)
    {
    }

    template <typename T, template <typename U> class Descriptor, int direction, int orientation>
    PoissonLocalBoundaryDynamics<T, Descriptor, direction, orientation>::
        PoissonLocalBoundaryDynamics(HierarchicUnserializer &unserializer)
        : StoreDensityDynamics<T, Descriptor>(0, false)
    {
        this->unserialize(unserializer);
    }

    /// Clone the object, based on its dynamic type
    template <typename T, template <typename U> class Descriptor, int direction, int orientation>
    PoissonLocalBoundaryDynamics<T, Descriptor, direction, orientation> *
    PoissonLocalBoundaryDynamics<T, Descriptor, direction, orientation>::clone() const
    {
        return new PoissonLocalBoundaryDynamics<T, Descriptor, direction, orientation>(*this);
    }

    // /// Return a unique ID for this class.
    template <typename T, template <typename U> class Descriptor, int direction, int orientation>
    int PoissonLocalBoundaryDynamics<T, Descriptor, direction, orientation>::getId()
    {
        return id;
    }

    template <typename T, template <typename U> class Descriptor, int direction, int orientation>
    void PoissonLocalBoundaryDynamics<T, Descriptor, direction, orientation>::completePopulations(Cell<T, Descriptor> &cell) const
    // virtual void completePopulations(Cell<T, Descriptor> &cell, Box2D domain, BlockLattice2D<T, Descriptor> &lattice) const
    {
        DensityClosurePoisson<T, Descriptor, direction, orientation>(cell, *this);

        // T rho = this->computeDensity(cell);
        // // T rhoBar = Descriptor<T>::rhoBar(rho);
        // T rhoBar = (rho);

        // Array<T, Descriptor<T>::d> jEq(0.0, 0.0);

        // for (plint iX = domain.x0; iX <= domain.x1; ++iX)
        // {
        //     for (plint iY = domain.y0; iY <= domain.y1; ++iY)
        //     {

        //         plint iX_prev = iX + ((direction == 0) ? (-orientation) : 0);
        //         plint iY_prev = iY + ((direction == 1) ? (-orientation) : 0);

        //         T phi1 = lattice.get(iX_prev, iY_prev).computeDensity();

        //         Cell<T, Descriptor> cell1 = lattice.get(iX_prev, iY_prev);

        //         for (plint i = 0; i < Descriptor<T>::q; ++i)
        //         {
        //             cell[i] = lattice.get(iX, iY).computeEquilibrium(i, rhoBar, jEq, T()) //
        //                       + cell1[i] - cell1.computeEquilibrium(i, phi1, jEq, T());
        //         }
        //     }
        // }
        // for (plint i = 0; i < Descriptor<T>::q; ++i)
        // {
        //    cell[i] = cell.computeEquilibrium(i, rhoBar, jEq, T());
        // }
    }

    // =============== 2D corners ===================//

    template <typename T, template <typename U> class Descriptor,
              int xNormal, int yNormal>
    int PoissonLocalCornerDynamics2D<T, Descriptor, xNormal, yNormal>::id =
        meta::registerGeneralDynamics<T, Descriptor, PoissonLocalBoundaryDynamics<T, Descriptor, xNormal, yNormal>>(
            std::string("Boundary_PoissonChai08Corner") + util::val2str(xNormal) +
            std::string("_") + util::val2str(yNormal));

    template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal>
    PoissonLocalCornerDynamics2D<T, Descriptor, xNormal, yNormal>::PoissonLocalCornerDynamics2D(
        Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision)
        : StoreDensityDynamics<T, Descriptor>(baseDynamics, automaticPrepareCollision)
    {
    }

    template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal>
    PoissonLocalCornerDynamics2D<T, Descriptor, xNormal, yNormal>::PoissonLocalCornerDynamics2D(
        HierarchicUnserializer &unserializer)
        : StoreDensityDynamics<T, Descriptor>(0, false)
    {
        this->unserialize(unserializer);
    }

    template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal>
    PoissonLocalCornerDynamics2D<T, Descriptor, xNormal, yNormal> *
    PoissonLocalCornerDynamics2D<T, Descriptor, xNormal, yNormal>::clone() const
    {
        return new PoissonLocalCornerDynamics2D<T, Descriptor, xNormal, yNormal>(*this);
    }

    template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal>
    int PoissonLocalCornerDynamics2D<T, Descriptor, xNormal, yNormal>::getId() const
    {
        return id;
    }

    template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal>
    void PoissonLocalCornerDynamics2D<T, Descriptor, xNormal, yNormal>::completePopulations(Cell<T, Descriptor> &cell) const
    {
        typedef Descriptor<T> D;
        typedef PoissonChai08DynamicsTemplatesImpl<T, Descriptor<T>> adTempl;

        T rho = this->computeDensity(cell);
        // I need to get Missing information on the corners !!!!
        std::vector<plint> unknownIndexes = utilAdvDiff::subIndexOutgoing2DonCorners<D, xNormal, yNormal>();
        // here I know all missing and non missing f_i

        // The collision procedure for D2Q5 and D3Q7 lattice is the same ...
        // Given the rule f_i_neq = -f_opposite(i)_neq
        // I have the right number of equations for the number of unknowns using these lattices

        for (pluint iPop = 0; iPop < unknownIndexes.size(); ++iPop)
        {
            cell[unknownIndexes[iPop]] =
                adTempl::bgk_ma1_equilibrium(unknownIndexes[iPop], rho)    //
                - (cell[indexTemplates::opposite<D>(unknownIndexes[iPop])] //
                   - adTempl::bgk_ma1_equilibrium(indexTemplates::opposite<D>(unknownIndexes[iPop]), rho));
        }
    }

    // =============== 3D edges ===================//

    template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
    int PoissonEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>::id =
        meta::registerGeneralDynamics<T, Descriptor, PoissonEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>>(
            std::string("Boundary_PoissonEdge") + util::val2str(plane) +
            std::string("_") + util::val2str(normal1) + std::string("_") + util::val2str(normal2));

    template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
    PoissonEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>::PoissonEdgeDynamics3D(
        Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision)
        : StoreDensityDynamics<T, Descriptor>(baseDynamics, automaticPrepareCollision)
    {
    }

    template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
    PoissonEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>::PoissonEdgeDynamics3D(
        HierarchicUnserializer &unserializer)
        : StoreDensityDynamics<T, Descriptor>(0, false)
    {
        this->unserialize(unserializer);
    }

    template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
    PoissonEdgeDynamics3D<T, Descriptor, plane, normal1, normal2> *
    PoissonEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>::clone() const
    {
        return new PoissonEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>(*this);
    }

    template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
    int PoissonEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>::getId() const
    {
        return id;
    }

    template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
    void PoissonEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>::completePopulations(Cell<T, Descriptor> &cell) const
    {
        typedef Descriptor<T> D;
        typedef PoissonChai08DynamicsTemplatesImpl<T, Descriptor<T>> adTempl;

        T rho = this->computeDensity(cell);
        // I need to get Missing information on the corners !!!!
        std::vector<plint> unknownIndexes = utilAdvDiff::subIndexOutgoing3DonEdges<D, plane, normal1, normal2>();

        // here I know all missing and non missing f_i

        // The collision procedure for D2Q5 and D3Q7 lattice is the same ...
        // Given the rule f_i_neq = -f_opposite(i)_neq
        // I have the right number of equations for the number of unknowns using these lattices

        for (pluint iPop = 0; iPop < unknownIndexes.size(); ++iPop)
        {
            cell[unknownIndexes[iPop]] =
                adTempl::bgk_ma1_equilibrium(unknownIndexes[iPop], rho)    //
                - (cell[indexTemplates::opposite<D>(unknownIndexes[iPop])] //
                   - adTempl::bgk_ma1_equilibrium(indexTemplates::opposite<D>(unknownIndexes[iPop]), rho));
        }
    }

    // =============== 3D corners ===================//

    template <typename T, template <typename U> class Descriptor,
              int xNormal, int yNormal, int zNormal>
    int PoissonCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal>::id =
        meta::registerGeneralDynamics<T, Descriptor, PoissonCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal>>(
            std::string("Boundary_PoissonCorner") + util::val2str(xNormal) +
            std::string("_") + util::val2str(yNormal) + std::string("_") + util::val2str(zNormal));

    template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
    PoissonCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal>::PoissonCornerDynamics3D(
        Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision)
        : StoreDensityDynamics<T, Descriptor>(baseDynamics, automaticPrepareCollision)
    {
    }

    template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
    PoissonCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal>::PoissonCornerDynamics3D(
        HierarchicUnserializer &unserializer)
        : StoreDensityDynamics<T, Descriptor>(0, false)
    {
        this->unserialize(unserializer);
    }

    template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
    PoissonCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal> *
    PoissonCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal>::clone() const
    {
        return new PoissonCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal>(*this);
    }

    template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
    int PoissonCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal>::getId() const
    {
        return id;
    }

    template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
    void PoissonCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal>::completePopulations(Cell<T, Descriptor> &cell) const
    {
        typedef Descriptor<T> D;
        typedef PoissonChai08DynamicsTemplatesImpl<T, Descriptor<T>> adTempl;

        T rho = this->computeDensity(cell);
        // I need to get Missing information on the corners !!!!
        std::vector<plint> unknownIndexes = utilAdvDiff::subIndexOutgoing3DonCorners<D, xNormal, yNormal, zNormal>();

        // here I know all missing and non missing f_i

        // The collision procedure for D2Q5 and D3Q7 lattice is the same ...
        // Given the rule f_i_neq = -f_opposite(i)_neq
        // I have the right number of equations for the number of unknowns using these lattices

        for (pluint iPop = 0; iPop < unknownIndexes.size(); ++iPop)
        {
            cell[unknownIndexes[iPop]] =
                adTempl::bgk_ma1_equilibrium(unknownIndexes[iPop], rho)    //
                - (cell[indexTemplates::opposite<D>(unknownIndexes[iPop])] //
                   - adTempl::bgk_ma1_equilibrium(indexTemplates::opposite<D>(unknownIndexes[iPop]), rho));
        }
    }

}
#endif
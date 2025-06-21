#ifndef POISSONCHAIBOUNDARIES_H
#define POISSONCHAIBOUNDARIES_H

#include "poissonChaiDynamics.h"

#include "../../src/core/globalDefs.h"
#include "../../src/complexDynamics/advectionDiffusionDynamics.h"
#include "../../src/boundaryCondition/boundaryDynamics.h"

namespace plb
{
    template <typename T, template <typename U> class Descriptor, int direction, int orientation>
    void DensityClosurePoisson(Cell<T, Descriptor> &cell, Dynamics<T, Descriptor> const &dynamics)
    {
        typedef Descriptor<T> D;

        T rho = dynamics.computeDensity(cell);
        T rhoBar = (rho);

        plint missingNormal = 0;
        std::vector<plint> missingDiagonal = indexTemplates::subIndexOutgoing<D, direction, orientation>();
        std::vector<plint> knownIndexes = indexTemplates::remainingIndexes<D>(missingDiagonal);
        // here I know all missing and non missing f_i
        for (pluint iPop = 0; iPop < missingDiagonal.size(); ++iPop)
        {
            plint numOfNonNullComp = 0;
            for (int iDim = 0; iDim < D::d; ++iDim)
                numOfNonNullComp += abs(D::c[missingDiagonal[iPop]][iDim]);

            if (numOfNonNullComp == 1)
            {
                missingNormal = missingDiagonal[iPop];
                missingDiagonal.erase(missingDiagonal.begin() + iPop);
                break;
            }
        }

        T sum = T();
        for (pluint iPop = 1; iPop < knownIndexes.size(); ++iPop)
        {
            sum += cell[knownIndexes[iPop]];
        }
        cell[missingNormal] = rhoBar * (T(1.) - D::t[0]) - sum;
    }

    /// Advection-diffusion dynamics on flat boundaries
    template <typename T, template <typename U> class Descriptor, int direction, int orientation>
    class PoissonLocalBoundaryDynamics : public StoreDensityDynamics<T, Descriptor>
    {
    public:
        /// Constructor
        PoissonLocalBoundaryDynamics(Dynamics<T, Descriptor> *baseDynamics,
                                     bool automaticPrepareCollision_ = true);
        // : StoreDensityDynamics<T, Descriptor>(baseDynamics, automaticPrepareCollision_) {}

        PoissonLocalBoundaryDynamics(HierarchicUnserializer &unserializer);
        // : StoreDensityDynamics<T, Descriptor>(0, false) { this->unserialize(unserializer); }

        /// Clone the object, based on its dynamic type
        virtual PoissonLocalBoundaryDynamics<T, Descriptor, direction, orientation> *clone() const;
        // {
        //     return new PoissonLocalBoundaryDynamics<T, Descriptor, direction, orientation>(*this);
        // }
        // /// Return a unique ID for this class.
        virtual int getId();
        virtual void completePopulations(Cell<T, Descriptor> &cell) const;

    private:
        static int id;
    };

    /// Poisson  dynamics on 2D corners
    template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal>
    class PoissonLocalCornerDynamics2D : public StoreDensityDynamics<T, Descriptor>
    {
    public:
        /// Constructor
        PoissonLocalCornerDynamics2D(Dynamics<T, Descriptor> *baseDynamics,
                                     bool automaticPrepareCollision_ = true);
        PoissonLocalCornerDynamics2D(HierarchicUnserializer &unserializer);

        /// Return a unique ID for this class.
        virtual int getId() const;
        /// Clone the object on its dynamic type.
        virtual PoissonLocalCornerDynamics2D<T, Descriptor, xNormal, yNormal> *clone() const;

        /// Execute completion scheme before base collision
        virtual void completePopulations(Cell<T, Descriptor> &cell) const;

    private:
        static int id;
    };

    /// Poisson dynamics on 3D edges
    template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
    class PoissonEdgeDynamics3D : public StoreDensityDynamics<T, Descriptor>
    {
    public:
        /// Constructor
        PoissonEdgeDynamics3D(Dynamics<T, Descriptor> *baseDynamics,
                              bool automaticPrepareCollision_ = true);
        PoissonEdgeDynamics3D(HierarchicUnserializer &unserializer);

        /// Return a unique ID for this class.
        virtual int getId() const;

        /// Clone the object, based on its dynamic type
        virtual PoissonEdgeDynamics3D<T, Descriptor, plane, normal1, normal2> *clone() const;

        /// Execute completion scheme before base collision
        virtual void completePopulations(Cell<T, Descriptor> &cell) const;

    private:
        static int id;
    };

    /// Poisson -dynamics on 3D corners
    template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
    class PoissonCornerDynamics3D : public StoreDensityDynamics<T, Descriptor>
    {
    public:
        /// Constructor
        PoissonCornerDynamics3D(Dynamics<T, Descriptor> *baseDynamics,
                                           bool automaticPrepareCollision_ = true);
        PoissonCornerDynamics3D(HierarchicUnserializer &unserializer);

        /// Return a unique ID for this class.
        virtual int getId() const;

        /// Clone the object on its dynamic type.
        virtual PoissonCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal> *clone() const;

        /// Execute completion scheme before base collision
        virtual void completePopulations(Cell<T, Descriptor> &cell) const;

    private:
        static int id;
    };
}
#endif
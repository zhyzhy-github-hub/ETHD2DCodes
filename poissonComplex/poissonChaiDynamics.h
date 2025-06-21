#ifndef POISSONCHAIDYNAMICS_H
#define POISSONCHAIDYNAMICS_H

#include "../../src/core/array.h"
#include "../../src/core/globalDefs.h"
#include "../../src/core/dynamics.h"

#include "poissonChaiLattices.h"

namespace plb
{
    template <typename T, template <typename U> class Descriptor>
    class PoissonChai08WithSourceBGKdynamics : public AdvectionDiffusionWithSourceBGKdynamics<T, Descriptor>
    // class PoissonChai08WithSourceBGKdynamics : public AdvectionDiffusionDynamics<T, Descriptor>
    {
    public:
        /// Constructor
        PoissonChai08WithSourceBGKdynamics(T omega_);

        PoissonChai08WithSourceBGKdynamics(HierarchicUnserializer &unserializer);
        /// Clone the object on its dynamic type.
        virtual PoissonChai08WithSourceBGKdynamics<T, Descriptor> *clone() const;
        /// Return a unique ID for this class.
        virtual int getId() const;
        /// Collision step
        virtual void collide(Cell<T, Descriptor> &cell,
                             BlockStatistics &statistics);
        /// Compute equilibrium distribution function
        virtual T computeEquilibrium(plint iPop, T rhoBar,
                                     Array<T, Descriptor<T>::d> const &j,
                                     T jSqr, T thetaBar);

        T getTau() { return T(1.0) / omega; }
        T getDiffusivity() { return Descriptor<T>::cs2 * (0.5 - this->getTau()); }

        virtual T computeDensity(Cell<T, Descriptor> const &cell) const;

        virtual void computeVelocity(Cell<T, Descriptor> const &cell,
                                     Array<T, Descriptor<T>::d> &q) const
        {
            computeNegativeGrad(cell, q);
        }

        void computeNegativeGrad(Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &grad) const
        {
            for (int iD = 0; iD < Descriptor<T>::d; ++iD)
            {
                grad[iD] = cell[0] * Descriptor<T>::c[0][iD];
            }
            for (plint iPop = 1; iPop < Descriptor<T>::q; ++iPop)
            {
                for (int iD = 0; iD < Descriptor<T>::d; ++iD)
                {
                    grad[iD] += cell[iPop] * Descriptor<T>::c[iPop][iD];
                }
            }
            T temp = omega * Descriptor<T>::invCs2;
            for (int iD = 0; iD < Descriptor<T>::d; ++iD)
            {
                grad[iD] *= temp;
            }
        }

    private:
        static int id;
        T omega;
        // T getDiffusivity() { return 0.5 * (0.5 - this->getTau()); }
        // T D;
    };

}
#endif
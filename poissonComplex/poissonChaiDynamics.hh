#ifndef POISSONCHAIDYNAMICS_HH
#define POISSONCHAIDYNAMICS_HH

#include "../../src/core/latticeStatistics.h"
#include "../../src/core/dynamicsIdentifiers.h"
#include "../../src/latticeBoltzmann/momentTemplates.h"
#include "../../src/latticeBoltzmann/offEquilibriumAdvectionDiffusionTemplates.h"

#include "poissonChaiDynamics.h"

namespace plb
{
    // This structure forwards the calls to the appropriate helper class
    template <typename T, template <typename U> class Descriptor>
    struct PoissonChai08MomentTemplates
    {
        static void get_rhoBar(Array<T, Descriptor<T>::q> const &f, T &rhoBar)
        {
            // rhoBar = momentTemplatesImpl<T, Descriptor<T>>::get_rhoBar(f);
            T phi = 0;
            for (plint i = 1; i < Descriptor<T>::q; ++i)
            {
                phi += f[i];
            }
            rhoBar = (phi) / (1.0 - Descriptor<T>::t[0]);
        }
    };

    /// All helper functions are inside this structure
    template <typename T, class Descriptor>
    struct PoissonChai08DynamicsTemplatesImpl
    {
        static T bgk_ma1_equilibrium(plint iPop, T rhoBar)
        {
            T feq{0.};
            if (iPop == 0)
            {
                feq = (Descriptor::t[iPop] - 1.0) * rhoBar;
                // feq = (Descriptor::t[iPop]) * rhoBar;
                // feq = -rhoBar;
            }
            else
            {
                feq = Descriptor::t[iPop] * rhoBar;
            }
            return feq;
        }

        static T no_corr_bgk_collision(
            Array<T, Descriptor::q> &f, T rhoBar, T omega, T sourceTerm)
        {
            // T invRho = Descriptor::invRho(rhoBar);
            // const T jSqr = <T, Descriptor::d>::normSqr(jEq);

            for (plint iPop = 0; iPop < Descriptor::q; ++iPop)
            {
                f[iPop] *= (T)1.0 - omega;
                f[iPop] += omega * PoissonChai08DynamicsTemplatesImpl<T, Descriptor>::bgk_ma1_equilibrium(
                                       iPop, rhoBar);
                // sourceTerm = 0.1;
                f[iPop] += Descriptor::t[iPop] * sourceTerm * Descriptor::cs2 * (0.5 - T(1.) / omega);

                // f[iPop] = PoissonChai08DynamicsTemplatesImpl<T, Descriptor>::bgk_ma1_equilibrium(iPop, rhoBar);
            }
            return rhoBar;
            // return T(0.0);
        }
    };

    /// This structure forwards the calls to the appropriate helper class
    template <typename T, template <typename U> class Descriptor>
    struct PoissonChai08DynamicsTemplates
    {
        static T bgk_ma1_equilibrium(plint iPop, T rhoBar)
        {
            return PoissonChai08DynamicsTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::bgk_ma1_equilibrium(iPop, rhoBar);
        }

        static T no_corr_bgk_collision(Cell<T, Descriptor> &cell, T rhoBar, T omega, T sourceTerm)
        {
            return PoissonChai08DynamicsTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::
                no_corr_bgk_collision(cell.getRawPopulations(), rhoBar, omega, sourceTerm);
        }
    };

    /* *************** Class AdvectionDiffusionWithSourceBGKdynamics *************** */

    template <typename T, template <typename U> class Descriptor>
    int PoissonChai08WithSourceBGKdynamics<T, Descriptor>::id =
        meta::registerGeneralDynamics<T, Descriptor, PoissonChai08WithSourceBGKdynamics<T, Descriptor>>("PoissonChai08WithSource_BGK");

    /** \param omega_ relaxation parameter, related to the dynamic diffusivity
 */
    template <typename T, template <typename U> class Descriptor>
    PoissonChai08WithSourceBGKdynamics<T, Descriptor>::PoissonChai08WithSourceBGKdynamics(T omega_)
        : AdvectionDiffusionWithSourceBGKdynamics<T, Descriptor>(omega_), omega(omega_)
    {
        // AdvectionDiffusionWithSourceBGKdynamics<T, Descriptor>(omega);
    }

    template <typename T, template <typename U> class Descriptor>
    PoissonChai08WithSourceBGKdynamics<T, Descriptor>::PoissonChai08WithSourceBGKdynamics(HierarchicUnserializer &unserializer)
        : AdvectionDiffusionWithSourceBGKdynamics<T, Descriptor>(T())
    {
        this->unserialize(unserializer);
    }

    /// Clone the object on its dynamic type.
    template <typename T, template <typename U> class Descriptor>
    PoissonChai08WithSourceBGKdynamics<T, Descriptor> *PoissonChai08WithSourceBGKdynamics<T, Descriptor>::clone() const
    {
        return new PoissonChai08WithSourceBGKdynamics<T, Descriptor>(*this);
    }

    template <typename T, template <typename U> class Descriptor>
    int PoissonChai08WithSourceBGKdynamics<T, Descriptor>::getId() const
    {
        return id;
    }

    /// Collision step
    template <typename T, template <typename U> class Descriptor>
    void PoissonChai08WithSourceBGKdynamics<T, Descriptor>::collide(
        Cell<T, Descriptor> &cell, BlockStatistics &statistics)
    {
        // T phi;
        // PoissonChai08MomentTemplates<T, Descriptor>::get_rhoBar(cell.getRawPopulations(), phi);
        // phi = (phi - cell[0]) / (1.0 - Descriptor<T>::t[0]);
        T phi = computeDensity(cell);
        // phi = 1.0;

        T sourceTerm = *cell.getExternal(Descriptor<T>::ExternalField::scalarBeginsAt);

        T uSqr = PoissonChai08DynamicsTemplates<T, Descriptor>::
            no_corr_bgk_collision(cell, phi, this->getOmega(), sourceTerm);

        if (cell.takesStatistics())
        {
            gatherStatistics(statistics, phi, uSqr);
        }
    }

    /// Compute equilibrium distribution function

    template <typename T, template <typename U> class Descriptor>
    T PoissonChai08WithSourceBGKdynamics<T, Descriptor>::computeEquilibrium(plint iPop, T rhoBar,
                                                                            Array<T, Descriptor<T>::d> const &j,
                                                                            T jSqr, T thetaBar)
    {
        return PoissonChai08DynamicsTemplates<T, Descriptor>::
            bgk_ma1_equilibrium(iPop, rhoBar);
    }

    template <typename T, template <typename U> class Descriptor>
    T PoissonChai08WithSourceBGKdynamics<T, Descriptor>::computeDensity(Cell<T, Descriptor> const &cell) const
    {
        T phi = 0;
        for (plint i = 1; i < Descriptor<T>::q; ++i)
        {
            phi += cell[i];
        }

        phi = (phi) / (1.0 - Descriptor<T>::t[0]);
        // return 0.4;
        return phi;
    }
}
#endif
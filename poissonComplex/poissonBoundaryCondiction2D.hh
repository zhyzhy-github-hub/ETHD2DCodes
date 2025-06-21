#ifndef poissonBoundaryCondiction2D_HH
#define poissonBoundaryCondiction2D_HH

#include "../../src/complexDynamics/advectionDiffusionBoundaryCondition2D.h"
#include "poissonBoundaryCondiction2D.h"
namespace plb
{
    template <typename T, template <typename U> class Descriptor>
    OnLatticeAdvectionDiffusionBoundaryCondition2D<T, Descriptor> *
    createLocalPoissonBoundaryCondition2D()
    {
        return new PoissonBoundaryConditionInstantiator2D<
            T, Descriptor, AdvectionDiffusionBoundaryManager2D<T, Descriptor>>();
    }

    template <typename T, template <typename U> class Descriptor>
    OnLatticeAdvectionDiffusionBoundaryCondition2D<T, Descriptor> *
    createLocalPoissonPoissonLocalBoundaryCondition2D()
    {
        return new PoissonLocalBoundaryConditionInstantiator2D<
            T, Descriptor, PoissonLocalBoundaryManager2D<T, Descriptor>>();
    }

    // template <typename T, template <typename U> class Descriptor, int direction, int orientation>
    // int PoissonLocalBoundaryDynamics<T, Descriptor, direction, orientation>::getId()
    // {
    //     return id;
    // }

    template <typename T, template <typename U> class Descriptor, class BoundaryManager>
    template <int direction, int orientation>
    void PoissonLocalBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
        addTemperatureBoundary(Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
    {
        PLB_PRECONDITION(domain.x0 == domain.x1 || domain.y0 == domain.y1);

        setCompositeDynamics(
            lattice, domain,
            BoundaryManager::template getPoissonLocalBoundaryDynamics<direction, orientation>(new NoDynamics<T, Descriptor>));

        // If the boundary condition has a non-local component, instantiate a corresponding data processor.
        BoxProcessingFunctional2D_L<T, Descriptor> *functional = BoundaryManager ::template getPoissonLocalBoundaryProcessor<direction, orientation>();
        if (functional)
        {
            integrateProcessingFunctional(functional, domain, lattice);
        }
    }

    template <typename T, template <typename U> class Descriptor, class BoundaryManager>
    template <int direction, int orientation>
    void PoissonLocalBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
        addTemperatureBoundary(Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice)
    {
        PLB_PRECONDITION(domain.x0 == domain.x1 || domain.y0 == domain.y1);

        setCompositeDynamics(
            lattice, domain,
            BoundaryManager::template getPoissonLocalBoundaryDynamics<direction, orientation>(new NoDynamics<T, Descriptor>));

        // If the boundary condition has a non-local component, instantiate a corresponding data processor.
        BoxProcessingFunctional2D_L<T, Descriptor> *functional = BoundaryManager::template getPoissonLocalBoundaryProcessor<direction, orientation>();
        if (functional)
        {
            integrateProcessingFunctional(functional, domain, lattice);
        }
    }


    template <typename T, template <typename U> class Descriptor, class BoundaryManager>
    template <int xNormal, int yNormal>
    void PoissonLocalBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
        PoissonLocalCorner(plint x, plint y, BlockLattice2D<T, Descriptor> &lattice)
    {
        setCompositeDynamics(
            lattice, Box2D(x, x, y, y),
            BoundaryManager::template getPoissonLocalCornerDynamics<xNormal, yNormal>(new NoDynamics<T, Descriptor>));

        // If the boundary condition has a non-local component, instantiate a corresponding data processor.
        BoxProcessingFunctional2D_L<T, Descriptor> *functional = BoundaryManager::template getPoissonLocalCornerProcessor<xNormal, yNormal>();
        if (functional)
        {
            integrateProcessingFunctional(functional, Box2D(x, x, y, y), lattice);
        }
    }

    template <typename T, template <typename U> class Descriptor, class BoundaryManager>
    template <int xNormal, int yNormal>
    void PoissonLocalBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
        PoissonLocalCorner(plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice)
    {
        setCompositeDynamics(
            lattice, Box2D(x, x, y, y),
            BoundaryManager::template getPoissonLocalCornerDynamics<xNormal, yNormal>(new NoDynamics<T, Descriptor>));

        // If the boundary condition has a non-local component, instantiate a corresponding data processor.
        BoxProcessingFunctional2D_L<T, Descriptor> *functional = BoundaryManager::template getPoissonLocalCornerProcessor<xNormal, yNormal>();
        if (functional)
        {
            integrateProcessingFunctional(functional, Box2D(x, x, y, y), lattice);
        }
    }
}

#endif
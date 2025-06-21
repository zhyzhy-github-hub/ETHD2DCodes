#ifndef poissonBoundaryCondiction3D_HH
#define poissonBoundaryCondiction3D_HH

#include "../../src/complexDynamics/advectionDiffusionBoundaryCondition3D.h"
#include "../../src/complexDynamics/adiabaticBoundaryProcessor3D.h"
#include "poissonBoundaryCondiction3D.h"
namespace plb
{

    ////////// AdvectionDiffusionBoundaryManager3D /////////////////////////////////////////

    template <typename T, template <typename U> class Descriptor>
    template <int direction, int orientation>
    BoundaryCompositeDynamics<T, Descriptor> *PoissonBoundaryManager3D<T, Descriptor>::
        getPoissonBoundaryDynamics(Dynamics<T, Descriptor> *baseDynamics)
    {
        return new PoissonLocalBoundaryDynamics<T, Descriptor, direction, orientation>(baseDynamics);
    }

    template <typename T, template <typename U> class Descriptor>
    template <int direction, int orientation>
    BoxProcessingFunctional3D_L<T, Descriptor> *PoissonBoundaryManager3D<T, Descriptor>::
        getPoissonBoundaryProcessor(Box3D domain)
    {
        return 0;
    }

    template <typename T, template <typename U> class Descriptor>
    template <int plane, int normal1, int normal2>
    BoundaryCompositeDynamics<T, Descriptor> *PoissonBoundaryManager3D<T, Descriptor>::
        getPoissonEdgeDynamics(Dynamics<T, Descriptor> *baseDynamics)
    {
        return new PoissonEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>(baseDynamics);
    }

    template <typename T, template <typename U> class Descriptor>
    template <int plane, int normal1, int normal2>
    BoxProcessingFunctional3D_L<T, Descriptor> *PoissonBoundaryManager3D<T, Descriptor>::
        getPoissonEdgeProcessor(Box3D domain)
    {
        return 0;
    }

    template <typename T, template <typename U> class Descriptor>
    template <int xNormal, int yNormal, int zNormal>
    BoundaryCompositeDynamics<T, Descriptor> *PoissonBoundaryManager3D<T, Descriptor>::
        getPoissonCornerDynamics(Dynamics<T, Descriptor> *baseDynamics)
    {
        return new PoissonCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal>(baseDynamics);
    }

    template <typename T, template <typename U> class Descriptor>
    template <int xNormal, int yNormal, int zNormal>
    BoxProcessingFunctional3D_L<T, Descriptor> *PoissonBoundaryManager3D<T, Descriptor>::
        getPoissonCornerProcessor(plint x, plint y, plint z)
    {
        return 0;
    }

    ///////// class PoissonBoundaryConditionInstantiator3D ////////////////////////

    template <typename T, template <typename U> class Descriptor, class BoundaryManager>
    PoissonChai08BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
        PoissonChai08BoundaryConditionInstantiator3D()
    {
    }

    template <typename T, template <typename U> class Descriptor, class BoundaryManager>
    template <int direction, int orientation>
    void PoissonChai08BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
        addPoissonBoundary(Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
                           boundary::BcType bcType)
    {
        PLB_PRECONDITION(domain.x0 == domain.x1 || domain.y0 == domain.y1 || domain.z0 == domain.z1);

        setCompositeDynamics(
            lattice, domain,
            BoundaryManager::template getPoissonBoundaryDynamics<direction, orientation>(new NoDynamics<T, Descriptor>));

        // In case an outflow condition is used, start by instantiating a data processor which copies
        //   all velocity values from the previous lattice cell.
        if (bcType == boundary::neumann)
        {
            integrateProcessingFunctional(
                new FlatAdiabaticBoundaryFunctional3D<T, Descriptor, direction, orientation>,
                domain, lattice);
        }

        BoxProcessingFunctional3D_L<T, Descriptor> *functional = BoundaryManager::template getPoissonBoundaryProcessor<direction, orientation>(domain);
        if (functional)
        {
            integrateProcessingFunctional(functional, domain, lattice);
        }
    }

    template <typename T, template <typename U> class Descriptor, class BoundaryManager>
    template <int direction, int orientation>
    void PoissonChai08BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
        addPoissonBoundary(Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
                           boundary::BcType bcType)
    {
        PLB_PRECONDITION(domain.x0 == domain.x1 || domain.y0 == domain.y1 || domain.z0 == domain.z1);

        setCompositeDynamics(
            lattice, domain,
            BoundaryManager::template getPoissonBoundaryDynamics<direction, orientation>(new NoDynamics<T, Descriptor>));

        BoxProcessingFunctional3D_L<T, Descriptor> *functional = BoundaryManager::template getPoissonBoundaryProcessor<direction, orientation>(domain);
        if (functional)
        {
            integrateProcessingFunctional(functional, domain, lattice);
        }
    }

    template <typename T, template <typename U> class Descriptor, class BoundaryManager>
    template <int plane, int normal1, int normal2>
    void PoissonChai08BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
        addPoissonEdge(Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
                       boundary::BcType bcType)
    {
        PLB_PRECONDITION(
            (domain.x0 == domain.x1 && domain.y0 == domain.y1) ||
            (domain.x0 == domain.x1 && domain.z0 == domain.z1) ||
            (domain.y0 == domain.y1 && domain.z0 == domain.z1));

        setCompositeDynamics(
            lattice, domain,
            BoundaryManager::template getPoissonEdgeDynamics<plane, normal1, normal2>(new NoDynamics<T, Descriptor>));

        BoxProcessingFunctional3D_L<T, Descriptor> *functional = BoundaryManager::template getPoissonEdgeProcessor<plane, normal1, normal2>(domain);
        if (functional)
        {
            integrateProcessingFunctional(functional, domain, lattice);
        }
    }

    template <typename T, template <typename U> class Descriptor, class BoundaryManager>
    template <int xNormal, int yNormal, int zNormal>
    void PoissonChai08BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
        addPoissonCorner(plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice,
                         boundary::BcType bcType)
    {
        setCompositeDynamics(
            lattice, Box3D(x, x, y, y, z, z),
            BoundaryManager::template getPoissonCornerDynamics<xNormal, yNormal, zNormal>(new NoDynamics<T, Descriptor>));

        BoxProcessingFunctional3D_L<T, Descriptor> *functional = BoundaryManager::template getPoissonCornerProcessor<xNormal, yNormal, zNormal>(x, y, z);
        if (functional)
        {
            integrateProcessingFunctional(functional, Box3D(x, x, y, y, z, z), lattice);
        }
    }

    template <typename T, template <typename U> class Descriptor, class BoundaryManager>
    template <int plane, int normal1, int normal2>
    void PoissonChai08BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
        addPoissonEdge(Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
                       boundary::BcType bcType)
    {
        PLB_PRECONDITION(
            (domain.x0 == domain.x1 && domain.y0 == domain.y1) ||
            (domain.x0 == domain.x1 && domain.z0 == domain.z1) ||
            (domain.y0 == domain.y1 && domain.z0 == domain.z1));

        setCompositeDynamics(
            lattice, domain,
            BoundaryManager::template getPoissonEdgeDynamics<plane, normal1, normal2>(new NoDynamics<T, Descriptor>));

        BoxProcessingFunctional3D_L<T, Descriptor> *functional = BoundaryManager::template getPoissonEdgeProcessor<plane, normal1, normal2>(domain);
        if (functional)
        {
            integrateProcessingFunctional(functional, domain, lattice);
        }
    }

    template <typename T, template <typename U> class Descriptor, class BoundaryManager>
    template <int xNormal, int yNormal, int zNormal>
    void PoissonChai08BoundaryConditionInstantiator3D<T, Descriptor, BoundaryManager>::
        addPoissonCorner(plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
                         boundary::BcType bcType)
    {
        setCompositeDynamics(
            lattice, Box3D(x, x, y, y, z, z),
            BoundaryManager::template getPoissonCornerDynamics<xNormal, yNormal, zNormal>(new NoDynamics<T, Descriptor>));

        BoxProcessingFunctional3D_L<T, Descriptor> *functional = BoundaryManager::template getPoissonCornerProcessor<xNormal, yNormal, zNormal>(x, y, z);
        if (functional)
        {
            integrateProcessingFunctional(functional, Box3D(x, x, y, y, z, z), lattice);
        }
    }

    template <typename T, template <typename U> class Descriptor>
    OnLatticeAdvectionDiffusionBoundaryCondition3D<T, Descriptor> *
    createLocalPoissonPoissonLocalBoundaryCondition3D()
    {
        return new PoissonChai08BoundaryConditionInstantiator3D<
            T, Descriptor, PoissonBoundaryManager3D<T, Descriptor>>();
    }

}
#endif
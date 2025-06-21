#ifndef poissonBoundaryCondiction3D_H
#define poissonBoundaryCondiction3D_H

#include "../../src/complexDynamics/advectionDiffusionBoundaryCondition3D.h"
// #include "../../src/complexDynamics/advectionDiffusionBoundaryInstantiator3D.h"

#include "poissonChaiDynamics.h"
#include "poissonChaiBoundaries.h"

namespace plb
{
    template <typename T, template <typename U> class Descriptor, class BoundaryManager>
    class PoissonChai08BoundaryConditionInstantiator3D : public OnLatticeAdvectionDiffusionBoundaryCondition3D<T, Descriptor>
    {
    public:
        PoissonChai08BoundaryConditionInstantiator3D();

        void addTemperatureBoundary0N(
            Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonBoundary<0, -1>(domain, lattice, bcType); }

        void addTemperatureBoundary0P(
            Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType = boundary::dirichlet) { addPoissonBoundary<0, 1>(domain, lattice, bcType); }

        void addTemperatureBoundary1N(
            Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonBoundary<1, -1>(domain, lattice, bcType); }

        void addTemperatureBoundary1P(
            Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonBoundary<1, 1>(domain, lattice, bcType); }

        void addTemperatureBoundary2N(
            Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonBoundary<2, -1>(domain, lattice, bcType); }

        void addTemperatureBoundary2P(
            Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonBoundary<2, 1>(domain, lattice, bcType); }

        void addTemperatureEdge0NN(
            Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<0, -1, -1>(domain, lattice, bcType); }

        void addTemperatureEdge0NP(
            Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<0, -1, 1>(domain, lattice, bcType); }

        void addTemperatureEdge0PN(
            Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<0, 1, -1>(domain, lattice, bcType); }

        void addTemperatureEdge0PP(
            Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<0, 1, 1>(domain, lattice, bcType); }

        void addTemperatureEdge1NN(
            Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<1, -1, -1>(domain, lattice, bcType); }

        void addTemperatureEdge1NP(
            Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<1, -1, 1>(domain, lattice, bcType); }

        void addTemperatureEdge1PN(
            Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<1, 1, -1>(domain, lattice, bcType); }

        void addTemperatureEdge1PP(
            Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<1, 1, 1>(domain, lattice, bcType); }

        void addTemperatureEdge2NN(
            Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<2, -1, -1>(domain, lattice, bcType); }

        void addTemperatureEdge2NP(
            Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<2, -1, 1>(domain, lattice, bcType); }

        void addTemperatureEdge2PN(
            Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<2, 1, -1>(domain, lattice, bcType); }

        void addTemperatureEdge2PP(
            Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<2, 1, 1>(domain, lattice, bcType); }

        void addTemperatureCornerNNN(
            plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonCorner<-1, -1, -1>(x, y, z, lattice, bcType); }

        void addTemperatureCornerNNP(
            plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonCorner<-1, -1, 1>(x, y, z, lattice, bcType); }

        void addTemperatureCornerNPN(
            plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonCorner<-1, 1, -1>(x, y, z, lattice, bcType); }

        void addTemperatureCornerNPP(
            plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonCorner<-1, 1, 1>(x, y, z, lattice, bcType); }
        void addTemperatureCornerPNN(
            plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonCorner<1, -1, -1>(x, y, z, lattice, bcType); }

        void addTemperatureCornerPNP(
            plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonCorner<1, -1, 1>(x, y, z, lattice, bcType); }

        void addTemperatureCornerPPN(
            plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonCorner<1, 1, -1>(x, y, z, lattice, bcType); }

        void addTemperatureCornerPPP(
            plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonCorner<1, 1, 1>(x, y, z, lattice, bcType); }

        void addTemperatureBoundary0N(
            Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonBoundary<0, -1>(domain, lattice, bcType); }

        void addTemperatureBoundary0P(
            Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonBoundary<0, 1>(domain, lattice, bcType); }

        void addTemperatureBoundary1N(
            Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonBoundary<1, -1>(domain, lattice, bcType); };

        void addTemperatureBoundary1P(
            Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonBoundary<1, 1>(domain, lattice, bcType); }

        void addTemperatureBoundary2N(
            Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonBoundary<2, -1>(domain, lattice, bcType); }

        void addTemperatureBoundary2P(
            Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonBoundary<2, 1>(domain, lattice, bcType); }

        void addTemperatureEdge0NN(
            Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<0, -1, -1>(domain, lattice, bcType); }

        void addTemperatureEdge0NP(
            Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<0, -1, 1>(domain, lattice, bcType); }

        void addTemperatureEdge0PN(
            Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<0, 1, -1>(domain, lattice, bcType); }

        void addTemperatureEdge0PP(
            Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<0, 1, 1>(domain, lattice, bcType); }

        void addTemperatureEdge1NN(
            Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<1, -1, -1>(domain, lattice, bcType); }

        void addTemperatureEdge1NP(
            Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<1, -1, 1>(domain, lattice, bcType); }

        void addTemperatureEdge1PN(
            Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<1, 1, -1>(domain, lattice, bcType); }

        void addTemperatureEdge1PP(
            Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<1, 1, 1>(domain, lattice, bcType); }

        void addTemperatureEdge2NN(
            Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<2, -1, -1>(domain, lattice, bcType); }

        void addTemperatureEdge2NP(
            Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<2, -1, 1>(domain, lattice, bcType); }

        void addTemperatureEdge2PN(
            Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<2, 1, -1>(domain, lattice, bcType); }

        void addTemperatureEdge2PP(
            Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonEdge<2, 1, 1>(domain, lattice, bcType); }

        void addTemperatureCornerNNN(
            plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonCorner<-1, -1, -1>(x, y, z, lattice, bcType); }

        void addTemperatureCornerNNP(
            plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonCorner<-1, -1, 1>(x, y, z, lattice, bcType); }

        void addTemperatureCornerNPN(
            plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonCorner<-1, 1, -1>(x, y, z, lattice, bcType); }

        void addTemperatureCornerNPP(
            plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonCorner<-1, 1, 1>(x, y, z, lattice, bcType); }

        void addTemperatureCornerPNN(
            plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonCorner<1, -1, -1>(x, y, z, lattice, bcType); }

        void addTemperatureCornerPNP(
            plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonCorner<1, -1, 1>(x, y, z, lattice, bcType); }

        void addTemperatureCornerPPN(
            plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonCorner<1, 1, -1>(x, y, z, lattice, bcType); }

        void addTemperatureCornerPPP(
            plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType) { addPoissonCorner<1, 1, 1>(x, y, z, lattice, bcType); };

    private:
        template <int direction, int orientation>
        void addPoissonBoundary(
            Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType);
        template <int plane, int normal1, int normal2>
        void addPoissonEdge(
            Box3D domain, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType);
        template <int normalX, int normalY, int normalZ>
        void addPoissonCorner(
            plint x, plint y, plint z, BlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType);

        template <int direction, int orientation>
        void addPoissonBoundary(
            Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType);
        template <int plane, int normal1, int normal2>
        void addPoissonEdge(
            Box3D domain, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType);
        template <int normalX, int normalY, int normalZ>
        void addPoissonCorner(
            plint x, plint y, plint z, MultiBlockLattice3D<T, Descriptor> &lattice,
            boundary::BcType bcType);
    };

    template <typename T, template <typename U> class Descriptor>
    class PoissonBoundaryManager3D
    {
    public:
        template <int direction, int orientation>
        static BoundaryCompositeDynamics<T, Descriptor> *
        getPoissonBoundaryDynamics(Dynamics<T, Descriptor> *baseDynamics);

        template <int direction, int orientation>
        static BoxProcessingFunctional3D_L<T, Descriptor> *
        getPoissonBoundaryProcessor(Box3D domain);

        template <int plane, int normal1, int normal2>
        static BoundaryCompositeDynamics<T, Descriptor> *
        getPoissonEdgeDynamics(Dynamics<T, Descriptor> *baseDynamics);

        template <int plane, int normal1, int normal2>
        static BoxProcessingFunctional3D_L<T, Descriptor> *
        getPoissonEdgeProcessor(Box3D domain);

        template <int xNormal, int yNormal, int zNormal>
        static BoundaryCompositeDynamics<T, Descriptor> *
        getPoissonCornerDynamics(Dynamics<T, Descriptor> *baseDynamics);
        template <int xNormal, int yNormal, int zNormal>
        static BoxProcessingFunctional3D_L<T, Descriptor> *
        getPoissonCornerProcessor(plint x, plint y, plint z);
    };

    template <typename T, template <typename U> class Descriptor>
    OnLatticeAdvectionDiffusionBoundaryCondition3D<T, Descriptor> *
    createLocalPoissonPoissonLocalBoundaryCondition3D();

}
#endif
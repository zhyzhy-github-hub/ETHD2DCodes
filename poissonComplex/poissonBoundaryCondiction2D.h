#ifndef poissonBoundaryCondiction2D_H
#define poissonBoundaryCondiction2D_H

#include "../../src/complexDynamics/advectionDiffusionBoundaryCondition2D.h"
#include "../../src/complexDynamics/advectionDiffusionBoundaryInstantiator2D.h"

#include "poissonChaiDynamics.h"
#include "poissonChaiBoundaries.h"

namespace plb
{
    template <typename T, template <typename U> class Descriptor>
    class PoissonLocalBoundaryManager2D
    {
    public:
        template <int direction, int orientation>
        static BoundaryCompositeDynamics<T, Descriptor> *
        getPoissonLocalBoundaryDynamics(Dynamics<T, Descriptor> *baseDynamics)
        {
            return new PoissonLocalBoundaryDynamics<T, Descriptor, direction, orientation>(baseDynamics);
        }

        template <int direction, int orientation>
        static BoxProcessingFunctional2D_L<T, Descriptor> *
        getPoissonLocalBoundaryProcessor() {return 0;}

        template <int xNormal, int yNormal>
        static BoundaryCompositeDynamics<T, Descriptor> *
        getPoissonLocalCornerDynamics(Dynamics<T, Descriptor> *baseDynamics)
        {
            return new PoissonLocalCornerDynamics2D<T, Descriptor, xNormal, yNormal>(baseDynamics);
        }

        template <int xNormal, int yNormal>
        static BoxProcessingFunctional2D_L<T, Descriptor> *
        getPoissonLocalCornerProcessor() {return 0;}
    };

    //================ Zou-He like advectionDiffusionBoundaryManager2D ==========//

    template <typename T, template <typename U> class Descriptor, class BoundaryManager>
    class PoissonBoundaryConditionInstantiator2D
        : public AdvectionDiffusionBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>
    {
    };

    template <typename T, template <typename U> class Descriptor, class BoundaryManager>
    class PoissonLocalBoundaryConditionInstantiator2D
        : public OnLatticeAdvectionDiffusionBoundaryCondition2D<T, Descriptor>
    {
    public:
        PoissonLocalBoundaryConditionInstantiator2D(){};

        void addTemperatureBoundary0N(Box2D domain, BlockLattice2D<T, Descriptor> &lattice) { addTemperatureBoundary<0, -1>(domain, lattice); }
        void addTemperatureBoundary0P(Box2D domain, BlockLattice2D<T, Descriptor> &lattice) { addTemperatureBoundary<0, 1>(domain, lattice); }
        void addTemperatureBoundary1N(Box2D domain, BlockLattice2D<T, Descriptor> &lattice) { addTemperatureBoundary<1, -1>(domain, lattice); }
        void addTemperatureBoundary1P(Box2D domain, BlockLattice2D<T, Descriptor> &lattice) { addTemperatureBoundary<1, 1>(domain, lattice); }

        void addTemperatureCornerNN(plint x, plint y, BlockLattice2D<T, Descriptor> &lattice) { PoissonLocalCorner<-1, -1>(x, y, lattice); }
        void addTemperatureCornerNP(plint x, plint y, BlockLattice2D<T, Descriptor> &lattice) { PoissonLocalCorner<-1, 1>(x, y, lattice); }
        void addTemperatureCornerPN(plint x, plint y, BlockLattice2D<T, Descriptor> &lattice) { PoissonLocalCorner<1, -1>(x, y, lattice); }
        void addTemperatureCornerPP(plint x, plint y, BlockLattice2D<T, Descriptor> &lattice) { PoissonLocalCorner<1, 1>(x, y, lattice); }

        void addTemperatureBoundary0N(Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice) { addTemperatureBoundary<0, -1>(domain, lattice); }
        void addTemperatureBoundary0P(Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice) { addTemperatureBoundary<0, 1>(domain, lattice); }
        void addTemperatureBoundary1N(Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice) { addTemperatureBoundary<1, -1>(domain, lattice); }
        void addTemperatureBoundary1P(Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice) { addTemperatureBoundary<1, 1>(domain, lattice); }

        void addTemperatureCornerNN(plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice) { PoissonLocalCorner<-1, -1>(x, y, lattice); }
        void addTemperatureCornerNP(plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice) { PoissonLocalCorner<-1, 1>(x, y, lattice); }
        void addTemperatureCornerPN(plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice) { PoissonLocalCorner<1, -1>(x, y, lattice); }
        void addTemperatureCornerPP(plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice) { PoissonLocalCorner<1, 1>(x, y, lattice); }

    private:
        template <int direction, int orientation>
        void addTemperatureBoundary(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);

        template <int normalX, int normalY>
        void PoissonLocalCorner(plint x, plint y, BlockLattice2D<T, Descriptor> &lattice); // {}

        template <int direction, int orientation>
        void addTemperatureBoundary(Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice);

        template <int normalX, int normalY>
        void PoissonLocalCorner(plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice); //{}
    };

    
    template <typename T, template <typename U> class Descriptor>
    OnLatticeAdvectionDiffusionBoundaryCondition2D<T, Descriptor> *
    createLocalPoissonBoundaryCondition2D();

    template <typename T, template <typename U> class Descriptor>
    OnLatticeAdvectionDiffusionBoundaryCondition2D<T, Descriptor> *
    createLocalPoissonPoissonLocalBoundaryCondition2D();
}

#endif
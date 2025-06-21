#ifndef POISSONCHAI_LATTICES_H
#define POISSONCHAI_LATTICES_H

#include "../../src/core/globalDefs.h"
#include "../../src/latticeBoltzmann/externalFields.h"
#include "../../src/latticeBoltzmann/roundOffPolicy.h"
#include "../../src/latticeBoltzmann/nearestNeighborLattices3D.h"
#include <vector>

namespace plb
{
    namespace descriptors
    {

        struct scalarDescriptor
        {
            static const int numScalars = 1;
            static const int numSpecies = 1;
            static const int scalarBeginsAt = 0;
            static const int sizeOfScalar = 1;
        };

        struct scalarDescriptorBase
        {
            typedef scalarDescriptor ExternalField;
        };
    }

    /// Descriptors for the 2D and 3D lattices.
    /** \warning Attention: The lattice directions must always be ordered in
 * such a way that c[i] = -c[i+(q-1)/2] for i=1..(q-1)/2, and c[0] = 0 must
 * be the rest velocity. Furthermore, the velocities c[i] for i=1..(q-1)/2
 * must verify
 *  - in 2D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *  - in 3D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *                       || (c[i][0]==0 && c[i][1]==0 && c[i][2]<0)
 * Otherwise some of the code will work erroneously, because the
 * aformentioned relations are taken as given to enable a few
 * optimizations.
*/
    namespace descriptors
    {
        /// D2Q5 lattice
        template <typename T>
        struct D2Q5ConstantsPoisson
        {
            enum
            {
                d = 2,
                q = 5
            };                            ///< number of dimensions/distr. functions
            static const T invD;          ///< 1 / (number of dimensions)
            static const int vicinity;    ///< size of neighborhood
            static const int c[q][d];     ///< lattice directions
            static const int cNormSqr[q]; ///< norm-square of the vector c
            static const T t[q];          ///< lattice weights
            static const T cs2;           ///< lattice constant cs2 (in BGK, this is the square-speed-of-sound)
            static const T invCs2;        ///< 1 / cs2
        };

        template <typename T>
        struct D2Q5DescriptorPoissonBase
            : public D2Q5ConstantsPoisson<T>,
              public DefaultRoundOffPolicy<T>
        {
            typedef D2Q5DescriptorPoissonBase<T> BaseDescriptor;
            enum
            {
                numPop = D2Q5ConstantsPoisson<T>::q
            };
        };

        ///  D2Q5 lattice
        template <typename T>
        struct poissonChai08D2Q5Descriptor
            : public D2Q5DescriptorPoissonBase<T>,
              public Velocity2dDescriptorBase
        {
            static const char name[];
        };

        template <typename T>
        struct poissonChai08WithSourceD2Q5Descriptor
            : public D2Q5DescriptorPoissonBase<T>,
              public VelocityAndScalar2dBase
        {
            static const char name[];
        };

        template <typename T>
        struct poissonChai08WithScalarSourceD2Q5Descriptor
            : public D2Q5DescriptorPoissonBase<T>,
              public VelocityAndScalar2dBase
        {
            static const char name[];
        };

    }
    namespace descriptors
    {
        /// D3Q7 lattice
        template <typename T>
        struct D3Q7ConstantsPoisson
        {
            enum
            {
                d = 3,
                q = 7
            };                            ///< number of dimensions/distr. functions
            static const T invD;          ///< 1 / (number of dimensions)
            static const int vicinity;    ///< size of neighborhood
            static const int c[q][d];     ///< lattice directions
            static const int cNormSqr[q]; ///< norm-square of the vector c
            static const T t[q];          ///< lattice weights
            static const T cs2;           ///< lattice constant cs2 (in BGK, this is the square-speed-of-sound)
            static const T invCs2;        ///< 1 / cs2
        };

        template <typename T>
        struct D3Q7PoissonDescriptorBase
            : public D3Q7ConstantsPoisson<T>,
              public DefaultRoundOffPolicy<T>
        {
            typedef D3Q7PoissonDescriptorBase<T> BaseDescriptor;
            enum
            {
                numPop = D3Q7ConstantsPoisson<T>::q
            };
        };

        template <typename T>
        struct PoissonD3Q7Descriptor
            : public D3Q7PoissonDescriptorBase<T>,
              public Velocity3dBase
        {
            static const char name[];
        };

        template <typename T>
        struct PoissonWithSourceD3Q7Descriptor
            : public D3Q7PoissonDescriptorBase<T>,
              public VelocityAndScalar3dBase
        {
            static const char name[];
        };

    }
}
#endif
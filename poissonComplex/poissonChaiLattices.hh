#ifndef POISSONCHAI_LATTICES_HH
#define POISSONCHAI_LATTICES_HH

#include "poissonChaiLattices.h"

namespace plb
{
    namespace descriptors
    {
        // Poisson Chai 2008 AMM D2Q5 //////////////////////////////////////////////

        template <typename T>
        const T D2Q5ConstantsPoisson<T>::invD = (T)1 / (T)d;

        template <typename T>
        const int D2Q5ConstantsPoisson<T>::vicinity = 1;

        template <typename T>
        const int D2Q5ConstantsPoisson<T>::c
            [D2Q5ConstantsPoisson<T>::q][D2Q5ConstantsPoisson<T>::d] =
                {
                    {0, 0},
                    {-1, 0},
                    {0, -1},
                    {1, 0},
                    {0, 1}};

        template <typename T>
        const int D2Q5ConstantsPoisson<T>::cNormSqr[D2Q5ConstantsPoisson<T>::q] =
            {0, 1, 1, 1, 1};

        template <typename T>
        const T D2Q5ConstantsPoisson<T>::t[D2Q5ConstantsPoisson<T>::q] =
            {T(0),
             T(1) / T(4.0), T(1) / T(4.0),
             T(1) / T(4.0), T(1) / T(4.0)};

        template <typename T>
        const T D2Q5ConstantsPoisson<T>::cs2 = (T)1. / (T)2;

        template <typename T>
        const T D2Q5ConstantsPoisson<T>::invCs2 = (T)2;

        template <typename T>
        const char poissonChai08D2Q5Descriptor<T>::name[] = "poissonChai08D2Q5";

        template <typename T>
        const char poissonChai08WithSourceD2Q5Descriptor<T>::name[] = "poissonChai08WithSourceD2Q5";

        template <typename T>
        const char poissonChai08WithScalarSourceD2Q5Descriptor<T>::name[] = "poissonChai08WithScalarSourceD2Q5";
    }

    namespace descriptors
    {
        // AdvectionDiffusion D3Q7 ////////////////////////////////////////////////////

        template <typename T>
        const T D3Q7ConstantsPoisson<T>::invD = (T)1 / (T)d;

        template <typename T>
        const int D3Q7ConstantsPoisson<T>::vicinity = 1;

        template <typename T>
        const int D3Q7ConstantsPoisson<T>::c
            [D3Q7ConstantsPoisson<T>::q][D3Q7ConstantsPoisson<T>::d] =
                {
                    {0, 0, 0},
                    {-1, 0, 0},
                    {0, -1, 0},
                    {0, 0, -1},
                    {1, 0, 0},
                    {0, 1, 0},
                    {0, 0, 1},
        };

        template <typename T>
        const int D3Q7ConstantsPoisson<T>::cNormSqr[D3Q7ConstantsPoisson<T>::q] =
            {0, 1, 1, 1, 1, 1, 1};

        template <typename T>
        const T D3Q7ConstantsPoisson<T>::t[D3Q7ConstantsPoisson<T>::q] =
            {T(0),
             T(1) / T(6.0), T(1) / T(6.0), T(1) / T(6.0),
             T(1) / T(6.0), T(1) / T(6.0), T(1) / T(6.0)};
        // {(T)1 - (T)3 / invCs2,
        //  (T)1 / (invCs2 * (T)2), (T)1 / (invCs2 * (T)2), (T)1 / (invCs2 * (T)2),
        //  (T)1 / (invCs2 * (T)2), (T)1 / (invCs2 * (T)2), (T)1 / (invCs2 * (T)2)};

        template <typename T>
        // const T D3Q7ConstantsPoisson<T>::cs2 = (T)1.0 / (T)4;
        const T D3Q7ConstantsPoisson<T>::cs2 = (T)1.0 / (T)3;

        template <typename T>
        // const T D3Q7ConstantsPoisson<T>::invCs2 = (T)4;
        const T D3Q7ConstantsPoisson<T>::invCs2 = (T)3;

        template <typename T>
        const char PoissonD3Q7Descriptor<T>::name[] = "PoissonD3Q7Descriptor";

        template <typename T>
        const char PoissonWithSourceD3Q7Descriptor<T>::name[] = "PoissonWithSourceD3Q7Descriptor";
    }
}

#endif
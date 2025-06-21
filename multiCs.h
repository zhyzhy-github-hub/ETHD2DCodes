#ifndef MULTI_CS_H
#define MULTI_CS_H

#include "../../src/core/units.h"
#include "../../src/core/array.h"
#include "../../src/core/globalDefs.h"
#include "../../src/core/dynamics.h"
#include "../../src/latticeBoltzmann/nearestNeighborLattices2D.h"
#include "../../src/latticeBoltzmann/mrtLattices.h"
#include "../../src/latticeBoltzmann/advectionDiffusionLattices.h"

#include <string>
#include <fstream>

#define ccDefine 5

#define lesModel

// typedef double T;

namespace plb
{
    namespace descriptors
    {
        struct advSourceVelocityTauD2dDescriptor
        {
            static const int numScalars = 3;
            static const int numSpecies = 2;
            static const int velocityBeginsAt = 0;
            static const int sizeOfVelocity = 2;
            static const int tauDBeginsAt = 2;
            static const int sizeOftauD = 1;
        };
        struct advSourceVelocityTauD2dDescriptorBase
        {
            typedef advSourceVelocityTauD2dDescriptor ExternalField;
        };

        /// AD D2Q5 lattice
        template <typename T>
        struct AdvectionDiffusionMultiLesD2Q5Descriptor
            : public D2Q5DescriptorBase<T>,
              public advSourceVelocityTauD2dDescriptorBase
        {
            static const char name[];
        };
        template <typename T>
        const char AdvectionDiffusionMultiLesD2Q5Descriptor<T>::name[] = "AdvectionDiffusionMultiLesD2Q5";
        // ------------------------------------------------//

        struct ForceVelocity2dDescriptor
        {
            static const int numScalars = 4;
            static const int numSpecies = 2;
            static const int velocityBeginsAt = 0;
            static const int sizeOfVelocity = 2;
            static const int forceBeginsAt = 2;
            static const int sizeOfForce = 2;
        };

        struct ForceVelocity2dDescriptorBase
        {
            typedef ForceVelocity2dDescriptor ExternalField;
        };
        // ------------------------------------------------//

        template <typename T>
        struct ForceVelocitydD2Q9Descriptor
            : public D2Q9DescriptorBase<T>,
              public ForceVelocity2dDescriptorBase
        {
            static const char name[];
        };

        template <typename T>
        const char ForceVelocitydD2Q9Descriptor<T>::name[] = "ForceVelocitydD2Q9Descriptor";

        template <typename T>
        struct ForceVelocitydMRTD2Q9Descriptor
            : public MRTD2Q9DescriptorBase<T>,
              public ForceVelocity2dDescriptorBase
        {
            static const char name[];
        };

        template <typename T>
        const char ForceVelocitydMRTD2Q9Descriptor<T>::name[] = "ForceVelocitydMRTD2Q9Descriptor";

        struct ForceVelocity2dTauLesDescriptor
        {
            static const int numScalars = 5;
            static const int numSpecies = 3;
            static const int velocityBeginsAt = 0;
            static const int sizeOfVelocity = 2;
            static const int forceBeginsAt = 2;
            static const int sizeOfForce = 2;
            static const int tauTBeginsAt = 4;
            static const int sizeOftauT = 1;
        };

        struct ForceVelocity2dTauLesDescriptorBase
        {
            typedef ForceVelocity2dTauLesDescriptor ExternalField;
        };
        // ------------------------------------------------//

        template <typename T>
        struct ForceVelocitydD2Q9TauLesDescriptor
            : public D2Q9DescriptorBase<T>,
              public ForceVelocity2dTauLesDescriptorBase
        {
            static const char name[];
        };

        template <typename T>
        const char ForceVelocitydD2Q9TauLesDescriptor<T>::name[] = "ForceVelocitydD2Q9TauLesDescriptor";

        struct ForceVelocity2dForceTauLesDescriptor
        {
            static const int numScalars = 6;
            static const int numSpecies = 4;
            static const int velocityBeginsAt = 0;
            static const int sizeOfVelocity = 2;
            static const int forceBeginsAt = 2;
            static const int sizeOfForce = 2;
            static const int tauTBeginsAt = 4;
            static const int sizeOftauT = 1;
            static const int forceTauBeginsAt = 5;
            static const int sizeOfForceTau = 1;
        };

        struct ForceVelocity2dForceTauLesDescriptorBase
        {
            typedef ForceVelocity2dForceTauLesDescriptor ExternalField;
        };
        // ------------------------------------------------//

        template <typename T>
        struct ForceVelocitydD2Q9ForceTauLesDescriptor
            : public D2Q9DescriptorBase<T>,
              public ForceVelocity2dForceTauLesDescriptorBase
        {
            static const char name[];
        };

        template <typename T>
        const char ForceVelocitydD2Q9ForceTauLesDescriptor<T>::name[] = "ForceVelocitydD2Q9ForceTauLesDescriptor";
    }

    // ----------------------------------------------------------------------------------------

    template <typename T, template <typename U> class Descriptor>
    class BoxTauLesFunctional2D : public BoxProcessingFunctional2D_LS<T, Descriptor, T>
    {
    public:
        virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
                             ScalarField2D<T> &scalarField)
        {
            Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
            for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            {
                for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                {
                    scalarField.get(iX + offset.x, iY + offset.y) = *(lattice.get(iX, iY).getExternal(Descriptor<T>::ExternalField::tauTBeginsAt));
                }
            }
        }
        virtual BoxTauLesFunctional2D<T, Descriptor> *clone() const
        {
            return new BoxTauLesFunctional2D<T, Descriptor>(*this);
        }
        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
        {
            modified[0] = modif::nothing;
            modified[1] = modif::staticVariables;
        }
        virtual BlockDomain::DomainT appliesTo() const
        {
            return BlockDomain::bulk;
        }
    };

    template <typename T, template <typename U> class Descriptor>
    void computeTauLes(MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &density, Box2D domain)
    {
        applyProcessingFunctional(
            new BoxTauLesFunctional2D<T, Descriptor>, domain, lattice, density);
    }

    template <typename T, template <typename U> class Descriptor>
    std::unique_ptr<MultiScalarField2D<T>> computeTauLes(MultiBlockLattice2D<T, Descriptor> &lattice, Box2D domain)
    {
        std::unique_ptr<MultiScalarField2D<T>> density =
            generateMultiScalarField<T>(lattice, domain);

        computeTauLes(lattice, *density, domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
        return density;
    }

    template <typename T, template <typename U> class Descriptor>
    std::unique_ptr<MultiScalarField2D<T>> computeTauLes(MultiBlockLattice2D<T, Descriptor> &lattice)
    {
        return computeTauLes(lattice, lattice.getBoundingBox());
    }

    // ----------------------------------------------------------------------------------------
    /// Implementation of the MRT collision step
    template <typename T, template <typename U> class Descriptor>
    class MyGuoD2Q9ExternalForceMRTdynamics01MultiCs : public ExternalForceDynamics<T, Descriptor>
    {
    public:
        /* *************** Construction / Destruction ************************ */
        MyGuoD2Q9ExternalForceMRTdynamics01MultiCs(T omega_)
            : ExternalForceDynamics<T, Descriptor>(omega_) //, cCs(1.0)
        {
        }
        MyGuoD2Q9ExternalForceMRTdynamics01MultiCs(T omega_, T cCs_)
            : ExternalForceDynamics<T, Descriptor>(omega_), cCs(cCs_)
        {
        }

        MyGuoD2Q9ExternalForceMRTdynamics01MultiCs(HierarchicUnserializer &unserializer)
            : ExternalForceDynamics<T, Descriptor>(T())
        {
            this->unserialize(unserializer);
        }

        /// Clone the object on its dynamic type.
        virtual MyGuoD2Q9ExternalForceMRTdynamics01MultiCs<T, Descriptor> *clone() const
        {
            return new MyGuoD2Q9ExternalForceMRTdynamics01MultiCs<T, Descriptor>(*this);
        }

        /// Return a unique ID for this class.
        virtual int getId() const
        {
            return id;
        }

        /* *************** Collision and Equilibrium ************************* */

        virtual void collide(Cell<T, Descriptor> &cell,
                             BlockStatistics &statistics_)
        {
            T rhoBar = momentTemplates<T, Descriptor>::get_rhoBar(cell);
            const T rhoFull = Descriptor<T>::fullRho(rhoBar);
            Array<T, Descriptor<T>::d> u;
            // this->computeVelocity(cell, u);

            Array<T, Descriptor<T>::d> j;
            momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
            // force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));
            for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
            {
                u[iD] = j[iD] / rhoFull * cCs; // + force[iD] / cCs * 0.5 * invRho;
            }

            Array<T, Descriptor<T>::d> force;
            force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));

            const T rho = rhoBar;
            T ux = u[0];
            T uy = u[1];
            T fx = force[0];
            T fy = force[1];

            T *vel = cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt);

            const T f0 = cell[0];
            const T f1 = cell[1];
            const T f2 = cell[2];
            const T f3 = cell[3];
            const T f4 = cell[4];
            const T f5 = cell[5];
            const T f6 = cell[6];
            const T f7 = cell[7];
            const T f8 = cell[8];

            const T m0 = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
            const T m1 = -4 * f0 + 2 * f1 - f2 + 2 * f3 - f4 + 2 * f5 - f6 + 2 * f7 - f8;
            const T m2 = 4 * f0 + f1 - 2 * f2 + f3 - 2 * f4 + f5 - 2 * f6 + f7 - 2 * f8;
            const T m3 = -f1 - f2 - f3 + f5 + f6 + f7;
            const T m4 = -f1 + 2 * f2 - f3 + f5 - 2 * f6 + f7;
            const T m5 = f1 - f3 - f4 - f5 + f7 + f8;
            const T m6 = f1 - f3 + 2 * f4 - f5 + f7 - 2 * f8;
            const T m7 = f2 - f4 + f6 - f8;
            const T m8 = -f1 + f3 - f5 + f7;

            const T meq0 = 1.0 * rhoBar;
            const T meq1 = (3.0 * rhoFull * ux * ux + 3.0 * rhoFull * uy * uy - 2.0 * rhoBar * cCs * cCs) / (cCs * cCs);
            const T meq2 = (-3.0 * rhoFull * ux * ux - 3.0 * rhoFull * uy * uy + rhoBar * cCs * cCs) / (cCs * cCs);
            const T meq3 = 1.0 * rhoFull * ux / cCs;
            const T meq4 = -1.0 * rhoFull * ux / cCs;
            const T meq5 = 1.0 * rhoFull * uy / cCs;
            const T meq6 = -1.0 * rhoFull * uy / cCs;
            const T meq7 = 1.0 * rhoFull * (ux * ux - uy * uy) / (cCs * cCs);
            const T meq8 = 1.0 * rhoFull * ux * uy / (cCs * cCs);

            const T mF0 = 0;
            const T mF1 = (6.0 * fx * ux + 6.0 * fy * uy) / (cCs * cCs);
            const T mF2 = (-6.0 * fx * ux - 6.0 * fy * uy) / (cCs * cCs);
            const T mF3 = 1.0 * fx / cCs;
            const T mF4 = -1.0 * fx / cCs;
            const T mF5 = 1.0 * fy / cCs;
            const T mF6 = -1.0 * fy / cCs;
            const T mF7 = (2.0 * fx * ux - 2.0 * fy * uy) / (cCs * cCs);
            const T mF8 = (1.0 * fx * uy + 1.0 * fy * ux) / (cCs * cCs);

            const T omega = ExternalForceDynamics<T, Descriptor>::getOmega();
            const T tau = 1. / omega;
            const T omg1 = 8 * (2 * tau - 1) / (8 * tau - 1);

            // pcout << " c = " << cCs << "\n";

            T s0 = 1;

            T s3 = 1;
            T s5 = 1;
            T s1 = 0.2;

            T s2 = 0.2;
            T s4 = 1.2;
            T s6 = 1.2;
            T s7 = omega;
            T s8 = omega;

            T mcp0 = m0 * (1.0 - s0) + s0 * meq0; //+ (1 - s0 / 2) * mF0;
            T mcp1 = m1 * (1.0 - s1) + s1 * meq1; //+ (1 - s1 / 2) * mF1;
            T mcp2 = m2 * (1.0 - s2) + s2 * meq2; //+ (1 - s2 / 2) * mF2;
            T mcp3 = m3 * (1.0 - s3) + s3 * meq3 + fx / cCs / cCs;
            T mcp4 = m4 * (1.0 - s4) + s4 * meq4; //+ (1 - s4 / 2) * mF4;
            T mcp5 = m5 * (1.0 - s5) + s5 * meq5 + fy / cCs / cCs;
            T mcp6 = m6 * (1.0 - s6) + s6 * meq6; // + (1 - s6 / 2) * mF6;
            T mcp7 = m7 * (1.0 - s7) + s7 * meq7; // + (1 - s7 / 2) * mF7;
            T mcp8 = m8 * (1.0 - s8) + s8 * meq8; // + (1 - s8 / 2) * mF8;

            cell[0] = mcp0 / 9 - mcp1 / 9 + mcp2 / 9;
            cell[1] = mcp0 / 9 + mcp1 / 18 + mcp2 / 36 - mcp3 / 6 - mcp4 / 12 + mcp5 / 6 + mcp6 / 12 - mcp8 / 4;
            cell[2] = mcp0 / 9 - mcp1 / 36 - mcp2 / 18 - mcp3 / 6 + mcp4 / 6 + mcp7 / 4;
            cell[3] = mcp0 / 9 + mcp1 / 18 + mcp2 / 36 - mcp3 / 6 - mcp4 / 12 - mcp5 / 6 - mcp6 / 12 + mcp8 / 4;
            cell[4] = mcp0 / 9 - mcp1 / 36 - mcp2 / 18 - mcp5 / 6 + mcp6 / 6 - mcp7 / 4;
            cell[5] = mcp0 / 9 + mcp1 / 18 + mcp2 / 36 + mcp3 / 6 + mcp4 / 12 - mcp5 / 6 - mcp6 / 12 - mcp8 / 4;
            cell[6] = mcp0 / 9 - mcp1 / 36 - mcp2 / 18 + mcp3 / 6 - mcp4 / 6 + mcp7 / 4;
            cell[7] = mcp0 / 9 + mcp1 / 18 + mcp2 / 36 + mcp3 / 6 + mcp4 / 12 + mcp5 / 6 + mcp6 / 12 + mcp8 / 4;
            cell[8] = mcp0 / 9 - mcp1 / 36 - mcp2 / 18 + mcp5 / 6 - mcp6 / 6 - mcp7 / 4;

            ux += fx / 2 / cCs / rhoFull;
            uy += fy / 2 / cCs / rhoFull;
            vel[0] = ux;
            vel[1] = uy;
        }

        virtual void computeVelocity(
            Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const
        {

            const T *vel = cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt);
            // pcout << "hello" << std::endl;
            u[0] = vel[0];
            u[1] = vel[1];
        }

        virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j,
                                     T jSqr, T thetaBar = T()) const
        {
            // pcout << "hello equ" << std::endl;
            T invRho = Descriptor<T>::invRho(rhoBar);
            T invCs2_01 = Descriptor<T>::invCs2 / (cCs * cCs);

            T c_j = Descriptor<T>::c[iPop][0] * j[0] * cCs;
            for (int iD = 1; iD < Descriptor<T>::d; ++iD)
            {
                c_j += Descriptor<T>::c[iPop][iD] * j[iD] * cCs;
            }
            return Descriptor<T>::t[iPop] * (rhoBar + invCs2_01 * c_j +
                                             invCs2_01 / (T)2 * invRho * (invCs2_01 * c_j * c_j - jSqr));
        }

        T getCCs() const { return cCs; }

    private:
        static int id;
        const T cCs = ccDefine;
    };

    template <typename T, template <typename U> class Descriptor>
    int MyGuoD2Q9ExternalForceMRTdynamics01MultiCs<T, Descriptor>::id =
        meta::registerGeneralDynamics<T, Descriptor, MyGuoD2Q9ExternalForceMRTdynamics01MultiCs<T, Descriptor>>("MyGuoD2Q9ExternalForceMRTdynamics01MultiCs");

    /// Implementation of the MRT collision step
    template <typename T, template <typename U> class Descriptor>
    class MyGuoD2Q9ExternalForceLesMRTdynamics01MultiCs : public ExternalForceDynamics<T, Descriptor>
    {
    public:
        /* *************** Construction / Destruction ************************ */
        MyGuoD2Q9ExternalForceLesMRTdynamics01MultiCs(T omega_)
            : ExternalForceDynamics<T, Descriptor>(omega_) //, cCs(1.0)
        {
        }
        MyGuoD2Q9ExternalForceLesMRTdynamics01MultiCs(T omega_, T cCs_)
            : ExternalForceDynamics<T, Descriptor>(omega_), cCs(cCs_)
        {
        }

        MyGuoD2Q9ExternalForceLesMRTdynamics01MultiCs(HierarchicUnserializer &unserializer)
            : ExternalForceDynamics<T, Descriptor>(T())
        {
            this->unserialize(unserializer);
        }

        /// Clone the object on its dynamic type.
        virtual MyGuoD2Q9ExternalForceLesMRTdynamics01MultiCs<T, Descriptor> *clone() const
        {
            return new MyGuoD2Q9ExternalForceLesMRTdynamics01MultiCs<T, Descriptor>(*this);
        }

        /// Return a unique ID for this class.
        virtual int getId() const
        {
            return id;
        }

        virtual T computeTemperature(Cell<T, Descriptor> const &cell) const
        {
            const T *s8_0 = cell.getExternal(Descriptor<T>::ExternalField::tauTBeginsAt);
            // return s8_0[0];

            const T omega = ExternalForceDynamics<T, Descriptor>::getOmega();
            // return Descriptor<T>::cs2 * (1. / omega - 0.5) + *s8_0;

            return *s8_0;
        }

        void compute_PiNeq_This(Cell<T, Descriptor> const &f, T rhoBar,
                                Array<T, Descriptor<T>::d> const &j, T jSqr,
                                Array<T, SymmetricTensorImpl<T, Descriptor<T>::d>::n> &PiNeq)
        {
            int iPi = 0;
            for (int iAlpha = 0; iAlpha < Descriptor<T>::d; ++iAlpha)
            {
                int iDiagonal = iPi;
                for (int iBeta = iAlpha; iBeta < Descriptor<T>::d; ++iBeta)
                {
                    T fneq0 = f[0] - computeEquilibrium(0, rhoBar, j, jSqr);
                    PiNeq[iPi] = Descriptor<T>::c[0][iAlpha] * Descriptor<T>::c[0][iBeta] * fneq0;

                    for (plint iPop = 1; iPop < Descriptor<T>::q; ++iPop)
                    {
                        T fneqi = f[iPop] - computeEquilibrium(iPop, rhoBar, j, jSqr);
                        PiNeq[iPi] += Descriptor<T>::c[iPop][iAlpha] * Descriptor<T>::c[iPop][iBeta] * fneqi;
                    }
                    PiNeq[iPi] *= (cCs * cCs);

                    ++iPi;
                }
            }
        }

        virtual void collide(Cell<T, Descriptor> &cell,
                             BlockStatistics &statistics_)
        {
            // ---------------------------------------------------------计算密度
            T rhoBar = momentTemplates<T, Descriptor>::get_rhoBar(cell);
            const T rhoFull = Descriptor<T>::fullRho(rhoBar);

            // ---------------------------------------------------------计算速度
            Array<T, Descriptor<T>::d> u;
            Array<T, Descriptor<T>::d> j;
            momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
            for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
            {
                u[iD] = j[iD] / rhoFull * cCs; // + force[iD] / cCs * 0.5 * invRho;
            }

            // ---------------------------------------------------------计算外力
            Array<T, Descriptor<T>::d> force;
            force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));

            // ---------------------------------------------------------计算速度和外力
            T ux = u[0];
            T uy = u[1];
            T fx = force[0];
            T fy = force[1];

            T *vel = cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt);

            // ---------------------------------------------------------碰撞
            const T f0 = cell[0];
            const T f1 = cell[1];
            const T f2 = cell[2];
            const T f3 = cell[3];
            const T f4 = cell[4];
            const T f5 = cell[5];
            const T f6 = cell[6];
            const T f7 = cell[7];
            const T f8 = cell[8];

            const T m0 = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
            const T m1 = -4 * f0 + 2 * f1 - f2 + 2 * f3 - f4 + 2 * f5 - f6 + 2 * f7 - f8;
            const T m2 = 4 * f0 + f1 - 2 * f2 + f3 - 2 * f4 + f5 - 2 * f6 + f7 - 2 * f8;
            const T m3 = -f1 - f2 - f3 + f5 + f6 + f7;
            const T m4 = -f1 + 2 * f2 - f3 + f5 - 2 * f6 + f7;
            const T m5 = f1 - f3 - f4 - f5 + f7 + f8;
            const T m6 = f1 - f3 + 2 * f4 - f5 + f7 - 2 * f8;
            const T m7 = f2 - f4 + f6 - f8;
            const T m8 = -f1 + f3 - f5 + f7;

            // ---------------------------------------------------------计算平衡态矩
            const T meq0 = 1.0 * rhoBar;
            const T meq1 = (3.0 * rhoFull * ux * ux + 3.0 * rhoFull * uy * uy - 2.0 * rhoBar * cCs * cCs) / (cCs * cCs);
            const T meq2 = (-3.0 * rhoFull * ux * ux - 3.0 * rhoFull * uy * uy + rhoBar * cCs * cCs) / (cCs * cCs);
            const T meq3 = 1.0 * rhoFull * ux / cCs;
            const T meq4 = -1.0 * rhoFull * ux / cCs;
            const T meq5 = 1.0 * rhoFull * uy / cCs;
            const T meq6 = -1.0 * rhoFull * uy / cCs;
            const T meq7 = 1.0 * rhoFull * (ux * ux - uy * uy) / (cCs * cCs);
            const T meq8 = 1.0 * rhoFull * ux * uy / (cCs * cCs);

            const T omega = ExternalForceDynamics<T, Descriptor>::getOmega();
            const T tau = 1. / omega;

            // ---------------------------------------------------------取出湍涡粘度相关的粘性
            T *nuTur0 = cell.getExternal(Descriptor<T>::ExternalField::tauTBeginsAt);
            T nuTur = nuTur0[0];
            // T tauLast = tau + nuTur / (Descriptor<T>::cs2);
            T tauLast = tau + nuTur / (Descriptor<T>::cs2 * cCs);
            T s8Last = 1. / tauLast;

            T s0 = 1;
            T s1 = 1.63;
            T s2 = 1.14;
            T s3 = 1;
            T s4 = 1.92;
            T s5 = 1;
            T s6 = 1.92;

            T s7 = s8Last;
            T s8 = s8Last;

#ifdef lesModel
            // 这是之前幂律流体做的东西
            Array<T, SymmetricTensor<T, Descriptor>::n> stress;

            compute_PiNeq_This(cell, rhoBar, j, 0, stress);
            T D2 = (pow(stress[0], 2) + pow(stress[2], 2) + 2 * pow(stress[1], 2) + 1e-20);

            T para = 1. / (2 * rhoFull * Descriptor<T>::cs2 * cCs * tauLast);

            T gamma = sqrt(2 * D2 * para * para);
            T SsSqrt = gamma;

            T nuT = pow((0.17 * 1.0), 2) * SsSqrt;
            // T nuT = pow((0.17 * 1.414213562), 1) * SsSqrt;
            // T tauTM05 = nuT / (Descriptor<T>::cs2);
            T tauTM05 = nuT / (Descriptor<T>::cs2 * cCs);

            T tau_u = tauTM05 + tau;
            nuTur0[0] = nuT;
#else
            T tau_u = tau;
            nuTur0[0] = 0;
            // pcout << "------------" << std::endl;
#endif

            s7 = 1. / tau_u;
            s8 = 1. / tau_u;

            // sMagin = 1.0 / (1.0 / (12 * (tau_u - 0.5)) + 0.5);
            // s1 = sMagin;
            // s2 = sMagin;
            // s4 = sMagin;
            // s6 = sMagin;

            // if (1. / s8 <= 0.5 || 1. / s7 <= 0.5)
            // {
            //     pcout << "error S7 s8 022 ~~~!!!!!!!!!!!\n";
            // }

            T mcp0 = m0 * (1.0 - s0) + s0 * meq0; //+ (1 - s0 / 2) * mF0;
            T mcp1 = m1 * (1.0 - s1) + s1 * meq1; //+ (1 - s1 / 2) * mF1;
            T mcp2 = m2 * (1.0 - s2) + s2 * meq2; //+ (1 - s2 / 2) * mF2;
            T mcp3 = m3 * (1.0 - s3) + s3 * meq3 + fx / cCs / cCs;
            T mcp4 = m4 * (1.0 - s4) + s4 * meq4; //+ (1 - s4 / 2) * mF4;
            T mcp5 = m5 * (1.0 - s5) + s5 * meq5 + fy / cCs / cCs;
            T mcp6 = m6 * (1.0 - s6) + s6 * meq6; // + (1 - s6 / 2) * mF6;
            T mcp7 = m7 * (1.0 - s7) + s7 * meq7; // + (1 - s7 / 2) * mF7;
            T mcp8 = m8 * (1.0 - s8) + s8 * meq8; // + (1 - s8 / 2) * mF8;

            cell[0] = mcp0 / 9 - mcp1 / 9 + mcp2 / 9;
            cell[1] = mcp0 / 9 + mcp1 / 18 + mcp2 / 36 - mcp3 / 6 - mcp4 / 12 + mcp5 / 6 + mcp6 / 12 - mcp8 / 4;
            cell[2] = mcp0 / 9 - mcp1 / 36 - mcp2 / 18 - mcp3 / 6 + mcp4 / 6 + mcp7 / 4;
            cell[3] = mcp0 / 9 + mcp1 / 18 + mcp2 / 36 - mcp3 / 6 - mcp4 / 12 - mcp5 / 6 - mcp6 / 12 + mcp8 / 4;
            cell[4] = mcp0 / 9 - mcp1 / 36 - mcp2 / 18 - mcp5 / 6 + mcp6 / 6 - mcp7 / 4;
            cell[5] = mcp0 / 9 + mcp1 / 18 + mcp2 / 36 + mcp3 / 6 + mcp4 / 12 - mcp5 / 6 - mcp6 / 12 - mcp8 / 4;
            cell[6] = mcp0 / 9 - mcp1 / 36 - mcp2 / 18 + mcp3 / 6 - mcp4 / 6 + mcp7 / 4;
            cell[7] = mcp0 / 9 + mcp1 / 18 + mcp2 / 36 + mcp3 / 6 + mcp4 / 12 + mcp5 / 6 + mcp6 / 12 + mcp8 / 4;
            cell[8] = mcp0 / 9 - mcp1 / 36 - mcp2 / 18 + mcp5 / 6 - mcp6 / 6 - mcp7 / 4;

            ux += fx / 2 / cCs / rhoFull;
            uy += fy / 2 / cCs / rhoFull;
            vel[0] = ux;
            vel[1] = uy;
        }

        virtual void computeVelocity(
            Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const
        {

            const T *vel = cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt);
            // pcout << "hello" << std::endl;
            u[0] = vel[0];
            u[1] = vel[1];
        }

        virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j,
                                     T jSqr, T thetaBar = T()) const
        {
            // pcout << "hello equ" << std::endl;
            T invRho = Descriptor<T>::invRho(rhoBar);
            T invCs2_01 = Descriptor<T>::invCs2 / (cCs * cCs);

            T c_j = Descriptor<T>::c[iPop][0] * j[0] * cCs;
            for (int iD = 1; iD < Descriptor<T>::d; ++iD)
            {
                c_j += Descriptor<T>::c[iPop][iD] * j[iD] * cCs;
            }
            return Descriptor<T>::t[iPop] * (rhoBar + invCs2_01 * c_j +
                                             invCs2_01 / (T)2 * invRho * (invCs2_01 * c_j * c_j - jSqr));
        }

        T getCCs() const { return cCs; }

    private:
        static int id;
        const T cCs = ccDefine;
    };
    template <typename T, template <typename U> class Descriptor>
    int MyGuoD2Q9ExternalForceLesMRTdynamics01MultiCs<T, Descriptor>::id =
        meta::registerGeneralDynamics<T, Descriptor, MyGuoD2Q9ExternalForceLesMRTdynamics01MultiCs<T, Descriptor>>("MyGuoD2Q9ExternalForceLesMRTdynamics01MultiCs");

    /// Implementation of the MRT collision step
    template <typename T, template <typename U> class Descriptor>
    class MyGuoD2Q9ExternalForceLesMRTdynamics02_direct_MultiCs : public ExternalForceDynamics<T, Descriptor>
    {
    public:
        /* *************** Construction / Destruction ************************ */
        MyGuoD2Q9ExternalForceLesMRTdynamics02_direct_MultiCs(T omega_)
            : ExternalForceDynamics<T, Descriptor>(omega_) //, cCs(1.0)
        {
        }
        MyGuoD2Q9ExternalForceLesMRTdynamics02_direct_MultiCs(T omega_, T cCs_)
            : ExternalForceDynamics<T, Descriptor>(omega_), cCs(cCs_)
        {
        }

        MyGuoD2Q9ExternalForceLesMRTdynamics02_direct_MultiCs(HierarchicUnserializer &unserializer)
            : ExternalForceDynamics<T, Descriptor>(T())
        {
            this->unserialize(unserializer);
        }

        /// Clone the object on its dynamic type.
        virtual MyGuoD2Q9ExternalForceLesMRTdynamics02_direct_MultiCs<T, Descriptor> *clone() const
        {
            return new MyGuoD2Q9ExternalForceLesMRTdynamics02_direct_MultiCs<T, Descriptor>(*this);
        }

        /// Return a unique ID for this class.
        virtual int getId() const
        {
            return id;
        }

        virtual T computeTemperature(Cell<T, Descriptor> const &cell) const
        {
            const T *s8_0 = cell.getExternal(Descriptor<T>::ExternalField::tauTBeginsAt);
            // return s8_0[0];

            const T omega = ExternalForceDynamics<T, Descriptor>::getOmega();
            // return Descriptor<T>::cs2 * (1. / omega - 0.5) + *s8_0;

            return *s8_0;
        }

        void compute_PiNeq_This(Cell<T, Descriptor> const &f, T rhoBar,
                                Array<T, Descriptor<T>::d> const &j, T jSqr,
                                Array<T, SymmetricTensorImpl<T, Descriptor<T>::d>::n> &PiNeq)
        {
            int iPi = 0;
            for (int iAlpha = 0; iAlpha < Descriptor<T>::d; ++iAlpha)
            {
                int iDiagonal = iPi;
                for (int iBeta = iAlpha; iBeta < Descriptor<T>::d; ++iBeta)
                {
                    T fneq0 = f[0] - computeEquilibrium(0, rhoBar, j, jSqr);
                    PiNeq[iPi] = Descriptor<T>::c[0][iAlpha] * Descriptor<T>::c[0][iBeta] * fneq0;

                    for (plint iPop = 1; iPop < Descriptor<T>::q; ++iPop)
                    {
                        T fneqi = f[iPop] - computeEquilibrium(iPop, rhoBar, j, jSqr);
                        PiNeq[iPi] += Descriptor<T>::c[iPop][iAlpha] * Descriptor<T>::c[iPop][iBeta] * fneqi;
                    }
                    PiNeq[iPi] *= (cCs * cCs);

                    ++iPi;
                }
            }
        }

        virtual void collide(Cell<T, Descriptor> &cell,
                             BlockStatistics &statistics_)
        {
            // ---------------------------------------------------------计算密度
            T rhoBar = momentTemplates<T, Descriptor>::get_rhoBar(cell);
            const T rhoFull = Descriptor<T>::fullRho(rhoBar);

            // ---------------------------------compute velocity
            Array<T, Descriptor<T>::d> u;
            Array<T, Descriptor<T>::d> j;
            momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
            for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
            {
                u[iD] = j[iD] / rhoFull * cCs; // + force[iD] / cCs * 0.5 * invRho;
            }

            // --------------------------------external force
            Array<T, Descriptor<T>::d> force;
            force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));

            // -------------------------------------velocity and force
            T ux = u[0];
            T uy = u[1];
            T fx = force[0];
            T fy = force[1];

            T *vel = cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt);

            // -------------------------------------------collision
            const T f0 = cell[0];
            const T f1 = cell[1];
            const T f2 = cell[2];
            const T f3 = cell[3];
            const T f4 = cell[4];
            const T f5 = cell[5];
            const T f6 = cell[6];
            const T f7 = cell[7];
            const T f8 = cell[8];

            const T m0 = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
            const T m1 = -4 * f0 + 2 * f1 - f2 + 2 * f3 - f4 + 2 * f5 - f6 + 2 * f7 - f8;
            const T m2 = 4 * f0 + f1 - 2 * f2 + f3 - 2 * f4 + f5 - 2 * f6 + f7 - 2 * f8;
            const T m3 = -f1 - f2 - f3 + f5 + f6 + f7;
            const T m4 = -f1 + 2 * f2 - f3 + f5 - 2 * f6 + f7;
            const T m5 = f1 - f3 - f4 - f5 + f7 + f8;
            const T m6 = f1 - f3 + 2 * f4 - f5 + f7 - 2 * f8;
            const T m7 = f2 - f4 + f6 - f8;
            const T m8 = -f1 + f3 - f5 + f7;

            // ----------------------------------Equilibrium moment-----------
            const T meq0 = 1.0 * rhoBar;
            const T meq1 = (3.0 * rhoFull * ux * ux + 3.0 * rhoFull * uy * uy - 2.0 * rhoBar * cCs * cCs) / (cCs * cCs);
            const T meq2 = (-3.0 * rhoFull * ux * ux - 3.0 * rhoFull * uy * uy + rhoBar * cCs * cCs) / (cCs * cCs);
            const T meq3 = 1.0 * rhoFull * ux / cCs;
            const T meq4 = -1.0 * rhoFull * ux / cCs;
            const T meq5 = 1.0 * rhoFull * uy / cCs;
            const T meq6 = -1.0 * rhoFull * uy / cCs;
            const T meq7 = 1.0 * rhoFull * (ux * ux - uy * uy) / (cCs * cCs);
            const T meq8 = 1.0 * rhoFull * ux * uy / (cCs * cCs);

            const T omega = ExternalForceDynamics<T, Descriptor>::getOmega();
            const T tau = 1. / omega;

            T s0 = 1;
            T s1 = 1.63;
            T s2 = 1.14;
            T s3 = 1;
            T s4 = 1.92;
            T s5 = 1;
            T s6 = 1.92;

            // --------------------viscisity of LES-------------------------------------
            T *nuTur0 = cell.getExternal(Descriptor<T>::ExternalField::tauTBeginsAt);
            T nuTur = nuTur0[0];
            // T tauLast = tau + nuTur / (Descriptor<T>::cs2);
            T tauLast = tau + nuTur / (Descriptor<T>::cs2 * cCs);

            T s7 = 1. / tauLast;
            T s8 = 1. / tauLast;

            T s8Last = 1. / tauLast;

            T sMagin = 1.0 / (1.0 / (12 * (tauLast - 0.5)) + 0.5);

            //    { (T)1, (T)1.63, (T)1.14, (T)1, (T)1.92, (T)1, (T)1.92, T(), T() };

#ifdef lesModel
            T mneq0 = m0 - meq0;
            T mneq1 = m1 - meq1;
            T mneq7 = m7 - meq7;
            T mneq8 = m8 - meq8;

            T Sxx = 2.0 / 3.0 * mneq0 * s0 + 1.0 / 6.0 * mneq1 * s1 + 1.0 / 2.0 * mneq7 * s7;
            T Sxy = mneq8 * s8;
            T Syx = mneq8 * s8;
            T Syy = 2.0 / 3.0 * mneq0 * s0 + 1.0 / 6.0 * mneq1 * s1 - 1.0 / 2.0 * mneq7 * s7;

            T Snorm2 = Sxx * Sxx + Sxy * Sxy //
                       + Syx * Syx + Syy * Syy + 1e-16;

            T SpreFactor = -1.0 * cCs / 2.0 / rhoFull / (Descriptor<T>::cs2);
            T SsSqrt = sqrt(2) * sqrt(Snorm2) * sqrt(SpreFactor * SpreFactor);

            T nuT = pow((0.17 * sqrt(1.)), 2) * SsSqrt;

            // T tauTM05 = nuT / (Descriptor<T>::cs2);
            T tauTM05 = nuT / (Descriptor<T>::cs2 * cCs);

            T tau_u = tauTM05 + tau;
            nuTur0[0] = nuT;
#else
            T tau_u = tau;
            nuTur0[0] = 0;
            // pcout << "------------" << std::endl;
#endif

            s7 = 1. / tau_u;
            s8 = 1. / tau_u;

            sMagin = 1.0 / (1.0 / (12 * (tau_u - 0.5)) + 0.5);
            s1 = sMagin;
            s2 = sMagin;
            s4 = sMagin;
            s6 = sMagin;

            T mcp0 = m0 * (1.0 - s0) + s0 * meq0; //+ (1 - s0 / 2) * mF0;
            T mcp1 = m1 * (1.0 - s1) + s1 * meq1; //+ (1 - s1 / 2) * mF1;
            T mcp2 = m2 * (1.0 - s2) + s2 * meq2; //+ (1 - s2 / 2) * mF2;
            T mcp3 = m3 * (1.0 - s3) + s3 * meq3 + fx / cCs / cCs;
            T mcp4 = m4 * (1.0 - s4) + s4 * meq4; //+ (1 - s4 / 2) * mF4;
            T mcp5 = m5 * (1.0 - s5) + s5 * meq5 + fy / cCs / cCs;
            T mcp6 = m6 * (1.0 - s6) + s6 * meq6; // + (1 - s6 / 2) * mF6;
            T mcp7 = m7 * (1.0 - s7) + s7 * meq7; // + (1 - s7 / 2) * mF7;
            T mcp8 = m8 * (1.0 - s8) + s8 * meq8; // + (1 - s8 / 2) * mF8;

            cell[0] = mcp0 / 9 - mcp1 / 9 + mcp2 / 9;
            cell[1] = mcp0 / 9 + mcp1 / 18 + mcp2 / 36 - mcp3 / 6 - mcp4 / 12 + mcp5 / 6 + mcp6 / 12 - mcp8 / 4;
            cell[2] = mcp0 / 9 - mcp1 / 36 - mcp2 / 18 - mcp3 / 6 + mcp4 / 6 + mcp7 / 4;
            cell[3] = mcp0 / 9 + mcp1 / 18 + mcp2 / 36 - mcp3 / 6 - mcp4 / 12 - mcp5 / 6 - mcp6 / 12 + mcp8 / 4;
            cell[4] = mcp0 / 9 - mcp1 / 36 - mcp2 / 18 - mcp5 / 6 + mcp6 / 6 - mcp7 / 4;
            cell[5] = mcp0 / 9 + mcp1 / 18 + mcp2 / 36 + mcp3 / 6 + mcp4 / 12 - mcp5 / 6 - mcp6 / 12 - mcp8 / 4;
            cell[6] = mcp0 / 9 - mcp1 / 36 - mcp2 / 18 + mcp3 / 6 - mcp4 / 6 + mcp7 / 4;
            cell[7] = mcp0 / 9 + mcp1 / 18 + mcp2 / 36 + mcp3 / 6 + mcp4 / 12 + mcp5 / 6 + mcp6 / 12 + mcp8 / 4;
            cell[8] = mcp0 / 9 - mcp1 / 36 - mcp2 / 18 + mcp5 / 6 - mcp6 / 6 - mcp7 / 4;

            ux += fx / 2 / cCs / rhoFull;
            uy += fy / 2 / cCs / rhoFull;
            vel[0] = ux;
            vel[1] = uy;
        }

        virtual void computeVelocity(
            Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const
        {

            const T *vel = cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt);
            // pcout << "hello" << std::endl;
            u[0] = vel[0];
            u[1] = vel[1];
        }

        virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j,
                                     T jSqr, T thetaBar = T()) const
        {
            // pcout << "hello equ" << std::endl;
            T invRho = Descriptor<T>::invRho(rhoBar);
            T invCs2_01 = Descriptor<T>::invCs2 / (cCs * cCs);

            T c_j = Descriptor<T>::c[iPop][0] * j[0] * cCs;
            for (int iD = 1; iD < Descriptor<T>::d; ++iD)
            {
                c_j += Descriptor<T>::c[iPop][iD] * j[iD] * cCs;
            }
            return Descriptor<T>::t[iPop] * (rhoBar + invCs2_01 * c_j +
                                             invCs2_01 / (T)2 * invRho * (invCs2_01 * c_j * c_j - jSqr));
        }

        T getCCs() const { return cCs; }

    private:
        static int id;
        const T cCs = ccDefine;
    };
    template <typename T, template <typename U> class Descriptor>
    int MyGuoD2Q9ExternalForceLesMRTdynamics02_direct_MultiCs<T, Descriptor>::id =
        meta::registerGeneralDynamics<T, Descriptor, MyGuoD2Q9ExternalForceLesMRTdynamics02_direct_MultiCs<T, Descriptor>>("MyGuoD2Q9ExternalForceLesMRTdynamics02_direct_MultiCs");

    /// Implementation of the MRT force LES collision step
    template <typename T, template <typename U> class Descriptor>
    class MyGuoD2Q9ExternalForceLesForceGradMRTdynamics01MultiCs : public ExternalForceDynamics<T, Descriptor>
    {
    public:
        /* *************** Construction / Destruction ************************ */
        MyGuoD2Q9ExternalForceLesForceGradMRTdynamics01MultiCs(T omega_)
            : ExternalForceDynamics<T, Descriptor>(omega_) //, cCs(1.0)
        {
        }
        MyGuoD2Q9ExternalForceLesForceGradMRTdynamics01MultiCs(T omega_, T cCs_)
            : ExternalForceDynamics<T, Descriptor>(omega_), cCs(cCs_)
        {
        }

        MyGuoD2Q9ExternalForceLesForceGradMRTdynamics01MultiCs(HierarchicUnserializer &unserializer)
            : ExternalForceDynamics<T, Descriptor>(T())
        {
            this->unserialize(unserializer);
        }

        /// Clone the object on its dynamic type.
        virtual MyGuoD2Q9ExternalForceLesForceGradMRTdynamics01MultiCs<T, Descriptor> *clone() const
        {
            return new MyGuoD2Q9ExternalForceLesForceGradMRTdynamics01MultiCs<T, Descriptor>(*this);
        }

        /// Return a unique ID for this class.
        virtual int getId() const
        {
            return id;
        }

        virtual T computeTemperature(Cell<T, Descriptor> const &cell) const
        {
            const T *s8_0 = cell.getExternal(Descriptor<T>::ExternalField::tauTBeginsAt);
            return s8_0[0];
            // const T *s8_0 = cell.getExternal(Descriptor<T>::ExternalField::tauTBeginsAt);
            // T s8 = s8_0[0];
            // return (1. / s8) * (Descriptor<T>::invCs2 / cCs);
        }

        virtual void collide(Cell<T, Descriptor> &cell,
                             BlockStatistics &statistics_)
        {
            T rhoBar = momentTemplates<T, Descriptor>::get_rhoBar(cell);
            const T rhoFull = Descriptor<T>::fullRho(rhoBar);
            Array<T, Descriptor<T>::d> u;
            // this->computeVelocity(cell, u);

            Array<T, Descriptor<T>::d> j;
            momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
            // force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));
            for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
            {
                u[iD] = j[iD] / rhoFull * cCs; // + force[iD] / cCs * 0.5 * invRho;
            }

            Array<T, Descriptor<T>::d> force;
            force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));

            const T rho = rhoBar;
            T ux = u[0];
            T uy = u[1];
            T fx = force[0];
            T fy = force[1];

            T *vel = cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt);

            const T f0 = cell[0];
            const T f1 = cell[1];
            const T f2 = cell[2];
            const T f3 = cell[3];
            const T f4 = cell[4];
            const T f5 = cell[5];
            const T f6 = cell[6];
            const T f7 = cell[7];
            const T f8 = cell[8];

            const T m0 = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
            const T m1 = -4 * f0 + 2 * f1 - f2 + 2 * f3 - f4 + 2 * f5 - f6 + 2 * f7 - f8;
            const T m2 = 4 * f0 + f1 - 2 * f2 + f3 - 2 * f4 + f5 - 2 * f6 + f7 - 2 * f8;
            const T m3 = -f1 - f2 - f3 + f5 + f6 + f7;
            const T m4 = -f1 + 2 * f2 - f3 + f5 - 2 * f6 + f7;
            const T m5 = f1 - f3 - f4 - f5 + f7 + f8;
            const T m6 = f1 - f3 + 2 * f4 - f5 + f7 - 2 * f8;
            const T m7 = f2 - f4 + f6 - f8;
            const T m8 = -f1 + f3 - f5 + f7;

            const T meq0 = 1.0 * rhoBar;
            const T meq1 = (3.0 * rhoFull * ux * ux + 3.0 * rhoFull * uy * uy - 2.0 * rhoBar * cCs * cCs) / (cCs * cCs);
            const T meq2 = (-3.0 * rhoFull * ux * ux - 3.0 * rhoFull * uy * uy + rhoBar * cCs * cCs) / (cCs * cCs);
            const T meq3 = 1.0 * rhoFull * ux / cCs;
            const T meq4 = -1.0 * rhoFull * ux / cCs;
            const T meq5 = 1.0 * rhoFull * uy / cCs;
            const T meq6 = -1.0 * rhoFull * uy / cCs;
            const T meq7 = 1.0 * rhoFull * (ux * ux - uy * uy) / (cCs * cCs);
            const T meq8 = 1.0 * rhoFull * ux * uy / (cCs * cCs);

            const T omega = ExternalForceDynamics<T, Descriptor>::getOmega();
            const T tau = 1. / omega;

            T s0 = 1;

            T s3 = 1;
            T s5 = 1;

            T *nuTur0 = cell.getExternal(Descriptor<T>::ExternalField::tauTBeginsAt);
            T nuTur = nuTur0[0];
            T tauLast = tau + nuTur / (Descriptor<T>::cs2 * cCs);
            T s8Last = 1. / tauLast;

            T sMagin = 1.0 / (1.0 / (12 * (tauLast - 0.5)) + 0.5);
            T s1 = sMagin;
            T s2 = sMagin;
            T s4 = sMagin;
            T s6 = sMagin;

            T s7 = s8Last;
            T s8 = s8Last;

            if (1. / s8 <= 0.5 || 1. / s7 <= 0.5)
            {
                pcout << "error S7 s8 01 ~~~!!!!!!!!!!!\n";
            }

            // double preFactor = -c * c / 2 / rhoFull / Cs2 / dt;
            T preFactor = -3 * cCs / 2 / rhoFull;
            T mneq0 = m0 - meq0;
            T mneq1 = m1 - meq1;
            T mneq7 = m7 - meq7;
            T mneq8 = m8 - meq8;

            // 这是之前幂律流体做的东西
            Array<T, SymmetricTensor<T, Descriptor>::n> stress;
            this->computePiNeq(cell, stress);
            T D2 = (pow(stress[0] * cCs, 2) + pow(stress[1] * cCs, 2) + 2 * pow(stress[2] * cCs, 2) + 1e-20);
            T para = s8Last / ((rhoBar + 1) * Descriptor<T>::cs2);
            para *= (0.5 * para);
            T gamma = sqrt(2 * D2 * para);

            T Sxx = 4 * mneq0 * s0 + mneq1 * s1 + 3 * mneq7 * s7;
            T Syy = 4 * mneq0 * s0 + mneq1 * s1 - 3 * mneq7 * s7;
            T Sxy = mneq8 * s8;
            Sxx *= (preFactor / 6.0);
            Syy *= (preFactor / 6.0);
            Sxy *= preFactor;

            double SsSqrt = sqrt(2) * sqrt(Sxx * Sxx + Syy * Syy + 2 * Sxy * Sxy + 1e-10);
            SsSqrt = gamma;

            double nuT = pow((0.18 * 1.414), 2) * SsSqrt;
            // double nuT = pow((0.16 * 1.414), 2) * SsSqrt;
            double tauTM05 = nuT / (Descriptor<T>::cs2 * cCs);
            double tau_u = tauTM05 + tau;
            s7 = 1. / tau_u;
            s8 = 1. / tau_u;

            sMagin = 1.0 / (1.0 / (12 * (tau_u - 0.5)) + 0.5);
            s1 = sMagin;
            s2 = sMagin;
            s4 = sMagin;
            s6 = sMagin;

            // nuT = 0;
            nuTur0[0] = nuT;

            if (1. / s8 <= 0.5 || 1. / s7 <= 0.5)
            {
                pcout << "error S7 s8 022 ~~~!!!!!!!!!!!\n";
            }

            T mcp0 = m0 * (1.0 - s0) + s0 * meq0; //+ (1 - s0 / 2) * mF0;
            T mcp1 = m1 * (1.0 - s1) + s1 * meq1; //+ (1 - s1 / 2) * mF1;
            T mcp2 = m2 * (1.0 - s2) + s2 * meq2; //+ (1 - s2 / 2) * mF2;
            T mcp3 = m3 * (1.0 - s3) + s3 * meq3 + fx / cCs / cCs;
            T mcp4 = m4 * (1.0 - s4) + s4 * meq4; //+ (1 - s4 / 2) * mF4;
            T mcp5 = m5 * (1.0 - s5) + s5 * meq5 + fy / cCs / cCs;
            T mcp6 = m6 * (1.0 - s6) + s6 * meq6; // + (1 - s6 / 2) * mF6;
            T mcp7 = m7 * (1.0 - s7) + s7 * meq7; // + (1 - s7 / 2) * mF7;
            T mcp8 = m8 * (1.0 - s8) + s8 * meq8; // + (1 - s8 / 2) * mF8;

            cell[0] = mcp0 / 9 - mcp1 / 9 + mcp2 / 9;
            cell[1] = mcp0 / 9 + mcp1 / 18 + mcp2 / 36 - mcp3 / 6 - mcp4 / 12 + mcp5 / 6 + mcp6 / 12 - mcp8 / 4;
            cell[2] = mcp0 / 9 - mcp1 / 36 - mcp2 / 18 - mcp3 / 6 + mcp4 / 6 + mcp7 / 4;
            cell[3] = mcp0 / 9 + mcp1 / 18 + mcp2 / 36 - mcp3 / 6 - mcp4 / 12 - mcp5 / 6 - mcp6 / 12 + mcp8 / 4;
            cell[4] = mcp0 / 9 - mcp1 / 36 - mcp2 / 18 - mcp5 / 6 + mcp6 / 6 - mcp7 / 4;
            cell[5] = mcp0 / 9 + mcp1 / 18 + mcp2 / 36 + mcp3 / 6 + mcp4 / 12 - mcp5 / 6 - mcp6 / 12 - mcp8 / 4;
            cell[6] = mcp0 / 9 - mcp1 / 36 - mcp2 / 18 + mcp3 / 6 - mcp4 / 6 + mcp7 / 4;
            cell[7] = mcp0 / 9 + mcp1 / 18 + mcp2 / 36 + mcp3 / 6 + mcp4 / 12 + mcp5 / 6 + mcp6 / 12 + mcp8 / 4;
            cell[8] = mcp0 / 9 - mcp1 / 36 - mcp2 / 18 + mcp5 / 6 - mcp6 / 6 - mcp7 / 4;

            ux += fx / 2 / cCs / rhoFull;
            uy += fy / 2 / cCs / rhoFull;
            vel[0] = ux;
            vel[1] = uy;
        }

        virtual void computeVelocity(
            Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const
        {

            const T *vel = cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt);
            // pcout << "hello" << std::endl;
            u[0] = vel[0];
            u[1] = vel[1];
        }

        virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j,
                                     T jSqr, T thetaBar = T()) const
        {
            // pcout << "hello equ" << std::endl;
            T invRho = Descriptor<T>::invRho(rhoBar);
            T invCs2_01 = Descriptor<T>::invCs2 / (cCs * cCs);

            T c_j = Descriptor<T>::c[iPop][0] * j[0] * cCs;
            for (int iD = 1; iD < Descriptor<T>::d; ++iD)
            {
                c_j += Descriptor<T>::c[iPop][iD] * j[iD] * cCs;
            }
            return Descriptor<T>::t[iPop] * (rhoBar + invCs2_01 * c_j +
                                             invCs2_01 / (T)2 * invRho * (invCs2_01 * c_j * c_j - jSqr));
        }

        T getCCs() const { return cCs; }

    private:
        static int id;
        const T cCs = ccDefine;
    };
    template <typename T, template <typename U> class Descriptor>
    int MyGuoD2Q9ExternalForceLesForceGradMRTdynamics01MultiCs<T, Descriptor>::id =
        meta::registerGeneralDynamics<T, Descriptor, MyGuoD2Q9ExternalForceLesForceGradMRTdynamics01MultiCs<T, Descriptor>>("MyGuoD2Q9ExternalForceLesForceGradMRTdynamics01MultiCs");

    // -------------------------------------------------------------------
    template <typename T, template <typename U> class Descriptor>
    class AdvectionDiffusionWithMultiNSLesMRT2Ddynamics : public AdvectionDiffusionDynamics<T, Descriptor>
    {
    public:
        /// Constructor
        /// Constructor
        AdvectionDiffusionWithMultiNSLesMRT2Ddynamics(T omega_)
            : AdvectionDiffusionDynamics<T, Descriptor>(omega_)
        {
        }
        AdvectionDiffusionWithMultiNSLesMRT2Ddynamics(HierarchicUnserializer &unserializer)
            : AdvectionDiffusionDynamics<T, Descriptor>(T())
        {
            this->unserialize(unserializer);
        }

        virtual T computeTemperature(Cell<T, Descriptor> const &cell) const
        {
            return *(cell.getExternal(Descriptor<T>::ExternalField::tauDBeginsAt));
        }

        /// Clone the object on its dynamic type.
        virtual AdvectionDiffusionWithMultiNSLesMRT2Ddynamics<T, Descriptor> *clone() const
        {
            return new AdvectionDiffusionWithMultiNSLesMRT2Ddynamics<T, Descriptor>(*this);
        }
        /// Return a unique ID for this class.
        virtual int getId() const
        {
            return id;
        }
        /// Collision step
        virtual void collide(Cell<T, Descriptor> &cell,
                             BlockStatistics &statistics)
        {
            T rhoBar;
            Array<T, Descriptor<T>::d> jEq;
            advectionDiffusionMomentTemplates<T, Descriptor>::get_rhoBar_jEq(cell, rhoBar, jEq);

            T invRho = Descriptor<T>::invRho(rhoBar);
            const T jSqr = VectorTemplateImpl<T, Descriptor<T>::d>::normSqr(jEq);

            T f0 = cell[0];
            T f1 = cell[3];
            T f2 = cell[4];
            T f3 = cell[1];
            T f4 = cell[2];

            T m0 = f0 + f1 + f2 + f3 + f4;
            T m1 = f1 - f3;
            T m2 = f2 - f4;
            T m3 = f1 - f2 + f3 - f4;
            T m4 = -4 * f0 + f1 + f2 + f3 + f4;

            T w0 = Descriptor<T>::t[0];
            T wc = Descriptor<T>::t[1];
            T cs2 = Descriptor<T>::cs2;

            T meq0 = rhoBar * (w0 + 4 * wc);
            T meq1 = 2 * jEq[0] * wc / cs2;
            T meq2 = 2 * jEq[1] * wc / cs2;
            T meq3 = 0;
            T meq4 = 4 * rhoBar * (-w0 + wc);

            T omg = AdvectionDiffusionDynamics<T, Descriptor>::getOmega();
            T tau = 1.0 / omg;
            const T *DTur = cell.getExternal(Descriptor<T>::ExternalField::tauDBeginsAt);
            T tauAll = tau + DTur[0] / cs2;
            T s0 = 1.0;
            T s1 = 1. / tauAll;
            T s2 = 1. / tauAll;
            // T s1 = omg;
            // T s2 = omg;
            T s3 = 1.0 / (1.0 / (12 * (tauAll - 0.5)) + 0.5);
            T s4 = s3;

            T mCollPost0 = m0 * (1.0 - s0) + s0 * meq0;
            T mCollPost1 = m1 * (1.0 - s1) + s1 * meq1;
            T mCollPost2 = m2 * (1.0 - s2) + s2 * meq2;
            T mCollPost3 = m3 * (1.0 - s3) + s3 * meq3;
            T mCollPost4 = m4 * (1.0 - s4) + s4 * meq4;

            cell[0] = mCollPost0 / 5. - mCollPost4 / 5.;
            cell[3] = mCollPost0 / 5. + mCollPost1 / 2. + mCollPost3 / 4. + mCollPost4 / 20.;
            cell[4] = mCollPost0 / 5. + mCollPost2 / 2. - mCollPost3 / 4. + mCollPost4 / 20.;
            cell[1] = mCollPost0 / 5. - mCollPost1 / 2. + mCollPost3 / 4. + mCollPost4 / 20.;
            cell[2] = mCollPost0 / 5. - mCollPost2 / 2. - mCollPost3 / 4. + mCollPost4 / 20.;
        }

        virtual void collideExternal(
            Cell<T, Descriptor> &cell, T rhoBar,
            Array<T, Descriptor<T>::d> const &jEq, T thetaBar, BlockStatistics &stat)
        {
        }

        virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &jEq,
                                     T jSqr, T thetaBar = T()) const
        {
            T c_j = Descriptor<T>::c[iPop][0] * jEq[0];
            for (int iD = 1; iD < Descriptor<T>::d; ++iD)
            {
                c_j += Descriptor<T>::c[iPop][iD] * jEq[iD];
            }
            return Descriptor<T>::t[iPop] * (rhoBar + Descriptor<T>::invCs2 * c_j);
        }

    private:
        // T s0, s1, s2, s3, s4;
        // static int id;
        static int id;
    };
    template <typename T, template <typename U> class Descriptor>
    int AdvectionDiffusionWithMultiNSLesMRT2Ddynamics<T, Descriptor>::id =
        meta::registerGeneralDynamics<T, Descriptor, AdvectionDiffusionWithMultiNSLesMRT2Ddynamics<T, Descriptor>>("AdvectionDiffusionWithMultiNSLesMRT2Ddynamics");

    template <typename T, template <typename U> class Descriptor>
    class AdvectionDiffusionLesBGKdynamics : public AdvectionDiffusionDynamics<T, Descriptor>
    {
    public:
        /// Constructor
        AdvectionDiffusionLesBGKdynamics(T omega_)
            : AdvectionDiffusionDynamics<T, Descriptor>(omega_)
        {
        }

        AdvectionDiffusionLesBGKdynamics(HierarchicUnserializer &unserializer)
            : AdvectionDiffusionDynamics<T, Descriptor>(T())
        {
            this->unserialize(unserializer);
        }
        /// Clone the object on its dynamic type.
        virtual AdvectionDiffusionLesBGKdynamics<T, Descriptor> *clone() const
        {
            return new AdvectionDiffusionLesBGKdynamics<T, Descriptor>(*this);
        }

        /// Return a unique ID for this class.
        virtual int getId() const
        {
            return id;
        }

        virtual T computeTemperature(Cell<T, Descriptor> const &cell) const
        {
            // return *(cell.getExternal(Descriptor<T>::ExternalField::tauDBeginsAt));
            return *(cell.getExternal(Descriptor<T>::ExternalField::tauDBeginsAt)) + 0;
            // return *(cell.getExternal(Descriptor<T>::ExternalField::tauDBeginsAt)) + Descriptor<T>::cs2 * (1. / this->getOmega() - 0.5);
            // return Descriptor<T>::cs2 * (1. / this->getOmega() - 0.5);
        }

        /// Collision step
        virtual void collide(Cell<T, Descriptor> &cell,
                             BlockStatistics &statistics)
        {
            T rhoBar;
            Array<T, Descriptor<T>::d> jEq;
            advectionDiffusionMomentTemplates<T, Descriptor>::get_rhoBar_jEq(cell, rhoBar, jEq);

            T cs2 = Descriptor<T>::cs2;
            const T *DTur = cell.getExternal(Descriptor<T>::ExternalField::tauDBeginsAt);
            T tauAll = 1. / this->getOmega() + DTur[0] / cs2;

            T uSqr = advectionDiffusionDynamicsTemplates<T, Descriptor>::
                no_corr_bgk_collision(cell, rhoBar, jEq, 1. / tauAll);
            // no_corr_bgk_collision(cell, rhoBar, jEq, this->getOmega());

            if (cell.takesStatistics())
            {
                gatherStatistics(statistics, rhoBar, uSqr);
            }
        }
        /// Implementation of the collision step, with imposed macroscopic variables
        /// The arguments:
        /// - rhoBar: the "rhoBar" version of the scalar rho.
        /// - jEq: the equilibrium part of the second-order moment. jEq = u*rho, where u is the external convective term.

        /// Compute equilibrium distribution function
        virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &jEq,
                                     T jSqr, T thetaBar = T()) const
        {
            return advectionDiffusionDynamicsTemplates<T, Descriptor>::
                bgk_ma1_equilibrium(iPop, rhoBar, jEq);
        }

    private:
        static int id;
    };
    template <typename T, template <typename U> class Descriptor>
    int AdvectionDiffusionLesBGKdynamics<T, Descriptor>::id =
        meta::registerGeneralDynamics<T, Descriptor, AdvectionDiffusionLesBGKdynamics<T, Descriptor>>("AdvectionDiffusionLesBGKdynamics");

    // ------------------------------------------------------------------------template< typename T,
    template <typename T,
              template <typename U1> class FluidDescriptor,
              template <typename U2> class TemperatureDescriptor>
    class TurbulencePrandalProcessor2D : public BoxProcessingFunctional2D_LL<T, FluidDescriptor, T, TemperatureDescriptor>
    {
    public:
        TurbulencePrandalProcessor2D(T c_, T cCs_) : Pr_tur(c_), cCs(cCs_)
        {
        }
        virtual void process(Box2D domain,
                             BlockLattice2D<T, FluidDescriptor> &fluid,
                             BlockLattice2D<T, TemperatureDescriptor> &temperature)
        {
            typedef FluidDescriptor<T> D;
            enum
            {
                velOffset = TemperatureDescriptor<T>::ExternalField::tauDBeginsAt,
                forceOffset = FluidDescriptor<T>::ExternalField::tauTBeginsAt
            };
            Dot2D offset = computeRelativeDisplacement(fluid, temperature);

            for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            {
                for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                {
                    // Velocity coupling
                    T *sAD = temperature.get(iX + offset.x, iY + offset.y).getExternal(velOffset);
                    T *sNS = fluid.get(iX, iY).getExternal(forceOffset);

                    sAD[0] = sNS[0] / Pr_tur;
                }
            }
        }
        virtual TurbulencePrandalProcessor2D<T, FluidDescriptor, TemperatureDescriptor> *clone() const
        {
            return new TurbulencePrandalProcessor2D<T, FluidDescriptor, TemperatureDescriptor>(*this);
        }
        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
        {
            modified[0] = modif::staticVariables;
            modified[1] = modif::staticVariables;
        }

    private:
        T Pr_tur, cCs;
    };

    //-------------------------bounce back -------------------------
    //-------------------------bounce back -------------------------
    template <typename T, template <typename U> class Descriptor>
    class myBounceBack : public BounceBack<T, Descriptor>
    {
    public:
        /* *************** Construction / Destruction ******************************* */

        /// You may fix a fictitious density value on bounce-back nodes via the constructor.
        myBounceBack(T rho_ = 0) : rho(rho_)
        {
            pcout << "my bounce back\n";
        }
        myBounceBack(HierarchicUnserializer &unserializer)
            : rho(0), cCs(1)
        {
            this->unserialize(unserializer);
        }
        /// Clone the object on its dynamic type.
        virtual myBounceBack<T, Descriptor> *clone() const
        {
            return new myBounceBack<T, Descriptor>(*this);
        }
        /// Return a unique ID for this class.
        virtual int getId() const
        {
            return id;
        }

        /// Serialize the dynamics object.
        virtual void serialize(HierarchicSerializer &serializer) const
        {
            Dynamics<T, Descriptor>::serialize(serializer);
            serializer.addValue(rho);
            // serializer.addValue(cCs);
        }

        /// Un-Serialize the dynamics object.
        virtual void unserialize(HierarchicUnserializer &unserializer)
        {
            PLB_PRECONDITION(unserializer.getId() == this->getId());
            Dynamics<T, Descriptor>::unserialize(unserializer);
            rho = unserializer.readValue<T>();
            // cCs = unserializer.readValue<T>();
        }

        virtual void collide(
            Cell<T, Descriptor> &cell,
            BlockStatistics &statistics)
        {
            for (plint iPop = 1; iPop <= Descriptor<T>::q / 2; ++iPop)
            {
                std::swap(cell[iPop], cell[iPop + Descriptor<T>::q / 2]);
            }
        }
        T getcCs() const { return cCs; }

    private:
        T rho;
        T cCs;

    private:
        static int id;
    };

    template <typename T, template <typename U> class Descriptor>
    int myBounceBack<T, Descriptor>::id =
        meta::registerGeneralDynamics<T, Descriptor, myBounceBack<T, Descriptor>>("myBounceBack");
}

#endif

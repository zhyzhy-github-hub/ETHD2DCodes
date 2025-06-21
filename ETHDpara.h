#ifndef PARARBC_H
#define PARARBC_H

#include "../../src/core/units.h"
#include "../../src/core/array.h"
#include "../../src/core/globalDefs.h"
#include "../../src/core/dynamics.h"
#include <string>
#include <fstream>

// typedef double T;

namespace plb
{
    /// A useful class for the conversion between dimensionless and lattice units.
    template <typename T,
              template <typename NSU> class nsDescriptor,
              template <typename ADELE> class adElecDescriptor,
              template <typename ADTEM> class adTemDescriptor>
    class MyETHDParameters2D_Nv_Sc
    {

    public:
        /// Constructor
        /** \param T   Electric Rayleigh number
         *  \param M   character velocity
         *  \param C   EHD C number
         *  \param Sc  Sc number,
         *  \param Ma  Mach number
         *  \param Ra  Thercal Rayleigh number
         *  \param Pr  Prandtl number, thermal
         *  \param lx  length x
         *  \param ly  length y
         *  \param resolution  grids in y dircetion
         *  \param deltaPhi voltage
         *  \param rho0  density, rho0
         */
        MyETHDParameters2D_Nv_Sc(T cc_, T Tehd_, T Mehd_, T Cehd_, T Sc_,
                                 T Nv_, T Ra_, T Pr_, T deltaPhi_, T deltaT_,
                                 T lx_, T ly_, T resolution_,
                                 T rho0_ = T(1.0))
            : cc(cc_), Tehd(Tehd_), Mehd(Mehd_), Cehd(Cehd_), Sc(Sc_),
              Nv(Nv_), Ra(Ra_), Pr(Pr_), deltaPhi(deltaPhi_), deltaT(deltaT_),
              lx(lx_), ly(ly_), resolution(resolution_),
              rho0(rho0_)
        {
        }

        /// x-length in dimensionless units
        T getLx() const { return lx; }
        /// y-length in dimensionless units
        T getLy() const { return ly; }

        T getTETHD() const { return Tehd; }
        T getMETHD() const { return Mehd; }
        T getCETHD() const { return Cehd; }
        T getScETHD() const { return Sc; }
        T getPrETHD() const { return Pr; }
        T getRaETHD() const { return Ra; }
        T getDeltaPhi() const { return deltaPhi; }
        T getResolution() const { return resolution; }
        T getRho0() const { return rho0; }
        T getDeltaT() const { return deltaT; }

        /// lattice spacing in dimensionless units
        T getDeltaX() const { return (T)1 / getResolution(); }
        /// conversion from dimensionless to lattice units for space coordinate
        plint nCell(T l) const { return (plint)(l / getDeltaX() + (T)0.5); }
        /// number of lattice cells in x-direction
        plint getNx() const { return nCell(lx) + 1; }
        /// number of lattice cells in y-direction
        plint getNy() const { return nCell(ly) + 1; }

        T getLatticeNu() const { return Nv * getLatticeCC(); }
        T getLatticeCC() const { return cc; }

        T getGBeta() const { return getRaETHD() * getLatticeNu() * getLatticeAlpha() / getDeltaT() / pow(getResolution(), 3); }
        T getLatticeAlpha() const { return getLatticeNu() / getPrETHD(); }
        T getLatticeK() const { return getLatticeNu() * getTETHD() / pow(getMETHD(), 2) / getDeltaPhi(); }
        T getMa() const { return getDeltaPhi() * getLatticeK() / getResolution() / (getLatticeCC() * sqrt(nsDescriptor<T>::cs2)); }
        T getLatticeEpsilon() const { return getRho0() * pow(getLatticeK() * getMETHD(), 2); }

        T getLatticeD() const { return getLatticeNu() / getScETHD(); }
        // T getLatticeD() const { return getDeltaPhi() * getLatticeK() / getScETHD(); }

        T getQ0() const { return getCETHD() * getLatticeEpsilon() * getDeltaPhi() / (pow(getResolution(), 2)); }
        T getReElec() const { return getTETHD() / (getMETHD() * getMETHD()); }
        T getReRBC() const { return sqrt(getRaETHD() / getPrETHD()); }

        T getTauNSE() const { return nsDescriptor<T>::invCs2 / getLatticeCC() * getLatticeNu() + (T)0.5; }

        T getTauElecADE() const { return adElecDescriptor<T>::invCs2 * getLatticeD() + (T)0.5; }
        T getTauTemADE() const { return adTemDescriptor<T>::invCs2 * getLatticeAlpha() + (T)0.5; }

        T getOmegaNSE() const { return 1.0 / getTauNSE(); }
        T getOmegaElecADE() const { return 1.0 / getTauElecADE(); }
        T getOmegaTemADE() const { return 1.0 / getTauTemADE(); }

        T getDimlessTime() const { return getLatticeNu() / pow(getResolution(), 2); }
        T getDimlessVelocity() const { return getResolution() / getLatticeK() / getDeltaPhi(); }
        ///****************if RBC_ELE
        // T getDimlessVelocity() const { return getResolution() / getLatticeNu(); }
        T getDimlessCurrent() const { return getDimlessVelocity() / getQ0(); }
        T getDimlessHeatFlux() const { return getDimlessVelocity() / getDeltaT(); }
        // T getDimlessCurrent() const { return getResolution() / getLatticeNu() / getQ0(); }

    private:
        T Tehd, Mehd, Cehd, Sc, Nv, Ra, Pr, lx, ly, deltaPhi, deltaT, resolution, rho0;
        T cc;
    };

    template <typename T,
              template <typename NSU> class nsDescriptor,
              template <typename ADELE> class adElecDescriptor,
              template <typename ADTEM> class adTemDescriptor>
    void writeLogFile(
        MyETHDParameters2D_Nv_Sc<T, nsDescriptor, adElecDescriptor, adTemDescriptor> const &parameters,
        std::string const &title)
    {
        std::string fullName = global::directories().getLogOutDir() + "plbLogEHD.dat";
        std::ofstream ofile(fullName.c_str());
        ofile << title << "\n\n";
        ofile << "Reynolds number Elec:           Re = " << parameters.getReElec() << "\n";
        ofile << "Reynolds number RBC:            Re = " << parameters.getReRBC() << "\n";
        ofile << "T number:                        T = " << parameters.getTETHD() << "\n";
        ofile << "M number:                        M = " << parameters.getMETHD() << "\n";
        ofile << "C number:                        C = " << parameters.getCETHD() << "\n";
        ofile << "Sc number:                      Sc = " << parameters.getScETHD() << "\n";
        ofile << "Ra number:                      Ra = " << parameters.getRaETHD() << "\n";
        ofile << "Pr number:                      Pr = " << parameters.getPrETHD() << "\n";
        ofile << "g x beta:                    gBeta = " << parameters.getGBeta() << "\n";

        ofile << "\n";
        ofile << "thermal ade diffusivity:     Alpha = " << parameters.getLatticeAlpha() << "\n";
        ofile << "Kinematic viscosity:            Nu = " << parameters.getLatticeNu() << "\n";
        ofile << "charge ade diffusivity:          D = " << parameters.getLatticeD() << "\n";
        ofile << "Solvent omega:           omega_nse = " << parameters.getOmegaNSE() << "\n";
        ofile << "electric omega:     elec omega_ade = " << parameters.getOmegaElecADE() << "\n";
        ofile << "Temperature omega:   te  omega_ade = " << parameters.getOmegaTemADE() << "\n";

        ofile << "Solvent tau:               tau_nse = " << parameters.getTauNSE() << "\n";
        ofile << "electric tau:         elec tau_ade = " << parameters.getTauElecADE() << "\n";
        ofile << "Temperature tau:       tem tau_ade = " << parameters.getTauTemADE() << "\n";

        ofile << "              cc :              cc = " << parameters.getLatticeCC() << "\n";
        ofile << "              Ma :              Ma = " << parameters.getMa() << "\n";
        ofile << "               N :               N = " << parameters.getResolution() << "\n";
        ofile << "              lx :              lx = " << parameters.getLx() << "\n";
        ofile << "              ly :              ly = " << parameters.getLy() << "\n";
        ofile << "              Nx :              Nx = " << parameters.getNx() << "\n";
        ofile << "              Ny :              Ny = " << parameters.getNy() << "\n";

        ofile << "       lattice K :       lattice K = " << parameters.getLatticeK() << "\n";
        ofile << " lattice epsilon : lattice epsilon = " << parameters.getLatticeEpsilon() << "\n";
        ofile << "              q0 :              q0 = " << parameters.getQ0() << "\n";
        ofile << "       delta phi :       delta phi = " << parameters.getDeltaPhi() << "\n";
        ofile << "       delta Tem :       delta Tem = " << parameters.getDeltaT() << "\n";

        ofile << "     getDimlessTime : = " << parameters.getDimlessTime() << "\n";
        ofile << " getDimlessVelocity : = " << parameters.getDimlessVelocity() << "\n";
        ofile << "  getDimlessCurrent : = " << parameters.getDimlessCurrent() << "\n";
        ofile << "  getDimlessHeatFlux: = " << parameters.getDimlessHeatFlux() << "\n";
    }

    /// A useful class for the conversion between dimensionless and lattice units.
    template <typename T,
              template <typename NSU> class nsDescriptor,
              template <typename ADELE> class adElecDescriptor,
              template <typename ADTEM> class adTemDescriptor>
    class MyETHDParameters2D_Nv
    {

    public:
        /// Constructor
        /** \param T   Electric Rayleigh number
         *  \param M   character velocity
         *  \param C   EHD C number
         *  \param Fe  Fe number,
         *  \param Ma  Mach number
         *  \param Ra  Thercal Rayleigh number
         *  \param Pr  Prandtl number, thermal
         *  \param lx  length x
         *  \param ly  length y
         *  \param resolution  grids in y dircetion
         *  \param deltaPhi voltage
         *  \param rho0  density, rho0
         */
        MyETHDParameters2D_Nv(T Tehd_, T Mehd_, T Cehd_, T Fe_,
                              T Nv_, T Ra_, T Pr_, T deltaPhi_, T deltaT_,
                              T lx_, T ly_, T resolution_,
                              T rho0_ = T(1.0))
            : Tehd(Tehd_), Mehd(Mehd_), Cehd(Cehd_), Fe(Fe_),
              Nv(Nv_), Ra(Ra_), Pr(Pr_), deltaPhi(deltaPhi_), deltaT(deltaT_),
              lx(lx_), ly(ly_), resolution(resolution_),
              rho0(rho0_)
        {
        }

        /// x-length in dimensionless units
        T getLx() const { return lx; }
        /// y-length in dimensionless units
        T getLy() const { return ly; }

        T getTETHD() const { return Tehd; }
        T getMETHD() const { return Mehd; }
        T getCETHD() const { return Cehd; }
        T getFeETHD() const { return Fe; }
        T getPrETHD() const { return Pr; }
        T getRaETHD() const { return Ra; }

        T getDeltaPhi() const { return deltaPhi; }
        /// resolution (a lattice of size 1 has getN()+1 cells)
        T getResolution() const { return resolution; }

        T getRho0() const { return rho0; }
        T getDeltaT() const { return deltaT; }

        /// lattice spacing in dimensionless units
        T getDeltaX() const { return (T)1 / getResolution(); }
        /// conversion from dimensionless to lattice units for space coordinate
        plint nCell(T l) const { return (plint)(l / getDeltaX() + (T)0.5); }
        /// number of lattice cells in x-direction
        plint getNx() const { return nCell(lx) + 1; }
        /// number of lattice cells in y-direction
        plint getNy() const { return nCell(ly) + 1; }

        T getLatticeNu() const { return Nv; }

        T getGBeta() const
        {
            return getRaETHD() * getLatticeNu() * getLatticeAlpha() / getDeltaT() / pow(getResolution(), 3);
        }

        T getLatticeAlpha() const
        {
            // return getMa() * sqrt(nsDescriptor<T>::cs2) * getResolution() / getPrETHD();
            return getLatticeNu() / getPrETHD();
        }
        T getLatticeK() const
        {
            return getLatticeNu() * getTETHD() / pow(getMETHD(), 2) / getDeltaPhi();
        }

        T getMa() const
        {
            return getDeltaPhi() * getLatticeK() / getResolution() / sqrt(nsDescriptor<T>::cs2);
        }
        T getLatticeEpsilon() const
        {
            return getRho0() * pow(getLatticeK() * getMETHD(), 2);
        }

        T getLatticeD() const
        {
            return getDeltaPhi() * getLatticeK() / getFeETHD();
        }
        T getQ0() const
        {
            return getCETHD() * getLatticeEpsilon() * getDeltaPhi() / (pow(getResolution(), 2));
        }

        T getReElec() const
        {
            return getTETHD() / (getMETHD() * getMETHD());
        }

        T getReRBC() const
        {
            return 0;
        }

        T getTauNSE() const { return nsDescriptor<T>::invCs2 * getLatticeNu() + (T)0.5; }
        T getTauElecADE() const { return adElecDescriptor<T>::invCs2 * getLatticeD() + (T)0.5; }
        T getTauTemADE() const { return adTemDescriptor<T>::invCs2 * getLatticeAlpha() + (T)0.5; }
        T getOmegaNSE() const { return 1.0 / getTauNSE(); }
        T getOmegaElecADE() const { return 1.0 / getTauElecADE(); }
        T getOmegaTemADE() const { return 1.0 / getTauTemADE(); }

        T getDimlessTime() const { return getLatticeNu() / pow(getResolution(), 2); }
        T getDimlessVelocity() const { return getResolution() / getLatticeK() / getDeltaPhi(); }
        ///****************if RBC_ELE
        // T getDimlessVelocity() const { return getResolution() / getLatticeNu(); }
        T getDimlessCurrent() const { return getDimlessVelocity() / getQ0(); }
        // T getDimlessCurrent() const { return getResolution() / getLatticeNu() / getQ0(); }

    private:
        T Tehd, Mehd, Cehd, Fe, Nv, Ra, Pr, lx, ly, deltaPhi, deltaT, resolution, rho0;
    };
    template <typename T,
              template <typename NSU> class nsDescriptor,
              template <typename ADELE> class adElecDescriptor,
              template <typename ADTEM> class adTemDescriptor>
    void writeLogFile(
        MyETHDParameters2D_Nv<T, nsDescriptor, adElecDescriptor, adTemDescriptor> const &parameters,
        std::string const &title)
    {
        std::string fullName = global::directories().getLogOutDir() + "plbLogEHD.dat";
        std::ofstream ofile(fullName.c_str());
        ofile << title << "\n\n";
        ofile << "Reynolds number Elec:           Re = " << parameters.getReElec() << "\n";
        ofile << "Reynolds number RBC:            Re = " << parameters.getReRBC() << "\n";
        ofile << "T number:                        T = " << parameters.getTETHD() << "\n";
        ofile << "M number:                        M = " << parameters.getMETHD() << "\n";
        ofile << "C number:                        C = " << parameters.getCETHD() << "\n";
        ofile << "Fe number:                      Fe = " << parameters.getFeETHD() << "\n";
        ofile << "Ra number:                      Ra = " << parameters.getRaETHD() << "\n";
        ofile << "Pr number:                      Pr = " << parameters.getPrETHD() << "\n";
        ofile << "g x beta:                    gBeta = " << parameters.getGBeta() << "\n";

        ofile << "\n";
        ofile << "thermal ade diffusivity:     Alpha = " << parameters.getLatticeAlpha() << "\n";
        ofile << "Kinematic viscosity:            Nu = " << parameters.getLatticeNu() << "\n";
        ofile << "charge ade diffusivity:          D = " << parameters.getLatticeD() << "\n";
        ofile << "Solvent omega:           omega_nse = " << parameters.getOmegaNSE() << "\n";
        ofile << "electric omega:     elec omega_ade = " << parameters.getOmegaElecADE() << "\n";
        ofile << "Temperature omega:   te  omega_ade = " << parameters.getOmegaTemADE() << "\n";

        ofile << "Solvent tau:               tau_nse = " << parameters.getTauNSE() << "\n";
        ofile << "electric tau:         elec tau_ade = " << parameters.getTauElecADE() << "\n";
        ofile << "Temperature tau:       tem tau_ade = " << parameters.getTauTemADE() << "\n";

        ofile << "              Ma :              Ma = " << parameters.getMa() << "\n";
        ofile << "               N :               N = " << parameters.getResolution() << "\n";
        ofile << "              lx :              lx = " << parameters.getLx() << "\n";
        ofile << "              ly :              ly = " << parameters.getLy() << "\n";
        ofile << "              Nx :              Nx = " << parameters.getNx() << "\n";
        ofile << "              Ny :              Ny = " << parameters.getNy() << "\n";

        ofile << "       lattice K :       lattice K = " << parameters.getLatticeK() << "\n";
        ofile << " lattice epsilon : lattice epsilon = " << parameters.getLatticeEpsilon() << "\n";
        ofile << "              q0 :              q0 = " << parameters.getQ0() << "\n";
        ofile << "       delta phi :       delta phi = " << parameters.getDeltaPhi() << "\n";
        ofile << "       delta Tem :       delta Tem = " << parameters.getDeltaT() << "\n";

        ofile << "     getDimlessTime : = " << parameters.getDimlessTime() << "\n";
        ofile << " getDimlessVelocity : = " << parameters.getDimlessVelocity() << "\n";
        ofile << "  getDimlessCurrent : = " << parameters.getDimlessCurrent() << "\n";
    }

    /// A useful class for the conversion between dimensionless and lattice units.
    template <typename T,
              template <typename NSU> class nsDescriptor,
              template <typename ADELE> class adElecDescriptor,
              template <typename ADTEM> class adTemDescriptor>
    class MyETHDParameters2D_Ma
    {

    public:
        /// Constructor
        /** \param T   Electric Rayleigh number
         *  \param M   character velocity
         *  \param C   EHD C number
         *  \param Fe  Fe number,
         *  \param Ma  Mach number
         *  \param Ra  Thercal Rayleigh number
         *  \param Pr  Prandtl number, thermal
         *  \param lx  length x
         *  \param ly  length y
         *  \param resolution  grids in y dircetion
         *  \param deltaPhi voltage
         *  \param rho0  density, rho0
         */
        MyETHDParameters2D_Ma(T Tehd_, T Mehd_, T Cehd_, T Fe_,
                              T Ma_, T Ra_, T Pr_, T deltaPhi_, T deltaT_,
                              T lx_, T ly_, T resolution_,
                              T rho0_ = T(1.0))
            : Tehd(Tehd_), Mehd(Mehd_), Cehd(Cehd_), Fe(Fe_),
              Ma(Ma_), Ra(Ra_), Pr(Pr_), deltaPhi(deltaPhi_), deltaT(deltaT_),
              lx(lx_), ly(ly_), resolution(resolution_),
              rho0(rho0_)
        {
        }

        /// x-length in dimensionless units
        T getLx() const { return lx; }
        /// y-length in dimensionless units
        T getLy() const { return ly; }

        T getTETHD() const { return Tehd; }
        T getMETHD() const { return Mehd; }
        T getCETHD() const { return Cehd; }
        T getFeETHD() const { return Fe; }
        T getPrETHD() const { return Pr; }
        T getRaETHD() const { return Ra; }

        T getDeltaPhi() const { return deltaPhi; }
        T getMa() const { return Ma; }
        /// resolution (a lattice of size 1 has getN()+1 cells)
        T getResolution() const { return resolution; }

        T getRho0() const { return rho0; }
        T getDeltaT() const { return deltaT; }

        /// lattice spacing in dimensionless units
        T getDeltaX() const { return (T)1 / getResolution(); }
        /// conversion from dimensionless to lattice units for space coordinate
        plint nCell(T l) const { return (plint)(l / getDeltaX() + (T)0.5); }
        /// number of lattice cells in x-direction
        plint getNx() const { return nCell(lx) + 1; }
        /// number of lattice cells in y-direction
        plint getNy() const { return nCell(ly) + 1; }

        T getLatticeNu() const
        {
            return getMa() * sqrt(nsDescriptor<T>::cs2) * getResolution();
        }

        T getGBeta() const
        {
            return getRaETHD() * getLatticeNu() * getLatticeAlpha() / getDeltaT() / pow(getResolution(), 3);
        }

        T getLatticeAlpha() const
        {
            // return getMa() * sqrt(nsDescriptor<T>::cs2) * getResolution() / getPrETHD();
            return getLatticeNu() / getPrETHD();
        }

        T getLatticeK() const
        {
            return getLatticeNu() * getTETHD() / pow(getMETHD(), 2) / getDeltaPhi();
        }

        T getLatticeEpsilon() const
        {
            return getRho0() * pow(getLatticeK() * getMETHD(), 2);
        }

        T getLatticeD() const
        {
            return getDeltaPhi() * getLatticeK() / getFeETHD();
        }
        T getQ0() const
        {
            return getCETHD() * getLatticeEpsilon() * getDeltaPhi() / (pow(getResolution(), 2));
        }

        T getReElec() const
        {
            return getTETHD() / (getMETHD() * getMETHD());
        }

        T getReRBC() const
        {
            return 0;
        }

        T getTauNSE() const { return nsDescriptor<T>::invCs2 * getLatticeNu() + (T)0.5; }
        T getTauElecADE() const { return adElecDescriptor<T>::invCs2 * getLatticeD() + (T)0.5; }
        T getTauTemADE() const { return adTemDescriptor<T>::invCs2 * getLatticeAlpha() + (T)0.5; }
        T getOmegaNSE() const { return 1.0 / getTauNSE(); }
        T getOmegaElecADE() const { return 1.0 / getTauElecADE(); }
        T getOmegaTemADE() const { return 1.0 / getTauTemADE(); }

        T getDimlessTime() const { return getLatticeNu() / pow(getResolution(), 2); }
        T getDimlessVelocity() const { return getResolution() / getLatticeK() / getDeltaPhi(); }
        ///****************if RBC_ELE
        // T getDimlessVelocity() const { return getResolution() / getLatticeNu(); }
        T getDimlessCurrent() const { return getDimlessVelocity() / getQ0(); }
        // T getDimlessCurrent() const { return getResolution() / getLatticeNu() / getQ0(); }

    private:
        T Tehd, Mehd, Cehd, Fe, Ma, Ra, Pr, lx, ly, deltaPhi, deltaT, resolution, rho0;
    };

    template <typename T,
              template <typename NSU> class nsDescriptor,
              template <typename ADELE> class adElecDescriptor,
              template <typename ADTEM> class adTemDescriptor>
    void writeLogFile(
        MyETHDParameters2D_Ma<T, nsDescriptor, adElecDescriptor, adTemDescriptor> const &parameters,
        std::string const &title)
    {
        std::string fullName = global::directories().getLogOutDir() + "plbLogEHD.dat";
        std::ofstream ofile(fullName.c_str());
        ofile << title << "\n\n";
        ofile << "Reynolds number Elec:           Re = " << parameters.getReElec() << "\n";
        ofile << "Reynolds number RBC:            Re = " << parameters.getReRBC() << "\n";
        ofile << "T number:                        T = " << parameters.getTETHD() << "\n";
        ofile << "M number:                        M = " << parameters.getMETHD() << "\n";
        ofile << "C number:                        C = " << parameters.getCETHD() << "\n";
        ofile << "Fe number:                      Fe = " << parameters.getFeETHD() << "\n";
        ofile << "Ra number:                      Ra = " << parameters.getRaETHD() << "\n";
        ofile << "Pr number:                      Pr = " << parameters.getPrETHD() << "\n";
        ofile << "g x beta:                    gBeta = " << parameters.getGBeta() << "\n";

        ofile << "\n";
        ofile << "thermal ade diffusivity:     Alpha = " << parameters.getLatticeAlpha() << "\n";
        ofile << "Kinematic viscosity:            Nu = " << parameters.getLatticeNu() << "\n";
        ofile << "charge ade diffusivity:          D = " << parameters.getLatticeD() << "\n";
        ofile << "Solvent omega:           omega_nse = " << parameters.getOmegaNSE() << "\n";
        ofile << "electric omega:     elec omega_ade = " << parameters.getOmegaElecADE() << "\n";
        ofile << "Temperature omega:   te  omega_ade = " << parameters.getOmegaTemADE() << "\n";

        ofile << "Solvent tau:               tau_nse = " << parameters.getTauNSE() << "\n";
        ofile << "electric tau:         elec tau_ade = " << parameters.getTauElecADE() << "\n";
        ofile << "Temperature tau:       tem tau_ade = " << parameters.getTauTemADE() << "\n";

        ofile << "              Ma :              Ma = " << parameters.getMa() << "\n";
        ofile << "               N :               N = " << parameters.getResolution() << "\n";
        ofile << "              lx :              lx = " << parameters.getLx() << "\n";
        ofile << "              ly :              ly = " << parameters.getLy() << "\n";
        ofile << "              Nx :              Nx = " << parameters.getNx() << "\n";
        ofile << "              Ny :              Ny = " << parameters.getNy() << "\n";

        ofile << "       lattice K :       lattice K = " << parameters.getLatticeK() << "\n";
        ofile << " lattice epsilon : lattice epsilon = " << parameters.getLatticeEpsilon() << "\n";
        ofile << "              q0 :              q0 = " << parameters.getQ0() << "\n";
        ofile << "       delta phi :       delta phi = " << parameters.getDeltaPhi() << "\n";
        ofile << "       delta Tem :       delta Tem = " << parameters.getDeltaT() << "\n";

        ofile << "     getDimlessTime : = " << parameters.getDimlessTime() << "\n";
        ofile << " getDimlessVelocity : = " << parameters.getDimlessVelocity() << "\n";
        ofile << "  getDimlessCurrent : = " << parameters.getDimlessCurrent() << "\n";
    }

    /// Initialization of the force field.
    template <typename T, typename T2, int nDim>
    class getComponentVector
        : public BoxProcessingFunctional2D_ST<T, T2, nDim>
    {
    public:
        getComponentVector(plint dir_) : dir(dir_) {}

        virtual void process(Box2D domain, ScalarField2D<T> &scalar,
                             TensorField2D<T2, nDim> &vector)
        {
            for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            {
                for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                {
                    scalar.get(iX, iY) = vector.get(iX, iY)[dir];
                }
            }
        }
        virtual getComponentVector<T, T2, nDim> *clone() const
        {
            return new getComponentVector<T, T2, nDim>(*this);
        }
        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
        {
            modified[0] = modif::staticVariables;
        }
        virtual BlockDomain::DomainT appliesTo() const
        {
            return BlockDomain::bulkAndEnvelope;
        }

    private:
        plint dir;
    };

    /// Initialization of the force field.
    template <typename T, typename T2, int nDim>
    class getElectricForce
        : public BoxProcessingFunctional2D_ST<T, T2, nDim>
    {
    public:
        getElectricForce(T para_) : para(para_) {}

        virtual void process(Box2D domain, ScalarField2D<T> &scalar,
                             TensorField2D<T2, nDim> &vector)
        {
            for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            {
                for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                {
                    for (plint i = 0; i < nDim; ++i)
                    {
                        // vector.get(iX, iY)[i] *= scalar.get(iX, iY);
                        vector.get(iX, iY)[i] = vector.get(iX, iY)[i] * scalar.get(iX, iY);
                        // vector.get(iX, iY)[i] = 0.0;
                    }
                }
            }
        }
        virtual getElectricForce<T, T2, nDim> *clone() const
        {
            return new getElectricForce<T, T2, nDim>(*this);
        }
        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
        {
            modified[0] = modif::staticVariables;
        }
        virtual BlockDomain::DomainT appliesTo() const
        {
            return BlockDomain::bulkAndEnvelope;
        }

    private:
        T para;
    };

    /// Initialization of the poisson field.
    template <typename T,
              template <typename NS> class nsDescriptor,
              typename T2>
    class getNSSuorceGradLes
        : public BoxProcessingFunctional2D_LS<T, nsDescriptor, T>
    {
    public:
        getNSSuorceGradLes(T C_) : C(C_) {}
        virtual void process(Box2D domain, BlockLattice2D<T, nsDescriptor> &lattice,
                             ScalarField2D<T> &field)
        {
            enum
            {
                nsOffset = nsDescriptor<T>::ExternalField::forceTauBeginsAt
            };
            Dot2D offset = computeRelativeDisplacement(lattice, field);
            for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            {
                for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                {

                    T *s = lattice.get(iX, iY).getExternal(nsOffset);
                    *s = -C * field.get(iX, iY);
                    // *s = lambda * lambda * phi;
                }
            }
        }
        virtual getNSSuorceGradLes<T, nsDescriptor, T> *clone() const
        {
            return new getNSSuorceGradLes<T, nsDescriptor, T>(*this);
        }
        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
        {
            modified[0] = modif::staticVariables;
        }
        virtual BlockDomain::DomainT appliesTo() const
        {
            return BlockDomain::bulkAndEnvelope;
        }

    private:
        T C;
    };

    /// Initialization of the force field.
    template <typename T,
              template <typename NSE> class nseDescriptor,
              typename T2, int nDim>
    class getForceTerm
        : public BoxProcessingFunctional2D_LT<T, nseDescriptor, T, nDim>
    {
    public:
        getForceTerm(T para_) : para(para_) {}

        virtual void process(Box2D domain, BlockLattice2D<T, nseDescriptor> &lattice,
                             TensorField2D<T, nDim> &field)
        {
            enum
            {
                forceOffset = nseDescriptor<T>::ExternalField::forceBeginsAt
            };
            // Dot2D offset = computeRelativeDisplacement(lattice, field);
            for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            {
                for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                {
                    Array<T, nseDescriptor<T>::d> velField = field.get(iX, iY);
                    T vel[] = {para * velField[0], para * velField[1]};
                    lattice.get(iX, iY).setExternalField(forceOffset, 2, vel);
                }
            }
        }
        virtual getForceTerm<T, nseDescriptor, T, nDim> *clone() const
        {
            return new getForceTerm<T, nseDescriptor, T, nDim>(*this);
        }
        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
        {
            modified[0] = modif::staticVariables;
        }
        virtual BlockDomain::DomainT appliesTo() const
        {
            return BlockDomain::bulkAndEnvelope;
        }

    private:
        T para;
    };

    /// Initialization of the temperature field.
    template <typename T,
              template <typename ADE> class adeDescriptor,
              typename T2, int nDim>
    class getAdvectionTerm
        : public BoxProcessingFunctional2D_LT<T, adeDescriptor, T, nDim>
    {
    public:
        getAdvectionTerm(T K_) : K(K_) {}
        virtual void process(Box2D domain, BlockLattice2D<T, adeDescriptor> &lattice,
                             TensorField2D<T, nDim> &field)
        {
            enum
            {
                velOffset = adeDescriptor<T>::ExternalField::velocityBeginsAt
            };

            Dot2D offset = computeRelativeDisplacement(lattice, field);
            // Dot2D absoluteOffset = lattice.getLocation();
            for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            {
                for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                {
                    // // Velocity coupling
                    // T *u = lattice.get(iX + offset.x, iY + offset.y).getExternal(velOffset);
                    // Array<T, adeDescriptor<T>::d> velField(0.0, 0.0);
                    // // // vel = field.get(iX + offset.x, iY + offset.y);
                    // // vel.to_cArray(u);
                    // Array<T, adeDescriptor<T>::d> velField(0.0, 0.01);
                    Array<T, adeDescriptor<T>::d> velField = field.get(iX, iY);
                    T vel[] = {K * velField[0], K * velField[1]};
                    lattice.get(iX, iY).setExternalField(velOffset, 2, vel);
                }
            }
        }

        virtual getAdvectionTerm<T, adeDescriptor, T, nDim> *clone() const
        {
            return new getAdvectionTerm<T, adeDescriptor, T, nDim>(*this);
        }

        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
        {
            modified[0] = modif::staticVariables;
        }

        virtual BlockDomain::DomainT appliesTo() const
        {
            return BlockDomain::bulkAndEnvelope;
        }

    private:
        T K;
    };

    /// Initialization of the poisson field.
    template <typename T,
              template <typename POI> class poiDescriptor,
              typename T2>
    class getPoissonSuorce
        : public BoxProcessingFunctional2D_LS<T, poiDescriptor, T>
    {
    public:
        getPoissonSuorce(T C_) : C(C_) {}
        virtual void process(Box2D domain, BlockLattice2D<T, poiDescriptor> &lattice,
                             ScalarField2D<T> &field)
        {
            enum
            {
                // velOffset = poiDescriptor<T>::ExternalField::velocityBeginsAt
                PoiOffset = poiDescriptor<T>::ExternalField::scalarBeginsAt
            };
            Dot2D offset = computeRelativeDisplacement(lattice, field);
            // Dot2D absoluteOffset = lattice.getLocation();
            for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            {
                for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                {
                    // T *velField = field.get(iX, iY);
                    // lattice.get(iX, iY).setExternalField(PoiOffset, 1, velField);

                    T *s = lattice.get(iX, iY).getExternal(PoiOffset);
                    *s = -C * field.get(iX, iY);
                    // *s = lambda * lambda * phi;
                }
            }
        }
        virtual getPoissonSuorce<T, poiDescriptor, T> *clone() const
        {
            return new getPoissonSuorce<T, poiDescriptor, T>(*this);
        }
        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
        {
            modified[0] = modif::staticVariables;
        }
        virtual BlockDomain::DomainT appliesTo() const
        {
            return BlockDomain::bulkAndEnvelope;
        }

    private:
        T C;
    };

    /// Initialization of the temperature field.
    template <typename T, template <typename POI> class poissonDescriptor>
    struct InitPoissonProcessor2D : public BoxProcessingFunctional2D_L<T, poissonDescriptor>
    {
        InitPoissonProcessor2D(plint nx, plint ny, T deltaPhi_) : NX(nx), NY(ny), deltaPhi(deltaPhi_)
        {
        }
        virtual void process(Box2D domain, BlockLattice2D<T, poissonDescriptor> &lattice)
        {
            Dot2D absoluteOffset = lattice.getLocation();
            double a = 1.4882;
            double b = 5.539e-3;
            double c = 1.0005;
            for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            {
                for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                {
                    plint absoluteY = absoluteOffset.y + iY;
                    // T phi = deltaPhi - deltaPhi / (T)(NY - 1) * (T)absoluteY;
                    double yStar = absoluteY / T(NY - 1);
                    T ini{0.0};
                    ini = -2.0 / 3 * a * pow(yStar + b, 1.5) + c;
                    Array<T, poissonDescriptor<T>::d> jEq(0.0, 0.0);

                    ini = ini * deltaPhi;
                    lattice.get(iX, iY).defineDensity(ini);
                    iniCellAtEquilibriumPoisson(lattice.get(iX, iY), ini);

                    if (absoluteY == 0)
                    {
                        T ini{deltaPhi};
                        lattice.get(iX, iY).defineDensity(ini);
                        iniCellAtEquilibriumPoisson(lattice.get(iX, iY), ini);
                    }
                    if (absoluteY == NY - 1)
                    {
                        T ini{0};
                        lattice.get(iX, iY).defineDensity(ini);
                        iniCellAtEquilibriumPoisson(lattice.get(iX, iY), ini);
                    }
                }
            }
        }

        virtual InitPoissonProcessor2D<T, poissonDescriptor> *clone() const
        {
            return new InitPoissonProcessor2D<T, poissonDescriptor>(*this);
        }

        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
        {
            modified[0] = modif::staticVariables;
        }

        virtual BlockDomain::DomainT appliesTo() const
        {
            return BlockDomain::bulkAndEnvelope;
        }

    private:
        plint NX;
        plint NY;
        T deltaPhi;
    };

    /// Initialization of the temperature field.
    template <typename T, template <typename POI> class adeDescriptor>
    struct InitADEProcessor2D : public BoxProcessingFunctional2D_L<T, adeDescriptor>
    {
        InitADEProcessor2D(plint nx, plint ny, T q0_) : NX(nx), NY(ny), q0(q0_)
        {
        }

        virtual void process(Box2D domain, BlockLattice2D<T, adeDescriptor> &lattice)
        {
            Dot2D absoluteOffset = lattice.getLocation();
            double a = 1.4882;
            double b = 5.539e-3;

            for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            {
                for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                {
                    plint absoluteY = absoluteOffset.y + iY;
                    double yStar = absoluteY / T(NY);
                    T ini{0.0};
                    ini = a / 2 / sqrt(yStar + b) / 10 * q0;

                    Array<T, adeDescriptor<T>::d> jEq(0.0, 0.0);
                    lattice.get(iX, iY).defineDensity(ini);

                    iniCellAtEquilibrium(lattice.get(iX, iY), ini, jEq);

                    if (absoluteY == 0)
                    {
                        T ini{q0};
                        lattice.get(iX, iY).defineDensity(ini);
                        iniCellAtEquilibrium(lattice.get(iX, iY), ini, jEq);
                    }
                }
            }
        }

        virtual InitADEProcessor2D<T, adeDescriptor> *clone() const
        {
            return new InitADEProcessor2D<T, adeDescriptor>(*this);
        }

        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
        {
            modified[0] = modif::staticVariables;
        }

        virtual BlockDomain::DomainT appliesTo() const
        {
            return BlockDomain::bulkAndEnvelope;
        }

    private:
        plint NX;
        plint NY;
        T q0;
    };

    /// Initialization of the velocity field.
    template <typename T, template <typename NSE> class nsDescriptor>
    struct InitNSEProcessor2D : public BoxProcessingFunctional2D_L<T, nsDescriptor>
    {
        InitNSEProcessor2D(plint nx, plint ny)
            : NX(nx), NY(ny)
        {
        }

        virtual void process(Box2D domain, BlockLattice2D<T, nsDescriptor> &lattice)
        {
            Dot2D absoluteOffset = lattice.getLocation();

            double pi = 3.1415926859;
            double kx = 2 * pi / (NX - 1);
            double ky = 2 * pi / (NY - 1);
            double a = 1e-4;

            for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            {
                for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                {
                    plint absoluteX = absoluteOffset.x + iX;
                    plint absoluteY = absoluteOffset.y + iY;
                    double xStar = absoluteX / T(NY - 1);
                    double yStar = absoluteY / T(NY - 1);

                    T uIni = -a * sin(ky * absoluteY) * cos(kx * absoluteX);
                    T vIni = a * cos(ky * absoluteY) * sin(kx * absoluteX);

                    Array<T, nsDescriptor<T>::d> veloEq(uIni, vIni);
                    T rho0{1.};

                    lattice.get(iX, iY).defineVelocity(veloEq);
                    lattice.get(iX, iY).defineDensity(rho0);

                    iniCellAtEquilibrium(lattice.get(iX, iY), rho0, veloEq);
                }
            }
        }

        virtual InitNSEProcessor2D<T, nsDescriptor> *clone() const
        {
            return new InitNSEProcessor2D<T, nsDescriptor>(*this);
        }

        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
        {
            modified[0] = modif::staticVariables;
        }

        virtual BlockDomain::DomainT appliesTo() const
        {
            return BlockDomain::bulkAndEnvelope;
        }

    private:
        plint NX;
        plint NY;
    };

    /// Initialization of the temperature field.
    template <typename T, template <typename POI> class adeDescriptor>
    struct InitADETemProcessor2D : public BoxProcessingFunctional2D_L<T, adeDescriptor>
    {
        InitADETemProcessor2D(plint nx, plint ny, T temBottom_, T temTop_)
            : NX(nx), NY(ny), temBottom(temBottom_), temTop(temTop_)
        {
        }

        virtual void process(Box2D domain, BlockLattice2D<T, adeDescriptor> &lattice)
        {
            Dot2D absoluteOffset = lattice.getLocation();

            T deltaT = temBottom - temTop;

            for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            {
                for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                {
                    plint absoluteY = absoluteOffset.y + iY;
                    double yStar = absoluteY / T(NY);
                    T ini{0.0};
                    // ini = a / 2 / sqrt(yStar + b) / 10 * q0;
                    ini = temBottom - deltaT / (T)(NY - 1) * (T)absoluteY;

                    Array<T, adeDescriptor<T>::d> jEq(0.0, 0.0);
                    lattice.get(iX, iY).defineDensity(ini);

                    iniCellAtEquilibrium(lattice.get(iX, iY), ini, jEq);

                    if (absoluteY == 0)
                    {
                        T ini{temBottom};
                        lattice.get(iX, iY).defineDensity(ini);
                        iniCellAtEquilibrium(lattice.get(iX, iY), ini, jEq);
                    }
                }
            }
        }

        virtual InitADETemProcessor2D<T, adeDescriptor> *clone() const
        {
            return new InitADETemProcessor2D<T, adeDescriptor>(*this);
        }

        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
        {
            modified[0] = modif::staticVariables;
        }

        virtual BlockDomain::DomainT appliesTo() const
        {
            return BlockDomain::bulkAndEnvelope;
        }

    private:
        plint NX;
        plint NY;
        T temBottom;
        T temTop;
    };

    /**
     * Class for the coupling between a Navier-Stokes (NS) lattice and an
     * Advection-Diffusion (AD) lattice in the boussinesq approximation.
     */
    template <typename T1, template <typename U> class TemperatureDescriptor,
              typename T2, int nDim>
    class MyBoussinesqThermalProcessor2D
        : public BoxProcessingFunctional2D_LT<T1, TemperatureDescriptor, T2, nDim>
    {
    public:
        MyBoussinesqThermalProcessor2D(T1 gravity_, T1 T0_, T1 deltaTemp_,
                                       Array<T2, nDim> dir_)
            : gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_), dir(dir_)
        {
            // We normalize the direction of the force vector.
            T2 normDir = std::sqrt(normSqr(dir));
            for (pluint iD = 0; iD < nDim; ++iD)
            {
                dir[iD] /= normDir;
            }
        }

        T2 normSqr(const Array<T2, nDim> u1)
        {
            T2 result = T2();
            for (int iD = 0; iD < nDim; ++iD)
            {
                result += u1[iD] * u1[iD];
            }
            return result;
        }

        virtual void process(Box2D domain, BlockLattice2D<T1, TemperatureDescriptor> &temperature,
                             TensorField2D<T2, nDim> &vector)
        {
            Dot2D offset = computeRelativeDisplacement(vector, temperature);
            for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            {
                for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                {
                    // Computation of the Boussinesq force
                    // Temperature is the order-0 moment of the advection-diffusion lattice.
                    //   You can compute it with the method computeDensity().
                    T1 localTemperature = temperature.get(iX + offset.x, iY + offset.y).computeDensity();
                    // T1 localTemperature = temperature.get(iX, iY).computeDensity();
                    const T1 diffT = localTemperature - T0;

                    // T1 force[] = {0, -0.000241667 * diffT};
                    T1 force[] = {dir[0] * gravity * diffT, dir[1] * gravity * diffT};

                    // for (pluint iD = 0; iD < nDim; ++iD)
                    // {
                    // force[iD] = -dir[iD] * diffT * gravity;
                    // force[iD] = -dir[iD] * gravity;
                    // }

                    vector.get(iX, iY).from_cArray(force);
                    // Array<T1, nDim> vecTemp = {dir[iD] * gravity, -dir[iD] * gravity};
                    // vector.get(iX, iY) = vecTemp;
                }
            }
        }
        virtual MyBoussinesqThermalProcessor2D<T1, TemperatureDescriptor, T2, nDim> *clone() const
        {
            return new MyBoussinesqThermalProcessor2D<T1, TemperatureDescriptor, T2, nDim>(*this);
        }
        virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
        {
            modified[0] = modif::staticVariables;
            modified[1] = modif::staticVariables;
        }
        virtual BlockDomain::DomainT appliesTo() const
        {
            return BlockDomain::bulkAndEnvelope;
        }

    private:
        T1 gravity, T0, deltaTemp;
        Array<T2, nDim> dir;
    };

    /// MRT2D Advection-Diffusion dynamics
    /** This approach contains a slight error in the diffusion
     *  term.
     */

    template <typename T, template <typename U> class Descriptor>
    class AdvectionDiffusionMRT2Ddynamics : public AdvectionDiffusionDynamics<T, Descriptor>
    {
    public:
        /// Constructor
        AdvectionDiffusionMRT2Ddynamics(T omega_)
            : AdvectionDiffusionDynamics<T, Descriptor>(omega_)
        // AdvectionDiffusionMRT2Ddynamics()
        {
            // omega = omega_;
            // // omega = 0.7;
            // T tau = 1.0 / omega_;
            // s0 = 1.0;
            // s3 = s4 = omega_;
            // s1 = s2 = 1.0 / (1.0 / (6 * (tau - 0.5)) + 0.5);
        }
        AdvectionDiffusionMRT2Ddynamics(HierarchicUnserializer &unserializer)
            : AdvectionDiffusionDynamics<T, Descriptor>(T())
        {
            this->unserialize(unserializer);
        }
        /// Clone the object on its dynamic type.
        virtual AdvectionDiffusionMRT2Ddynamics<T, Descriptor> *clone() const
        {
            return new AdvectionDiffusionMRT2Ddynamics<T, Descriptor>(*this);
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

            // T feq0 = computeEquilibrium(0, rhoBar, jEq, 0, 0);
            // T feq1 = computeEquilibrium(3, rhoBar, jEq, 0, 0);
            // T feq2 = computeEquilibrium(4, rhoBar, jEq, 0, 0);
            // T feq3 = computeEquilibrium(1, rhoBar, jEq, 0, 0);
            // T feq4 = computeEquilibrium(2, rhoBar, jEq, 0, 0);

            T w0 = Descriptor<T>::t[0];
            T wc = Descriptor<T>::t[1];
            T cs2 = Descriptor<T>::cs2;

            T meq0 = rhoBar * (w0 + 4 * wc);
            T meq1 = 2 * jEq[0] * wc / cs2;
            T meq2 = 2 * jEq[1] * wc / cs2;
            T meq3 = 0;
            T meq4 = 4 * rhoBar * (-w0 + wc);

            // T meq0 = feq0 + feq1 + feq2 + feq3 + feq4;
            // T meq1 = feq1 - feq3;
            // T meq2 = feq2 - feq4;
            // T meq3 = feq1 - feq2 + feq3 - feq4;
            // T meq4 = -4 * feq0 + feq1 + feq2 + feq3 + feq4;

            T omg = AdvectionDiffusionDynamics<T, Descriptor>::getOmega();
            T tau = 1.0 / omg;
            T s0 = 1.0;
            T s1 = omg;
            T s2 = omg;
            T s3 = 1.0 / (1.0 / (6 * (tau - 0.5)) + 0.5);
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
    int AdvectionDiffusionMRT2Ddynamics<T, Descriptor>::id =
        meta::registerGeneralDynamics<T, Descriptor, AdvectionDiffusionMRT2Ddynamics<T, Descriptor>>("AdvectionDiffusion_MRT2D");

    /// Implementation of the MRT collision step
    template <typename T, template <typename U> class Descriptor>
    class MyGuoD2Q9ExternalForceMRTdynamics01 : public ExternalForceDynamics<T, Descriptor>
    {
    public:
        /* *************** Construction / Destruction ************************ */
        MyGuoD2Q9ExternalForceMRTdynamics01(T omega_)
            : ExternalForceDynamics<T, Descriptor>(omega_)
        {
        }

        MyGuoD2Q9ExternalForceMRTdynamics01(HierarchicUnserializer &unserializer)
            : ExternalForceDynamics<T, Descriptor>(T())
        {
            this->unserialize(unserializer);
        }

        /// Clone the object on its dynamic type.
        virtual MyGuoD2Q9ExternalForceMRTdynamics01<T, Descriptor> *clone() const
        {
            return new MyGuoD2Q9ExternalForceMRTdynamics01<T, Descriptor>(*this);
        }

        /// Return a unique ID for this class.
        virtual int getId() const
        {
            return id;
        }

        /* *************** Collision and Equilibrium ************************* */

        /// Implementation of the collision step
        virtual void collide(Cell<T, Descriptor> &cell,
                             BlockStatistics &statistics_)
        {
            T rhoBar = momentTemplates<T, Descriptor>::get_rhoBar(cell);
            Array<T, Descriptor<T>::d> u;
            this->computeVelocity(cell, u);
            Array<T, Descriptor<T>::d> force;
            force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));

            const T rhoFull = Descriptor<T>::fullRho(rhoBar);
            const T rho = rhoBar;
            const T ux = u[0];
            const T uy = u[1];
            const T fx = force[0];
            const T fy = force[1];

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

            // momentsEq[0] = rhoBar;
            // momentsEq[1] = (T)3*jSqr*invRho-2*rhoBar;
            // momentsEq[2] = -(T)3*jSqr*invRho+rhoBar;
            // momentsEq[3] = j[0];
            // momentsEq[4] = -j[0];
            // momentsEq[5] = j[1];
            // momentsEq[6] = -j[1];
            // momentsEq[7] = (j[0]*j[0]-j[1]*j[1])*invRho;
            // momentsEq[8] = j[1]*j[0]*invRho;

            const T meq0 = 1.0 * rhoBar;
            const T meq1 = 3.0 * rhoFull * ux * ux + 3.0 * rhoFull * uy * uy - 2.0 * rhoBar;
            const T meq2 = -3.0 * rhoFull * ux * ux - 3.0 * rhoFull * uy * uy + rhoBar;
            const T meq3 = 1.0 * rhoFull * ux;
            const T meq4 = -1.0 * rhoFull * ux;
            const T meq5 = 1.0 * rhoFull * uy;
            const T meq6 = -1.0 * rhoFull * uy;
            const T meq7 = 1.0 * rhoFull * (ux * ux - uy * uy);
            const T meq8 = 1.0 * rhoFull * ux * uy;

            const T mF0 = 0;
            const T mF1 = 6.0 * fx * ux + 6.0 * fy * uy;
            const T mF2 = -6.0 * fx * ux - 6.0 * fy * uy;
            const T mF3 = 1.0 * fx;
            const T mF4 = -1.0 * fx;
            const T mF5 = 1.0 * fy;
            const T mF6 = -1.0 * fy;
            const T mF7 = 2.0 * fx * ux - 2.0 * fy * uy;
            const T mF8 = 1.0 * fx * uy + 1.0 * fy * ux;

            const T omega = ExternalForceDynamics<T, Descriptor>::getOmega();
            const T tau = 1. / omega;
            const T omg1 = 8 * (2 * tau - 1) / (8 * tau - 1);

            T s0 = 0;

            T s3 = omega;
            T s5 = omega;
            T s1 = omega;

            T s2 = omega;
            T s4 = omg1;
            T s6 = omg1;
            T s7 = omega;
            T s8 = omega;

            T mcp0 = m0 * (1.0 - s0) + s0 * meq0 + (1 - s0 / 2) * mF0;
            T mcp1 = m1 * (1.0 - s1) + s1 * meq1 + (1 - s1 / 2) * mF1;
            T mcp2 = m2 * (1.0 - s2) + s2 * meq2 + (1 - s2 / 2) * mF2;
            T mcp3 = m3 * (1.0 - s3) + s3 * meq3 + (1 - s3 / 2) * mF3;
            T mcp4 = m4 * (1.0 - s4) + s4 * meq4 + (1 - s4 / 2) * mF4;
            T mcp5 = m5 * (1.0 - s5) + s5 * meq5 + (1 - s5 / 2) * mF5;
            T mcp6 = m6 * (1.0 - s6) + s6 * meq6 + (1 - s6 / 2) * mF6;
            T mcp7 = m7 * (1.0 - s7) + s7 * meq7 + (1 - s7 / 2) * mF7;
            T mcp8 = m8 * (1.0 - s8) + s8 * meq8 + (1 - s8 / 2) * mF8;

            cell[0] = mcp0 / 9 - mcp1 / 9 + mcp2 / 9;
            cell[1] = mcp0 / 9 + mcp1 / 18 + mcp2 / 36 - mcp3 / 6 - mcp4 / 12 + mcp5 / 6 + mcp6 / 12 - mcp8 / 4;
            cell[2] = mcp0 / 9 - mcp1 / 36 - mcp2 / 18 - mcp3 / 6 + mcp4 / 6 + mcp7 / 4;
            cell[3] = mcp0 / 9 + mcp1 / 18 + mcp2 / 36 - mcp3 / 6 - mcp4 / 12 - mcp5 / 6 - mcp6 / 12 + mcp8 / 4;
            cell[4] = mcp0 / 9 - mcp1 / 36 - mcp2 / 18 - mcp5 / 6 + mcp6 / 6 - mcp7 / 4;
            cell[5] = mcp0 / 9 + mcp1 / 18 + mcp2 / 36 + mcp3 / 6 + mcp4 / 12 - mcp5 / 6 - mcp6 / 12 - mcp8 / 4;
            cell[6] = mcp0 / 9 - mcp1 / 36 - mcp2 / 18 + mcp3 / 6 - mcp4 / 6 + mcp7 / 4;
            cell[7] = mcp0 / 9 + mcp1 / 18 + mcp2 / 36 + mcp3 / 6 + mcp4 / 12 + mcp5 / 6 + mcp6 / 12 + mcp8 / 4;
            cell[8] = mcp0 / 9 - mcp1 / 36 - mcp2 / 18 + mcp5 / 6 - mcp6 / 6 - mcp7 / 4;
        }

        // virtual void collideExternal(
        //     Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar, BlockStatistics &stat)
        // {
        // }

        virtual void computeVelocity(
            Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const
        {
            T rhoBar;
            Array<T, Descriptor<T>::d> force, j;
            momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
            force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));

            T invRho = Descriptor<T>::invRho(rhoBar);
            T rhoFull = Descriptor<T>::fullRho(rhoBar);
            for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
            {
                u[iD] = j[iD] * invRho + force[iD] / (T)2 * invRho;
            }
        }

        virtual void computeVelocityExternal(
            Cell<T, Descriptor> const &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j,
            Array<T, Descriptor<T>::d> &u) const
        {
            Array<T, Descriptor<T>::d> force;
            force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));

            T invRho = Descriptor<T>::invRho(rhoBar);
            for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
            {
                u[iD] = j[iD] * invRho + force[iD] / (T)2 * invRho;
            }
        }

        /// Compute equilibrium distribution function
        virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j,
                                     T jSqr, T thetaBar = T()) const
        {
            T invRho = Descriptor<T>::invRho(rhoBar);
            return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
        }

    private:
        static int id;
    };

    template <typename T, template <typename U> class Descriptor>
    int MyGuoD2Q9ExternalForceMRTdynamics01<T, Descriptor>::id =
        meta::registerGeneralDynamics<T, Descriptor, MyGuoD2Q9ExternalForceMRTdynamics01<T, Descriptor>>("MyGuoD2Q9_ExternalForceMRTdynamics01");

}

#endif
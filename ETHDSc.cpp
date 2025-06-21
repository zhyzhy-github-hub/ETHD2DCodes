
#include "palabos2D.h"
#include "palabos2D.hh"
#include "poissonComplex/poisson.h"
#include "poissonComplex/poisson.hh"
#include "NeeBc/neeBC.h"
#include "ETHDpara.h"
#include "multiCs.h"

#include <cstdlib>
#include <iostream>

using namespace plb;
using namespace std;

typedef double T;

#define NSDESCRIPTOR descriptors::ForceVelocitydD2Q9ForceTauLesDescriptor
// #define NSDYNAMICS MyGuoD2Q9ExternalForceLesMRTdynamics01MultiCs
#define NSDYNAMICS MyGuoD2Q9ExternalForceLesMRTdynamics02_direct_MultiCs

#define ADESCRIPTOR descriptors::AdvectionDiffusionMultiLesD2Q5Descriptor
// #define ADESCRIPTOR descriptors::AdvectionDiffusionD2Q5Descriptor
// #define ADYNAMICS AdvectionDiffusionWithMultiNSLesMRT2Ddynamics
#define ADYNAMICS AdvectionDiffusionLesBGKdynamics

#define POISCRIPTOR descriptors::poissonChai08WithScalarSourceD2Q5Descriptor
#define POIDYNAMICS PoissonChai08WithSourceBGKdynamics

// 是否定义周期边界，如果有这个，那么就是周期边界条件

// #define PERIOD_COMP

// #define IF_TEST
#define lesModel

const T cc = ccDefine;

plint continueSimulation = 0;

plint outStatistics = 20;
plint vtkOutNum = outStatistics * 2000 / cc;
plint continueIter = outStatistics * 3000;
plint maxIter = 1000000;

/// Poisson 方程内迭代
/// @brief
/// @param poiLattice Poisson lattice
/// @param phiOld the phi field at last time step
/// @param phiSub phi field at this time step
/// @param chargeOld charge density at source term
/// @param C parameter C
/// @return

T poissonIter(MultiBlockLattice2D<T, POISCRIPTOR> &poiLattice,
              MultiScalarField2D<T> phiOld,
              MultiScalarField2D<T> phiSub,
              MultiScalarField2D<T> chargeOld,
              T C)
{

    plint poissonInt = 0;
    T phiMax = 1.0;

    applyProcessingFunctional(
        new getPoissonSuorce<T, POISCRIPTOR, T>(C),
        poiLattice.getBoundingBox(),
        poiLattice, chargeOld);

    for (poissonInt = 0; poissonInt < 1000000; ++poissonInt)
    {
        phiOld = (*extractSubDomain(*computeDensity(poiLattice), phiOld.getBoundingBox()));

        poiLattice.collideAndStream();

        phiSub = *computeAbsoluteValue(*subtract(phiOld, *computeDensity(poiLattice)));
        phiMax = computeMax(phiSub);
        // pcout << "Iter num = " << poissonInt << ", phiMax ... " << phiMax << endl;

        if (phiMax < 1e-6)
        {
            // pcout << "Iter num = " << poissonInt << ", phiMax ... " << phiMax << endl;
            break;
        }
    }
    return phiMax;
}

/// @brief Iteration of Nernst Planck equations
/// @param lattice NP equs
/// @param boundaryCondition
/// @param nx
/// @param ny
/// @param q0 The charge density at the bottom plate
/// @param initial
void adeChargeSetup(
    MultiBlockLattice2D<T, ADESCRIPTOR> &lattice,
    OnLatticeAdvectionDiffusionBoundaryCondition2D<T, ADESCRIPTOR> &boundaryCondition,
    plint nx, plint ny, T q0, plint initial = 0)
{
    Box2D top(0, nx - 1, ny - 1, ny - 1);
    Box2D bottom(0, nx - 1, 0, 0);

    // 下板是Dirichlet边界条件，定值电荷密度

    integrateProcessingFunctional(new NeeFlatDirichletBoundaryFunctional2D<T, ADESCRIPTOR, 1, -1>(q0), bottom, lattice);

    // 上班零梯度

    integrateProcessingFunctional(new NeeFlatAdiabaticBoundaryFunctional2D<T, ADESCRIPTOR, 1, 1>, top, lattice);

#ifndef PERIOD_COMP
    Box2D left(0, 0, 1, ny - 2);
    Box2D right(nx - 1, nx - 1, 1, ny - 2);

    // 左右绝缘

    integrateProcessingFunctional(new NeeFlatAdiabaticBoundaryFunctional2D<T, ADESCRIPTOR, 0, 1>, right, lattice);
    integrateProcessingFunctional(new NeeFlatAdiabaticBoundaryFunctional2D<T, ADESCRIPTOR, 0, -1>, left, lattice);
#endif

    //--------------------------------------------------------------------------------
    if (!initial)
    {
        Array<T, POISCRIPTOR<T>::d>
            jEq(0.0, 0.0);
        initializeAtEquilibrium(lattice, lattice.getBoundingBox(), (T)0.);

        applyProcessingFunctional(
            new InitADEProcessor2D<T, ADESCRIPTOR>(nx, ny, q0),
            lattice.getBoundingBox(), lattice);
    }
    lattice.initialize();
}

/// @brief
/// @param lattice temperature lattice
/// @param boundaryCondition
/// @param nx
/// @param ny
/// @param temBottom the temperature at below plate
/// @param temTop temperature at upper plate
/// @param initial

void adeTemperatureSetup(
    MultiBlockLattice2D<T, ADESCRIPTOR> &lattice,
    OnLatticeAdvectionDiffusionBoundaryCondition2D<T, ADESCRIPTOR> &boundaryCondition,
    plint nx, plint ny, T temBottom, T temTop, plint initial = 0)
{
    Box2D top(0, nx - 1, ny - 1, ny - 1);
    Box2D bottom(0, nx - 1, 0, 0);

    integrateProcessingFunctional(new NeeFlatDirichletBoundaryFunctional2D<T, ADESCRIPTOR, 1, -1>(temBottom), bottom, lattice);
    integrateProcessingFunctional(new NeeFlatDirichletBoundaryFunctional2D<T, ADESCRIPTOR, 1, 1>(temTop), top, lattice);

#ifndef PERIOD_COMP
    Box2D left(0, 0, 1, ny - 2);
    Box2D right(nx - 1, nx - 1, 1, ny - 2);
    integrateProcessingFunctional(new NeeFlatAdiabaticBoundaryFunctional2D<T, ADESCRIPTOR, 0, 1>, right, lattice);
    integrateProcessingFunctional(new NeeFlatAdiabaticBoundaryFunctional2D<T, ADESCRIPTOR, 0, -1>, left, lattice);
#endif

    //--------------------------------------------------------------------------------
    if (!initial)
    {
        Array<T, POISCRIPTOR<T>::d>
            jEq(0.0, 0.0);
        initializeAtEquilibrium(lattice, lattice.getBoundingBox(), (T)0.);

        applyProcessingFunctional(
            new InitADETemProcessor2D<T, ADESCRIPTOR>(nx, ny, temBottom, temTop),
            lattice.getBoundingBox(), lattice);
    }
    // lattice.initialize();
}

void nseSetup(
    MultiBlockLattice2D<T, NSDESCRIPTOR> &lattice,
    OnLatticeBoundaryCondition2D<T, NSDESCRIPTOR> &boundaryCondition,
    plint nx, plint ny, plint initial = 0)
{
#ifndef PERIOD_COMP

    Box2D top(1, nx - 2, ny - 1, ny - 1);
    Box2D bot(1, nx - 2, 0, 0);
    Box2D lef(0, 0, 1, ny - 2);
    Box2D rig(nx - 1, nx - 1, 1, ny - 2);

    Box2D LB(0, 0, 0, 0);
    Box2D LT(0, 0, ny - 1, ny - 1);
    Box2D RB(nx - 1, nx - 1, 0, 0);
    Box2D RT(nx - 1, nx - 1, ny - 1, ny - 1);

    defineDynamics(lattice, LB, new BounceBack<T, NSDESCRIPTOR>);
    defineDynamics(lattice, LT, new BounceBack<T, NSDESCRIPTOR>);
    defineDynamics(lattice, RB, new BounceBack<T, NSDESCRIPTOR>);
    defineDynamics(lattice, RT, new BounceBack<T, NSDESCRIPTOR>);
    defineDynamics(lattice, top, new BounceBack<T, NSDESCRIPTOR>);
    defineDynamics(lattice, bot, new BounceBack<T, NSDESCRIPTOR>);
    defineDynamics(lattice, lef, new BounceBack<T, NSDESCRIPTOR>);
    defineDynamics(lattice, rig, new BounceBack<T, NSDESCRIPTOR>);
#else

    Box2D top(0, nx - 1, ny - 1, ny - 1);
    Box2D bot(0, nx - 1, 0, 0);
    defineDynamics(lattice, top, new BounceBack<T, NSDESCRIPTOR>);
    defineDynamics(lattice, bot, new BounceBack<T, NSDESCRIPTOR>);

#endif

    if (!initial)
    {
        initializeAtEquilibrium(lattice, lattice.getBoundingBox(), (T)1., Array<T, 2>((T)0., (T)0.));
    }
    lattice.initialize();
}

/// @brief
/// @param lattice
/// @param boundaryCondition
/// @param nx
/// @param ny
/// @param deltaPhi
/// @param initial

void poissonSetup(
    MultiBlockLattice2D<T, POISCRIPTOR> &lattice,
    OnLatticeAdvectionDiffusionBoundaryCondition2D<T, POISCRIPTOR> &boundaryCondition,
    plint nx, plint ny, T deltaPhi, plint initial = 0)
{

    Box2D top(0, nx - 1, ny - 1, ny - 1);
    Box2D bottom(0, nx - 1, 0, 0);
    boundaryCondition.addTemperatureBoundary1N(bottom, lattice);
    boundaryCondition.addTemperatureBoundary1P(top, lattice);

#ifndef PERIOD_COMP
    Box2D left(0, 0, 1, ny - 2);
    Box2D right(nx - 1, nx - 1, 1, ny - 2);
    integrateProcessingFunctional(new NeeFlatAdiabaticBoundaryFunctional2D<T, POISCRIPTOR, 0, 1>, right, lattice);
    integrateProcessingFunctional(new NeeFlatAdiabaticBoundaryFunctional2D<T, POISCRIPTOR, 0, -1>, left, lattice);
#endif

    //--------------------------------------------------------------------------------
    if (!initial)
    {
        Array<T, POISCRIPTOR<T>::d>
            jEq(0.0, 0.0);

        initializeAtEquilibrium(lattice, lattice.getBoundingBox(), (T)0.);

        applyProcessingFunctional(
            new InitPoissonProcessor2D<T, POISCRIPTOR>(nx, ny, deltaPhi),
            lattice.getBoundingBox(), lattice);
    }
    lattice.initialize();
}

void computeDissipation(plint i,
                        MultiBlockLattice2D<T, NSDESCRIPTOR> &nsLattice,
                        MultiBlockLattice2D<T, ADESCRIPTOR> &adLattice,
                        MultiBlockLattice2D<T, ADESCRIPTOR> &adTemLattice,
                        MultiBlockLattice2D<T, POISCRIPTOR> &poiLattice,
                        const T nx, const T ny,
                        const MyETHDParameters2D_Nv_Sc<T, NSDESCRIPTOR, ADESCRIPTOR, ADESCRIPTOR> parameters,
                        string prefix)
{
    pcout << "i = " << i << ", computing dissipation ..." << endl;
    plint iT = i;

    Box2D domain = nsLattice.getBoundingBox();

    MultiScalarField2D<T> uxDomain = (*computeVelocityComponent(nsLattice, domain, 0));
    MultiScalarField2D<T> uyDomain = (*computeVelocityComponent(nsLattice, domain, 1));
    MultiScalarField2D<T> dUxDx(nx, ny);
    MultiScalarField2D<T> dUxDy(nx, ny);
    MultiScalarField2D<T> dUyDx(nx, ny);
    MultiScalarField2D<T> dUyDy(nx, ny);

    MultiScalarField2D<T> nuTur = (*computeTemperature(nsLattice, domain));
    MultiScalarField2D<T> nuEff = *add(nuTur, parameters.getLatticeNu());
    MultiScalarField2D<T> density = (*computeDensity(nsLattice, domain));
    MultiScalarField2D<T> nuEffRho = *multiply(density, nuEff);

    applyProcessingFunctional(new BoxXderivativeFunctional2D<T>(), uxDomain.getBoundingBox(), uxDomain, dUxDx, 1);
    applyProcessingFunctional(new BoxYderivativeFunctional2D<T>(), uxDomain.getBoundingBox(), uxDomain, dUxDy, 1);
    applyProcessingFunctional(new BoxXderivativeFunctional2D<T>(), uyDomain.getBoundingBox(), uyDomain, dUyDx, 1);
    applyProcessingFunctional(new BoxYderivativeFunctional2D<T>(), uyDomain.getBoundingBox(), uyDomain, dUyDy, 1);

    T dUxDx2 = computeAverage(*multiply(*multiply(dUxDx, dUxDx, dUxDx.getBoundingBox()), nuEff));
    T dUxDy2 = computeAverage(*multiply(*multiply(dUxDy, dUxDy, dUxDy.getBoundingBox()), nuEff));
    T dUyDx2 = computeAverage(*multiply(*multiply(dUyDx, dUyDx, dUyDx.getBoundingBox()), nuEff));
    T dUyDy2 = computeAverage(*multiply(*multiply(dUyDy, dUyDy, dUyDy.getBoundingBox()), nuEff));

    T dUxDyDUyDx = computeAverage(*multiply(*multiply(dUxDy, dUyDx, dUxDy.getBoundingBox()), nuEff));
    T dUxDxDUyDy = computeAverage(*multiply(*multiply(dUxDx, dUyDy, dUxDx.getBoundingBox()), nuEff));

    T nuEffAve = computeAverage(nuEff);
    T nuEffRhoAve = computeAverage(nuEffRho);

    T kiniteEnergyAve = computeAverageEnergy(nsLattice);

    plb_ofstream dissFile;
    std::ostringstream name;
    name << prefix
         << "diss_" << cc << ".dat";
    // name << "T_" << std::setfill('0') << std::setw(6) << parameters.getTETHD()
    //      << "_M_" << std::setw(3) << parameters.getMETHD()
    //      << "_C_" << std::setw(3) << parameters.getMETHD()
    //      << "_ny_" << std::setw(3) << ny << "_nx_" << std::setw(3) << nx
    //      << "diss_" << cc << ".dat";

    if (i == 0)
    {
        dissFile.open(name.str().c_str());
    }
    else
    {
        dissFile.open(name.str().c_str(), std::ostream::app);
    }

    T Cnn = nx * ny / (nx - 0) / (ny - 0);
    T Cepsilon = 1.;
    T Cenergy = 1.;
    T tdd = i * parameters.getDimlessTime();

    dissFile << i
             << " " << tdd
             << " " << nuEffAve
             << " " << nuEffRhoAve
             << " " << dUxDx2 / Cepsilon * Cnn
             << " " << dUxDy2 / Cepsilon * Cnn
             << " " << dUyDx2 / Cepsilon * Cnn
             << " " << dUyDy2 / Cepsilon * Cnn
             << " " << dUxDyDUyDx / Cepsilon * Cnn
             << " " << dUxDxDUyDy / Cepsilon * Cnn
             << " " << (dUxDx2 + dUxDy2 //
                        + dUyDx2 + dUyDy2) /
                           Cepsilon * Cnn
             << " " << kiniteEnergyAve / Cenergy * Cnn
             << "\n ";
    dissFile.close();
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmpLES2D/");

    global::timer("simTime").start();

    plint NyIn;
    T tauPhi = 1.5;
    T TEthd;
    T RaEthd;
    T PrEthd;
    T Ma;
    T lx = 1.220;
    T deltaPhi = 1.0;
    T Ttop = -0.5;
    T TBottom = 0.5;
    T MEthd = 60;
    T CEthd = 10;
    T ScEthd = 10000;
    bool ifNse;

    pcout << "Ny     :" << NyIn << endl;
    pcout << "T      :" << TEthd << endl;
    pcout << "Ra      :" << RaEthd << endl;
    pcout << "Pr      :" << PrEthd << endl;
    pcout << "Ma     :" << Ma << endl;
    pcout << "tauPhi :" << tauPhi << endl;
    pcout << "lx :" << lx << endl;
    pcout << "temTop :" << Ttop << endl;
    pcout << "temBottom :" << TBottom << endl;
    pcout << "ifNse  :" << ifNse << endl;

    try
    {
        global::argv(1).read(NyIn);
        global::argv(2).read(TEthd);
        global::argv(3).read(ScEthd);
        global::argv(4).read(CEthd);
        global::argv(5).read(MEthd);
        global::argv(6).read(RaEthd);
        global::argv(7).read(PrEthd);
        global::argv(8).read(Ma);
        global::argv(9).read(ifNse);
        global::argv(10).read(lx);
        global::argv(11).read(continueSimulation);
    }
    catch (PlbIOException &exception)
    {
        pcout << exception.what() << endl;
        pcout << "The structure of the input parameters should be : "
              << (string)global::argv(0) << " NyIn, T, continueSimulation" << endl;
        // Exit the program, because wrong input data is a fatal error.
        exit(1);
    }

    plint H = NyIn;
    T ly = 1.0;

    // -------------------设置参数

    MyETHDParameters2D_Nv_Sc<T, NSDESCRIPTOR, ADESCRIPTOR, ADESCRIPTOR>
        parameters(cc, TEthd, MEthd, CEthd, ScEthd, Ma,
                   RaEthd, PrEthd, deltaPhi, TBottom - Ttop,
                   lx, ly, H, 1.0);

    writeLogFile(parameters, "palabosEHD.log");

    plint nx = parameters.getNx();
    plint ny = parameters.getNy();

    T nsOmega = parameters.getOmegaNSE();
    T adElecOmega = parameters.getOmegaElecADE();
    T adTemOmega = parameters.getOmegaTemADE();

    // 电荷弛豫时间

    T poiOmega = 1.0 / tauPhi;

    pcout << "fluid  nsOmega / tau =  " << nsOmega << "/" << 1.0 / nsOmega << endl;
    pcout << "adElecOmega /tau =  " << adElecOmega << " / " << 1.0 / adElecOmega << endl;
    pcout << "adTemOmega /tau =  " << adTemOmega << " / " << 1.0 / adTemOmega << endl;
    pcout << "Ma               = " << parameters.getMa() << endl;
    pcout << "nu/ny           =  " << parameters.getLatticeNu() << " " << ny << endl;
    // ------------------------------------------------------------------
    // 定义NS方程格子

    MultiBlockLattice2D<T, NSDESCRIPTOR> nsLattice(
        nx, ny, new NSDYNAMICS<T, NSDESCRIPTOR>(nsOmega));
    // Use periodic boundary conditions.
    nsLattice.periodicity().toggleAll(true);

    // ------------------------------------------------------------------
    // 定义电场格子

    MultiBlockLattice2D<T, ADESCRIPTOR> adLattice(
        nx, ny, new ADYNAMICS<T, ADESCRIPTOR>(adElecOmega));
    // Use periodic boundary conditions.
    adLattice.periodicity().toggleAll(true);

    // ------------------------------------------------------------------
    // 定义温度方程格子

    MultiBlockLattice2D<T, ADESCRIPTOR> adTemLattice(
        // nx, ny, new ADYNAMICS<T, ADESCRIPTOR>(adTemOmega, parameters.getQ0() / 2, 0.18));
        nx, ny, new ADYNAMICS<T, ADESCRIPTOR>(adTemOmega));
    adTemLattice.periodicity().toggleAll(true);

    // ------------------------------------------------------------------
    // 定义Poisson 格子

    MultiBlockLattice2D<T, POISCRIPTOR> poiLattice(
        nx, ny, new POIDYNAMICS<T, POISCRIPTOR>(poiOmega));
    // Use periodic boundary conditions.
    poiLattice.periodicity().toggleAll(true);

    // ------------------------------------------------------------------
    // 定义NS边界条件，实际上都用了反弹

    OnLatticeBoundaryCondition2D<T, NSDESCRIPTOR>
        *nsBoundaryCondition =
            createLocalBoundaryCondition2D<T, NSDESCRIPTOR>();
    // 电场边界条件

    OnLatticeAdvectionDiffusionBoundaryCondition2D<T, ADESCRIPTOR>
        *adBoundaryCondition =
            createLocalAdvectionDiffusionBoundaryCondition2D<T, ADESCRIPTOR>();
    // 温度场边界条件

    OnLatticeAdvectionDiffusionBoundaryCondition2D<T, ADESCRIPTOR>
        *adTemBoundaryCondition =
            createLocalAdvectionDiffusionBoundaryCondition2D<T, ADESCRIPTOR>();

    // Poisson边界条件

    OnLatticeAdvectionDiffusionBoundaryCondition2D<T, POISCRIPTOR>
        *poiBoundaryCondition =
            createLocalPoissonPoissonLocalBoundaryCondition2D<T, POISCRIPTOR>();

    // ------------------------------------------------------------------
    // 续点参数，只保存电场，温度和流动

    std::vector<MultiBlock2D *> checkpointBlocks;
    checkpointBlocks.push_back(&adLattice);
    checkpointBlocks.push_back(&adTemLattice);
    checkpointBlocks.push_back(&nsLattice);

    const string continueFileName = "continue.xml";
    plint iter = 0;

    // 初始化

    poissonSetup(poiLattice, *poiBoundaryCondition, nx, ny, parameters.getDeltaPhi());
    adeChargeSetup(adLattice, *adBoundaryCondition, nx, ny, parameters.getQ0());
    adeTemperatureSetup(adTemLattice, *adTemBoundaryCondition, nx, ny, TBottom, Ttop);
    nseSetup(nsLattice, *nsBoundaryCondition, nx, ny);

    // 判断是否从上一次的断点 Restart 计算

    if (continueSimulation)
    {
        // pcout << "Reading state of the simulation from file: " << endl;
        loadBinaryBlock(adLattice, "checkpoint_AD.dat");
        loadBinaryBlock(adTemLattice, "checkpoint_TemAD.dat");
        loadBinaryBlock(nsLattice, "checkpoint_NS.dat");
        pcout << "Reading state of the simulation from file: "
              << " " << std::endl;
        loadState(checkpointBlocks, iter, true, continueFileName);
        nsLattice.resetTime(iter);
        adTemLattice.resetTime(iter);
        adLattice.resetTime(iter);
        pcout << "End reading state of the simulation from file: " << std::endl;
        poissonSetup(poiLattice, *poiBoundaryCondition, nx, ny, parameters.getDeltaPhi(), 0);
    }
    else
    {
        pcout << "Start new simulation" << endl;
        poissonSetup(poiLattice, *poiBoundaryCondition, nx, ny, parameters.getDeltaPhi());
        adeChargeSetup(adLattice, *adBoundaryCondition, nx, ny, parameters.getQ0());
        nseSetup(nsLattice, *nsBoundaryCondition, nx, ny);
        adeTemperatureSetup(adTemLattice, *adTemBoundaryCondition, nx, ny, TBottom, Ttop);
    }

    plb_ofstream energyFile;
    std::stringstream prefix;

    prefix << "T_" << std::setfill('0') << std::setw(6) << parameters.getTETHD()
           << "_M_" << std::setw(3) << parameters.getMETHD()
           << "_C_" << std::setw(3) << parameters.getMETHD()
           << "_Sc_" << std::setw(4) << parameters.getScETHD()
           << "_Pr_" << std::setw(3) << parameters.getPrETHD()
           << "_Ra_" << std::setw(8) << parameters.getRaETHD()
           << "_ny_" << std::setw(3) << ny << "_nx_" << std::setw(3) << nx;

    // 根据格子声速的不同来命名输出量的名称
    std::ostringstream name;
    name << prefix.str() << "_average_" << cc << ".dat";

    if (continueSimulation)
    {
        energyFile.open((name).str().c_str(), std::ostream::app);
    }
    else
    {
        energyFile.open(name.str().c_str());
    }

    plint processorLevel = 1;

    MultiScalarField2D<T> chargeOld(nx, ny, (T)0.01);
    MultiScalarField2D<T> temOld(nx, ny, (T)0.01);
    MultiScalarField2D<T> phiOld(nx, ny);
    MultiScalarField2D<T> phiSub(nx, ny);
    MultiTensorField2D<T, 2> electricField(nx, ny);

    MultiTensorField2D<T, 2> temperatureForceField(nx, ny);

    MultiTensorField2D<T, 2> adeConvectionTerm(nx, ny);
    MultiTensorField2D<T, 2> nseForce(nx, ny);
    MultiTensorField2D<T, 2> nseVelocity(nx, ny);
    plint poissonInt = 0;

    chargeOld = (*extractSubDomain(*computeDensity(adLattice), chargeOld.getBoundingBox()));

    Array<T, NSDESCRIPTOR<T>::d> forceTemOrientation((T)0.0, (T)1.0);
    T averageTemIni = (Ttop + TBottom) / 2.0;
    T kiniteEnergyAveold;

    const T prTur = 0.8;
    const T ScTur = 0.4;
    ///****************** get gradient of charge
    MultiScalarField2D<T> gradChargeX(nx, ny);
    MultiScalarField2D<T> gradChargeY(nx, ny);
    MultiScalarField2D<T> gradXPhi(nx, ny);
    MultiScalarField2D<T> gradYPhi(nx, ny);
    MultiScalarField2D<T> gradForceLes(nx, ny);
    MultiScalarField2D<T> phiTemp{*computeDensity(poiLattice)};

    for (plint i = iter; i < maxIter; ++i)
    {

        // ----------------iteration of Poisson LBE------------
        T phiMax = poissonIter(poiLattice, phiOld, phiSub,
                               chargeOld, 1.0 / parameters.getLatticeEpsilon());
        // ---------compute electric field ---------------
        electricField = (*extractSubDomain(*computeVelocity(poiLattice), electricField.getBoundingBox()));
        // ------- compute velocity of NS Equ
        nseVelocity = (*extractSubDomain(*computeVelocity(nsLattice), nseVelocity.getBoundingBox()));
        // the advection terms of NP equation
        adeConvectionTerm = (*add(*multiply(parameters.getLatticeK(), electricField), nseVelocity));
        // --------add the advection term into the LBE of NP equ--------
        applyProcessingFunctional(
            new getAdvectionTerm<T, ADESCRIPTOR, T, 2>(1.0),
            adLattice.getBoundingBox(),
            adLattice, adeConvectionTerm);

#ifdef lesModel
        // Eddy viscosity coefficient of temperature equ
        applyProcessingFunctional(
            new TurbulencePrandalProcessor2D<T, NSDESCRIPTOR, ADESCRIPTOR>(prTur, ccDefine),
            nsLattice.getBoundingBox(),
            nsLattice, adTemLattice);

        // Eddy viscosity coefficient of np equ
        applyProcessingFunctional(
            new TurbulencePrandalProcessor2D<T, NSDESCRIPTOR, ADESCRIPTOR>(ScTur, ccDefine),
            nsLattice.getBoundingBox(),
            nsLattice, adLattice);
#endif

        // collision of advection diffusion equation
        adLattice.collideAndStream();

        // pcout << "End poi Iter, phiMax ... " << phiMax << endl;
        // get the charge
        chargeOld = (*extractSubDomain(*computeDensity(adLattice), chargeOld.getBoundingBox()));
        // electric field force, Fe = qE
        applyProcessingFunctional(
            new getElectricForce<T, T, 2>(1.0),
            chargeOld.getBoundingBox(),
            chargeOld, electricField);

        //  ------------------------------
        // add the advection term into the temperature equation
        applyProcessingFunctional(
            new getAdvectionTerm<T, ADESCRIPTOR, T, 2>(1.0),
            adTemLattice.getBoundingBox(),
            adTemLattice, nseVelocity);

        // collision of temperature LBE
        adTemLattice.collideAndStream();

        // pcout << "Iteration= " << i << " " <<averageTemIni<< " " <<parameters.getDeltaT() << endl;
        //  BoussinesqThermal hypothesis
        applyProcessingFunctional(
            new MyBoussinesqThermalProcessor2D<T, ADESCRIPTOR, T, 2> //
            (parameters.getGBeta(), averageTemIni, parameters.getDeltaT(), forceTemOrientation),
            temperatureForceField.getBoundingBox(), adTemLattice, temperatureForceField);

        // add electric field force and bouyance
        nseForce = *add(electricField, temperatureForceField);
        // add external force into NS equ
        applyProcessingFunctional(
            new getForceTerm<T, NSDESCRIPTOR, T, 2>(1.0),
            nsLattice.getBoundingBox(),
            nsLattice, nseForce);
        // get the electric potential
        phiTemp = *computeDensity(poiLattice);

        //*************************** charge gradient******************************************************
        applyProcessingFunctional(
            new BoxXderivativeFunctional2D<T>(),
            gradChargeX.getBoundingBox(), chargeOld, gradChargeX, 1);
        applyProcessingFunctional(
            new BoxYderivativeFunctional2D<T>(),
            gradChargeY.getBoundingBox(), chargeOld, gradChargeY, 1);
        applyProcessingFunctional(
            new BoxXderivativeFunctional2D<T>(),
            gradXPhi.getBoundingBox(), phiTemp, gradXPhi, 1);
        applyProcessingFunctional(
            new BoxYderivativeFunctional2D<T>(),
            gradYPhi.getBoundingBox(), phiTemp, gradYPhi, 1);
        gradChargeX = *multiply(gradChargeX, gradXPhi);
        gradChargeY = *multiply(gradChargeY, gradYPhi);
        // ----------------------(nabla q) * (nabla phi)
        gradForceLes = *add(gradChargeX, gradChargeY);

        // source term of NS equation
        applyProcessingFunctional(
            new getNSSuorceGradLes<T, NSDESCRIPTOR, T>(1.0),
            nsLattice.getBoundingBox(),
            nsLattice, gradForceLes);

        // if flow?
        if (ifNse)
        {
            for (int i = 0; i < plint(cc); ++i)
            {
                nsLattice.collideAndStream();
            }
        }

        // some statistics quantity
        if (i % outStatistics == 0)
        {
            ///****************** get negative electric field strength

            applyProcessingFunctional(
                new BoxXderivativeFunctional2D<T>(),
                gradXPhi.getBoundingBox(), phiTemp, gradXPhi, 1);
            applyProcessingFunctional(
                new BoxYderivativeFunctional2D<T>(),
                gradYPhi.getBoundingBox(), phiTemp, gradYPhi, 1);

            ///****************** get gradient of charge
            MultiScalarField2D<T> gradChargeX(nx, ny);
            MultiScalarField2D<T> gradChargeY(nx, ny);
            applyProcessingFunctional(
                new BoxXderivativeFunctional2D<T>(),
                gradChargeX.getBoundingBox(), chargeOld, gradChargeX, 1);
            applyProcessingFunctional(
                new BoxYderivativeFunctional2D<T>(),
                gradChargeY.getBoundingBox(), chargeOld, gradChargeY, 1);

            ///************ get current of every component
            MultiScalarField2D<T> JDiffX(nx, ny);
            MultiScalarField2D<T> JDiffY(nx, ny);
            MultiScalarField2D<T> JVeloX(nx, ny);
            MultiScalarField2D<T> JVeloY(nx, ny);
            MultiScalarField2D<T> JVeloX0(nx, ny);
            MultiScalarField2D<T> JVeloY0(nx, ny);
            MultiScalarField2D<T> JElecX(nx, ny);
            MultiScalarField2D<T> JElecY(nx, ny);
            MultiScalarField2D<T> TVeloY(nx, ny);
            MultiScalarField2D<T> TDiffY(nx, ny);

            ///************ get current of electric --- Electromigration
            T kElec = parameters.getLatticeK();
            JElecX = (*multiply(*multiply(-kElec, gradXPhi), chargeOld));
            JElecY = (*multiply(*multiply(-kElec, gradYPhi), chargeOld));

            ///************ get current of diffusivity
            JDiffX = (*multiply(-parameters.getLatticeD(), gradChargeX));
            JDiffY = (*multiply(-parameters.getLatticeD(), gradChargeY));

            ///************ get current of fluid motion
            applyProcessingFunctional(
                new getComponentVector<T, T, 2>(0),
                JVeloX.getBoundingBox(), JVeloX0, *computeVelocity(nsLattice));
            applyProcessingFunctional(
                new getComponentVector<T, T, 2>(1),
                JVeloY.getBoundingBox(), JVeloY0, *computeVelocity(nsLattice));

            JVeloX = (*multiply(JVeloX0, chargeOld));
            JVeloY = (*multiply(JVeloY0, chargeOld));

            Box2D hisCenter(nx / 4, nx / 4, ny / 4, ny / 4);

            // average electric current
            T JDiffYAve = computeAverage(JDiffY) * parameters.getDimlessCurrent();
            T JVeloYAve = computeAverage(JVeloY) * parameters.getDimlessCurrent();
            T JElecYAve = computeAverage(JElecY) * parameters.getDimlessCurrent();
            // 计算总和电流

            // T JDiffYAve = computeSum(JDiffY) * parameters.getDimlessCurrent();
            // T JVeloYAve = computeSum(JVeloY) * parameters.getDimlessCurrent();
            // T JElecYAve = computeSum(JElecY) * parameters.getDimlessCurrent();

            T J = JDiffYAve + JVeloYAve + JElecYAve;

            // a reference electric Nusselt number
            T NuEle = J / (JDiffYAve + JElecYAve);

            MultiTensorField2D<T, 2> velCenter = (*extractSubDomain(*computeVelocity(nsLattice), hisCenter));

            // variable of one point
            T uCenter = computeAverage(*computeVelocityNorm(nsLattice, hisCenter));
            T uxCenter = computeAverage(*computeVelocityComponent(nsLattice, hisCenter, 0));
            T uyCenter = computeAverage(*computeVelocityComponent(nsLattice, hisCenter, 1));
            T veloMax = computeMax(*computeVelocityNorm(nsLattice));

            plint yDirection = 1;

            // thermal Nusselt number
            T nusselt = computeNusseltNumber(
                nsLattice, adTemLattice, nsLattice.getBoundingBox(),
                yDirection, parameters.getDeltaX(),
                parameters.getLatticeAlpha(), parameters.getDeltaT());

            MultiScalarField2D<T> gradTemY(nx, ny);
            MultiScalarField2D<T> tem{*computeDensity(adTemLattice)};

            // bottom and top side
            Box2D sliceBottom(0, nx - 1, 0, 0);
            Box2D sliceTop(0, nx - 1, ny - 1, ny - 1);

            applyProcessingFunctional(
                new BoxYderivativeFunctional2D<T>(),
                gradTemY.getBoundingBox(), tem, gradTemY, 1);
            // gradient of temperature
            MultiScalarField2D<T> temGradBottom(*extractSubDomain(gradTemY, sliceBottom));
            MultiScalarField2D<T> temGradTop(*extractSubDomain(gradTemY, sliceTop));
            // compure thermal nusselt number
            T nuAve = -computeAverage(gradTemY) * parameters.getResolution();
            T nuBottom = -computeAverage(temGradBottom) * parameters.getResolution();
            T nuTop = -computeAverage(temGradTop) * parameters.getResolution();

            TDiffY = (*multiply(-parameters.getLatticeAlpha(), gradTemY));
            TVeloY = (*multiply(JVeloY0, tem));
            T TheatFluxV = computeAverage(TVeloY) * parameters.getDimlessHeatFlux();
            T TheatFluxD = computeAverage(TDiffY) * parameters.getDimlessHeatFlux();

            T kiniteEnergyAve = computeAverageEnergy(nsLattice) * parameters.getDimlessVelocity() * parameters.getDimlessVelocity();

            T rhoAll = computeSum(*computeDensity(nsLattice));

            pcout << "Iteration= " << i << " " << i * parameters.getDimlessTime()
                  << "; maxV = " << veloMax
                  << " " << veloMax * parameters.getDimlessVelocity()
                  << " " << uCenter * parameters.getDimlessVelocity()
                  << " " << uxCenter * parameters.getDimlessVelocity()
                  << " " << uyCenter * parameters.getDimlessVelocity()
                  //   << " " << vCenter * parameters.getDimlessVelocity()
                  << " " << JDiffYAve
                  << " " << JElecYAve
                  << " " << JVeloYAve
                  << " " << J
                  << " " << NuEle
                  << " " << nusselt
                  << " " << nuAve
                  << " " << nuBottom
                  << " " << nuTop
                  << " " << kiniteEnergyAve
                  << " " << rhoAll
                  << std::endl;

            energyFile << i << " " << i * parameters.getDimlessTime()
                       << std::setprecision(16)
                       << " " << veloMax //<< " " << veloMax * parameters.getDimlessVelocity()
                       << " " << veloMax * parameters.getDimlessVelocity()
                       << " " << uCenter * parameters.getDimlessVelocity()
                       << " " << uxCenter * parameters.getDimlessVelocity()
                       << " " << uyCenter * parameters.getDimlessVelocity()
                       //    << " " << vCenter * parameters.getDimlessVelocity()
                       << " " << JDiffYAve
                       << " " << JElecYAve
                       << " " << JVeloYAve
                       << " " << J
                       << " " << NuEle
                       << " " << nusselt
                       << " " << nuAve
                       << " " << nuBottom
                       << " " << nuTop
                       << " " << kiniteEnergyAve
                       << " " << TheatFluxV
                       << " " << TheatFluxD
                       << " " << TheatFluxD + TheatFluxV
                       << std::endl;
            computeDissipation(i, nsLattice, adLattice, adTemLattice, poiLattice, nx, ny, parameters, prefix.str());

            if (i % vtkOutNum == 0)
            {
#ifdef IF_TEST
                MultiScalarField2D<T> ux = *computeVelocityComponent(nsLattice, 0);
                MultiScalarField2D<T> uy = *computeVelocityComponent(nsLattice, 1);

                MultiScalarField2D<T> dUxDx(nx, ny);
                MultiScalarField2D<T> dUxDy(nx, ny);
                MultiScalarField2D<T> dUyDx(nx, ny);
                MultiScalarField2D<T> dUyDy(nx, ny);

                applyProcessingFunctional(
                    new BoxXderivativeFunctional2D<T>(),
                    ux.getBoundingBox(), ux, dUxDx, 1);
                applyProcessingFunctional(
                    new BoxXderivativeFunctional2D<T>(),
                    uy.getBoundingBox(), uy, dUyDx, 1);
                applyProcessingFunctional(
                    new BoxYderivativeFunctional2D<T>(),
                    ux.getBoundingBox(), ux, dUxDy, 1);
                applyProcessingFunctional(
                    new BoxYderivativeFunctional2D<T>(),
                    uy.getBoundingBox(), uy, dUyDy, 1);

                // ux = multiply
                ux = *add((*multiply(dUxDx, dUxDx)), (*multiply(dUyDy, dUyDy)));
                uy = *add(dUxDy, dUyDx);
                uy = *multiply(uy, uy);
                uy = *multiply(0.5, uy);

                ux = *add(ux, uy);
                ux = *multiply(2., ux);
                ux = *computeSqrt(ux);
                ux = *multiply(pow((0.18 * 1.414), 2), ux);
#endif
                pcout << "i = " << i << ", Saving Gif ..." << endl;
                plint iT = i;
                VtkImageOutput2D<T> vtkOut(createFileName(prefix.str() + "vtk", iT, 8), 1);
                vtkOut.writeData<float>(*computeDensity(poiLattice), "phi", (T)1);
                vtkOut.writeData<float>(*computeDensity(nsLattice), "rho", 1);
                vtkOut.writeData<float>(*computeTemperature(nsLattice), "nuLestem", 1);
                vtkOut.writeData<float>(*computeTemperature(adLattice), "DElecLestem", 1);
                vtkOut.writeData<float>(*computeTemperature(adTemLattice), "TemlecLestem", 1);
                // vtkOut.writeData<float>(ux, "stressSquare", 1);

                // vtkOut.writeData<float>(*computeTauLes(nsLattice), "nuLes", 1);
                vtkOut.writeData<float>(*computeDensity(adLattice), "charge1", 1);
                vtkOut.writeData<float>(*computeDensity(adTemLattice), "temperature", 1);
                vtkOut.writeData<2, float>(*computeVelocity(nsLattice), "velocity1", 1);
                vtkOut.writeData<2, float>(*computeVelocity(poiLattice), "EIndependent", 1);
#ifdef IF_TEST
                // vtkOut.writeData<2, float>(temperatureForceField, "temperatureForce", 1);
                // vtkOut.writeData<float>(JElecX, "JElecX", 1);
                // vtkOut.writeData<float>(JElecY, "JElecY", 1);
                // vtkOut.writeData<float>(JDiffX, "JDiffX", 1);
                // vtkOut.writeData<float>(JDiffY, "JDiffY", 1);
                // vtkOut.writeData<float>(JVeloX, "JVeloX", 1);
                // vtkOut.writeData<float>(JVeloY, "JVeloY", 1);
#endif
            }
        }

        // whether save the restart computation
        if (i % continueIter == 0)
        {
            const int PADD = 8;
            const string continueFile = "continue.xml";
            const string checkpointFile = "checkpoint_";
            pcout << "Saving the state of the simulation at iteration: " << i << endl;
            saveState(checkpointBlocks, i, true, continueFile, checkpointFile, PADD);
            pcout << std::endl;

            pcout << "Saving the state of the simulation at iteration: " << i << std::endl;
            saveBinaryBlock(adLattice, "checkpoint_AD.dat");
            saveBinaryBlock(poiLattice, "checkpoint_TemAD.dat");
            saveBinaryBlock(nsLattice, "checkpoint_NS.dat");

            pcout << "i = " << i << "Saving continue Simulation ..." << endl;
        }
    }

    delete poiBoundaryCondition;
    delete nsBoundaryCondition;
    delete adBoundaryCondition;
    delete adTemBoundaryCondition;
}

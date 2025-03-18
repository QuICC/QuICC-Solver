/**
 * @file TestUtils.cpp
 * @brief Source of the Test utils for IO validation
 */

// System includes
//
#include <catch2/catch.hpp>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/TestSuite/Framework/Io/TestUtils.hpp"
#include "QuICC/TestSuite/Framework/Io/TestBackend.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Generator/States/SphereExactScalarState.hpp"
#include "QuICC/Generator/States/SphereExactVectorState.hpp"
#include "QuICC/Io/Variable/SphereDipolarityWriter.hpp"
#include "QuICC/SpectralKernels/Typedefs.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/PhysicalNames/Magnetic.hpp"
#include "QuICC/Model/IModelBackend.hpp"
#include "QuICC/Generator/StateGenerator.hpp"
#include "QuICC/Io/Variable/SphereAngularMomentumWriter.hpp"
#include "QuICC/Io/Variable/SphereDipolarityWriter.hpp"
#include "QuICC/Io/Variable/SphereMaxAbsoluteFieldValueWriter.hpp"
#include "QuICC/Io/Variable/SphereNusseltWriter.hpp"
#include "QuICC/Io/Variable/SphereScalarEnergyWriter.hpp"
#include "QuICC/Io/Variable/SphereScalarLSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereScalarMSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereScalarMeanWriter.hpp"
#include "QuICC/Io/Variable/SphereScalarModeSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereScalarNSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereScalarRSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolEnergyWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolLSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolMSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolNSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolRSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolModeSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolEnstrophyWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolEnstrophyLSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolEnstrophyMSpectrumWriter.hpp"

#include "QuICC/Io/Variable/ShellNusseltWriter.hpp"
#include "QuICC/Io/Variable/ShellScalarEnergyWriter.hpp"
#include "QuICC/Io/Variable/ShellScalarLSpectrumWriter.hpp"
#include "QuICC/Io/Variable/ShellScalarMSpectrumWriter.hpp"
#include "QuICC/Io/Variable/ShellTorPolEnergyWriter.hpp"
#include "QuICC/Io/Variable/ShellTorPolLSpectrumWriter.hpp"
#include "QuICC/Io/Variable/ShellTorPolMSpectrumWriter.hpp"
#include "QuICC/Io/Variable/ShellTorPolEnstrophyWriter.hpp"
#include "QuICC/NonDimensional/RRatio.hpp"
#include "QuICC/NonDimensional/Lower1d.hpp"
#include "QuICC/NonDimensional/Upper1d.hpp"
#include "TestSuite/Io.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace Io {

std::shared_ptr<StateGenerator> createRunner(std::shared_ptr<SpatialScheme::ISpatialScheme> spScheme, TestParameters& test)
{
   // Create simulation
   auto spRunner = std::make_shared<StateGenerator>();

   // Create list of field ID strings for boundary conditions
   std::vector<std::string> bcNames = {};

   // Create list of nondimensional ID strings for physical parameters
   std::vector<std::string> ndNames = {};
   if(spScheme->tag() == "SLFl")
   {
      ndNames.push_back(NonDimensional::RRatio().tag());
   }

   // Get model configuration tags
   std::map<std::string, std::map<std::string, int>> modelCfg = {};

   // Geometry setup
   std::vector<bool> isPeriodic = {false, false, false};

   // Add configuration file and parameters
   spRunner->setConfiguration(spScheme->dimension(), spScheme->tag(), isPeriodic, bcNames, ndNames, modelCfg);

   // Initialise simulation
   std::map<std::string,MHDFloat> cfg;
   std::set<SpatialScheme::Feature> features;
   spRunner->getConfig(cfg, features, "IO test");
   spScheme->enable(features);

   std::map<std::string, MHDFloat> extra;
   if(cfg.count(NonDimensional::RRatio().tag()) > 0)
   {
      auto rratio = cfg.at(NonDimensional::RRatio().tag());
      extra.emplace(NonDimensional::Lower1d().tag(), rratio / (1.0 - rratio));
      extra.emplace(NonDimensional::Upper1d().tag(), 1.0 / (1.0 - rratio));
   }
   spRunner->updateConfig(extra);

   spRunner->initBase();

   // Initialise resolution
   spRunner->initResolution(spScheme);

   // Create the reference states
   auto spEq = createStates(spRunner);

   // Add sphere ascii files
   if(spScheme->tag() == "WLFl" || spScheme->tag() == "WLFm")
   {
      addSphereFiles(spRunner);
   }
   // Add shell ascii files
   else if(spScheme->tag() == "SLFl" || spScheme->tag() == "SLFm")
   {
      addShellFiles(spRunner);
   }
   else
   {
      throw std::logic_error("No ASCII files have been setup for " + spScheme->tag() + " scheme");
   }

   // Set the boundary conditions
   SharedSimulationBoundary spBcs = spRunner->createBoundary();

   // Initialise the simulation
   spRunner->init(spBcs);

   test.spRes = spEq->spRes();

   return spRunner;
}

Equations::SharedIEquation createStates(std::shared_ptr<StateGenerator> spRunner)
{
   std::shared_ptr<Model::IModelBackend> spBackend = std::make_shared<TestBackend>(spRunner->ss().tag());

   // Add state generation equations
   Equations::SharedSphereExactScalarState spScalar;
   Equations::SharedSphereExactVectorState spVector;

   Spectral::Kernel::Complex3DMapType tSH;
   std::pair<Spectral::Kernel::Complex3DMapType::iterator, bool> ptSH;

   spScalar =
      spRunner->addEquation<Equations::SphereExactScalarState>(spBackend);
   spScalar->setIdentity(PhysicalNames::Temperature::id());
   tSH.clear();
   ptSH = tSH.insert(
      std::make_pair(std::make_pair(3, 3), std::map<int, MHDComplex>()));
   ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0, 2.0)));
   spScalar->setSpectralModes(tSH);

   spVector =
      spRunner->addEquation<Equations::SphereExactVectorState>(spBackend);
   spVector->setIdentity(PhysicalNames::Velocity::id());
   // Toroidal
   tSH.clear();
   ptSH = tSH.insert(
         std::make_pair(std::make_pair(1, 1), std::map<int, MHDComplex>()));
   ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
   spVector->setSpectralModes(FieldComponents::Spectral::TOR, tSH);
   // Poloidal
   tSH.clear();
   ptSH = tSH.insert(
         std::make_pair(std::make_pair(2, 0), std::map<int, MHDComplex>()));
   ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
   spVector->setSpectralModes(FieldComponents::Spectral::POL, tSH);

   spVector =
      spRunner->addEquation<Equations::SphereExactVectorState>(spBackend);
   spVector->setIdentity(PhysicalNames::Magnetic::id());
   // Toroidal
   tSH.clear();
   ptSH = tSH.insert(
         std::make_pair(std::make_pair(1, 1), std::map<int, MHDComplex>()));
   ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
   spVector->setSpectralModes(FieldComponents::Spectral::TOR, tSH);
   // Poloidal
   tSH.clear();
   ptSH = tSH.insert(
         std::make_pair(std::make_pair(2, 0), std::map<int, MHDComplex>()));
   ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
   spVector->setSpectralModes(FieldComponents::Spectral::POL, tSH);

   return spVector;
}

void addSphereFiles(std::shared_ptr<StateGenerator> spRunner)
{
   // Add ASCII output files
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::SphereAngularMomentumWriter>("", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Velocity::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::SphereDipolarityWriter>("magnetic", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Magnetic::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::SphereMaxAbsoluteFieldValueWriter>("velocity", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Velocity::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::SphereNusseltWriter>("", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Temperature::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::SphereScalarEnergyWriter>("temperature", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Temperature::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::SphereScalarLSpectrumWriter>("temperature", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Temperature::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::SphereScalarMSpectrumWriter>("temperature", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Temperature::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::SphereScalarMeanWriter>("temperature", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Temperature::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::SphereScalarModeSpectrumWriter>("temperature", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Temperature::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::SphereScalarNSpectrumWriter>("temperature", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Temperature::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
   auto spFile = std::make_shared<QuICC::Io::Variable::SphereScalarRSpectrumWriter>("temperature", spRunner->ss().tag());
   spFile->expect(PhysicalNames::Temperature::id());
   spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::SphereTorPolEnergyWriter>("velocity", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Velocity::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::SphereTorPolLSpectrumWriter>("velocity", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Velocity::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::SphereTorPolMSpectrumWriter>("velocity", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Velocity::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::SphereTorPolNSpectrumWriter>("velocity", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Velocity::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::SphereTorPolRSpectrumWriter>("velocity", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Velocity::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::SphereTorPolModeSpectrumWriter>("velocity", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Velocity::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::SphereTorPolEnstrophyWriter>("velocity", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Velocity::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::SphereTorPolEnstrophyLSpectrumWriter>("velocity", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Velocity::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::SphereTorPolEnstrophyMSpectrumWriter>("velocity", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Velocity::id());
      spRunner->addAsciiOutputFile(spFile);
   }
}

void addShellFiles(std::shared_ptr<StateGenerator> spRunner)
{
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::ShellNusseltWriter>("", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Temperature::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::ShellScalarEnergyWriter>("temperature", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Temperature::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::ShellScalarLSpectrumWriter>("temperature", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Temperature::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::ShellScalarMSpectrumWriter>("temperature", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Temperature::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::ShellTorPolEnergyWriter>("velocity", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Velocity::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::ShellTorPolLSpectrumWriter>("velocity", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Velocity::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::ShellTorPolMSpectrumWriter>("velocity", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Velocity::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::ShellTorPolEnstrophyWriter>("velocity", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Velocity::id());
      spRunner->addAsciiOutputFile(spFile);
   }
#if 0
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::ShellTorPolEnstrophyLSpectrumWriter>("velocity", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Velocity::id());
      spRunner->addAsciiOutputFile(spFile);
   }
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::ShellTorPolEnstrophyMSpectrumWriter>("velocity", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Velocity::id());
      spRunner->addAsciiOutputFile(spFile);
   }
#endif
}

void checkFiles(std::shared_ptr<SpatialScheme::ISpatialScheme> spScheme, const TestParameters& test)
{
   // Add sphere ascii files
   if(spScheme->tag() == "WLFl" || spScheme->tag() == "WLFm")
   {
      checkSphereFiles(test);
   }
   // Add shell ascii files
   else if(spScheme->tag() == "SLFl" || spScheme->tag() == "SLFm")
   {
      checkShellFiles(test);
   }
   else
   {
      throw std::logic_error("No ASCII files have been setup for " + spScheme->tag() + " scheme");
   }
}

void checkSphereFiles(const TestParameters& test)
{
   const std::string& datadir = test.datadir;
   const std::string& refdir = test.refdir;

   const auto& nN = test.spRes->sim().dim(QuICC::Dimensions::Simulation::SIM1D, QuICC::Dimensions::Space::SPECTRAL);
   const auto& nL = test.spRes->sim().dim(QuICC::Dimensions::Simulation::SIM2D, QuICC::Dimensions::Space::SPECTRAL);
   const auto& nM = test.spRes->sim().dim(QuICC::Dimensions::Simulation::SIM3D, QuICC::Dimensions::Space::SPECTRAL);
   const auto& nR = test.spRes->sim().dim(QuICC::Dimensions::Simulation::SIM1D, QuICC::Dimensions::Space::PHYSICAL);
   int nH = nL*(nL+1)/2;
   const auto& maxUlp = test.maxUlp;

   MHDFloat eps = std::numeric_limits<MHDFloat>::epsilon();
   MHDFloat tol = maxUlp*eps;

   auto compareRef = [=](std::string fname, int rows, int cols, int blocks = 1)
   {
      std::vector<Matrix> data;
      std::vector<Matrix> ref;
      for(int i = 0; i < blocks; i++)
      {
         data.emplace(data.end(), rows, cols);
         ref.emplace(ref.end(), rows, cols);
      }
      QuICC::TestSuite::readBlockData(data, datadir + fname);
      QuICC::TestSuite::readBlockData(ref, refdir + fname);

      INFO( "Checking " + fname );
      CHECK( data.size() == ref.size() );
      for(int k = 0; k < blocks; k++)
      {
         const auto& d = data.at(k);
         const auto& r = data.at(k);
         CHECK( d.rows() == r.rows() );
         CHECK( d.cols() == r.cols() );

         for(int i = 0; i < d.rows(); i++)
         {
            for(int j = 0; j < d.cols(); j++)
            {
               auto err = computeUlp(d(i,j), r(i,j), std::abs(r(i,j)), tol, eps);
               INFO( "i,j,k: " << i << "," << j << "," << k );
               INFO( "data: " << std::scientific << std::setprecision(16) << d(i,j) );
               INFO( "ref: " << std::scientific << std::setprecision(16) << r(i,j) );
               INFO( "max ulp: " << maxUlp);
               INFO( "measured ulp: " << std::get<1>(err) );
               CHECK( std::get<0>(err) );
            }
         }
      }
   };

   //
   {
      std::string fname = "angular_momentum.dat";
      int r = 1;
      int c = 5;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "magnetic_dipolarity.dat";
      int r = 1;
      int c = 5;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "velocity_max.dat";
      int r = 1;
      int c = 2;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "nusselt.dat";
      int r = 1;
      int c = 2;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "temperature_energy.dat";
      int r = 1;
      int c = 2;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "temperature_l_spectrum.dat";
      int r = nL;
      int c = 2;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "temperature_m_spectrum.dat";
      int r = nM;
      int c = 2;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "temperature_mean.dat";
      int r = 1;
      int c = 2;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "temperature_mode_spectrum.dat";
      int r = nH;
      int c = 3;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "temperature_n_spectrum.dat";
      int r = nN;
      int c = nL + 1;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "temperature_r_spectrum.dat";
      int r = nR;
      int c = nL + 1;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "velocity_energy.dat";
      int r = 1;
      int c = 4;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "velocity_l_spectrum.dat";
      int r = nL;
      int c = 4;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "velocity_m_spectrum.dat";
      int r = nM;
      int c = 4;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "velocity_mode_spectrum.dat";
      int r = nH;
      int c = 5;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "velocity_n_spectrum.dat";
      int r = nN;
      int c = nL + 1;
      int b = 3;
      compareRef(fname, r, c, b);
   }

   //
   {
      std::string fname = "velocity_r_spectrum.dat";
      int r = nR;
      int c = nL + 1;
      int b = 3;
      compareRef(fname, r, c, b);
   }

   //
   {
      std::string fname = "velocity_enstrophy.dat";
      int r = 1;
      int c = 4;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "velocity_enstrophy_l_spectrum.dat";
      int r = nL;
      int c = 4;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "velocity_enstrophy_m_spectrum.dat";
      int r = nM;
      int c = 4;
      compareRef(fname, r, c);
   }
}

void checkShellFiles(const TestParameters& test)
{
   std::string datadir = "./";
   std::string refdir = "./";

   const auto& nN = test.spRes->sim().dim(QuICC::Dimensions::Simulation::SIM1D, QuICC::Dimensions::Space::SPECTRAL);
   const auto& nL = test.spRes->sim().dim(QuICC::Dimensions::Simulation::SIM2D, QuICC::Dimensions::Space::SPECTRAL);
   const auto& nM = test.spRes->sim().dim(QuICC::Dimensions::Simulation::SIM3D, QuICC::Dimensions::Space::SPECTRAL);
   const auto& nR = test.spRes->sim().dim(QuICC::Dimensions::Simulation::SIM1D, QuICC::Dimensions::Space::PHYSICAL);
   int nH = nL*(nL+1)/2;
   const auto& maxUlp = test.maxUlp;

   std::cerr << "PARAMS: " << std::endl;
   std::cerr << " nN = " << nN << std::endl;
   std::cerr << " nL = " << nL << std::endl;
   std::cerr << " nM = " << nM << std::endl;
   std::cerr << " nR = " << nR << std::endl;
   std::cerr << " nR = " << maxUlp << std::endl;

   MHDFloat eps = std::numeric_limits<MHDFloat>::epsilon();
   MHDFloat tol = maxUlp*eps;

   auto compareRef = [=](std::string fname, int rows, int cols, int blocks = 1)
   {
      std::vector<Matrix> data;
      std::vector<Matrix> ref;
      for(int i = 0; i < blocks; i++)
      {
         data.emplace(data.end(), rows, cols);
         ref.emplace(ref.end(), rows, cols);
      }
      QuICC::TestSuite::readBlockData(data, datadir + fname);
      QuICC::TestSuite::readBlockData(ref, refdir + fname);

      INFO( "Checking " + fname );
      CHECK( data.size() == ref.size() );
      for(int k = 0; k < blocks; k++)
      {
         const auto& d = data.at(k);
         const auto& r = data.at(k);
         CHECK( d.rows() == r.rows() );
         CHECK( d.cols() == r.cols() );

         for(int i = 0; i < d.rows(); i++)
         {
            for(int j = 0; j < d.cols(); j++)
            {
               auto err = computeUlp(d(i,j), r(i,j), std::abs(r(i,j)), tol, eps);
               INFO( "i,j,k: " << i << "," << j << "," << k );
               INFO( "data: " << std::scientific << std::setprecision(16) << d(i,j) );
               INFO( "ref: " << std::scientific << std::setprecision(16) << r(i,j) );
               INFO( "max ulp: " << maxUlp);
               INFO( "measured ulp: " << std::get<1>(err) );
               CHECK( std::get<0>(err) );
            }
         }
      }
   };

   //
   {
      std::string fname = "nusselt.dat";
      int r = 1;
      int c = 2;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "temperature_energy.dat";
      int r = 1;
      int c = 2;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "temperature_l_spectrum.dat";
      int r = nL;
      int c = 2;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "temperature_m_spectrum.dat";
      int r = nM;
      int c = 2;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "velocity_energy.dat";
      int r = 1;
      int c = 4;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "velocity_l_spectrum.dat";
      int r = nL;
      int c = 4;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "velocity_m_spectrum.dat";
      int r = nM;
      int c = 4;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "velocity_enstrophy.dat";
      int r = 1;
      int c = 4;
      compareRef(fname, r, c);
   }

#if 0
   //
   {
      std::string fname = "velocity_enstrophy_l_spectrum.dat";
      int r = nL;
      int c = 4;
      compareRef(fname, r, c);
   }

   //
   {
      std::string fname = "velocity_enstrophy_m_spectrum.dat";
      int r = nM;
      int c = 4;
      compareRef(fname, r, c);
   }
#endif
}


ErrorType computeUlp(const MHDFloat data, const MHDFloat ref, MHDFloat refMod, const MHDFloat tol, const MHDFloat eps)
{
   bool isEqual = false;
   if(ref == 0.0)
   {
      refMod = 1.0;
   }

   auto diff = std::abs(data-ref);

   if(diff < tol)
   {
      isEqual = diff < (tol * refMod);
   }
   else
   {
      isEqual = (diff / refMod ) < tol;
   }

   auto ulp = diff / (refMod * eps);

   return std::make_tuple(isEqual, ulp, diff);
}

} // namespace Io
} // namespace Framework
} // namespace TestSuite
} // namespace QuICC

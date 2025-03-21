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
#include "QuICC/SpectralKernels/Typedefs.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/PhysicalNames/Magnetic.hpp"
#include "QuICC/Model/IModelBackend.hpp"
#include "QuICC/Generator/StateGenerator.hpp"
#include "QuICC/NonDimensional/RRatio.hpp"
#include "QuICC/NonDimensional/Heating.hpp"
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
   if(spScheme->tag() == "SLFl" || spScheme->tag() == "SLFm")
   {
      MHDFloat rratio = 0.35;
      MHDFloat heating = 0;

      extra.emplace(NonDimensional::RRatio().tag(), rratio);
      extra.emplace(NonDimensional::Heating().tag(), heating);
      extra.emplace(NonDimensional::Lower1d().tag(), rratio / (1.0 - rratio));
      extra.emplace(NonDimensional::Upper1d().tag(), 1.0 / (1.0 - rratio));
   }
   spRunner->updateConfig(extra);

   spRunner->initBase();

   // Initialise resolution
   spRunner->initResolution(spScheme);

   // Create the reference states
   auto spEq = createStates(spRunner);


   // Add ASCII output files
   for(auto&& f: test.files)
   {
      spRunner->addAsciiOutputFile(f);
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
   int maxN = 5;
   int maxL = 5;
   int maxM = 5;

   auto setScalarModes = [](auto& tSH, const int maxN, const int maxL, const int maxM, const int l0)
   {
      std::pair<Spectral::Kernel::Complex3DMapType::iterator, bool> ptSH;
      for(int m = 0; m <= maxM; m++)
      {
         for(int l = std::min(l0,m); l <= maxL; l++)
         {
            ptSH = tSH.insert(
               std::make_pair(std::make_pair(l, m), std::map<int, MHDComplex>()));
            for(int n = 0; n <= maxN; n++)
            {
               MHDComplex val(3.0/static_cast<MHDFloat>(n+1), 2.0/static_cast<MHDFloat>(n+1));
               if(m == 0)
               {
                  val = std::real(val);
               }
               ptSH.first->second.insert(std::make_pair(n, val));
            }
         }
      }
   };

   auto setTorModes = [](auto& tSH, const int maxN, const int maxL, const int maxM, const int l0)
   {
      std::pair<Spectral::Kernel::Complex3DMapType::iterator, bool> ptSH;
      for(int m = 0; m <= maxM; m++)
      {
         for(int l = std::min(l0,m); l <= maxL; l++)
         {
            ptSH = tSH.insert(
               std::make_pair(std::make_pair(l, m), std::map<int, MHDComplex>()));
            for(int n = 0; n <= maxN; n++)
            {
               MHDComplex val(1.0/static_cast<MHDFloat>(n+1), 2.0/static_cast<MHDFloat>(n+1));
               if(m == 0)
               {
                  val = std::real(val);
               }
               ptSH.first->second.insert(std::make_pair(n, val));
            }
         }
      }
   };

   auto setPolModes = [](auto& tSH, const int maxN, const int maxL, const int maxM, const int l0)
   {
      std::pair<Spectral::Kernel::Complex3DMapType::iterator, bool> ptSH;
      for(int m = 0; m <= maxM; m++)
      {
         for(int l = std::min(l0,m); l <= maxL; l++)
         {
            ptSH = tSH.insert(
               std::make_pair(std::make_pair(l, m), std::map<int, MHDComplex>()));
            for(int n = 0; n <= maxN; n++)
            {
               MHDComplex val(2.0/static_cast<MHDFloat>(n+1), -1.0/static_cast<MHDFloat>(n+1));
               if(m == 0)
               {
                  val = std::real(val);
               }
               ptSH.first->second.insert(std::make_pair(n, val));
            }
         }
      }
   };

   std::shared_ptr<Model::IModelBackend> spBackend = std::make_shared<TestBackend>(spRunner->ss().tag());

   // Add state generation equations
   Equations::SharedSphereExactScalarState spScalar;
   Equations::SharedSphereExactVectorState spVector;

   Spectral::Kernel::Complex3DMapType tSH;

   spScalar =
      spRunner->addEquation<Equations::SphereExactScalarState>(spBackend);
   spScalar->setIdentity(PhysicalNames::Temperature::id());
   tSH.clear();
   setScalarModes(tSH, maxN, maxL, maxM, 0);
   spScalar->setSpectralModes(tSH);

   spVector =
      spRunner->addEquation<Equations::SphereExactVectorState>(spBackend);
   spVector->setIdentity(PhysicalNames::Velocity::id());
   // Toroidal
   tSH.clear();
   setTorModes(tSH, maxN, maxL, maxM, 1);
   spVector->setSpectralModes(FieldComponents::Spectral::TOR, tSH);
   // Poloidal
   tSH.clear();
   setPolModes(tSH, maxN, maxL, maxM, 1);
   spVector->setSpectralModes(FieldComponents::Spectral::POL, tSH);

   spVector =
      spRunner->addEquation<Equations::SphereExactVectorState>(spBackend);
   spVector->setIdentity(PhysicalNames::Magnetic::id());
   // Toroidal
   tSH.clear();
   setTorModes(tSH, maxN, maxL, maxM, 1);
   spVector->setSpectralModes(FieldComponents::Spectral::TOR, tSH);
   // Poloidal
   tSH.clear();
   setPolModes(tSH, maxN, maxL, maxM, 1);
   spVector->setSpectralModes(FieldComponents::Spectral::POL, tSH);

   return spVector;
}

#if 0
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
#if 0
   //
   {
      auto spFile = std::make_shared<QuICC::Io::Variable::ShellTorPolEnstrophyWriter>("velocity", spRunner->ss().tag());
      spFile->expect(PhysicalNames::Velocity::id());
      spRunner->addAsciiOutputFile(spFile);
   }
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
#endif

void checkFiles(std::shared_ptr<SpatialScheme::ISpatialScheme> spScheme, const TestParameters& test)
{
   const std::string& datadir = test.datadir;
   const std::string& refdir = test.refdir;

   const auto& maxUlp = test.maxUlp;

   auto compareRef = [](std::string fname, int rows, int cols, int blocks, const std::string& datadir, const std::string& refdir, const int maxUlp)
   {
      MHDFloat eps = std::numeric_limits<MHDFloat>::epsilon();
      MHDFloat tol = maxUlp*eps;

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
      bool foundFiles = true;
      for(auto&& blk: data)
      {
         CHECK( blk.size() > 0 );
         if(blk.size() == 0)
         {
            foundFiles = false;
         }
      }
      for(auto&& blk: ref)
      {
         CHECK( blk.size() > 0 );
         if(blk.size() == 0)
         {
            foundFiles = false;
         }
      }
      if(foundFiles)
      {
         for(int k = 0; k < blocks; k++)
         {
            const auto& d = data.at(k);
            const auto& r = ref.at(k);
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
      }
   };

   std::vector<std::tuple<std::string,int,int,int>> fileList;

   // Add sphere ascii files
   if(spScheme->tag() == "WLFl" || spScheme->tag() == "WLFm")
   {
      fileList = checkSphereFiles(test);
   }
   // Add shell ascii files
   else if(spScheme->tag() == "SLFl" || spScheme->tag() == "SLFm")
   {
      fileList = checkShellFiles(test);
   }
   else
   {
      throw std::logic_error("No ASCII files have been setup for " + spScheme->tag() + " scheme");
   }

   // Check all files
   for(auto&& f: fileList)
   {
      std::string fname = std::get<0>(f);
      int r = std::get<1>(f);
      int c = std::get<2>(f);
      int b = std::get<3>(f);
      compareRef(fname, r, c, b, datadir, refdir, maxUlp);
   }
}

std::vector<std::tuple<std::string,int,int,int>> checkSphereFiles(const TestParameters& test)
{
   const auto& nN = test.spRes->sim().dim(QuICC::Dimensions::Simulation::SIM1D, QuICC::Dimensions::Space::SPECTRAL);
   const auto& nL = test.spRes->sim().dim(QuICC::Dimensions::Simulation::SIM2D, QuICC::Dimensions::Space::SPECTRAL);
   const auto& nM = test.spRes->sim().dim(QuICC::Dimensions::Simulation::SIM3D, QuICC::Dimensions::Space::SPECTRAL);
   const auto& nR = test.spRes->sim().dim(QuICC::Dimensions::Simulation::SIM1D, QuICC::Dimensions::Space::PHYSICAL);
   int nH = nL*(nL+1)/2;

   // List of files to check: fname, rows, cols, blocks
   std::vector<std::tuple<std::string,int,int,int>> fileList;
   for(auto&& f: test.files)
   {
      std::string n = f->filename();

      if(n == "velocityangular_momentum.dat")
      {
         fileList.emplace_back(n, 1, 5, 1);
      }
      else if(n == "magnetic_dipolarity.dat")
      {
         fileList.emplace_back(n, 1,  5, 1);
      }
      else if(n == "velocity_max.dat")
      {
         fileList.emplace_back(n, 1, 2, 1);
      }
      else if(n == "temperaturenusselt.dat")
      {
         fileList.emplace_back(n, 1, 2, 1);
      }
      else if(n == "temperature_energy.dat")
      {
         fileList.emplace_back(n, 1, 2, 1);
      }
      else if(n == "temperature_l_spectrum.dat")
      {
         fileList.emplace_back(n, nL, 2, 1);
      }
      else if(n == "temperature_m_spectrum.dat")
      {
         fileList.emplace_back(n, nM, 2, 1);
      }
      else if(n == "temperature_mean.dat")
      {
         fileList.emplace_back(n, 1, 2, 1);
      }
      else if(n == "temperature_mode_spectrum.dat")
      {
         fileList.emplace_back(n, nH, 3, 1);
      }
      else if(n == "temperature_n_spectrum.dat")
      {
         fileList.emplace_back(n, nN, nL + 1, 1);
      }
      else if(n == "temperature_r_spectrum.dat")
      {
         fileList.emplace_back(n, nR, nL + 1, 1);
      }
      else if(n == "velocity_energy.dat")
      {
         fileList.emplace_back(n, 1, 4, 1);
      }
      else if(n == "velocity_l_spectrum.dat")
      {
         fileList.emplace_back(n, nL, 4, 1);
      }
      else if(n == "velocity_m_spectrum.dat")
      {
         fileList.emplace_back(n, nM, 4, 1);
      }
      else if(n == "velocity_mode_spectrum.dat")
      {
         fileList.emplace_back(n, nH, 5, 1);
      }
      else if(n == "velocity_n_spectrum.dat")
      {
         fileList.emplace_back(n, nN, nL + 1, 3);
      }
      else if(n == "velocity_r_spectrum.dat")
      {
         fileList.emplace_back(n, nR, nL + 1, 3);
      }
      else if(n == "velocity_enstrophy.dat")
      {
         fileList.emplace_back(n, 1, 4, 1);
      }
      else if(n == "velocity_enstrophy_l_spectrum.dat")
      {
         fileList.emplace_back(n, nL, 4, 1);
      }
      else if(n == "velocity_enstrophy_m_spectrum.dat")
      {
         fileList.emplace_back(n, nM, 4, 1);
      }
      else
      {
         throw std::logic_error("Could not identify file " + n);
      }
   }

   return fileList;
}

std::vector<std::tuple<std::string,int,int,int>> checkShellFiles(const TestParameters& test)
{
   const auto& nN = test.spRes->sim().dim(QuICC::Dimensions::Simulation::SIM1D, QuICC::Dimensions::Space::SPECTRAL);
   const auto& nL = test.spRes->sim().dim(QuICC::Dimensions::Simulation::SIM2D, QuICC::Dimensions::Space::SPECTRAL);
   const auto& nM = test.spRes->sim().dim(QuICC::Dimensions::Simulation::SIM3D, QuICC::Dimensions::Space::SPECTRAL);
   const auto& nR = test.spRes->sim().dim(QuICC::Dimensions::Simulation::SIM1D, QuICC::Dimensions::Space::PHYSICAL);
   int nH = nL*(nL+1)/2;

   // List of files to check: fname, rows, cols, blocks
   std::vector<std::tuple<std::string,int,int,int>> fileList;
   fileList.emplace_back("nusselt.dat", 1, 2, 1);
   fileList.emplace_back("temperature_energy.dat", 1, 2, 1);
   fileList.emplace_back("temperature_l_spectrum.dat", nL, 2, 1);
   fileList.emplace_back("temperature_m_spectrum.dat", nM, 2, 1);
   fileList.emplace_back("velocity_energy.dat", 1, 4, 1);
   fileList.emplace_back("velocity_l_spectrum.dat", nL, 4, 1);
   fileList.emplace_back("velocity_m_spectrum.dat", nM, 4, 1);
#if 0
   fileList.emplace_back("velocity_enstrophy.dat", 1, 4, 1);
   fileList.emplace_back("velocity_enstrophy_l_spectrum.dat", nL, 4, 1);
   fileList.emplace_back("velocity_enstrophy_m_spectrum.dat", nM, 4, 1);
#endif

   return fileList;
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

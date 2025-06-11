/**
 * @file TestHelper.cpp
 * @brief Setup helpers for the TransformCoordinator
 */

// System includes
//
#include <catch2/catch.hpp>
#include <fstream>
#include <limits>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/NonDimensional/Lower1d.hpp"
#include "QuICC/NonDimensional/Upper1d.hpp"
#include "QuICC/SpatialScheme/Feature.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/SpatialScheme/Tools/SpectralTriangularSH.hpp"
#include "QuICC/TestSuite/Framework/TransformCoordinator/TestHelper.hpp"
#include "QuICC/TransformCoordinators/TransformCoordinatorTools.hpp"
#include "QuICC/TransformConfigurators/TransformTreeTools.hpp"
#include "QuICC/TransformConfigurators/TransformStepsFactory.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/Transform/Path/Empty.hpp"
#include "QuICC/Transform/Path/Scalar.hpp"
#include "QuICC/Transform/Path/ValueScalar.hpp"
#include "QuICC/Transform/Path/InsulatingScalar.hpp"
#include "QuICC/Transform/Path/TorPol.hpp"
#include "QuICC/Transform/Path/ValueTorPol.hpp"
#include "QuICC/Transform/Path/NoPenetrationTorPol.hpp"
#include "QuICC/Transform/Path/InsulatingTorPol.hpp"
#include "QuICC/Transform/Path/ScalarNl.hpp"
#include "QuICC/Transform/Path/ValueScalarNl.hpp"
#include "QuICC/Transform/Path/InsulatingScalarNl.hpp"
#include "QuICC/Transform/Path/CurlNl.hpp"
#include "QuICC/Transform/Path/CurlCurlNl.hpp"
#include "QuICC/Transform/Path/ValueCurlNl.hpp"
#include "QuICC/Transform/Path/ValueCurlCurlNl.hpp"
#include "QuICC/Transform/Path/InsulatingCurlNl.hpp"
#include "QuICC/Transform/Path/InsulatingCurlCurlNl.hpp"
#include "QuICC/TypeSelectors/ParallelSelector.hpp"
#include "QuICC/Variables/FieldRequirement.hpp"
#include "QuICC/Variables/RequirementTools.hpp"
#include "QuICC/PhysicalKernels/Passthrough.hpp"
#include "QuICC/TestSuite/Framework/TransformCoordinator/TestArgs.hpp"
#include "QuICC/TestSuite/Framework/TransformCoordinator/CurlKernel.hpp"
#include "QuICC/TestSuite/Framework/TransformCoordinator/InertiaKernel.hpp"
#include "QuICC/TestSuite/Framework/TransformCoordinator/IReference.hpp"
#include "TestSuite/Io.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace TCoord {

   ArrayI processCmdLine(Test& test)
   {
      // Set default arguments if required
      if(test.fbase != "")
      {
         Array meta;
         std::string path = test.fbase + "_meta.dat";
         readList(meta, path);

         args().dim1D = static_cast<int>(meta(0));
         args().dim2D = static_cast<int>(meta(1));
         args().dim3D = static_cast<int>(meta(2));

         if(args().algorithm == "")
         {
            args().algorithm = "tubular";
         }

         if(args().grouper == "")
         {
            args().grouper = "transform";
         }

         if(args().truncation == "")
         {
            args().truncation = "uniform";
         }
      }
      else if(args().useDefault)
      {
         args().dim1D = 15;
         args().dim2D = 31;
         args().dim3D = 31;

         if(args().algorithm == "")
         {
            args().algorithm = "tubular";
         }

         if(args().grouper == "")
         {
            args().grouper = "transform";
         }

         if(args().truncation == "")
         {
            args().truncation = "uniform";
         }

         if(args().truncation == "triangular")
         {
            SpatialScheme::Tools::SpectralTriangularSH t;

            args().dim1D = args().dim2D/2 + t.min() - 1;
         }
         else if(args().truncation == "trapezoidal")
         {
            SpatialScheme::Tools::SpectralTrapezoidalSH t;

            args().dim1D += args().dim2D/2 + t.min() - 1;
         }
      }
      else
      {
         if(args().algorithm == "")
         {
            args().algorithm = "serial";
         }

         if(args().grouper == "")
         {
            args().grouper = "equation";
         }

         if(args().truncation == "")
         {
            args().truncation = "uniform";
         }
      }

      if(args().dim1D == 0 || args().dim2D == 0 || args().dim3D == 0)
      {
         throw std::logic_error("Dimensions are not set properly");
      }

      // Set simulation truncation
      QuICC::ArrayI dim(3);
      dim << args().dim1D, args().dim2D, args().dim3D;

      INFO( "Input parameters" );
      INFO( "dim1D: " << dim(0) );
      INFO( "dim2D: " << dim(1) );
      INFO( "dim3D: " << dim(2) );

      if(args().params.size() == 0)
      {
         args().params.push_back(0);
      }

      return dim;
   }

   void initVariables(Test& test)
   {
      auto& ss = *test.spRes->sim().spSpatialScheme();

      VariableRequirement info;

      if(test.fieldId == Test::FieldId::SCALAR || test.fieldId == Test::FieldId::SCALAR_AND_TORPOL)
      {
         auto& tempReq = info.addField(PhysicalNames::Temperature::id(), FieldRequirement(true, ss.spectral(), ss.physical()));
         tempReq.enableSpectral();
         tempReq.enablePhysical();
      }

      if(test.fieldId == Test::FieldId::TOR ||
         test.fieldId == Test::FieldId::POL ||
         test.fieldId == Test::FieldId::TORPOL ||
         test.fieldId == Test::FieldId::SCALAR_AND_TORPOL)
      {
         auto& velReq = info.addField(PhysicalNames::Velocity::id(), FieldRequirement(false, ss.spectral(), ss.physical()));
         velReq.enableSpectral();
         velReq.enablePhysical();

         if(test.kernelId == Test::KernelId::INERTIA || test.kernelId == Test::KernelId::CURL)
         {
            velReq.enableCurl();
         }
      }

      RequirementTools::initVariables(test.scalars, test.vectors, info, test.spRes);
   }

   void initKernels(Test& test)
   {
      if(test.kernelId == Test::KernelId::PASSTHROUGH)
      {
         // Use passthrough kernel for scalars
         for(auto&& f: test.scalars)
         {
            auto spKernel = std::make_shared<Physical::Kernel::Passthrough>();
            spKernel->setField(f.first, f.second);
            test.kernels.emplace(f.first, spKernel);
         }

         // Use passthrough kernel for vectors
         for(auto&& f: test.vectors)
         {
            auto spKernel = std::make_shared<Physical::Kernel::Passthrough>();
            spKernel->setField(f.first, f.second);
            test.kernels.emplace(f.first, spKernel);
         }
      }
      else if(test.kernelId == Test::KernelId::INERTIA)
      {
         if(test.scalars.size() > 0)
         {
            throw std::logic_error("Inertial kernel not implemented for scalar");
         }

         // Use inertia kernel for vectors
         for(auto&& f: test.vectors)
         {
            auto spKernel = std::make_shared<Physical::Kernel::InertiaKernel>();
            spKernel->setVelocity(f.first, f.second);
            spKernel->init(1.0);
            test.kernels.emplace(f.first, spKernel);
         }
      }
      else if(test.kernelId == Test::KernelId::CURL)
      {
         if(test.scalars.size() > 0)
         {
            throw std::logic_error("Curl kernel not implemented for scalar");
         }

         // Use inertia kernel for vectors
         for(auto&& f: test.vectors)
         {
            auto spKernel = std::make_shared<Physical::Kernel::CurlKernel>();
            spKernel->setVelocity(f.first, f.second);
            spKernel->init(1.0);
            test.kernels.emplace(f.first, spKernel);
         }
      }
   }

   void initTrees(Test& test)
   {
      std::map<size_t, std::vector<Transform::TransformPath> > mt;

      auto spSteps = Transform::createTransformSteps(test.spRes->sim().spSpatialScheme());
      
      std::size_t sFwdPathId;
      std::size_t sBwdPathId;
      std::size_t torFwdPathId;
      std::size_t polFwdPathId;
      std::size_t vBwdPathId;
      bool isNonlinearPath;
      switch(test.pathId)
      {
         case Test::PathId::BFLOOP:
            isNonlinearPath = false;
            sFwdPathId = Transform::Path::Scalar::id();
            sBwdPathId = Transform::Path::Scalar::id();
            torFwdPathId = Transform::Path::TorPol::id();
            polFwdPathId = Transform::Path::TorPol::id();
            vBwdPathId = Transform::Path::TorPol::id();
            break;
         case Test::PathId::NLLOOP:
            isNonlinearPath = true;
            sFwdPathId = Transform::Path::ScalarNl::id();
            sBwdPathId = Transform::Path::Scalar::id();
            torFwdPathId = Transform::Path::CurlNl::id();
            polFwdPathId = Transform::Path::CurlCurlNl::id();
            vBwdPathId = Transform::Path::TorPol::id();
            break;
         default:
            throw std::logic_error("Test for this transform path id = " + std::to_string(static_cast<int>(test.pathId)) + " is not implemented");
      }

      for(auto&& f: test.scalars)
      {
         // Create forward scalar transform tree
         std::vector<Transform::ITransformSteps::PathId> comps = {{FieldComponents::Spectral::SCALAR,sFwdPathId}};
         if(isNonlinearPath)
         {
            auto t = spSteps->forwardNLScalar(comps);
            mt.insert(std::make_pair(f.first, t));
         }
         else
         {
            auto t = spSteps->forwardScalar(comps);
            mt.insert(std::make_pair(f.first, t));
         }
      }

      for(auto&& f: test.vectors)
      {
         // Create forward vector transform tree
         std::vector<Transform::ITransformSteps::PathId> comps = {{FieldComponents::Spectral::TOR,torFwdPathId},{FieldComponents::Spectral::POL,polFwdPathId}};
         if(isNonlinearPath)
         {
            auto t = spSteps->forwardNLVector(comps);
            mt.insert(std::make_pair(f.first, t));
         }
         else
         {
            auto t = spSteps->forwardVector(comps);
            mt.insert(std::make_pair(f.first, t));
         }
      }

      Transform::TransformTreeTools::generateTrees(test.fwdTree, mt, TransformDirection::FORWARD);

      // Create backward transform tree based on variables
      RequirementTools::buildBackwardTree(test.bwdTree, test.scalars, test.vectors, sBwdPathId, vBwdPathId);
   }

   void initCoordinator(Test& test, const Parallel::SplittingDescription& descr)
   {
      // Set grouper
      Parallel::setGrouper(descr, test.spFwdGrouper, test.spBwdGrouper);

      // Initialize coordinator
      std::map<std::size_t,NonDimensional::SharedINumber> runOptions;
      runOptions.emplace(NonDimensional::Lower1d::id(), std::make_shared<NonDimensional::Lower1d>(0.3));
      runOptions.emplace(NonDimensional::Upper1d::id(), std::make_shared<NonDimensional::Upper1d>(1.3));
      std::vector<ArrayI> packs;
      Transform::TransformCoordinatorTools::computePacks(packs, test.spFwdGrouper, test.spBwdGrouper, {{0, test.fwdTree}}, {{0, test.bwdTree}}, {0}, test.spRes);
      Transform::TransformCoordinatorTools::init(test.coord, test.spFwdGrouper, test.spBwdGrouper, packs, test.spRes, runOptions);
   }

   void setVariables(Test& test)
   {
      // Set unit spectrum for scalar fields
      for(auto&& f: test.scalars)
      {
         std::visit(
               [&](auto&& p)
               {
                  const auto& tRes = *test.spRes->cpu()->dim(Dimensions::Transform::SPECTRAL);
                  const auto& sRes = test.spRes->sim();
                  for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
                  {
                     auto k_ = tRes.idx<Dimensions::Data::DAT3D>(k);
                     for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
                     {
                        auto j_ = tRes.idx<Dimensions::Data::DAT2D>(j, k);
                        for(int i = 0; i < tRes.dim<Dimensions::Data::DATB1D>(j,k); i++)
                        {
                           auto val = test.spRef->inScalar(test, i, j_, k_);
                           p->rDom(0).rPerturbation().setPoint(val, i,j,k);
                        }
                     }
                  }
               },
            f.second);
      }

      // Set unit spectrum for vector fields
      for(auto&& f: test.vectors)
      {
         std::visit(
               [&](auto&& p)
               {
                  const auto& tRes = *test.spRes->cpu()->dim(Dimensions::Transform::SPECTRAL);
                  const auto& sRes = test.spRes->sim();
                  for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
                  {
                     auto k_ = tRes.idx<Dimensions::Data::DAT3D>(k);
                     for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
                     {
                        auto j_ = tRes.idx<Dimensions::Data::DAT2D>(j, k);
                        for(int i = 0; i < tRes.dim<Dimensions::Data::DATF1D>(j,k); i++)
                        {
                           // Toroidal component
                           if(test.fieldId == Test::FieldId::TOR ||
                              test.fieldId == Test::FieldId::TORPOL ||
                              test.fieldId == Test::FieldId::SCALAR_AND_TORPOL)
                           {
                              auto val = test.spRef->inTor(test, i, j_, k_);
                              p->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::TOR).setPoint(val, i,j,k);
                           }

                           // Poloidal component
                           if(test.fieldId == Test::FieldId::POL ||
                              test.fieldId == Test::FieldId::TORPOL ||
                              test.fieldId == Test::FieldId::SCALAR_AND_TORPOL)
                           {
                              auto val = test.spRef->inPol(test, i, j_, k_);
                              p->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::POL).setPoint(val, i,j,k);
                           }
                        }
                     }
                  }
                  std::cerr << sRes.dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL) << std::endl;
                  std::cerr << sRes.dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL) << std::endl;
                  std::cerr << sRes.dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL) << std::endl;
                  std::cerr << sRes.dim(Dimensions::Simulation::SIM1D, Dimensions::Space::PHYSICAL) << std::endl;
                  std::cerr << sRes.dim(Dimensions::Simulation::SIM2D, Dimensions::Space::PHYSICAL) << std::endl;
                  std::cerr << sRes.dim(Dimensions::Simulation::SIM3D, Dimensions::Space::PHYSICAL) << std::endl;
                  std::cerr << std::string(50, '&') << std::endl;
                  std::cerr << p->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::TOR).data() << std::endl;
                  std::cerr << std::string(50, '@') << std::endl;
                  std::cerr << p->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::POL).data() << std::endl;
                  std::cerr << std::string(50, 'q') << std::endl;
               },
            f.second);
      }
   }

   void scrambleVariables(Test& test)
   {
      // Set bad values for scalar fields
      for(auto&& f: test.scalars)
      {
         std::visit(
               [&](auto&& p)
               {
                  MHDComplex val = std::numeric_limits<MHDFloat>::max()/2.0;
                  const auto& tRes = *test.spRes->cpu()->dim(Dimensions::Transform::SPECTRAL);
                  const auto& sRes = test.spRes->sim();
                  for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
                  {
                     for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
                     {
                        for(int i = 0; i < sRes.dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL); i++)
                        {
                           p->rDom(0).rPerturbation().setPoint(val, i,j,k);
                        }
                     }
                  }
               },
            f.second);
      }

      // Set bad values for vector fields
      for(auto&& f: test.vectors)
      {
         std::visit(
               [&](auto&& p)
               {
                  MHDComplex val = std::numeric_limits<MHDFloat>::max()/2.0;
                  const auto& tRes = *test.spRes->cpu()->dim(Dimensions::Transform::SPECTRAL);
                  const auto& sRes = test.spRes->sim();
                  for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
                  {
                     for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
                     {
                        for(int i = 0; i < sRes.dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL); i++)
                        {
                           p->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::TOR).setPoint(val, i,j,k);
                           p->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::POL).setPoint(val, i,j,k);
                        }
                     }
                  }
               },
            f.second);
      }
   }

   void checkVariables(Test& test)
   {
      // Set unit spectrum for scalar fields
      for(auto&& f: test.scalars)
      {
         MHDFloat worstUlp = 0.0;
         std::visit(
               [&](auto&& p)
               {
                  const auto& tRes = *test.spRes->cpu()->dim(Dimensions::Transform::SPECTRAL);
                  const auto& sRes = test.spRes->sim();
                  for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
                  {
                     auto k_ = tRes.idx<Dimensions::Data::DAT3D>(k);
                     for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
                     {
                        auto j_ = tRes.idx<Dimensions::Data::DAT2D>(j, k);
                        for(int i = 0; i < tRes.dim<Dimensions::Data::DATF1D>(j,k); i++)
                        {
                           auto ref = test.spRef->refScalar(test, i, j_, k_);
                           INFO( "Checking Scalar" );
                           auto data = p->dom(0).perturbation().point(i,j,k);
                           auto err = computeUlp(data, ref, std::abs(ref), test.tolerance(), test.epsilon);
                           worstUlp = std::max(worstUlp, std::get<1>(err));
                           INFO( "i,j,k: " << i << "," << j_ << "," << k_ );
                           INFO( "data: " << std::scientific << std::setprecision(16) << data );
                           INFO( "ref: " << std::scientific << std::setprecision(16) << ref );
                           INFO( "max ulp: " << test.maxUlp);
                           INFO( "measured ulp: " << std::get<1>(err) );
                           CHECK( std::get<0>(err) );
                        }
                     }
                  }
               },
            f.second);
         std::cerr << "Max measured ULP for Scalar: " << worstUlp << std::endl;
      }

      // Set unit spectrum for vector fields
      for(auto&& f: test.vectors)
      {
         MHDFloat worstUlpTor = 0.0;
         MHDFloat worstUlpPol = 0.0;
         std::visit(
               [&](auto&& p)
               {
               std::cerr << std::string(80,'-') << std::endl;
               std::cerr << p->dom(0).perturbation().comp(FieldComponents::Spectral::TOR).data() << std::endl;
               std::cerr << std::string(80,'_') << std::endl;
               std::cerr << p->dom(0).perturbation().comp(FieldComponents::Spectral::POL).data() << std::endl;
               std::cerr << std::string(80,'=') << std::endl;
                  const auto& tRes = *test.spRes->cpu()->dim(Dimensions::Transform::SPECTRAL);
                  const auto& sRes = test.spRes->sim();
                  for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
                  {
                     auto k_ = tRes.idx<Dimensions::Data::DAT3D>(k);
                     for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
                     {
                        auto j_ = tRes.idx<Dimensions::Data::DAT2D>(j, k);
                        for(int i = 0; i < tRes.dim<Dimensions::Data::DATF1D>(j,k); i++)
                        {
                           auto checkFct = [&](auto comp, auto ref, auto scale, auto& worstUlp, const std::string name)
                           {
                              auto data = p->dom(0).perturbation().comp(comp).point(i,j,k);
                              auto err = computeUlp(data, ref, std::abs(scale), test.tolerance(), test.epsilon);
                              worstUlp = std::max(worstUlp, std::get<1>(err));
                              INFO( "Checking " + name + " component" );
                              INFO( "i,j,k: " << i << "," << j_ << "," << k_ );
                              INFO( "data: " << std::scientific << std::setprecision(16) << data );
                              INFO( "ref: " << std::scientific << std::setprecision(16) << ref );
                              INFO( "max ulp: " << test.maxUlp);
                              INFO( "measured ulp: " << std::get<1>(err) );
                              CHECK( std::get<0>(err) );
                           };

                           auto refTor = test.spRef->refTor(test, i, j_, k_);
                           auto scaleTor = refTor;
                           auto refPol = test.spRef->refPol(test, i, j_, k_);
                           auto scalePol = refPol;

                           if(refTor == 0.0)
                           {
                              scaleTor = 1.0;
                           }
                           else if(refPol == 0.0)
                           {
                              scalePol = 1.0;
                           }
                           checkFct(FieldComponents::Spectral::TOR, refTor, scaleTor, worstUlpTor, "Toroidal");
                           checkFct(FieldComponents::Spectral::POL, refPol, scalePol, worstUlpPol, "Poloidal");
                        }
                     }
                  }
               },
            f.second);
         std::cerr << "Max measured ULP for Toroidal: " << worstUlpTor << std::endl;
         std::cerr << "Max measured ULP for Poloidal: " << worstUlpPol << std::endl;
      }
   }

   void backward(Test& test)
   {
      test.coord.defineBwdTransforms(test.bwdTree);
      test.spBwdGrouper->transform(test.scalars, test.vectors, test.coord);
   }

   void nonlinearAndForward(Test& test)
   {
      test.coord.defineFwdTransforms(test.fwdTree);
      test.spFwdGrouper->transform(test.scalars, test.vectors, test.kernels, test.coord);
   }

   ErrorType computeUlp(const MHDComplex data, const MHDComplex ref, const MHDFloat refMod, const MHDFloat tol, const MHDFloat eps)
   {
      auto errReal = computeUlp(data.real(), ref.real(), refMod, tol, eps);
      auto errImag = computeUlp(data.real(), ref.real(), refMod, tol, eps);
      if(std::get<1>(errReal) > std::get<1>(errImag))
      {
         return errReal;
      }
      else
      {
         return errImag;
      }
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

}
}
}
}

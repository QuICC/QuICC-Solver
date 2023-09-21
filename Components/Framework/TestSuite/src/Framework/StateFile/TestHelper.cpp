/**
 * @file TestHelper.cpp
 * @brief Setup helpers for the StateFile tests
 */

// System includes
//
#include <catch2/catch.hpp>
#include <fstream>
#include <limits>

// Project includes
//
#include "QuICC/TestSuite/Framework/StateFile/TestHelper.hpp"
#include "QuICC/TestSuite/Framework/StateFile/TestArgs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/NonDimensional/Lower1d.hpp"
#include "QuICC/NonDimensional/Upper1d.hpp"
#include "QuICC/SpatialScheme/Feature.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/Transform/Path/Scalar.hpp"
#include "QuICC/Transform/Path/TorPol.hpp"
#include "QuICC/TypeSelectors/ParallelSelector.hpp"
#include "QuICC/Variables/FieldRequirement.hpp"
#include "QuICC/Variables/RequirementTools.hpp"
#include "QuICC/PhysicalKernels/Passthrough.hpp"
#include "QuICC/Io/Variable/StateFileReader.hpp"
#include "QuICC/Io/Variable/StateFileWriter.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace StateFile {

   Test::Test()
      : epsilon(std::numeric_limits<MHDFloat>::epsilon()), maxUlp(11), fileTag(""), refFileTag("")
   {
   }

   MHDFloat Test::tolerance() const
   {
      return this->maxUlp*this->epsilon;
   }

   void Test::configure(const int id)
   {
      switch(id)
      {
         case 0:
            this->fieldId = FieldId::SCALAR;
            this->spectrumId = SpectrumId::TAG;
            break;
         case 1:
            this->fieldId = FieldId::TOR;
            this->spectrumId = SpectrumId::TAG;
            break;
         case 2:
            this->fieldId = FieldId::POL;
            this->spectrumId = SpectrumId::TAG;
            break;
         case 3:
            this->fieldId = FieldId::TORPOL;
            this->spectrumId = SpectrumId::TAG;
            break;
         case 4:
            this->fieldId = FieldId::SCALAR_AND_TORPOL;
            this->spectrumId = SpectrumId::TAG;
            break;
         default:
            throw std::logic_error("Undefined test case was requested");
      }
   }

   ArrayI processCmdLine()
   {
      // Set default arguments if required
      if(args().useDefault)
      {
         args().dim1D = 63;
         args().dim2D = 127;
         args().dim3D = 127;

         if(args().algorithm == "")
         {
            args().algorithm = "tubular";
         }

         if(args().grouper == "")
         {
            args().grouper = "transform";
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
      CHECK( dim(0) > 0 );
      CHECK( dim(1) > 0 );
      CHECK( dim(2) > 0 );

      if(args().params.size() == 0)
      {
         std::cerr << "Setting zero" << std::endl;
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
      }

      RequirementTools::initVariables(test.scalars, test.vectors, info, test.spRes);
   }

   MHDComplex tagReference(const Test& test, const int i, const int j, const int k)
   {
      auto&& ss = test.spRes->sim().ss();
      if(ss.has(SpatialScheme::Feature::SphereGeometry) || ss.has(SpatialScheme::Feature::ShellGeometry))
      {
         if(ss.has(SpatialScheme::Feature::SpectralOrdering123))
         {
            return tagReferenceSH(i,j,k);
         }
         else
         {
            return tagReferenceSH(i,k,j);
         }
      }
      else if(ss.has(SpatialScheme::Feature::CartesianGeometry) && ss.has(SpatialScheme::Feature::FourierIndex23))
      {
         return tagReferenceFF(i,j,k);
      }
      else
      {
         throw std::logic_error("Unit spectrum for this geometry has not been implemented");
      }
   }

   MHDComplex tagReferenceSH(const int n, const int l, const int m)
   {
      MHDFloat rTag = static_cast<MHDFloat>(n) + static_cast<MHDFloat>(l)/1e4 + static_cast<MHDFloat>(m)/1e8;
      MHDFloat iTag = static_cast<MHDFloat>(l) + static_cast<MHDFloat>(m)/1e4 + static_cast<MHDFloat>(n)/1e8;
      MHDComplex ref(rTag,iTag);

      return ref;
   }

   MHDComplex tagReferenceFF(const int n, const int k1, const int k2)
   {
      MHDComplex ref(std::sqrt(2.0),-std::sqrt(2.0));

      if(k1 == 0 && k2 == 0)
      {
         ref.imag(0.0);
      }

      return ref;
   }

   void setVariables(Test& test)
   {
      MHDComplex (*refFct)(const Test& tet, int,int,int);
      if(test.spectrumId == Test::SpectrumId::TAG)
      {
         refFct = &tagReference;
      }
      else
      {
         throw std::logic_error("Reference type not implemented");
      }

      // Set unit spectrum for scalar fields
      for(auto&& f: test.scalars)
      {
         std::visit(
               [&](auto&& p)
               {
                  const auto& tRes = *test.spRes->cpu()->dim(Dimensions::Transform::SPECTRAL);
                  for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
                  {
                     auto k_ = tRes.idx<Dimensions::Data::DAT3D>(k);
                     for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
                     {
                        auto j_ = tRes.idx<Dimensions::Data::DAT2D>(j, k);
                        for(int i = 0; i < tRes.dim<Dimensions::Data::DATB1D>(j,k); i++)
                        {
                           auto val = (*refFct)(test, i, j_, k_);
                           p->rDom(0).rPerturbation().setPoint(val, i,j,k);
                        }

                        // Fill extra space with bad values for non-uniform
                        // truncation
                        for(int i = tRes.dim<Dimensions::Data::DATB1D>(j,k); i < p->rDom(0).rPerturbation().slice(k).rows(); i++)
                        {
                           auto val = (*refFct)(test, i, j_, k_);
                           val = std::numeric_limits<MHDFloat>::max();
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
                  for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
                  {
                     auto k_ = tRes.idx<Dimensions::Data::DAT3D>(k);
                     for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
                     {
                        auto iN = tRes.dim<Dimensions::Data::DATB1D>(j,k);
                        auto j_ = tRes.idx<Dimensions::Data::DAT2D>(j, k);
                        for(int i = 0; i < iN; i++)
                        {
                           auto val = (*refFct)(test, i, j_, k_);

                           // Toroidal component
                           if(test.fieldId == Test::FieldId::TOR ||
                              test.fieldId == Test::FieldId::TORPOL ||
                              test.fieldId == Test::FieldId::SCALAR_AND_TORPOL)
                           {
                              p->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::TOR).setPoint(val, i,j,k);
                           }

                           // Poloidal component
                           if(test.fieldId == Test::FieldId::POL ||
                              test.fieldId == Test::FieldId::TORPOL ||
                              test.fieldId == Test::FieldId::SCALAR_AND_TORPOL)
                           {
                              p->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::POL).setPoint(val, i,j,k);
                           }
                        }

                        // Fill extra space with bad values for non-uniform
                        // truncation
                        int sN = 0;
                        // Toroidal component
                        if(test.fieldId == Test::FieldId::TOR ||
                           test.fieldId == Test::FieldId::TORPOL ||
                           test.fieldId == Test::FieldId::SCALAR_AND_TORPOL)
                        {
                           sN = p->dom(0).perturbation().comp(FieldComponents::Spectral::TOR).slice(k).rows();
                        }

                        // Poloidal component
                        if(test.fieldId == Test::FieldId::POL ||
                           test.fieldId == Test::FieldId::TORPOL ||
                           test.fieldId == Test::FieldId::SCALAR_AND_TORPOL)
                        {
                           sN = p->dom(0).perturbation().comp(FieldComponents::Spectral::POL).slice(k).rows();
                        }
                        for(int i = iN; i < sN; i++)
                        {
                           auto val = (*refFct)(test, i, j_, k_);
                           val = std::numeric_limits<MHDFloat>::max();

                           // Toroidal component
                           if(test.fieldId == Test::FieldId::TOR ||
                              test.fieldId == Test::FieldId::TORPOL ||
                              test.fieldId == Test::FieldId::SCALAR_AND_TORPOL)
                           {
                              p->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::TOR).setPoint(val, i,j,k);
                           }

                           // Poloidal component
                           if(test.fieldId == Test::FieldId::POL ||
                              test.fieldId == Test::FieldId::TORPOL ||
                              test.fieldId == Test::FieldId::SCALAR_AND_TORPOL)
                           {
                              p->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::POL).setPoint(val, i,j,k);
                           }
                        }
                     }
                  }
               },
            f.second);
      }
   }

   void writeStateFile(Test& test)
   {
      auto& ss = *test.spRes->sim().spSpatialScheme();

      std::string postfix = "_" + ss.tag() + test.fileTag;
      auto spState = std::make_shared<Io::Variable::StateFileWriter>(postfix, ss.tag(), ss.has(SpatialScheme::Feature::RegularSpectrum));

      // Add scalars
      for(auto const& f: test.scalars)
      {
         spState->expect(f.first);
         spState->addScalar(f);
      }

      // Add vectors
      for(auto const& f: test.vectors)
      {
         spState->expect(f.first);
         spState->addVector(f);
      }

      if(!spState->isFull())
      {
         throw std::logic_error("Fields are missing in HDF5 file");
      }
      spState->init();
      spState->setSimTime(1.0, 1e-3);
      spState->write();
      spState->finalize();
   }

   void readStateFile(Test& test)
   {
      auto& ss = *test.spRes->sim().spSpatialScheme();

      std::string postfix = "_" + ss.tag() + test.fileTag + test.refFileTag;
      auto spState = std::make_shared<Io::Variable::StateFileReader>(postfix, ss.tag(), ss.has(SpatialScheme::Feature::RegularSpectrum));

      // Add scalars
      for(auto const& f: test.scalars)
      {
         spState->expect(f.first);
         spState->addScalar(f);
      }

      // Add vectors
      for(auto const& f: test.vectors)
      {
         spState->expect(f.first);
         spState->addVector(f);
      }

      if(!spState->isFull())
      {
         throw std::logic_error("Fields are missing in HDF5 file");
      }
      spState->init();
      spState->read();
      spState->finalize();
   }

   void checkVariables(Test& test)
   {
      MHDComplex (*refFct)(const Test& tet, int,int,int);
      if(test.spectrumId == Test::SpectrumId::TAG)
      {
         refFct = &tagReference;
      }
      else
      {
         throw std::logic_error("Reference type not implemented");
      }

      // Set unit spectrum for scalar fields
      for(auto&& f: test.scalars)
      {
         MHDFloat worstUlp = 0.0;
         std::visit(
               [&](auto&& p)
               {
                  const auto& tRes = *test.spRes->cpu()->dim(Dimensions::Transform::SPECTRAL);
                  for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
                  {
                     auto k_ = tRes.idx<Dimensions::Data::DAT3D>(k);
                     const int sN = p->dom(0).perturbation().slice(k).rows();
                     for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
                     {
                        auto j_ = tRes.idx<Dimensions::Data::DAT2D>(j, k);
                        const int iN = tRes.dim<Dimensions::Data::DATB1D>(j,k);
                        for(int i = 0; i < iN; i++)
                        {
                           auto ref = (*refFct)(test, i, j_, k_);
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

                        // Checking unused storage
                        for(int i = iN; i < sN; i++)
                        {
                           auto ref = 0.0;
                           INFO( "Checking Scalar" );
                           auto data = p->dom(0).perturbation().point(i,j,k);
                           auto err = computeUlp(data, ref, std::abs(ref), test.tolerance(), test.epsilon);
                           worstUlp = std::max(worstUlp, std::get<1>(err));
                           INFO( "Unused storage shoud be zero" );
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
                  const auto& tRes = *test.spRes->cpu()->dim(Dimensions::Transform::SPECTRAL);
                  for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
                  {
                     auto k_ = tRes.idx<Dimensions::Data::DAT3D>(k);
                     const int sN = p->dom(0).perturbation().comp(FieldComponents::Spectral::TOR).slice(k).rows();
                     for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
                     {
                        auto j_ = tRes.idx<Dimensions::Data::DAT2D>(j, k);
                        const int iN = tRes.dim<Dimensions::Data::DATB1D>(j,k);
                        for(int i = 0; i < iN; i++)
                        {
                           auto ref = (*refFct)(test, i, j_, k_);

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

                           auto refTor = ref;
                           auto scaleTor = ref;
                           auto refPol = ref;
                           auto scalePol = ref;

                           if(test.fieldId == Test::FieldId::TOR)
                           {
                              refPol = 0.0;
                              scalePol = 1.0;
                           }
                           else if(test.fieldId == Test::FieldId::POL)
                           {
                              refTor = 0.0;
                              scaleTor = 1.0;
                           }
                           checkFct(FieldComponents::Spectral::TOR, refTor, scaleTor, worstUlpTor, "Toroidal");
                           checkFct(FieldComponents::Spectral::POL, refPol, scalePol, worstUlpPol, "Poloidal");
                        }

                        // Checking unused storage
                        for(int i = iN; i < sN; i++)
                        {
                           auto ref = 0.0;
                           auto scale= 1.0;

                           auto checkFct = [&](auto comp, auto ref, auto scale, auto& worstUlp, const std::string name)
                           {
                              auto data = p->dom(0).perturbation().comp(comp).point(i,j,k);
                              auto err = computeUlp(data, ref, std::abs(scale), test.tolerance(), test.epsilon);
                              worstUlp = std::max(worstUlp, std::get<1>(err));
                              INFO( "Checking " + name + " component" );
                              INFO( "Unused storage shoud be zero" );
                              INFO( "i,j,k: " << i << "," << j_ << "," << k_ );
                              INFO( "data: " << std::scientific << std::setprecision(16) << data );
                              INFO( "ref: " << std::scientific << std::setprecision(16) << ref );
                              INFO( "max ulp: " << test.maxUlp);
                              INFO( "measured ulp: " << std::get<1>(err) );
                              CHECK( std::get<0>(err) );
                           };

                           checkFct(FieldComponents::Spectral::TOR, ref, scale, worstUlpTor, "Toroidal");
                           checkFct(FieldComponents::Spectral::POL, ref, scale, worstUlpPol, "Poloidal");
                        }
                     }
                  }
               },
            f.second);
         std::cerr << "Max measured ULP for Toroidal: " << worstUlpTor << std::endl;
         std::cerr << "Max measured ULP for Poloidal: " << worstUlpPol << std::endl;
      }
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

} // StateFile
} // Framework
} // TestSuite
} // QuICC

#include <algorithm>
#include <catch2/catch.hpp>
#include <cstdint>

//#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/ScalarFields/ScalarFieldSetup.hpp"
#include "details/Helper.hpp"
#include "QuICC/PhysicalOperators/SphericalSComponent.hpp"
#include "Profiler/Interface.hpp"
#include "Types/Math.hpp"

template <typename T>
   using sflat_t = QuICC::Datatypes::FlatScalarField<T>;
template <typename T>
   using sview_t = QuICC::Datatypes::ViewScalarField<T>;

using TestOp = QuICC::Physical::SphericalSComponent;

TEST_CASE("SphericalSComponent product validation FlatScalarField vs ViewScalarField", "[SphericalSComponentValidation]")
{
   std::vector<std::uint32_t> variants = details::validationVariants();
   for(std::uint32_t varBase: variants)
   {
      INFO( "variant = " << varBase );

      std::size_t dim3D;
      std::size_t dim2D;
      std::size_t dim1D;
      std::vector<std::vector<std::size_t>> idx2D;
      std::vector<std::size_t> idx3D;
      auto spSetup = details::createSetup(dim1D, dim2D, dim3D, idx2D, idx3D, details::SetupType::UniformUp, varBase + 0);

      auto spSFlatA = details::createScalarField<double, sflat_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 100);
      auto spVFlatA = details::createVectorField<double, sflat_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 1);

      auto spSViewA = details::createScalarField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 100);
      auto spVViewA = details::createVectorField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 1);

      QuICC::Array cosTheta(dim2D);
      QuICC::Array sinTheta(dim2D);
      for(std::uint32_t i = 0; i < dim2D; i++)
      {
         auto theta = 2.0*QuICC::Math::PI*static_cast<double>(i+1)/static_cast<double>(dim3D + 1);
         cosTheta(i) = std::cos(theta);
         sinTheta(i) = std::sin(theta);
      }

      details::IdxFunctor idxFunc(idx2D, idx3D);

      constexpr auto R = QuICC::FieldComponents::Physical::R;
      constexpr auto T = QuICC::FieldComponents::Physical::THETA;
      constexpr auto P = QuICC::FieldComponents::Physical::PHI;

      {
         INFO( "Check initialization" );
         details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         details::checkComputation(spSetup, spVFlatA->comp(R), spVViewA->comp(R));
         details::checkComputation(spSetup, spVFlatA->comp(T), spVViewA->comp(T));
         details::checkComputation(spSetup, spVFlatA->comp(P), spVViewA->comp(P));
      }

      std::vector<double> cs = {1.0, 3.0};
      for(auto&& c: cs)
      {
         INFO( "c = " << c );

         // Check set operation
         {
            INFO ( "set operation" );
            TestOp::set(*spSFlatA, idxFunc, cosTheta, sinTheta, *spVFlatA, c);
            TestOp::set(*spSViewA, idxFunc, cosTheta, sinTheta, *spVViewA, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         }

         // Check add operation
         {
            INFO ( "add operation" );
            TestOp::add(*spSFlatA, idxFunc, cosTheta, sinTheta, *spVFlatA, c);
            TestOp::add(*spSViewA, idxFunc, cosTheta, sinTheta, *spVViewA, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         }

         // Check sub operation
         {
            INFO ( "sub operation" );
            TestOp::sub(*spSFlatA, idxFunc, cosTheta, sinTheta, *spVFlatA, c);
            TestOp::sub(*spSViewA, idxFunc, cosTheta, sinTheta, *spVViewA, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         }
      }
   }
}

TEST_CASE("SphericalSComponent product timing", "[SphericalSComponentTiming]")
{
   const std::string testName = "SphericalSComponentTests";

   std::uint32_t itMax = 100;
   std::vector<double> cs = {1.0, 3.0};
   std::vector<std::uint32_t> variants = details::performanceVariants();

   std::string ctag = "";
   for(std::uint32_t varBase: variants)
   {
      std::size_t dim3D;
      std::size_t dim2D;
      std::size_t dim1D;
      std::vector<std::vector<std::size_t>> idx2D;
      std::vector<std::size_t> idx3D;
      auto spSetup = details::createSetup(dim1D, dim2D, dim3D, idx2D, idx3D, details::SetupType::UniformUp, varBase + 0);
      QuICC::Array cosTheta(dim2D);
      QuICC::Array sinTheta(dim2D);
      for(std::uint32_t i = 0; i < dim2D; i++)
      {
         auto theta = 2.0*QuICC::Math::PI*static_cast<double>(i+1)/static_cast<double>(dim3D + 1);
         cosTheta(i) = std::cos(theta);
         sinTheta(i) = std::sin(theta);
      }

      details::IdxFunctor idxFunc(idx2D, idx3D);

      for(auto&& c: cs)
      {
         if(c == 1.0)
         {
            ctag = "C";
         }
         else
         {
            ctag = "";
         }

         // Timing FlatScalarField
         {
            auto spScalarA = details::createScalarField<double, sflat_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 100);
            auto spVectorA = details::createVectorField<double, sflat_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 1);

            for(std::uint32_t it = 0; it < itMax; it++)
            {
               QuICC::Profiler::RegionStart<1>(testName + "::Flat" + ctag + "Set_" + std::to_string(varBase));
               TestOp::set(*spScalarA, idxFunc, cosTheta, sinTheta, *spVectorA, c);
               QuICC::Profiler::RegionStop<1>(testName + "::Flat" + ctag + "Set_" + std::to_string(varBase));

               QuICC::Profiler::RegionStart<1>(testName + "::Flat" + ctag + "Add_" + std::to_string(varBase));
               TestOp::add(*spScalarA, idxFunc, cosTheta, sinTheta, *spVectorA, c);
               QuICC::Profiler::RegionStop<1>(testName + "::Flat" + ctag + "Add_" + std::to_string(varBase));

               QuICC::Profiler::RegionStart<1>(testName + "::Flat" + ctag + "Sub_" + std::to_string(varBase));
               TestOp::sub(*spScalarA, idxFunc, cosTheta, sinTheta, *spVectorA, c);
               QuICC::Profiler::RegionStop<1>(testName + "::Flat" + ctag + "Sub_" + std::to_string(varBase));
            }
         }

         // Timing ViewScalarField
         {
            auto spScalarA = details::createScalarField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 100);
            auto spVectorA = details::createVectorField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 1);

            for(std::uint32_t it = 0; it < itMax; it++)
            {
               QuICC::Profiler::RegionStart<1>(testName + "::View" + ctag + "Set_" + std::to_string(varBase));
               TestOp::set(*spScalarA, idxFunc, cosTheta, sinTheta, *spVectorA, c);
               QuICC::Profiler::RegionStop<1>(testName + "::View" + ctag + "Set_" + std::to_string(varBase));

               QuICC::Profiler::RegionStart<1>(testName + "::View" + ctag + "Add_" + std::to_string(varBase));
               TestOp::add(*spScalarA, idxFunc, cosTheta, sinTheta, *spVectorA, c);
               QuICC::Profiler::RegionStop<1>(testName + "::View" + ctag + "Add_" + std::to_string(varBase));

               QuICC::Profiler::RegionStart<1>(testName + "::View" + ctag + "Sub_" + std::to_string(varBase));
               TestOp::sub(*spScalarA, idxFunc, cosTheta, sinTheta, *spVectorA, c);
               QuICC::Profiler::RegionStop<1>(testName + "::View" + ctag + "Sub_" + std::to_string(varBase));
            }
         }
      }
   }

   CHECK( true );
}

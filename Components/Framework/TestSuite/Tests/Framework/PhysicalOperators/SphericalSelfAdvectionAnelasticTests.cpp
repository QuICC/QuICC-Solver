#include <algorithm>
#include <catch2/catch.hpp>
#include <cstdint>

//#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/ScalarFields/ScalarFieldSetup.hpp"
#include "details/Helper.hpp"
#include "QuICC/PhysicalOperators/SphericalSelfAdvectionAnelastic.hpp"
#include "Profiler/Interface.hpp"
#include "Types/Math.hpp"

template <typename T>
   using sflat_t = QuICC::Datatypes::FlatScalarField<T>;
template <typename T>
   using sview_t = QuICC::Datatypes::ViewScalarField<T>;

using TestOp = QuICC::Physical::SphericalSelfAdvectionAnelastic;

TEST_CASE("SphericalSelfAdvectionAnelastic product validation FlatScalarField vs ViewScalarField", "[SphericalSelfAdvectionAnelasticValidation]")
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
      auto spVFlatB = details::createVectorField<double, sflat_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 2);

      auto spSViewA = details::createScalarField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 100);
      auto spVViewA = details::createVectorField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 1);
      auto spVViewB = details::createVectorField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 2);

      QuICC::Array rho(dim3D);
      for(std::uint32_t i = 0; i < dim3D; i++)
      {
         rho(i) = static_cast<double>(i+1)/static_cast<double>(dim3D + 1);
      }

      QuICC::Array dLogRho(dim3D);
      for(std::uint32_t i = 0; i < dim3D; i++)
      {
         dLogRho(i) = static_cast<double>(7*i+1)/static_cast<double>(dim3D + 1);
      }

      details::IdxFunctor idxFunc(idx2D, idx3D);

      constexpr auto R = QuICC::FieldComponents::Physical::R;
      constexpr auto T = QuICC::FieldComponents::Physical::THETA;
      constexpr auto P = QuICC::FieldComponents::Physical::PHI;

      details::checkComputation(spSetup, *spSFlatA, *spSViewA);
      details::checkComputation(spSetup, spVFlatA->comp(R), spVViewA->comp(R));
      details::checkComputation(spSetup, spVFlatA->comp(T), spVViewA->comp(T));
      details::checkComputation(spSetup, spVFlatA->comp(P), spVViewA->comp(P));

      std::vector<double> cs = {1.0, 3.0};
      for(auto&& c: cs)
      {
         INFO( "c = " << c );

         // Check set operation
         {
            INFO( "Set operation" );

            INFO( "R component" );
            TestOp::set(*spSFlatA, R, idxFunc, rho, dLogRho, *spVFlatA, *spVFlatB, c);
            TestOp::set(*spSViewA, R, idxFunc, rho, dLogRho, *spVViewA, *spVViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);

            INFO( "T component" );
            TestOp::set(*spSFlatA, T, idxFunc, rho, dLogRho, *spVFlatA, *spVFlatB, c);
            TestOp::set(*spSViewA, T, idxFunc, rho, dLogRho, *spVViewA, *spVViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);

            INFO( "P component" );
            TestOp::set(*spSFlatA, P, idxFunc, rho, dLogRho, *spVFlatA, *spVFlatB, c);
            TestOp::set(*spSViewA, P, idxFunc, rho, dLogRho, *spVViewA, *spVViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         }

         // Check add operation
         {
            INFO( "Add operation" );

            INFO( "R component" );
            TestOp::add(*spSFlatA, R, idxFunc, rho, dLogRho, *spVFlatA, *spVFlatB, c);
            TestOp::add(*spSViewA, R, idxFunc, rho, dLogRho, *spVViewA, *spVViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);

            INFO( "T component" );
            TestOp::add(*spSFlatA, T, idxFunc, rho, dLogRho, *spVFlatA, *spVFlatB, c);
            TestOp::add(*spSViewA, T, idxFunc, rho, dLogRho, *spVViewA, *spVViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);

            INFO( "P component" );
            TestOp::add(*spSFlatA, P, idxFunc, rho, dLogRho, *spVFlatA, *spVFlatB, c);
            TestOp::add(*spSViewA, P, idxFunc, rho, dLogRho, *spVViewA, *spVViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         }

         // Check sub operation
         {
            INFO( "Sub operation" );

            INFO( "R component" );
            TestOp::sub(*spSFlatA, R, idxFunc, rho, dLogRho, *spVFlatA, *spVFlatB, c);
            TestOp::sub(*spSViewA, R, idxFunc, rho, dLogRho, *spVViewA, *spVViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);

            INFO( "T component" );
            TestOp::sub(*spSFlatA, T, idxFunc, rho, dLogRho, *spVFlatA, *spVFlatB, c);
            TestOp::sub(*spSViewA, T, idxFunc, rho, dLogRho, *spVViewA, *spVViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);

            INFO( "P component" );
            TestOp::sub(*spSFlatA, P, idxFunc, rho, dLogRho, *spVFlatA, *spVFlatB, c);
            TestOp::sub(*spSViewA, P, idxFunc, rho, dLogRho, *spVViewA, *spVViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         }
      }
   }
}

TEST_CASE("SphericalSelfAdvectionAnelastic product timing", "[SphericalSelfAdvectionAnelasticTiming]")
{
   const std::string testName = "SphericalSelfAdvectionAnelasticTests";

   constexpr auto R = QuICC::FieldComponents::Physical::R;
   constexpr auto T = QuICC::FieldComponents::Physical::THETA;
   constexpr auto P = QuICC::FieldComponents::Physical::PHI;

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

      QuICC::Array rho(dim3D);
      for(std::uint32_t i = 0; i < dim3D; i++)
      {
         rho(i) = static_cast<double>(i+1)/static_cast<double>(dim3D + 1);
      }

      QuICC::Array dLogRho(dim3D);
      for(std::uint32_t i = 0; i < dim3D; i++)
      {
         dLogRho(i) = static_cast<double>(7*i+1)/static_cast<double>(dim3D + 1);
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
            auto spVectorB = details::createVectorField<double, sflat_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 2);

            for(std::uint32_t it = 0; it < itMax; it++)
            {
               QuICC::Profiler::RegionStart<1>(testName + "::Flat" + ctag + "Set_" + std::to_string(varBase));
               TestOp::set(*spScalarA, R, idxFunc, rho, dLogRho, *spVectorA, *spVectorB, c);
               TestOp::set(*spScalarA, T, idxFunc, rho, dLogRho, *spVectorA, *spVectorB, c);
               TestOp::set(*spScalarA, P, idxFunc, rho, dLogRho, *spVectorA, *spVectorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::Flat" + ctag + "Set_" + std::to_string(varBase));

               QuICC::Profiler::RegionStart<1>(testName + "::Flat" + ctag + "Add_" + std::to_string(varBase));
               TestOp::add(*spScalarA, R, idxFunc, rho, dLogRho, *spVectorA, *spVectorB, c);
               TestOp::add(*spScalarA, T, idxFunc, rho, dLogRho, *spVectorA, *spVectorB, c);
               TestOp::add(*spScalarA, P, idxFunc, rho, dLogRho, *spVectorA, *spVectorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::Flat" + ctag + "Add_" + std::to_string(varBase));

               QuICC::Profiler::RegionStart<1>(testName + "::Flat" + ctag + "Sub_" + std::to_string(varBase));
               TestOp::sub(*spScalarA, R, idxFunc, rho, dLogRho, *spVectorA, *spVectorB, c);
               TestOp::sub(*spScalarA, T, idxFunc, rho, dLogRho, *spVectorA, *spVectorB, c);
               TestOp::sub(*spScalarA, P, idxFunc, rho, dLogRho, *spVectorA, *spVectorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::Flat" + ctag + "Sub_" + std::to_string(varBase));
            }
         }

         // Timing ViewScalarField
         {
            auto spScalarA = details::createScalarField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 100);
            auto spVectorA = details::createVectorField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 1);
            auto spVectorB = details::createVectorField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 2);

            for(std::uint32_t it = 0; it < itMax; it++)
            {
               QuICC::Profiler::RegionStart<1>(testName + "::View" + ctag + "Set_" + std::to_string(varBase));
               TestOp::set(*spScalarA, R, idxFunc, rho, dLogRho, *spVectorA, *spVectorB, c);
               TestOp::set(*spScalarA, T, idxFunc, rho, dLogRho, *spVectorA, *spVectorB, c);
               TestOp::set(*spScalarA, P, idxFunc, rho, dLogRho, *spVectorA, *spVectorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::View" + ctag + "Set_" + std::to_string(varBase));

               QuICC::Profiler::RegionStart<1>(testName + "::View" + ctag + "Add_" + std::to_string(varBase));
               TestOp::add(*spScalarA, R, idxFunc, rho, dLogRho, *spVectorA, *spVectorB, c);
               TestOp::add(*spScalarA, T, idxFunc, rho, dLogRho, *spVectorA, *spVectorB, c);
               TestOp::add(*spScalarA, P, idxFunc, rho, dLogRho, *spVectorA, *spVectorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::View" + ctag + "Add_" + std::to_string(varBase));

               QuICC::Profiler::RegionStart<1>(testName + "::View" + ctag + "Sub_" + std::to_string(varBase));
               TestOp::sub(*spScalarA, R, idxFunc, rho, dLogRho, *spVectorA, *spVectorB, c);
               TestOp::sub(*spScalarA, T, idxFunc, rho, dLogRho, *spVectorA, *spVectorB, c);
               TestOp::sub(*spScalarA, P, idxFunc, rho, dLogRho, *spVectorA, *spVectorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::View" + ctag + "Sub_" + std::to_string(varBase));
            }
         }
      }
   }

   CHECK( true );
}

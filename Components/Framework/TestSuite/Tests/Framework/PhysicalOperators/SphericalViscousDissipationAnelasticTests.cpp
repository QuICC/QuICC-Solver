#include <algorithm>
#include <catch2/catch.hpp>
#include <cstdint>

//#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/ScalarFields/ScalarFieldSetup.hpp"
#include "details/Helper.hpp"
#include "QuICC/PhysicalOperators/SphericalViscousDissipationAnelastic.hpp"
#include "Profiler/Interface.hpp"
#include "Types/Math.hpp"

template <typename T>
   using sflat_t = QuICC::Datatypes::FlatScalarField<T>;
template <typename T>
   using sview_t = QuICC::Datatypes::ViewScalarField<T>;

using TestOp = QuICC::Physical::SphericalViscousDissipationAnelastic;

TEST_CASE("SphericalViscousDissipationAnelastic product validation FlatScalarField vs ViewScalarField", "[SphericalViscousDissipationAnelasticValidation]")
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

      auto spSFlatA = details::createScalarField<double, sflat_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 100 );
      auto spVFlatA = details::createVectorField<double, sflat_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 1);
      auto spTFlatB = details::createTensorField<double, sflat_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 2);

      auto spSViewA = details::createScalarField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 100);
      auto spVViewA = details::createVectorField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 1);
      auto spTViewB = details::createTensorField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 2);

      QuICC::Array rGrid(dim3D);
      for(std::uint32_t i = 0; i < dim3D; i++)
      {
         rGrid(i) = static_cast<double>(i+1)/static_cast<double>(dim3D + 1);
      }

      QuICC::Array nu(dim3D);
      for(std::uint32_t i = 0; i < dim3D; i++)
      {
         nu(i) = static_cast<double>(i+1)/static_cast<double>(dim3D + 1);
      }

      QuICC::Array temp(dim3D);
      for(std::uint32_t i = 0; i < dim3D; i++)
      {
         temp(i) = static_cast<double>(i+1)/static_cast<double>(dim3D + 1);
      }

      QuICC::Array rho(dim3D);
      for(std::uint32_t i = 0; i < dim3D; i++)
      {
         rho(i) = static_cast<double>(i+1)/static_cast<double>(dim3D + 1);
      }

      QuICC::Array dLogRho(dim3D);
      for(std::uint32_t i = 0; i < dim3D; i++)
      {
         dLogRho(i) = static_cast<double>(i+1)/static_cast<double>(dim3D + 1);
      }

      QuICC::Array tGrid(dim2D);
      for(std::uint32_t i = 0; i < dim2D; i++)
      {
         double theta = QuICC::Math::PI*static_cast<double>(i+1)/static_cast<double>(dim3D + 1);
         tGrid(i) = theta;
      }

      QuICC::Array pGrid(dim1D);
      for(std::uint32_t i = 0; i < dim1D; i++)
      {
         double phi = 2.0*QuICC::Math::PI*static_cast<double>(i+1)/static_cast<double>(dim3D + 1);
         pGrid(i) = phi;
      }

      details::IdxFunctor idxFunc(idx2D, idx3D);

      constexpr auto R = QuICC::FieldComponents::Physical::R;
      constexpr auto T = QuICC::FieldComponents::Physical::THETA;
      constexpr auto P = QuICC::FieldComponents::Physical::PHI;

      details::checkComputation(spSetup, *spSFlatA, *spSViewA);
      details::checkComputation(spSetup, spVFlatA->comp(R), spVViewA->comp(R));
      details::checkComputation(spSetup, spVFlatA->comp(T), spVViewA->comp(T));
      details::checkComputation(spSetup, spVFlatA->comp(P), spVViewA->comp(P));
      details::checkComputation(spSetup, spTFlatB->comp(R,R), spTViewB->comp(R,R));
      details::checkComputation(spSetup, spTFlatB->comp(T,T), spTViewB->comp(T,T));
      details::checkComputation(spSetup, spTFlatB->comp(P,P), spTViewB->comp(P,P));

      std::vector<double> cs = {1.0, 3.0};
      for(auto&& c: cs)
      {
         INFO( "c = " << c );

         // Check set operation
         {
            INFO( "Set operation" );

            TestOp::set(*spSFlatA, idxFunc, nu, temp, rho, dLogRho, *spVFlatA, *spTFlatB, c);
            TestOp::set(*spSViewA, idxFunc, nu, temp, rho, dLogRho, *spVViewA, *spTViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         }

         // Check add operation
         {
            INFO( "Add operation" );

            TestOp::add(*spSFlatA, idxFunc, nu, temp, rho, dLogRho, *spVFlatA, *spTFlatB, c);
            TestOp::add(*spSViewA, idxFunc, nu, temp, rho, dLogRho, *spVViewA, *spTViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         }

         // Check sub operation
         {
            INFO( "Sub operation" );

            TestOp::sub(*spSFlatA, idxFunc, nu, temp, rho, dLogRho, *spVFlatA, *spTFlatB, c);
            TestOp::sub(*spSViewA, idxFunc, nu, temp, rho, dLogRho, *spVViewA, *spTViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         }

#if 0
         // Check test operation
         {
            INFO( "Test operation" );

            TestOp::test(*spSFlatA, idxFunc, rGrid, nu, temp, rho, dLogRho, tGrid, pGrid, *spVFlatA, *spTFlatB, c);
            TestOp::test(*spSViewA, idxFunc, rGrid, nu, temp, rho, dLogRho, tGrid, pGrid, *spVViewA, *spTViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         }
#endif
      }
   }
}

TEST_CASE("SphericalViscousDissipationAnelastic product timing", "[SphericalViscousDissipationAnelasticTiming]")
{
   const std::string testName = "SphericalViscousDissipationAnelasticTests";

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

      QuICC::Array rGrid(dim3D);
      for(std::uint32_t i = 0; i < dim3D; i++)
      {
         rGrid(i) = static_cast<double>(i+1)/static_cast<double>(dim3D + 1);
      }

      QuICC::Array nu(dim3D);
      for(std::uint32_t i = 0; i < dim3D; i++)
      {
         nu(i) = static_cast<double>(i+1)/static_cast<double>(dim3D + 1);
      }

      QuICC::Array temp(dim3D);
      for(std::uint32_t i = 0; i < dim3D; i++)
      {
         temp(i) = static_cast<double>(i+1)/static_cast<double>(dim3D + 1);
      }

      QuICC::Array rho(dim3D);
      for(std::uint32_t i = 0; i < dim3D; i++)
      {
         rho(i) = static_cast<double>(i+1)/static_cast<double>(dim3D + 1);
      }

      QuICC::Array dLogRho(dim3D);
      for(std::uint32_t i = 0; i < dim3D; i++)
      {
         dLogRho(i) = static_cast<double>(i+1)/static_cast<double>(dim3D + 1);
      }

      QuICC::Array tGrid(dim2D);
      for(std::uint32_t i = 0; i < dim2D; i++)
      {
         double theta = QuICC::Math::PI*static_cast<double>(i+1)/static_cast<double>(dim3D + 1);
         tGrid(i) = theta;
      }

      QuICC::Array pGrid(dim1D);
      for(std::uint32_t i = 0; i < dim1D; i++)
      {
         double phi = 2.0*QuICC::Math::PI*static_cast<double>(i+1)/static_cast<double>(dim3D + 1);
         pGrid(i) = phi;
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
            auto spTensorB = details::createTensorField<double, sflat_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 1);

            for(std::uint32_t it = 0; it < itMax; it++)
            {
               QuICC::Profiler::RegionStart<1>(testName + "::Flat" + ctag + "Set_" + std::to_string(varBase));
               TestOp::set(*spScalarA, idxFunc, nu, temp, rho, dLogRho, *spVectorA, *spTensorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::Flat" + ctag + "Set_" + std::to_string(varBase));

               QuICC::Profiler::RegionStart<1>(testName + "::Flat" + ctag + "Add_" + std::to_string(varBase));
               TestOp::add(*spScalarA, idxFunc, nu, temp, rho, dLogRho, *spVectorA, *spTensorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::Flat" + ctag + "Add_" + std::to_string(varBase));

               QuICC::Profiler::RegionStart<1>(testName + "::Flat" + ctag + "Sub_" + std::to_string(varBase));
               TestOp::sub(*spScalarA, idxFunc, nu, temp, rho, dLogRho, *spVectorA, *spTensorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::Flat" + ctag + "Sub_" + std::to_string(varBase));

#if 0
               QuICC::Profiler::RegionStart<1>(testName + "::Flat" + ctag + "Sub_" + std::to_string(varBase));
               TestOp::test(*spScalarA, idxFunc, rGrid, nu, temp, rho, dLogRho, tGrid, pGrid, *spVectorA, *spTensorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::Flat" + ctag + "Sub_" + std::to_string(varBase));
#endif
            }
         }

         // Timing ViewScalarField
         {
            auto spScalarA = details::createScalarField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 100);
            auto spVectorA = details::createVectorField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 1);
            auto spTensorB = details::createTensorField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 1);

            for(std::uint32_t it = 0; it < itMax; it++)
            {
               QuICC::Profiler::RegionStart<1>(testName + "::View" + ctag + "Set_" + std::to_string(varBase));
               TestOp::set(*spScalarA, idxFunc, nu, temp, rho, dLogRho, *spVectorA, *spTensorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::View" + ctag + "Set_" + std::to_string(varBase));

               QuICC::Profiler::RegionStart<1>(testName + "::View" + ctag + "Add_" + std::to_string(varBase));
               TestOp::add(*spScalarA, idxFunc, nu, temp, rho, dLogRho, *spVectorA, *spTensorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::View" + ctag + "Add_" + std::to_string(varBase));

               QuICC::Profiler::RegionStart<1>(testName + "::View" + ctag + "Sub_" + std::to_string(varBase));
               TestOp::sub(*spScalarA, idxFunc, nu, temp, rho, dLogRho, *spVectorA, *spTensorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::View" + ctag + "Sub_" + std::to_string(varBase));

#if 0
               QuICC::Profiler::RegionStart<1>(testName + "::View" + ctag + "Sub_" + std::to_string(varBase));
               TestOp::test(*spScalarA, idxFunc, rGrid, nu, temp, rho, dLogRho, tGrid, pGrid, *spVectorA, *spTensorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::View" + ctag + "Sub_" + std::to_string(varBase));
#endif
            }
         }
      }
   }

   CHECK( true );
}

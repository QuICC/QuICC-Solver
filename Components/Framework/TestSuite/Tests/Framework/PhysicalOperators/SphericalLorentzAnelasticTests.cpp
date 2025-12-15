#include <algorithm>
#include <catch2/catch.hpp>
#include <cstdint>

//#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/ScalarFields/ScalarFieldSetup.hpp"
#include "details/Helper.hpp"
#include "QuICC/PhysicalOperators/SphericalLorentzAnelastic.hpp"
#include "Profiler/Interface.hpp"
#include "Types/Math.hpp"

template <typename T>
   using sflat_t = QuICC::Datatypes::FlatScalarField<T>;
template <typename T>
   using sview_t = QuICC::Datatypes::ViewScalarField<T>;

using TestOp = QuICC::Physical::SphericalLorentzAnelastic;

TEST_CASE("SphericalLorentzAnelastic product validation FlatScalarField vs ViewScalarField", "[SphericalLorentzAnelasticValidation]")
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
      auto spVFlatB = details::createVectorField<double, sflat_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 2);

      auto spSViewA = details::createScalarField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 100);
      auto spVViewA = details::createVectorField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 1);
      auto spVViewB = details::createVectorField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 2);

      QuICC::Array rGrid(dim3D);
      for(std::uint32_t i = 0; i < dim3D; i++)
      {
         rGrid(i) = static_cast<double>(i+1)/static_cast<double>(dim3D + 1);
      }

      QuICC::Array rho(dim3D);
      for(std::uint32_t i = 0; i < dim3D; i++)
      {
         rho(i) = static_cast<double>(i+1)/static_cast<double>(dim3D + 1);
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
      details::checkComputation(spSetup, spVFlatB->comp(R), spVViewB->comp(R));
      details::checkComputation(spSetup, spVFlatB->comp(T), spVViewB->comp(T));
      details::checkComputation(spSetup, spVFlatB->comp(P), spVViewB->comp(P));

      std::vector<double> cs = {1.0, 3.0};
      for(auto&& c: cs)
      {
         INFO( "c = " << c );

         // Check set operation
         {
            INFO( "Set operation" );

            INFO( "R component" );
            TestOp::set(*spSFlatA, R, idxFunc, rho, *spVFlatA, *spVFlatB, c);
            TestOp::set(*spSViewA, R, idxFunc, rho, *spVViewA, *spVViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);

            INFO( "T component" );
            TestOp::set(*spSFlatA, T, idxFunc, rho, *spVFlatA, *spVFlatB, c);
            TestOp::set(*spSViewA, T, idxFunc, rho, *spVViewA, *spVViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);

            INFO( "P component" );
            TestOp::set(*spSFlatA, P, idxFunc, rho, *spVFlatA, *spVFlatB, c);
            TestOp::set(*spSViewA, P, idxFunc, rho, *spVViewA, *spVViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         }

         // Check add operation
         {
            INFO( "Add operation" );

            INFO( "R component" );
            TestOp::add(*spSFlatA, R, idxFunc, rho, *spVFlatA, *spVFlatB, c);
            TestOp::add(*spSViewA, R, idxFunc, rho, *spVViewA, *spVViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);

            INFO( "T component" );
            TestOp::add(*spSFlatA, T, idxFunc, rho, *spVFlatA, *spVFlatB, c);
            TestOp::add(*spSViewA, T, idxFunc, rho, *spVViewA, *spVViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);

            INFO( "P component" );
            TestOp::add(*spSFlatA, P, idxFunc, rho, *spVFlatA, *spVFlatB, c);
            TestOp::add(*spSViewA, P, idxFunc, rho, *spVViewA, *spVViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         }

#if 0
         // Check sub operation
         {
            INFO( "Sub operation" );

            INFO( "R component" );
            TestOp::sub(*spSFlatA, R, idxFunc, rho, *spVFlatA, *spVFlatB, c);
            TestOp::sub(*spSViewA, R, idxFunc, rho, *spVViewA, *spVViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);

            INFO( "T component" );
            TestOp::sub(*spSFlatA, T, idxFunc, rho, *spVFlatA, *spVFlatB, c);
            TestOp::sub(*spSViewA, T, idxFunc, rho, *spVViewA, *spVViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);

            INFO( "P component" );
            TestOp::sub(*spSFlatA, P, idxFunc, rho, *spVFlatA, *spVFlatB, c);
            TestOp::sub(*spSViewA, P, idxFunc, rho, *spVViewA, *spVViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         }
#endif

#if 0
         // Check test operation
         {
            INFO( "Test operation" );

            INFO( "R component" );
            TestOp::test(*spSFlatA, R, idxFunc, rGrid, rho, tGrid, pGrid, *spVFlatA, *spVFlatB, c);
            TestOp::test(*spSViewA, R, idxFunc, rGrid, rho, tGrid, pGrid, *spVViewA, *spVViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);

            INFO( "T component" );
            TestOp::test(*spSFlatA, T, idxFunc, rGrid, rho, tGrid, pGrid, *spVFlatA, *spVFlatB, c);
            TestOp::test(*spSViewA, T, idxFunc, rGrid, rho, tGrid, pGrid, *spVViewA, *spVViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);

            INFO( "P component" );
            TestOp::test(*spSFlatA, P, idxFunc, rGrid, rho, tGrid, pGrid, *spVFlatA, *spVFlatB, c);
            TestOp::test(*spSViewA, P, idxFunc, rGrid, rho, tGrid, pGrid, *spVViewA, *spVViewB, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         }
#endif
      }
   }
}

TEST_CASE("SphericalLorentzAnelastic product timing", "[SphericalLorentzAnelasticTiming]")
{
   const std::string testName = "SphericalLorentzAnelasticTests";

   std::uint32_t itMax = 100;
   std::vector<double> cs = {1.0, 3.0};
   std::vector<std::uint32_t> variants = details::performanceVariants();

   constexpr auto R = QuICC::FieldComponents::Physical::R;
   constexpr auto T = QuICC::FieldComponents::Physical::THETA;
   constexpr auto P = QuICC::FieldComponents::Physical::PHI;

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

      QuICC::Array rho(dim3D);
      for(std::uint32_t i = 0; i < dim3D; i++)
      {
         rho(i) = static_cast<double>(i+1)/static_cast<double>(dim3D + 1);
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
            auto spVectorB = details::createVectorField<double, sflat_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 2);

            for(std::uint32_t it = 0; it < itMax; it++)
            {
               QuICC::Profiler::RegionStart<1>(testName + "::Flat" + ctag + "Set_" + std::to_string(varBase));
               TestOp::set(*spScalarA, R, idxFunc, rho, *spVectorA, *spVectorB, c);
               TestOp::set(*spScalarA, T, idxFunc, rho, *spVectorA, *spVectorB, c);
               TestOp::set(*spScalarA, P, idxFunc, rho, *spVectorA, *spVectorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::Flat" + ctag + "Set_" + std::to_string(varBase));

               QuICC::Profiler::RegionStart<1>(testName + "::Flat" + ctag + "Add_" + std::to_string(varBase));
               TestOp::add(*spScalarA, R, idxFunc, rho, *spVectorA, *spVectorB, c);
               TestOp::add(*spScalarA, T, idxFunc, rho, *spVectorA, *spVectorB, c);
               TestOp::add(*spScalarA, P, idxFunc, rho, *spVectorA, *spVectorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::Flat" + ctag + "Add_" + std::to_string(varBase));

#if 0
               QuICC::Profiler::RegionStart<1>(testName + "::Flat" + ctag + "Sub_" + std::to_string(varBase));
               TestOp::sub(*spScalarA, R, idxFunc, rho, *spVectorA, *spVectorB, c);
               TestOp::sub(*spScalarA, T, idxFunc, rho, *spVectorA, *spVectorB, c);
               TestOp::sub(*spScalarA, P, idxFunc, rho, *spVectorA, *spVectorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::Flat" + ctag + "Sub_" + std::to_string(varBase));
#endif

#if 0
               QuICC::Profiler::RegionStart<1>(testName + "::Flat" + ctag + "Sub_" + std::to_string(varBase));
               TestOp::test(*spScalarA, R, idxFunc, rGrid, rho, tGrid, pGrid, *spVectorA, *spVectorB, c);
               TestOp::test(*spScalarA, T, idxFunc, rGrid, rho, tGrid, pGrid, *spVectorA, *spVectorB, c);
               TestOp::test(*spScalarA, P, idxFunc, rGrid, rho, tGrid, pGrid, *spVectorA, *spVectorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::Flat" + ctag + "Sub_" + std::to_string(varBase));
#endif
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
               TestOp::set(*spScalarA, R, idxFunc, rho, *spVectorA, *spVectorB, c);
               TestOp::set(*spScalarA, T, idxFunc, rho, *spVectorA, *spVectorB, c);
               TestOp::set(*spScalarA, P, idxFunc, rho, *spVectorA, *spVectorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::View" + ctag + "Set_" + std::to_string(varBase));

               QuICC::Profiler::RegionStart<1>(testName + "::View" + ctag + "Add_" + std::to_string(varBase));
               TestOp::add(*spScalarA, R, idxFunc, rho, *spVectorA, *spVectorB, c);
               TestOp::add(*spScalarA, T, idxFunc, rho, *spVectorA, *spVectorB, c);
               TestOp::add(*spScalarA, P, idxFunc, rho, *spVectorA, *spVectorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::View" + ctag + "Add_" + std::to_string(varBase));

#if 0
               QuICC::Profiler::RegionStart<1>(testName + "::View" + ctag + "Sub_" + std::to_string(varBase));
               TestOp::sub(*spScalarA, R, idxFunc, rho, *spVectorA, *spVectorB, c);
               TestOp::sub(*spScalarA, T, idxFunc, rho, *spVectorA, *spVectorB, c);
               TestOp::sub(*spScalarA, P, idxFunc, rho, *spVectorA, *spVectorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::View" + ctag + "Sub_" + std::to_string(varBase));
#endif

#if 0
               QuICC::Profiler::RegionStart<1>(testName + "::View" + ctag + "Sub_" + std::to_string(varBase));
               TestOp::test(*spScalarA, R, idxFunc, rGrid, rho, tGrid, pGrid, *spVectorA, *spVectorB, c);
               TestOp::test(*spScalarA, T, idxFunc, rGrid, rho, tGrid, pGrid, *spVectorA, *spVectorB, c);
               TestOp::test(*spScalarA, P, idxFunc, rGrid, rho, tGrid, pGrid, *spVectorA, *spVectorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::View" + ctag + "Sub_" + std::to_string(varBase));
#endif
            }
         }
      }
   }

   CHECK( true );
}

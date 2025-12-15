#include <algorithm>
#include <catch2/catch.hpp>
#include <cstdint>

//#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/ScalarFields/ScalarFieldSetup.hpp"
#include "details/Helper.hpp"
#include "QuICC/PhysicalOperators/SphericalPoincare.hpp"
#include "Profiler/Interface.hpp"
#include "Types/Math.hpp"

template <typename T>
   using sflat_t = QuICC::Datatypes::FlatScalarField<T>;
template <typename T>
   using sview_t = QuICC::Datatypes::ViewScalarField<T>;

using TestOp = QuICC::Physical::SphericalPoincare;

TEST_CASE("SphericalPoincare product validation FlatScalarField vs ViewScalarField", "[SphericalPoincareValidation]")
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

      auto spSViewA = details::createScalarField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 100);

      QuICC::Array rGrid(dim3D);
      for(std::uint32_t i = 0; i < dim3D; i++)
      {
         rGrid(i) = static_cast<double>(i+1)/static_cast<double>(dim3D + 1);
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

      const double t = 0.13;

      details::IdxFunctor idxFunc(idx2D, idx3D);

      constexpr auto R = QuICC::FieldComponents::Physical::R;
      constexpr auto T = QuICC::FieldComponents::Physical::THETA;
      constexpr auto P = QuICC::FieldComponents::Physical::PHI;

      details::checkComputation(spSetup, *spSFlatA, *spSViewA);

      std::vector<double> cs = {1.0, 3.0};
      for(auto&& c: cs)
      {
         INFO( "c = " << c );

         // Check set operation
         {
            INFO( "Set operation" );

            INFO( "R component" );
            TestOp::set(*spSFlatA, R, idxFunc, rGrid, tGrid, pGrid, t, c);
            TestOp::set(*spSViewA, R, idxFunc, rGrid, tGrid, pGrid, t, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);

            INFO( "T component" );
            TestOp::set(*spSFlatA, T, idxFunc, rGrid, tGrid, pGrid, t, c);
            TestOp::set(*spSViewA, T, idxFunc, rGrid, tGrid, pGrid, t, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);

            INFO( "P component" );
            TestOp::set(*spSFlatA, P, idxFunc, rGrid, tGrid, pGrid, t, c);
            TestOp::set(*spSViewA, P, idxFunc, rGrid, tGrid, pGrid, t, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         }

         // Check add operation
         {
            INFO( "Add operation" );

            INFO( "R component" );
            TestOp::add(*spSFlatA, R, idxFunc, rGrid, tGrid, pGrid, t, c);
            TestOp::add(*spSViewA, R, idxFunc, rGrid, tGrid, pGrid, t, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);

            INFO( "T component" );
            TestOp::add(*spSFlatA, T, idxFunc, rGrid, tGrid, pGrid, t, c);
            TestOp::add(*spSViewA, T, idxFunc, rGrid, tGrid, pGrid, t, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);

            INFO( "P component" );
            TestOp::add(*spSFlatA, P, idxFunc, rGrid, tGrid, pGrid, t, c);
            TestOp::add(*spSViewA, P, idxFunc, rGrid, tGrid, pGrid, t, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         }

         // Check sub operation
         {
            INFO( "Sub operation" );

            INFO( "R component" );
            TestOp::sub(*spSFlatA, R, idxFunc, rGrid, tGrid, pGrid, t, c);
            TestOp::sub(*spSViewA, R, idxFunc, rGrid, tGrid, pGrid, t, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);

            INFO( "T component" );
            TestOp::sub(*spSFlatA, T, idxFunc, rGrid, tGrid, pGrid, t, c);
            TestOp::sub(*spSViewA, T, idxFunc, rGrid, tGrid, pGrid, t, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);

            INFO( "P component" );
            TestOp::sub(*spSFlatA, P, idxFunc, rGrid, tGrid, pGrid, t, c);
            TestOp::sub(*spSViewA, P, idxFunc, rGrid, tGrid, pGrid, t, c);
            details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         }
      }
   }
}

TEST_CASE("SphericalPoincare product timing", "[SphericalPoincareTiming]")
{
   const std::string testName = "SphericalPoincareTests";

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

      details::IdxFunctor idxFunc(idx2D, idx3D);

      QuICC::Array rGrid(dim3D);
      for(std::uint32_t i = 0; i < dim3D; i++)
      {
         rGrid(i) = static_cast<double>(i+1)/static_cast<double>(dim3D + 1);
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

      const double t = 0.13;

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

            for(std::uint32_t it = 0; it < itMax; it++)
            {
               QuICC::Profiler::RegionStart<1>(testName + "::Flat" + ctag + "Set_" + std::to_string(varBase));
               TestOp::set(*spScalarA, R, idxFunc, rGrid, tGrid, pGrid, t, c);
               TestOp::set(*spScalarA, T, idxFunc, rGrid, tGrid, pGrid, t, c);
               TestOp::set(*spScalarA, P, idxFunc, rGrid, tGrid, pGrid, t, c);
               QuICC::Profiler::RegionStop<1>(testName + "::Flat" + ctag + "Set_" + std::to_string(varBase));

               QuICC::Profiler::RegionStart<1>(testName + "::Flat" + ctag + "Add_" + std::to_string(varBase));
               TestOp::add(*spScalarA, R, idxFunc, rGrid, tGrid, pGrid, t, c);
               TestOp::add(*spScalarA, T, idxFunc, rGrid, tGrid, pGrid, t, c);
               TestOp::add(*spScalarA, P, idxFunc, rGrid, tGrid, pGrid, t, c);
               QuICC::Profiler::RegionStop<1>(testName + "::Flat" + ctag + "Add_" + std::to_string(varBase));

               QuICC::Profiler::RegionStart<1>(testName + "::Flat" + ctag + "Sub_" + std::to_string(varBase));
               TestOp::sub(*spScalarA, R, idxFunc, rGrid, tGrid, pGrid, t, c);
               TestOp::sub(*spScalarA, T, idxFunc, rGrid, tGrid, pGrid, t, c);
               TestOp::sub(*spScalarA, P, idxFunc, rGrid, tGrid, pGrid, t, c);
               QuICC::Profiler::RegionStop<1>(testName + "::Flat" + ctag + "Sub_" + std::to_string(varBase));
            }
         }

         // Timing ViewScalarField
         {
            auto spScalarA = details::createScalarField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 100);

            for(std::uint32_t it = 0; it < itMax; it++)
            {
               QuICC::Profiler::RegionStart<1>(testName + "::View" + ctag + "Set_" + std::to_string(varBase));
               TestOp::set(*spScalarA, R, idxFunc, rGrid, tGrid, pGrid, t, c);
               TestOp::set(*spScalarA, T, idxFunc, rGrid, tGrid, pGrid, t, c);
               TestOp::set(*spScalarA, P, idxFunc, rGrid, tGrid, pGrid, t, c);
               QuICC::Profiler::RegionStop<1>(testName + "::View" + ctag + "Set_" + std::to_string(varBase));

               QuICC::Profiler::RegionStart<1>(testName + "::View" + ctag + "Add_" + std::to_string(varBase));
               TestOp::add(*spScalarA, R, idxFunc, rGrid, tGrid, pGrid, t, c);
               TestOp::add(*spScalarA, T, idxFunc, rGrid, tGrid, pGrid, t, c);
               TestOp::add(*spScalarA, P, idxFunc, rGrid, tGrid, pGrid, t, c);
               QuICC::Profiler::RegionStop<1>(testName + "::View" + ctag + "Add_" + std::to_string(varBase));

               QuICC::Profiler::RegionStart<1>(testName + "::View" + ctag + "Sub_" + std::to_string(varBase));
               TestOp::sub(*spScalarA, R, idxFunc, rGrid, tGrid, pGrid, t, c);
               TestOp::sub(*spScalarA, T, idxFunc, rGrid, tGrid, pGrid, t, c);
               TestOp::sub(*spScalarA, P, idxFunc, rGrid, tGrid, pGrid, t, c);
               QuICC::Profiler::RegionStop<1>(testName + "::View" + ctag + "Sub_" + std::to_string(varBase));
            }
         }
      }
   }

   CHECK( true );
}

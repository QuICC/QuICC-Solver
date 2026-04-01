#include <algorithm>
#include <catch2/catch.hpp>
#include <cstdint>

//#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/ScalarFields/ScalarFieldSetup.hpp"
#include "details/Helper.hpp"
#include "QuICC/PhysicalOperators/Cross.hpp"
#include "Profiler/Interface.hpp"

template <typename T>
   using sflat_t = QuICC::Datatypes::FlatScalarField<T>;
template <typename T>
   using sview_t = QuICC::Datatypes::ViewScalarField<T>;

template <QuICC::FieldComponents::Physical::Id TFIRST, QuICC::FieldComponents::Physical::Id TSECOND>
   using TestOp = QuICC::Physical::Cross<TFIRST,TSECOND>;

TEST_CASE("Cross product validation FlatScalarField vs ViewScalarField", "[CrossValidation]")
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
      auto spVFlatA = details::createVectorField<double, sflat_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 0);
      auto spVFlatB = details::createVectorField<double, sflat_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 1);

      auto spSViewA = details::createScalarField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 100);
      auto spVViewA = details::createVectorField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 0);
      auto spVViewB = details::createVectorField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 1);

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
         TestOp<T,P>::set(*spSFlatA, *spVFlatA, *spVFlatB, c);
         TestOp<T,P>::set(*spSViewA, *spVViewA, *spVViewB, c);
         details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         TestOp<P,R>::set(*spSFlatA, *spVFlatA, *spVFlatB, c);
         TestOp<P,R>::set(*spSViewA, *spVViewA, *spVViewB, c);
         details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         TestOp<R,T>::set(*spSFlatA, *spVFlatA, *spVFlatB, c);
         TestOp<R,T>::set(*spSViewA, *spVViewA, *spVViewB, c);
         details::checkComputation(spSetup, *spSFlatA, *spSViewA);

         // Check add operation
         TestOp<T,P>::add(*spSFlatA, *spVFlatA, *spVFlatB, c);
         TestOp<T,P>::add(*spSViewA, *spVViewA, *spVViewB, c);
         details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         TestOp<P,R>::add(*spSFlatA, *spVFlatA, *spVFlatB, c);
         TestOp<P,R>::add(*spSViewA, *spVViewA, *spVViewB, c);
         details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         TestOp<R,T>::add(*spSFlatA, *spVFlatA, *spVFlatB, c);
         TestOp<R,T>::add(*spSViewA, *spVViewA, *spVViewB, c);
         details::checkComputation(spSetup, *spSFlatA, *spSViewA);

         // Check sub operation
         TestOp<T,P>::sub(*spSFlatA, *spVFlatA, *spVFlatB, c);
         TestOp<T,P>::sub(*spSViewA, *spVViewA, *spVViewB, c);
         details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         TestOp<P,R>::sub(*spSFlatA, *spVFlatA, *spVFlatB, c);
         TestOp<P,R>::sub(*spSViewA, *spVViewA, *spVViewB, c);
         details::checkComputation(spSetup, *spSFlatA, *spSViewA);
         TestOp<R,T>::sub(*spSFlatA, *spVFlatA, *spVFlatB, c);
         TestOp<R,T>::sub(*spSViewA, *spVViewA, *spVViewB, c);
         details::checkComputation(spSetup, *spSFlatA, *spSViewA);
      }
   }
}

TEST_CASE("Cross product timing", "[CrossTiming]")
{
   const std::string testName = "CrossTests";

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
            auto spVectorA = details::createVectorField<double, sflat_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 0);
            auto spVectorB = details::createVectorField<double, sflat_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 1);

            for(std::uint32_t it = 0; it < itMax; it++)
            {
               QuICC::Profiler::RegionStart<1>(testName + "::Flat" + ctag + "Set_" + std::to_string(varBase));
               TestOp<T,P>::set(*spScalarA, *spVectorA, *spVectorB, c);
               TestOp<P,R>::set(*spScalarA, *spVectorA, *spVectorB, c);
               TestOp<R,T>::set(*spScalarA, *spVectorA, *spVectorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::Flat" + ctag + "Set_" + std::to_string(varBase));

               QuICC::Profiler::RegionStart<1>(testName + "::Flat" + ctag + "Add_" + std::to_string(varBase));
               TestOp<T,P>::add(*spScalarA, *spVectorA, *spVectorB, c);
               TestOp<P,R>::add(*spScalarA, *spVectorA, *spVectorB, c);
               TestOp<R,T>::add(*spScalarA, *spVectorA, *spVectorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::Flat" + ctag + "Add_" + std::to_string(varBase));

               QuICC::Profiler::RegionStart<1>(testName + "::Flat" + ctag + "Sub_" + std::to_string(varBase));
               TestOp<T,P>::sub(*spScalarA, *spVectorA, *spVectorB, c);
               TestOp<P,R>::sub(*spScalarA, *spVectorA, *spVectorB, c);
               TestOp<R,T>::sub(*spScalarA, *spVectorA, *spVectorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::Flat" + ctag + "Sub_" + std::to_string(varBase));
            }
         }

         // Timing ViewScalarField
         {
            auto spScalarA = details::createScalarField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 100);
            auto spVectorA = details::createVectorField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 0);
            auto spVectorB = details::createVectorField<double, sview_t>(dim1D, dim3D, idx2D, idx3D, spSetup, varBase + 1);

            for(std::uint32_t it = 0; it < itMax; it++)
            {
               QuICC::Profiler::RegionStart<1>(testName + "::View" + ctag + "Set_" + std::to_string(varBase));
               TestOp<T,P>::set(*spScalarA, *spVectorA, *spVectorB, c);
               TestOp<P,R>::set(*spScalarA, *spVectorA, *spVectorB, c);
               TestOp<R,T>::set(*spScalarA, *spVectorA, *spVectorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::View" + ctag + "Set_" + std::to_string(varBase));

               QuICC::Profiler::RegionStart<1>(testName + "::View" + ctag + "Add_" + std::to_string(varBase));
               TestOp<T,P>::add(*spScalarA, *spVectorA, *spVectorB, c);
               TestOp<P,R>::add(*spScalarA, *spVectorA, *spVectorB, c);
               TestOp<R,T>::add(*spScalarA, *spVectorA, *spVectorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::View" + ctag + "Add_" + std::to_string(varBase));

               QuICC::Profiler::RegionStart<1>(testName + "::View" + ctag + "Sub_" + std::to_string(varBase));
               TestOp<T,P>::sub(*spScalarA, *spVectorA, *spVectorB, c);
               TestOp<P,R>::sub(*spScalarA, *spVectorA, *spVectorB, c);
               TestOp<R,T>::sub(*spScalarA, *spVectorA, *spVectorB, c);
               QuICC::Profiler::RegionStop<1>(testName + "::View" + ctag + "Sub_" + std::to_string(varBase));
            }
         }
      }
   }

   CHECK( true );
}

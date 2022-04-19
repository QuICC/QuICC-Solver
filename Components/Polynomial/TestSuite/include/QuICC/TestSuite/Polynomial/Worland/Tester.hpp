/**
 * @file Tester.hpp
 * @brief Tester for Worland polynomials
 */

#ifndef QUICC_TESTSUITE_POLYNOMIAL_WORLAND_TESTER_HPP
#define QUICC_TESTSUITE_POLYNOMIAL_WORLAND_TESTER_HPP


// Configuration includes
//
#include <catch2/catch.hpp>

// System includes
//
#include <string>
#include <set>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Precision.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandRule.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Reduce.hpp"
#include "QuICC/TestSuite/Polynomial/TesterBase.hpp"

namespace poly = ::QuICC::Polynomial;
namespace current = ::QuICC::Polynomial::Worland;

namespace QuICC {

namespace TestSuite {

namespace Polynomial {

namespace Worland {

   template <typename TOp> class Tester: public Polynomial::TesterBase<TOp>
   {
      public:
         /// Typedef for parameter type
         typedef typename Polynomial::TesterBase<TOp>::ParameterType ParameterType;

         /**
          * @brief Constructor
          */
         Tester(const std::string& fname, const bool keepData);

         /*
          * @brief Destructor
          */
         virtual ~Tester() = default;

      protected:
         /**
          * @brief Build filename extension with resolution information
          */
         virtual std::string resname(const int specN, const int physN, const ParameterType& param) const override;

         /**
          * @brief Test matrix operator
          */
         virtual Matrix buildOperator(const int specN, const int physN, const ParameterType& param, const TestType type) const override;

         /**
          * @brief Format the parameters
          */
         virtual std::string formatParameter(const ParameterType& param) const override;

      private:
         /**
          * @brief Append specific path
          */
         void appendPath();

         /**
          * @brief Test quadrature matrix operator
          */
         Matrix buildMatrix(const int specN, const int physN, const ParameterType& param, const bool isWeighted) const;

         /**
          * @brief Test on-the-fly inner product
          */
         Matrix buildInner(const int specN, const int physN, const ParameterType& param) const;

         /**
          * @brief Test on-the-fly outer product
          */
         Matrix buildOuter(const int specN, const int physN, const ParameterType& param) const;

         /**
          * @brief Test on-the-fly reduction
          */
         Matrix buildReduce(const int specN, const int physN, const ParameterType& param) const;
   };

   template <typename TOp> Tester<TOp>::Tester(const std::string& fname, const bool keepData)
      : Polynomial::TesterBase<TOp>(fname, keepData)
   {
      this->appendPath();
   }

   template <typename TOp> void Tester<TOp>::appendPath()
   {
      #if defined QUICC_WORLAND_TYPE_CHEBYSHEV
         this->mPath += "Worland/Chebyshev/";
      #elif defined QUICC_WORLAND_TYPE_LEGENDRE
         this->mPath += "Worland/Legendre/";
      #elif defined QUICC_WORLAND_TYPE_CYLENERGY
         this->mPath += "Worland/CylEnergy/";
      #elif defined QUICC_WORLAND_TYPE_SPHENERGY
         this->mPath += "Worland/SphEnergy/";
      #else
         this->mPath += "Worland/Unknown/";
      #endif
   }

   template <typename TOp> Matrix Tester<TOp>::buildOperator(const int specN, const int physN, const ParameterType& param, const TestType type) const
   {
      Matrix outData;
      switch(type)
      {
         case TestType::MATRIX:
            outData = this->buildMatrix(specN, physN, param, false);
            break;
         case TestType::WEIGHTED_MATRIX:
            outData = this->buildMatrix(specN, physN, param, true);
            break;
         case TestType::OTF_INNER:
            outData = this->buildInner(specN, physN, param);
            break;
         case TestType::OTF_OUTER:
            outData = this->buildOuter(specN, physN, param);
            break;
         case TestType::OTF_REDUCE:
            outData = this->buildReduce(specN, physN, param);
            break;
         default:
            throw std::logic_error("Test type not implemented");
            break;
      }

      return outData;
   }

   template <typename TOp> Matrix Tester<TOp>::buildMatrix(const int specN, const int physN, const ParameterType& param, const bool isWeighted) const
   {
      // Get harmonic degree from parameters
      int l = static_cast<int>(param.at(0));

      // Create quadrature
      internal::Array igrid;
      internal::Array iweights;
      poly::Quadrature::WorlandRule quad;
      quad.computeQuadrature(igrid, iweights, physN);

      TOp op;
      Matrix outData(physN, specN+1);

      if(isWeighted)
      {
         op.template compute<MHDFloat>(outData, specN+1, l, igrid, iweights, current::Evaluator::Set());
      } else
      {
         op.template compute<MHDFloat>(outData, specN+1, l, igrid, internal::Array(), current::Evaluator::Set());
      }

      return outData;
   }

   template <typename TOp> Matrix Tester<TOp>::buildOuter(const int specN, const int physN, const ParameterType& param) const
   {
      // Get harmonic degree from parameters
      int l = static_cast<int>(param.at(0));

      // Create quadrature
      internal::Array igrid;
      internal::Array iweights;
      poly::Quadrature::WorlandRule quad;
      quad.computeQuadrature(igrid, iweights, physN);

      TOp op;
      Matrix outData(physN, specN+1);
      Matrix inData = Matrix::Identity(specN+1, specN+1);

      op.template compute<MHDFloat>(outData, specN+1, l, igrid, internal::Array(), current::Evaluator::OuterProduct<MHDFloat>(inData));

      return outData;
   }

   template <typename TOp> Matrix Tester<TOp>::buildInner(const int specN, const int physN, const ParameterType& param) const
   {
      // Get harmonic degree from parameters
      int l = static_cast<int>(param.at(0));

      // Create quadrature
      internal::Array igrid;
      internal::Array iweights;
      poly::Quadrature::WorlandRule quad;
      quad.computeQuadrature(igrid, iweights, physN);

      TOp op;
      Matrix outData(specN+1, specN+1);
      Matrix inData = Matrix::Identity(physN, specN+1);
      op.template compute<MHDFloat>(inData, specN+1, l, igrid, iweights, current::Evaluator::Set());

      op.template compute<MHDFloat>(outData, specN+1, l, igrid, internal::Array(), current::Evaluator::InnerProduct<MHDFloat>(inData));

      Matrix outDataT = outData.transpose();

      return outDataT;
   }

   template <typename TOp> Matrix Tester<TOp>::buildReduce(const int specN, const int physN, const ParameterType& param) const
   {
      // Get harmonic degree from parameters
      int l = static_cast<int>(param.at(0));

      // Create quadrature
      internal::Array igrid;
      internal::Array iweights;
      poly::Quadrature::WorlandRule quad;
      quad.computeQuadrature(igrid, iweights, physN);

      TOp op;
      Matrix outData(specN+1, 1);

      op.template compute<MHDFloat>(outData, specN+1, l, igrid, internal::Array(), current::Evaluator::Reduce());

      return outData;
   }

   template <typename TOp> std::string Tester<TOp>::resname(const int specN, const int physN, const ParameterType& param) const
   {
      // Get harmonic degree from parameters
      int l = static_cast<int>(param.at(0));

      std::string s = "_l" + std::to_string(l) + "_n" + std::to_string(specN) + "_g" + std::to_string(physN);
      return s;
   }

   template <typename TOp> std::string Tester<TOp>::formatParameter(const ParameterType& param) const
   {
      // Get harmonic degree from parameters
      int l = static_cast<int>(param.at(0));

      std::stringstream ss;
      ss << "l: " << l;

      return ss.str();
   }

}
}
}
}

#endif //QUICC_TESTSUITE_POLYNOMIAL_WORLAND_TESTER_HPP

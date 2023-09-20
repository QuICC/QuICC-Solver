/**
 * @file Tester.hpp
 * @brief Tester for Legendre polynomials
 */

#ifndef QUICC_TESTSUITE_POLYNOMIAL_LEGENDRE_TESTER_HPP
#define QUICC_TESTSUITE_POLYNOMIAL_LEGENDRE_TESTER_HPP


// Configuration includes
//

// System includes
//
#include <catch2/catch.hpp>
#include <string>
#include <set>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Precision.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"
#include "QuICC/Polynomial/Legendre/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Legendre/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Legendre/Evaluator/OuterProduct.hpp"
#include "QuICC/Polynomial/Legendre/Evaluator/Reduce.hpp"
#include "QuICC/TestSuite/Polynomial/TesterBase.hpp"

namespace poly = ::QuICC::Polynomial;
namespace current = ::QuICC::Polynomial::Legendre;

namespace QuICC {

namespace TestSuite {

namespace Polynomial {

namespace Legendre {

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
      this->mPath += "Legendre/";
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
      // Create quadrature
      internal::Array igrid;
      internal::Array iweights;
      poly::Quadrature::LegendreRule quad;
      quad.computeQuadrature(igrid, iweights, physN);

      TOp op;
      Matrix outData(physN, specN);

      if(isWeighted)
      {
         op.template compute<MHDFloat>(outData, specN, igrid, iweights, current::Evaluator::Set());
      }
      else
      {
         op.template compute<MHDFloat>(outData, specN, igrid, internal::Array(), current::Evaluator::Set());
      }

      return outData;
   }

   template <typename TOp> Matrix Tester<TOp>::buildOuter(const int specN, const int physN, const ParameterType& param) const
   {
      // Create quadrature
      internal::Array igrid;
      internal::Array iweights;
      poly::Quadrature::LegendreRule quad;
      quad.computeQuadrature(igrid, iweights, physN);

      TOp op;
      Matrix outData(physN, specN);
      Matrix inData = Matrix::Identity(specN, specN);
      op.template compute<MHDFloat>(outData, specN, igrid, internal::Array(), current::Evaluator::OuterProduct<MHDFloat>(inData));

      return outData;
   }

   template <typename TOp> Matrix Tester<TOp>::buildInner(const int specN, const int physN, const ParameterType& param) const
   {
      // Create quadrature
      internal::Array igrid;
      internal::Array iweights;
      poly::Quadrature::LegendreRule quad;
      quad.computeQuadrature(igrid, iweights, physN);

      TOp op;
      Matrix outData(specN, specN);
      Matrix inData(physN, specN);
      op.template compute<MHDFloat>(inData, specN, igrid, internal::Array(), current::Evaluator::Set());
      op.template compute<MHDFloat>(outData, specN, igrid, iweights, current::Evaluator::InnerProduct<MHDFloat>(inData));

      Matrix outDataT = outData.transpose();

      return outDataT;
   }

   template <typename TOp> Matrix Tester<TOp>::buildReduce(const int specN, const int physN, const ParameterType& param) const
   {
      // Create quadrature
      internal::Array igrid;
      internal::Array iweights;
      poly::Quadrature::LegendreRule quad;
      quad.computeQuadrature(igrid, iweights, physN);

      TOp op;
      Matrix outData(specN, 1);
      op.template compute<MHDFloat>(outData, specN, igrid, internal::Array(), current::Evaluator::Reduce());

      return outData;
   }

   template <typename TOp> std::string Tester<TOp>::resname(const int specN, const int physN, const ParameterType& param) const
   {
      std::string s =  "_n" + std::to_string(specN) + "_g" + std::to_string(physN);

      return s;
   }

   template <typename TOp> std::string Tester<TOp>::formatParameter(const ParameterType& param) const
   {
      return "";
   }

}
}
}
}

#endif //QUICC_TESTSUITE_POLYNOMIAL_LEGENDRE_TESTER_HPP

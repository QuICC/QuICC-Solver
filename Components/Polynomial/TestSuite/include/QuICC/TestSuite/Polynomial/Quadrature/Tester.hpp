/**
 * @file Tester.hpp
 * @brief Tester for quadratures
 */

#ifndef QUICC_TESTSUITE_POLYNOMIAL_QUADRATURE_TESTER_HPP
#define QUICC_TESTSUITE_POLYNOMIAL_QUADRATURE_TESTER_HPP


// Configuration includes
//

// System includes
//
#include <catch2/catch.hpp>
#include <string>
#include <vector>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Internal/Casts.hpp"
#include "QuICC/TestSuite/Polynomial/TesterBase.hpp"

namespace QuICC {

namespace TestSuite {

namespace Polynomial {

namespace Quadrature {

   template <typename TOp, unsigned int N = 0> class Tester: public Polynomial::TesterBase<TOp>
   {
      public:
         /// Typedef for parameter type
         typedef typename Polynomial::TesterBase<TOp>::ParameterType ParameterType;
         typedef typename Polynomial::TesterBase<TOp>::ErrorType ErrorType;

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
          * @brief Build quadrature
          */
         virtual Matrix buildOperator(const int specN, const int physN, const ParameterType& param, const TestType type) const override;

         /**
          * @brief Check if reference is special value
          */
         virtual bool isSpecial(const int i, const int j, const MHDFloat data, const MHDFloat ref) const override;

         /**
          * @brief Check special values
          */
         ErrorType checkSpecial(const int i, const int j, const MHDFloat data, const Matrix& ref) const override;

         /**
          * @brief Format the parameters
          */
         virtual std::string formatParameter(const ParameterType& param) const override;

         template <typename Q = TOp> typename std::enable_if< N == 0, Q >::type createQuad(const ParameterType& param) const
         {
            TOp quad;

            return quad;
         };

         template <typename Q = TOp> typename std::enable_if< N == 2, Q >::type createQuad(const ParameterType& param) const
         {
            TOp quad(param.at(0), param.at(1));

            return quad;
         };

      private:
         /**
          * @brief Append specific path
          */
         void appendPath();

   };

   template <typename TOp, unsigned int N> Tester<TOp,N>::Tester(const std::string& fname, const bool keepData)
      : Polynomial::TesterBase<TOp>(fname, keepData)
   {
      this->appendPath();
   }

   template <typename TOp, unsigned int N> void Tester<TOp,N>::appendPath()
   {
      this->mPath += "Quadrature/";
   }

   template <typename TOp, unsigned int N> Matrix Tester<TOp,N>::buildOperator(const int specN, const int physN, const ParameterType& param, const TestType type) const
   {
      assert(type == TestType::QUADRATURE);

      // Create quadrature
      Internal::Array igrid;
      Internal::Array iweights;

      TOp quad = this->createQuad<>(param);
      quad.computeQuadrature(igrid, iweights, physN);

      Matrix outData(igrid.size(),2);
      outData.col(0) = igrid.cast<MHDFloat>();
      outData.col(1) = iweights.cast<MHDFloat>();

      return outData;
   }

   template <typename TOp, unsigned int N> std::string Tester<TOp,N>::resname(const int specN, const int physN, const ParameterType& param) const
   {
      std::string s = "";

      unsigned int ia = static_cast<unsigned int>('a');;
      for(std::size_t i = 0; i < param.size(); i++)
      {
         s += "_" + std::string(1,static_cast<char>(ia+i)) + std::to_string(param.at(i));
      }

      s += "_g" + std::to_string(physN);

      return s;
   }

   template <typename TOp, unsigned int N> std::string Tester<TOp,N>::formatParameter(const ParameterType& param) const
   {
      std::stringstream ss;
      unsigned int ia = static_cast<unsigned int>('a');;
      for(std::size_t i = 0; i < param.size(); i++)
      {
         ss << std::string(1,static_cast<char>(ia+i)) + ": " + std::to_string(param.at(i));
         if(i < param.size()-1)
         {
            ss << ", ";
         }
      }

      return ss.str();
   }

   template <typename TOp, unsigned int N> bool Tester<TOp,N>::isSpecial(const int i, const int j, const MHDFloat data, const MHDFloat ref) const
   {
      if(ref == 0.0 && data != ref)
      {
         return true;
      }
      else
      {
         return false;
      }
   }

   template <typename TOp, unsigned int N> typename Tester<TOp,N>::ErrorType Tester<TOp,N>::checkSpecial(const int i, const int j, const MHDFloat data, const Matrix& refData) const
   {
      // quadratures defined on O(1) interval
      MHDFloat ref = 1.0;
      return this->checkNormal(data, refData(i,j), ref);
   }

} // Quadrature
} // Polynomial
} // TestSuite
} // QuICC

#endif //QUICC_TESTSUITE_POLYNOMIAL_QUADRATURE_TESTER_HPP

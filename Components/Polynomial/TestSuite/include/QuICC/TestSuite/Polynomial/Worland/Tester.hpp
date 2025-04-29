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
#include <type_traits>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandRule.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
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
          *
          * @param fname      Filename
          * @param wtype      Worland type
          * @param keepData   Write computed data?
          */
         Tester(const std::string& fname, const std::string& wtype, const bool keepData);

         /*
          * @brief Destructor
          */
         virtual ~Tester() = default;

      protected:
         /**
          * @brief Build filename extension with resolution information
          *
          * @param specN   Number of spectral modes
          * @param physN   Number of physical grid points
          * @param param   Run parameters
          */
         virtual std::string resname(const int specN, const int physN, const ParameterType& param) const override;

         /**
          * @brief Test matrix operator
          *
          * @param specN   Number of spectral modes
          * @param physN   Number of physical grid points
          * @param param   Run parameters
          * @param type    Type of test
          */
         virtual Matrix buildOperator(const int specN, const int physN, const ParameterType& param, const TestType type) const override;

         /**
          * @brief Format the parameters
          *
          * @param param   Run parameters
          */
         virtual std::string formatParameter(const ParameterType& param) const override;

      private:
         /**
          * @brief Append specific path
          *
          * @param subdir  Sub-directory to append to path
          */
         void appendPath(const std::string& subdir);

         /**
          * @brief Test quadrature matrix operator
          *
          * @param specN      Number of spectral modes
          * @param physN      Number of physical grid points
          * @param param      Run parameters
          * @param isWeight   Quadrature weights included in matrix?
          */
         Matrix buildMatrix(const int specN, const int physN, const ParameterType& param, const bool isWeighted) const;

         /**
          * @brief Test on-the-fly inner product
          *
          * @param specN   Number of spectral modes
          * @param physN   Number of physical grid points
          * @param param   Run parameters
          */
         Matrix buildInner(const int specN, const int physN, const ParameterType& param) const;

         /**
          * @brief Test on-the-fly outer product
          *
          * @param specN   Number of spectral modes
          * @param physN   Number of physical grid points
          * @param param   Run parameters
          */
         Matrix buildOuter(const int specN, const int physN, const ParameterType& param) const;

         /**
          * @brief Test on-the-fly reduction
          *
          * @param specN   Number of spectral modes
          * @param physN   Number of physical grid points
          * @param param   Run parameters
          */
         Matrix buildReduce(const int specN, const int physN, const ParameterType& param) const;

         /**
          * @brief Make operator and compute grid and weights
          *
          * @param igrid      Array of grid points
          * @param iweights   Array of weights
          * @param physN      Number of physical grid points
          */
         std::shared_ptr<TOp> makeOp(Internal::Array& igrid, Internal::Array& iweights, const int physN) const;

         /**
          * @brief Worland type
          */
         std::string mWType;
   };

   template <typename TOp> Tester<TOp>::Tester(const std::string& fname, const std::string& wtype, const bool keepData)
      : Polynomial::TesterBase<TOp>(fname, keepData), mWType(wtype)
   {

      this->appendPath(this->mWType);
   }

   template <typename TOp> void Tester<TOp>::appendPath(const std::string& subdir)
   {
      this->mPath += "Worland/" + subdir + "/";
   }

   template <typename TOp> std::shared_ptr<TOp> Tester<TOp>::makeOp(Internal::Array& igrid, Internal::Array& iweights, const int physN) const
   {
      auto make = [&](auto& w)
      {
         typename std::remove_reference<decltype(w)>::type::Rule quad;
         quad.computeQuadrature(igrid, iweights, physN);
         auto pOp = std::make_shared<TOp>(w.ALPHA, w.DBETA);
         return pOp;
      };

      std::shared_ptr<TOp> pOp;


      if(this->mWType == "Chebyshev")
      {
         ::QuICC::Polynomial::Worland::worland_chebyshev_t wt;
         pOp = make(wt);
      }
      else if(this->mWType == "Legendre")
      {
         ::QuICC::Polynomial::Worland::worland_legendre_t wt;
         pOp = make(wt);
      }
      else if(this->mWType == "CylEnergy")
      {
         ::QuICC::Polynomial::Worland::worland_cylenergy_t wt;
         pOp = make(wt);
      }
      else if(this->mWType == "SphEnergy")
      {
         ::QuICC::Polynomial::Worland::worland_sphenergy_t wt;
         pOp = make(wt);
      }
      else
      {
         throw std::logic_error("Unknown Worland type");
      }

      return pOp;
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

      // Compute quadrature and make operator
      Internal::Array igrid;
      Internal::Array iweights;
      auto pOp = this->makeOp(igrid, iweights, physN);

      Matrix outData(physN, specN+1);

      if(isWeighted)
      {
         pOp->template compute<MHDFloat>(outData, specN+1, l, igrid, iweights, current::Evaluator::Set());
      } else
      {
         pOp->template compute<MHDFloat>(outData, specN+1, l, igrid, Internal::Array(), current::Evaluator::Set());
      }

      return outData;
   }

   template <typename TOp> Matrix Tester<TOp>::buildOuter(const int specN, const int physN, const ParameterType& param) const
   {
      // Get harmonic degree from parameters
      int l = static_cast<int>(param.at(0));

      // Create quadrature
      Internal::Array igrid;
      Internal::Array iweights;
      auto pOp = this->makeOp(igrid, iweights, physN);

      Matrix outData(physN, specN+1);
      Matrix inData = Matrix::Identity(specN+1, specN+1);

      pOp->template compute<MHDFloat>(outData, specN+1, l, igrid, Internal::Array(), current::Evaluator::OuterProduct<MHDFloat>(inData));

      return outData;
   }

   template <typename TOp> Matrix Tester<TOp>::buildInner(const int specN, const int physN, const ParameterType& param) const
   {
      // Get harmonic degree from parameters
      int l = static_cast<int>(param.at(0));

      // Create quadrature
      Internal::Array igrid;
      Internal::Array iweights;
      auto pOp = this->makeOp(igrid, iweights, physN);

      Matrix outData(specN+1, specN+1);
      Matrix inData = Matrix::Identity(physN, specN+1);
      pOp->template compute<MHDFloat>(inData, specN+1, l, igrid, iweights, current::Evaluator::Set());

      pOp->template compute<MHDFloat>(outData, specN+1, l, igrid, Internal::Array(), current::Evaluator::InnerProduct<MHDFloat>(inData));

      Matrix outDataT = outData.transpose();

      return outDataT;
   }

   template <typename TOp> Matrix Tester<TOp>::buildReduce(const int specN, const int physN, const ParameterType& param) const
   {
      // Get harmonic degree from parameters
      int l = static_cast<int>(param.at(0));

      // Create quadrature
      Internal::Array igrid;
      Internal::Array iweights;
      auto pOp = this->makeOp(igrid, iweights, physN);

      Matrix outData(specN+1, 1);

      pOp->template compute<MHDFloat>(outData, specN+1, l, igrid, Internal::Array(), current::Evaluator::Reduce());

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

/**
 * @file Tester.hpp
 * @brief Tester for Worland transforms
 */

#ifndef QUICC_TESTSUITE_DENSESM_WORLAND_TESTER_HPP
#define QUICC_TESTSUITE_DENSESM_WORLAND_TESTER_HPP

// System includes
//
#include <catch2/catch.hpp>
#include <string>
#include <set>
#include <sstream>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/GridPurpose.hpp"
#include "TestSuite/DenseSM/TesterBase.hpp"
#include "DenseSM/Worland/Geostrophic2Tor.hpp"
#include "DenseSM/Worland/PyGeostrophic2Tor.hpp"
#include "DenseSM/Worland/GeostrophicAngularMomentum.hpp"
#include "QuICC/Polynomial/Worland/WorlandBase.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"

namespace dsm = ::QuICC::DenseSM::Worland;

namespace QuICC {

namespace TestSuite {

namespace DenseSM {

namespace Worland {

   template <typename TOp, typename TOp2 = void> class Tester: public DenseSM::TesterBase<TOp>
   {
      public:
         /// Typedef for parameter type
         typedef typename DenseSM::TesterBase<TOp>::ParameterType ParameterType;

         /**
          * @brief Constructor
          */
         Tester(const std::string& fname, const bool keepData);

         /*
          * @brief Destructor
          */
         virtual ~Tester() = default;

      protected:
         /// Typedef ContentType from base
         typedef typename DenseSM::TesterBase<TOp>::ContentType ContentType;

         /**
          * @brief Build filename extension with resolution information
          */
         virtual std::string resname(const ParameterType& param) const override;

         /**
          * @brief Read real data from file
          */
         virtual void readFile(Matrix& data, const ParameterType& param, const TestType type, const ContentType ctype) const override;

         /**
          * @brief Read complex data from file
          */
         virtual void readFile(MatrixZ& data, const ParameterType& param, const TestType type, const ContentType ctype) const override;

         /**
          * @brief Test operator
          */
         virtual Matrix applyOperator(const ParameterType& param, const TestType type) const override;

         /**
          * @brief Format the parameters
          */
         virtual std::string formatParameter(const ParameterType& param) const override;

      private:
         /**
          * @brief Append specific path
          */
         void appendPath();
   };

   template <typename TOp, typename TOp2> Tester<TOp,TOp2>::Tester(const std::string& fname, const bool keepData)
      : DenseSM::TesterBase<TOp>(fname, keepData, false)
   {
      this->appendPath();
   }

   template <typename TOp, typename TOp2> void Tester<TOp,TOp2>::appendPath()
   {
      this->mPath += "Worland/";
   }

   template <typename TOp, typename TOp2> void Tester<TOp,TOp2>::readFile(Matrix& data, const ParameterType& param, const TestType type, const ContentType ctype) const
   {
      TesterBase<TOp>::basicReadFile(data, param, type, ctype);
   }

   template <typename TOp, typename TOp2> void Tester<TOp,TOp2>::readFile(MatrixZ& data, const ParameterType& param, const TestType type, const ContentType ctype) const
   {
      TesterBase<TOp>::basicReadFile(data, param, type, ctype);
   }

   template <typename TOp, typename TOp2> Matrix Tester<TOp,TOp2>::applyOperator(const ParameterType& param, const TestType type) const
   {
      Matrix outData;

      if constexpr(std::is_same_v<TOp, dsm::Geostrophic2Tor> || std::is_same_v<TOp, dsm::PyGeostrophic2Tor>)
      {
         auto a = ::QuICC::Polynomial::Worland::WorlandBase::ALPHA_CHEBYSHEV;
         auto b = ::QuICC::Polynomial::Worland::WorlandBase::DBETA_CHEBYSHEV;

         Array meta(0);
         std::string fullname = this->makeFilename(param, this->refRoot(), type, ContentType::META);
         readList(meta, fullname);

         int nN = meta(0) + 1;
         int maxnl = meta(1) + 1;
         QuICC::internal::MHDFloat ugAlpha = static_cast<QuICC::internal::MHDFloat>(meta(2));
         QuICC::internal::MHDFloat ugDBeta = static_cast<QuICC::internal::MHDFloat>(meta(3));

         int nR = int(((maxnl + 1) - (maxnl + 1) % 2) / 2 + 2) + 1;
         int maxNug = int(((maxnl - 2) - (maxnl - 2) % 2) / 2);
         int& nug = maxNug;
         QuICC::ArrayI nli;
         if (2 * nug + 1 > maxnl - 1)
         {
            throw std::logic_error("L truncation is not enough to capture all geostrophic modes required");
         }

         nli.resize(maxnl);
         nli.setConstant(0);

         for(int l = 0; l < maxnl; l++)
         {
            if(l % 2 == 0) nli(l) = -1;
         }

         for(int k = 0; k <= nug; k++)
         {
            nli(2*k+1) = nug - k;
         }

         std::vector<int> nIdx;
         for(int n = 0; n < maxNug+1; n++)
         {
            if(QuICC::QuICCEnv().id() == n%QuICCEnv().size())
            {
               nIdx.push_back(n);
            }
         }
         TOp op(nN, maxnl, nR, maxNug, nli, nIdx, ugAlpha, ugDBeta, a, b, 0);

         outData = op.mat();

         #if defined QUICC_MPI
            MPI_Allreduce(MPI_IN_PLACE, outData.data(), outData.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         #endif
      }
      else if constexpr(std::is_same_v<TOp, dsm::GeostrophicAngularMomentum>)
      {
         auto a = ::QuICC::Polynomial::Worland::WorlandBase::ALPHA_CHEBYSHEV;
         auto b = ::QuICC::Polynomial::Worland::WorlandBase::DBETA_CHEBYSHEV;

         Array meta(0);
         std::string fullname = this->makeFilename(param, this->refRoot(), type, ContentType::META);
         readList(meta, fullname);

         int nN = meta(0) + 1;
         int maxnl = meta(1) + 1;
         QuICC::internal::MHDFloat ugAlpha = static_cast<QuICC::internal::MHDFloat>(meta(2));
         QuICC::internal::MHDFloat ugDBeta = static_cast<QuICC::internal::MHDFloat>(meta(3));

         int nr = (maxnl - 3)/2 + 1;

         TOp op(ugAlpha, ugDBeta, nr, a, b, 0);

         outData = op.mat();
      }
      else
      {
         throw std::logic_error("Not implemented");
      }

      return outData;
   }

   template <typename TOp, typename TOp2> std::string Tester<TOp,TOp2>::resname(const ParameterType& param) const
   {
      auto id = param.at(0);

      std::stringstream ss;
      ss.precision(10);
      ss << "_id" << id;

      // Distributed data meta file
      if(param.size() == 3)
      {
         int np = param.at(1);
         ss << "_np" << np;
         int r = param.at(2);
         ss << "_r" << r;
         ss << "_stage0";
      }

      return ss.str();
   }

   template <typename TOp, typename TOp2> std::string Tester<TOp,TOp2>::formatParameter(const ParameterType& param) const
   {
      auto id = param.at(0);

      std::stringstream ss;
      ss << "id: " << id;

      if(param.size() == 3)
      {
         int np = param.at(1);
         ss << ", np: " << np;
         int r = param.at(2);
         ss << ", r: " << r;
      }

      return ss.str();
   }

} // Worland
} // DenseSM
} // TestSuite
} // QuICC

#endif //QUICC_TESTSUITE_DENSESM_WORLAND_TESTER_HPP

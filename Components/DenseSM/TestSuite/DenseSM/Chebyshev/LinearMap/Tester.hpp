/**
 * @file Tester.hpp
 * @brief Tester for Worland transforms
 */

#ifndef QUICC_TESTSUITE_DENSESM_CHEBYSHEV_LINEARMAP_TESTER_HPP
#define QUICC_TESTSUITE_DENSESM_CHEBYSHEV_LINEARMAP_TESTER_HPP

// System includes
//
#include <catch2/catch.hpp>
#include <string>
#include <sstream>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/ILinearMapOperator.hpp"
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"
#include "DenseSM/Chebyshev/LinearMap/ITripleHarmonicOperator.hpp"
#include "DenseSM/Chebyshev/LinearMap/IProjCrossOperator.hpp"
#include "TestSuite/DenseSM/TesterBase.hpp"
#include "TestSuite/DenseSM/Chebyshev/LinearMap/DipolarS1.hpp"
#include "TestSuite/DenseSM/Chebyshev/LinearMap/QuadrupolarS2.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/Insulating.hpp"

namespace dsm = ::QuICC::DenseSM::Chebyshev::LinearMap;

namespace QuICC {

namespace TestSuite {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

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
      this->mPath += "Chebyshev/LinearMap/";
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

      if constexpr(std::is_base_of_v<dsm::ITripleHarmonicOperator, TOp>)
      {
         Array meta(0);
         std::string fullname = this->makeFilename(param, this->refRoot(), type, ContentType::META);
         readList(meta, fullname);
         if(meta.size() != 11)
         {
            throw std::logic_error("Test meta data is wrong");
         }

         int nNr = meta(0) + 1;
         int nNc = meta(1) + 1;
         int lOut = meta(2);
         int mOut = meta(3);
         int lF = meta(4);
         int mF = meta(5);
         int lIn = meta(6);
         int mIn = meta(7);
         int fId = meta(8);
         auto lb = static_cast<QuICC::Internal::MHDFloat>(meta(9));
         auto ub = static_cast<QuICC::Internal::MHDFloat>(meta(10));

         std::shared_ptr<dsm::RadialTorPolFunction> pF;
         if (fId == 0)
         {
            pF = std::make_shared<DipolarS1>();
         }
         else if (fId == 1)
         {
            pF = std::make_shared<QuadrupolarS2>();
         }
         else
         {
            throw std::logic_error("Unknown forcing function ID");
         }

         TOp op(nNr, nNc, lOut, mOut, lF, mF, lIn, mIn, pF, lb, ub);

         outData = op.mat();
      }
      else if constexpr(std::is_base_of_v<dsm::IProjCrossOperator, TOp>)
      {
         Array meta(0);
         std::string fullname = this->makeFilename(param, this->refRoot(), type, ContentType::META);
         readList(meta, fullname);
         if(meta.size() != 13)
         {
            throw std::logic_error("Test meta data is wrong");
         }

         int nNr = meta(0) + 1;
         int nNc = meta(1) + 1;
         int p = meta(2);
         int lOut = meta(3);
         int mOut = meta(4);
         int lF = meta(5);
         int mF = meta(6);
         int lIn = meta(7);
         int mIn = meta(8);
         int fAId = meta(9);
         int fBId = meta(10);
         auto lb = static_cast<QuICC::Internal::MHDFloat>(meta(11));
         auto ub = static_cast<QuICC::Internal::MHDFloat>(meta(12));

         int lA, mA, lB, mB;
         std::shared_ptr<dsm::RadialTorPolFunction> pFa = nullptr;
         std::shared_ptr<dsm::RadialTorPolFunction> pFb = nullptr;
         if (fAId >= 0)
         {
            if (fAId == 0)
            {
               pFa = std::make_shared<DipolarS1>();
            }
            else if(fAId == 1)
            {
               pFa = std::make_shared<QuadrupolarS2>();
            }
            else
            {
               throw std::logic_error("Unknown forcing function");
            }
            lA = lF;
            mA = mF;
            lB = lIn;
            mB = mIn;
         }
         if (fBId >= 0)
         {
            if (fBId == 0)
            {
               pFb = std::make_shared<DipolarS1>();
            }
            else if(fBId == 1)
            {
               pFb = std::make_shared<QuadrupolarS2>();
            }
            else
            {
               throw std::logic_error("Unknown forcing function");
            }
            lA = lIn;
            mA = mIn;
            lB = lF;
            mB = mF;
         }

         TOp op(nNr, nNc, p, lOut, mOut, lA, mA, lB, mB, pFa, pFb, lb, ub);

         outData = op.mat();
      }
      else if constexpr(std::is_base_of_v<dsm::ILinearMapOperator, TOp>)
      {
         Array meta(0);
         std::string fullname = this->makeFilename(param, this->refRoot(), type, ContentType::META);
         readList(meta, fullname);
         if(meta.size() != 6)
         {
            throw std::logic_error("Test meta data is wrong");
         }

         int outRows = meta(0) + 1;
         int bc = static_cast<int>(meta(1));
         std::size_t bcId;
         int nN = meta(2) + 1;
         auto lb = static_cast<QuICC::Internal::MHDFloat>(meta(3));
         auto ub = static_cast<QuICC::Internal::MHDFloat>(meta(4));

         // Identify boundary condition
         if(bc == 0)
         {
            bcId = Bc::Name::FixedTemperature::id();
         }
         else if(bc == 1)
         {
            bcId = Bc::Name::FixedFlux::id();
         }
         else if(bc == 2)
         {
            bcId = Bc::Name::Insulating::id();
         }
         else
         {
            throw std::logic_error("Unknown boundary condition");
         }

         TOp op(outRows, bcId, nN, nN, lb, ub);

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

} // LinearMap
} // Chebyshev
} // DenseSM
} // TestSuite
} // QuICC

#endif //QUICC_TESTSUITE_DENSESM_CHEBYSHEV_LINEARMAP_TESTER_HPP

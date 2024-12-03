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
#include "DenseSM/Worland/DipolarS1.hpp"
#include "DenseSM/Worland/IWorlandOperator.hpp"
#include "DenseSM/Worland/RadialTorPolFunction.hpp"
#include "Types/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Enums/GridPurpose.hpp"
#include "TestSuite/DenseSM/TesterBase.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/Insulating.hpp"
#include "DenseSM/Worland/ITripleHarmonicOperator.hpp"

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

      if constexpr(std::is_base_of_v<dsm::ITripleHarmonicOperator, TOp>)
      {
         Array meta(0);
         std::string fullname = this->makeFilename(param, this->refRoot(), type, ContentType::META);
         readList(meta, fullname);
         if(meta.size() != 8)
         {
            throw std::logic_error("Test meta data is wrong");
         }

         int nNr = meta(0) + 1;
         int nNc = meta(1) + 1;
         int lOut = meta(2);
         int lF = meta(3);
         int lIn = meta(4);
         int fId = meta(5);
         auto alpha = static_cast<QuICC::Internal::MHDFloat>(meta(6));
         auto dBeta = static_cast<QuICC::Internal::MHDFloat>(meta(7));

         std::shared_ptr<dsm::RadialTorPolFunction> pF;
         if (fId == 0)
         {
            pF = std::make_shared<dsm::DipolarS1>();
         }

         TOp op(nNr, nNc, lOut, lF, lIn, pF, alpha, dBeta);

         outData = op.mat();

         #if defined QUICC_MPI
            MPI_Allreduce(MPI_IN_PLACE, outData.data(), outData.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         #endif
      }
      else if constexpr(std::is_base_of_v<dsm::IWorlandOperator, TOp>)
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
         auto a = static_cast<QuICC::Internal::MHDFloat>(meta(3));
         auto b = static_cast<QuICC::Internal::MHDFloat>(meta(4));
         auto l = static_cast<int>(meta(5));

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

         TOp op(outRows, bcId, nN, nN, a, b, l);

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

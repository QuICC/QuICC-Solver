/**
 * @file Tester.hpp
 * @brief Tester for Bessel DenseSM
 */

#ifndef QUICC_TESTSUITE_DENSESM_BESSEL_TESTER_HPP
#define QUICC_TESTSUITE_DENSESM_BESSEL_TESTER_HPP

// System includes
//
#include <catch2/catch.hpp>
#include <string>
#include <set>
#include <sstream>

// Project includes
//
#include "Types/BasicTypes.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/GridPurpose.hpp"
#include "TestSuite/DenseSM/TesterBase.hpp"
#include "DenseSM/IMatrixSMOperator.hpp"
#include "DenseSM/Bessel/CoriolisQm.hpp"
#include "DenseSM/Bessel/CoriolisQp.hpp"
#include "DenseSM/Bessel/Geostrophic2Tor.hpp"
#include "DenseSM/Bessel/GeostrophicAngularMomentum.hpp"
#include "DenseSM/Bessel/GeostrophicEnergy.hpp"
#include "DenseSM/Bessel/Tor2GridS.hpp"
#include "DenseSM/Bessel/Tor2Geostrophic.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/Insulating.hpp"

namespace dsm = ::QuICC::DenseSM;

namespace QuICC {

namespace TestSuite {

namespace DenseSM {

namespace Bessel {

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
      this->mPath += "Bessel/";
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

      if constexpr(std::is_same_v<dsm::Bessel::CoriolisQm, TOp> || std::is_same_v<dsm::Bessel::CoriolisQp, TOp>)
      {
         Array meta(0);
         std::string fullname = this->makeFilename(param, this->refRoot(), type, ContentType::META);
         readList(meta, fullname);
         assert(meta.size() == 4);

         Internal::MHDFloat outDNu = meta(0);
         Internal::MHDFloat inDNu = meta(1);
         int nN = meta(2) + 1;
         auto l = static_cast<int>(meta(3));

         TOp op(outDNu, inDNu, nN, nN, l);

         outData = op.mat();
      }
      else if constexpr(std::is_same_v<TOp, dsm::Bessel::Geostrophic2Tor>)
      {
         Array meta(0);
         std::string fullname = this->makeFilename(param, this->refRoot(), type, ContentType::META);
         readList(meta, fullname);
         if(meta.size() != 4)
         {
            throw std::logic_error("Test meta data is wrong");
         }

         int nN = meta(0) + 1;
         int nL = meta(1) + 1;
         auto torDNu = static_cast<QuICC::Internal::MHDFloat>(meta(2));
         auto sDNu = static_cast<QuICC::Internal::MHDFloat>(meta(3));

         int nCpu = QuICCEnv().size();
         TOp op(nN, nL, nCpu, sDNu, torDNu);

         outData = op.mat();

         #if defined QUICC_MPI
            MPI_Allreduce(MPI_IN_PLACE, outData.data(), outData.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         #endif
      }
      else if constexpr(std::is_same_v<TOp, dsm::Bessel::GeostrophicAngularMomentum>)
      {
         Array meta(0);
         std::string fullname = this->makeFilename(param, this->refRoot(), type, ContentType::META);
         readList(meta, fullname);
         if(meta.size() != 2)
         {
            throw std::logic_error("Test meta data is wrong");
         }

         int nN = meta(0) + 1;
         auto sDNu = static_cast<QuICC::Internal::MHDFloat>(meta(1));

         TOp op(nN, sDNu);

         outData = op.mat();
      }
      else if constexpr(std::is_same_v<TOp, dsm::Bessel::GeostrophicEnergy>)
      {
         Array meta(0);
         std::string fullname = this->makeFilename(param, this->refRoot(), type, ContentType::META);
         readList(meta, fullname);
         if(meta.size() != 3)
         {
            throw std::logic_error("Test meta data is wrong");
         }

         int nN = meta(0) + 1;
         int nL = meta(1) + 1;
         auto sDNu = static_cast<QuICC::Internal::MHDFloat>(meta(2));

         TOp op(nN, nL, sDNu);

         outData = op.mat();
      }
      else if constexpr(std::is_same_v<dsm::Bessel::Tor2Geostrophic, TOp>)
      {
         Array meta(0);
         std::string fullname = this->makeFilename(param, this->refRoot(), type, ContentType::META);
         readList(meta, fullname);
         if(meta.size() != 4)
         {
            throw std::logic_error("Test meta data is wrong");
         }

         int nN = meta(0) + 1;
         int nL = meta(1) + 1;
         auto torDNu = static_cast<QuICC::Internal::MHDFloat>(meta(2));
         auto sDNu = static_cast<QuICC::Internal::MHDFloat>(meta(3));

         int nCpu = QuICCEnv().size();
         TOp op(nN, nL, nCpu, sDNu, torDNu);

         outData = op.mat();
      }
      else if constexpr(std::is_same_v<dsm::Bessel::Tor2GridS, TOp>)
      {
         Array meta(0);
         std::string fullname = this->makeFilename(param, this->refRoot(), type, ContentType::META);
         readList(meta, fullname);
         if(meta.size() != 4)
         {
            throw std::logic_error("Test meta data is wrong");
         }

         int nN = meta(0) + 1;
         int nL = meta(1) + 1;
         auto torDNu = static_cast<QuICC::Internal::MHDFloat>(meta(2));
         auto sDNu = static_cast<QuICC::Internal::MHDFloat>(meta(3));

         int nCpu = QuICCEnv().size();
         TOp op(nN, nL, nCpu, sDNu, torDNu);

         outData = op.mat();
      }
      else if constexpr(std::is_base_of_v<dsm::IMatrixSMOperator, TOp>)
      {
         Array meta(0);
         std::string fullname = this->makeFilename(param, this->refRoot(), type, ContentType::META);
         readList(meta, fullname);
         assert(meta.size() == 4);

         int outRows = meta(0) + 1;
         int bc = static_cast<int>(meta(1));
         std::size_t bcId;
         int nN = meta(2) + 1;
         auto l = static_cast<int>(meta(3));

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

         TOp op(outRows, bcId, nN, nN, l);

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

} // Bessel
} // DenseSM
} // TestSuite
} // QuICC

#endif //QUICC_TESTSUITE_DENSESM_BESSEL_TESTER_HPP

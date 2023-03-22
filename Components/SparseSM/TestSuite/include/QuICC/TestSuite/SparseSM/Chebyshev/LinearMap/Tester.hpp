/**
 * @file Tester.hpp
 * @brief Tester for Worland transforms
 */

#ifndef QUICC_TESTSUITE_SPARSESM_CHEBYSHEV_LINEARMAP_TESTER_HPP
#define QUICC_TESTSUITE_SPARSESM_CHEBYSHEV_LINEARMAP_TESTER_HPP


// System includes
//
#include <catch2/catch.hpp>
#include <string>
#include <set>
#include <sstream>
#include <type_traits>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/GridPurpose.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/ISphericalOperator.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/ICondition.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/Operator.hpp"
#include "QuICC/TestSuite/SparseSM/TesterBase.hpp"

namespace QuICC {

namespace TestSuite {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   template <typename TOp> class Tester: public SparseSM::TesterBase<TOp>
   {
      public:
         /// Typedef for parameter type
         typedef typename SparseSM::TesterBase<TOp>::ParameterType ParameterType;

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
         typedef typename SparseSM::TesterBase<TOp>::ContentType ContentType;

         /**
          * @brief Build filename extension with resolution information
          */
         virtual std::string resname(const ParameterType& param) const override;

         /**
          * @brief Test operator
          */
         virtual void buildOperator(SparseMatrix& mat, const ParameterType& param, const TestType type) const override;

         /**
          * @brief Test operator
          */
         virtual void buildOperator(Matrix& mat, const ParameterType& param, const TestType type) const override;

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
          * @brief Read meta data
          */
         Array readMeta(const ParameterType& param, const TestType type) const;

   };

   template <typename TOp> Tester<TOp>::Tester(const std::string& fname, const bool keepData)
      : SparseSM::TesterBase<TOp>(fname, keepData)
   {
      this->appendPath();
   }

   template <typename TOp> void Tester<TOp>::appendPath()
   {
      this->mPath += "Chebyshev/LinearMap/";
   }

   template <typename TOp> void Tester<TOp>::buildOperator(SparseMatrix& mat, const ParameterType& param, const TestType type) const
   {
      Array meta = this->readMeta(param, type);

      if constexpr(std::is_base_of_v<QuICC::SparseSM::Chebyshev::LinearMap::ISphericalOperator, TOp>)
      {
         assert(meta.size() == 5);
         int rows = static_cast<int>(meta(0));
         int cols = static_cast<int>(meta(1));
         auto lower = meta(2);
         auto upper = meta(3);
         auto l = meta(4);

         TOp op(rows, cols, lower, upper, l);
         mat = op.mat();
      }
      else if constexpr(std::is_base_of_v<QuICC::SparseSM::Chebyshev::ILinearMapOperator, TOp>)
      {
         assert(meta.size() == 4);
         int rows = static_cast<int>(meta(0));
         int cols = static_cast<int>(meta(1));
         auto lower = meta(2);
         auto upper = meta(3);

         TOp op(rows, cols, lower, upper);
         mat = op.mat();
      }
      else
      {
         throw std::logic_error("Could not build sparse operator");
      }
   }

   template <typename TOp> void Tester<TOp>::buildOperator(Matrix& mat, const ParameterType& param, const TestType type) const
   {
      Array meta = this->readMeta(param, type);

      unsigned int KL, KU;
      if constexpr(std::is_base_of_v<QuICC::SparseSM::Chebyshev::LinearMap::ISphericalOperator, TOp>)
      {
         int rows = static_cast<int>(meta(0));
         int cols = static_cast<int>(meta(1));
         auto lower = meta(2);
         auto upper = meta(3);

         assert(meta.size() == 5);
         auto l = meta(4);
         TOp op(rows, cols, lower, upper, l);
         mat = op.banded(KL, KU);
      }
      else if constexpr(std::is_base_of_v<QuICC::SparseSM::Chebyshev::ILinearMapOperator, TOp>)
      {
         int rows = static_cast<int>(meta(0));
         int cols = static_cast<int>(meta(1));
         auto lower = meta(2);
         auto upper = meta(3);

         TOp op(rows, cols, lower, upper);
         mat = op.banded(KL, KU);
      }
      else if constexpr(std::is_base_of_v<QuICC::SparseSM::Chebyshev::LinearMap::Boundary::ICondition, TOp>)
      {
         auto lower = meta(0);
         auto upper = meta(1);
         int pos = static_cast<int>(meta(2));
         int maxN = static_cast<int>(meta(3));

         QuICC::SparseSM::Chebyshev::LinearMap::Boundary::Operator bcOp(1, maxN+1, lower, upper);

         if constexpr(std::is_constructible_v<TOp, internal::MHDFloat, internal::MHDFloat, typename TOp::Position, int>)
         {
            assert(meta.size() == 5);
            auto l = meta(4);
            bcOp.addRow<TOp>(static_cast<typename TOp::Position>(pos), l);
            mat = bcOp.mat();
         }
         else
         {
            bcOp.addRow<TOp>(static_cast<typename TOp::Position>(pos));
            mat = bcOp.mat();
         }
      }
      else
      {
         throw std::logic_error("Could not build dense operator");
      }
   }

   template <typename TOp> std::string Tester<TOp>::resname(const ParameterType& param) const
   {
      auto id = param.at(0);

      std::stringstream ss;
      ss.precision(10);
      ss << "_id" << id;

      return ss.str();
   }

   template <typename TOp> std::string Tester<TOp>::formatParameter(const ParameterType& param) const
   {
      auto id = param.at(0);

      std::stringstream ss;
      ss << "id: " << id;

      return ss.str();
   }

   template <typename TOp> Array Tester<TOp>::readMeta(const ParameterType& param, const TestType type) const
   {
      // Read metadata
      Array meta(0);
      std::string fullname = this->makeFilename(param, this->refRoot(), type, ContentType::META);
      readList(meta, fullname);

      return meta;
   }

}
}
}
}
}

#endif //QUICC_TESTSUITE_SPARSESM_CHEBYSHEV_LINEARMAP_TESTER_HPP

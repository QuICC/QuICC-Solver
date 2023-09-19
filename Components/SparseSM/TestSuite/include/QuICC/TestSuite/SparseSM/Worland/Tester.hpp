/**
 * @file Tester.hpp
 * @brief Tester for Worland transforms
 */

#ifndef QUICC_TESTSUITE_SPARSESM_WORLAND_TESTER_HPP
#define QUICC_TESTSUITE_SPARSESM_WORLAND_TESTER_HPP


// System includes
//
#include <catch2/catch.hpp>
#include <string>
#include <set>
#include <sstream>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/GridPurpose.hpp"
#include "QuICC/TestSuite/SparseSM/TesterBase.hpp"
#include "QuICC/SparseSM/Worland/Boundary/ICondition.hpp"
#include "QuICC/SparseSM/Worland/Boundary/Operator.hpp"

namespace QuICC {

namespace TestSuite {

namespace SparseSM {

namespace Worland {

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
      this->mPath += "Worland/";
   }

   template <typename TOp> void Tester<TOp>::buildOperator(SparseMatrix& mat, const ParameterType& param, const TestType type) const
   {
      Array meta = this->readMeta(param, type);

      if constexpr(std::is_base_of_v<QuICC::SparseSM::IWorlandOperator, TOp>)
      {
         int rows = static_cast<int>(meta(0));
         int cols = static_cast<int>(meta(1));
         auto alpha = meta(2);
         auto dbeta = meta(3);
         auto l = meta(4);
         int q = 0;
         // Check for truncation in meta data
         if(meta.size() > 5)
         {
            q = meta(5);
         }

         TOp op(rows, cols, alpha, dbeta, l, q);
         mat = op.mat();
      }
   }

   template <typename TOp> void Tester<TOp>::buildOperator(Matrix& mat, const ParameterType& param, const TestType type) const
   {
      Array meta = this->readMeta(param, type);

      if constexpr(std::is_base_of_v<QuICC::SparseSM::IWorlandOperator, TOp>)
      {
         int rows = static_cast<int>(meta(0));
         int cols = static_cast<int>(meta(1));
         auto alpha = meta(2);
         auto dbeta = meta(3);
         auto l = meta(4);
         int q = 0;
         // Check for truncation in meta data
         if(meta.size() > 5)
         {
            q = meta(5);
         }

         TOp op(rows, cols, alpha, dbeta, l, q);
         unsigned int KL, KU;
         mat = op.banded(KL, KU);
      }
      else if constexpr(std::is_base_of_v<QuICC::SparseSM::Worland::Boundary::ICondition, TOp>)
      {
         auto alpha = meta(0);
         auto dbeta = meta(1);
         auto l = meta(2);
         int maxN = static_cast<int>(meta(3));

         QuICC::SparseSM::Worland::Boundary::Operator bcOp(1, maxN+1, alpha, dbeta, l, true);

         bcOp.addRow<TOp>();
         mat = bcOp.mat();
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

#endif //QUICC_TESTSUITE_SPARSESM_WORLAND_TESTER_HPP

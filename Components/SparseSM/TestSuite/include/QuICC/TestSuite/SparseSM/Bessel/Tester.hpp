/**
 * @file Tester.hpp
 * @brief Tester for Bessel transforms
 */

#ifndef QUICC_TESTSUITE_SPARSESM_BESSEL_TESTER_HPP
#define QUICC_TESTSUITE_SPARSESM_BESSEL_TESTER_HPP


// System includes
//
#include <catch2/catch.hpp>
#include <set>
#include <sstream>
#include <string>

// Project includes
//
#include "QuICC/Enums/GridPurpose.hpp"
#include "QuICC/SparseSM/IBesselOperator.hpp"
#include "QuICC/TestSuite/SparseSM/TesterBase.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace TestSuite {

namespace SparseSM {

namespace Bessel {

template <typename TOp> class Tester : public SparseSM::TesterBase<TOp>
{
public:
   /// Typedef for parameter type
   typedef typename SparseSM::TesterBase<TOp>::ParameterType ParameterType;

   /**
    * @brief Constructor
    *
    * @param fname      Filename
    * @param wtype      Bessel type
    * @param keepData   Write computed data
    */
   Tester(const std::string& fname, const std::string& wtype,
      const bool keepData);

   /*
    * @brief Destructor
    */
   virtual ~Tester() = default;

protected:
   /// Typedef ContentType from base
   typedef typename SparseSM::TesterBase<TOp>::ContentType ContentType;

   /**
    * @brief Build filename extension with resolution information
    *
    * @param param Run parameters
    */
   virtual std::string resname(const ParameterType& param) const override;

   /**
    * @brief Test operator
    *
    * @param mat     Sparse matrix output
    * @param param   Run parameters
    * @param type    Type of test
    */
   virtual void buildOperator(SparseMatrix& mat, const ParameterType& param,
      const TestType type) const override;

   /**
    * @brief Test operator
    */
   virtual void buildOperator(Matrix& mat, const ParameterType& param,
      const TestType type) const override;

   /**
    * @brief Format the parameters
    *
    * @param param   Run parameters
    */
   virtual std::string formatParameter(
      const ParameterType& param) const override;

private:
   /**
    * @brief Append specific path
    *
    * @param subdir Sub-directory to append
    */
   void appendPath(const std::string& subdir);

   /**
    * @brief Read meta data
    *
    * @param param   Run parameters
    * @param type    Type of test
    */
   Array readMeta(const ParameterType& param, const TestType type) const;
};

template <typename TOp>
Tester<TOp>::Tester(const std::string& fname, const std::string& btype,
   const bool keepData) :
    SparseSM::TesterBase<TOp>(fname, keepData)
{
   this->appendPath(btype);
}

template <typename TOp> void Tester<TOp>::appendPath(const std::string& subdir)
{
   this->mPath += "Bessel/" + subdir + "/";
}

template <typename TOp>
void Tester<TOp>::buildOperator(SparseMatrix& mat, const ParameterType& param,
   const TestType type) const
{
   Array meta = this->readMeta(param, type);

   if constexpr (std::is_base_of_v<QuICC::SparseSM::IBesselOperator, TOp>)
   {
      int rows = static_cast<int>(meta(0));
      int cols = static_cast<int>(meta(1));
      auto btype = static_cast<QuICC::SparseSM::Bessel::BesselKind>(
         static_cast<int>(meta(2)));
      auto l = meta(3);

      TOp op(rows, cols, btype, l);
      mat = op.mat();
   }
}

template <typename TOp>
void Tester<TOp>::buildOperator(Matrix& mat, const ParameterType& param,
   const TestType type) const
{
   Array meta = this->readMeta(param, type);

   if constexpr (std::is_base_of_v<QuICC::SparseSM::IBesselOperator, TOp>)
   {
      int rows = static_cast<int>(meta(0));
      int cols = static_cast<int>(meta(1));
      auto btype = static_cast<QuICC::SparseSM::Bessel::BesselKind>(
         static_cast<int>(meta(2)));
      auto l = meta(3);

      TOp op(rows, cols, btype, l);
      unsigned int KL, KU;
      mat = op.banded(KL, KU);
   }
}

template <typename TOp>
std::string Tester<TOp>::resname(const ParameterType& param) const
{
   auto id = param.at(0);

   std::stringstream ss;
   ss.precision(10);
   ss << "_id" << id;

   return ss.str();
}

template <typename TOp>
std::string Tester<TOp>::formatParameter(const ParameterType& param) const
{
   auto id = param.at(0);

   std::stringstream ss;
   ss << "id: " << id;

   return ss.str();
}

template <typename TOp>
Array Tester<TOp>::readMeta(const ParameterType& param,
   const TestType type) const
{
   // Read metadata
   Array meta(0);
   std::string fullname =
      this->makeFilename(param, this->refRoot(), type, ContentType::META);
   readList(meta, fullname);

   return meta;
}

} // namespace Bessel
} // namespace SparseSM
} // namespace TestSuite
} // namespace QuICC

#endif // QUICC_TESTSUITE_SPARSESM_BESSEL_TESTER_HPP

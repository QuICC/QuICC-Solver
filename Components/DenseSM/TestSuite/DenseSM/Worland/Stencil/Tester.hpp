/**
 * @file Tester.hpp
 * @brief Tester for Worland transforms
 */

#ifndef QUICC_TESTSUITE_DENSESM_WORLAND_STENCIL_TESTER_HPP
#define QUICC_TESTSUITE_DENSESM_WORLAND_STENCIL_TESTER_HPP

// System includes
//
#include <catch2/catch.hpp>
#include <sstream>
#include <string>

// Project includes
//
#include "DenseSM/Worland/IWorlandOperator.hpp"
#include "DenseSM/Worland/Stencil/IStencilOperator.hpp"
#include "TestSuite/DenseSM/TesterBase.hpp"

namespace dsm = ::QuICC::DenseSM::Worland::Stencil;

namespace QuICC {

namespace TestSuite {

namespace DenseSM {

namespace Worland {

namespace Stencil {

template <typename TOp> class Tester : public DenseSM::TesterBase<TOp>
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
   virtual void readFile(Matrix& data, const ParameterType& param,
      const TestType type, const ContentType ctype) const override;

   /**
    * @brief Read complex data from file
    */
   virtual void readFile(MatrixZ& data, const ParameterType& param,
      const TestType type, const ContentType ctype) const override;

   /**
    * @brief Test operator
    */
   virtual Matrix applyOperator(const ParameterType& param,
      const TestType type) const override;

   /**
    * @brief Format the parameters
    */
   virtual std::string formatParameter(
      const ParameterType& param) const override;

private:
   /**
    * @brief Append specific path
    */
   void appendPath();
};

template <typename TOp>
Tester<TOp>::Tester(const std::string& fname, const bool keepData) :
    DenseSM::TesterBase<TOp>(fname, keepData, false)
{
   this->appendPath();
}

template <typename TOp> void Tester<TOp>::appendPath()
{
   this->mPath += "Worland/";
}

template <typename TOp>
void Tester<TOp>::readFile(Matrix& data, const ParameterType& param,
   const TestType type, const ContentType ctype) const
{
   TesterBase<TOp>::basicReadFile(data, param, type, ctype);
}

template <typename TOp>
void Tester<TOp>::readFile(MatrixZ& data, const ParameterType& param,
   const TestType type, const ContentType ctype) const
{
   TesterBase<TOp>::basicReadFile(data, param, type, ctype);
}

template <typename TOp>
Matrix Tester<TOp>::applyOperator(const ParameterType& param,
   const TestType type) const
{
   Matrix outData;

   if constexpr (std::is_base_of_v<dsm::IStencilOperator, TOp>)
   {
      Array meta(0);
      std::string fullname =
         this->makeFilename(param, this->refRoot(), type, ContentType::META);
      readList(meta, fullname);
      if (meta.size() != 7)
      {
         throw std::logic_error("Test meta data is wrong");
      }
      std::cerr << meta.transpose() << std::endl;

      int nNr = meta(0);
      int nNc = meta(1);
      auto a = static_cast<MHDFloat>(meta(2));
      auto b = static_cast<MHDFloat>(meta(3));
      auto l = static_cast<int>(meta(4));
      auto nId = static_cast<int>(meta(5));
      auto c = static_cast<int>(meta(6));

      TOp op(nNr, nNc, a, b, l, nId, c);

      outData = op.mat();
   }
   else
   {
      throw std::logic_error("Not implemented");
   }

   return outData;
}

template <typename TOp>
std::string Tester<TOp>::resname(const ParameterType& param) const
{
   auto id = param.at(0);

   std::stringstream ss;
   ss.precision(10);
   ss << "_id" << id;

   // Distributed data meta file
   if (param.size() == 3)
   {
      int np = param.at(1);
      ss << "_np" << np;
      int r = param.at(2);
      ss << "_r" << r;
      ss << "_stage0";
   }

   return ss.str();
}

template <typename TOp>
std::string Tester<TOp>::formatParameter(const ParameterType& param) const
{
   auto id = param.at(0);

   std::stringstream ss;
   ss << "id: " << id;

   if (param.size() == 3)
   {
      int np = param.at(1);
      ss << ", np: " << np;
      int r = param.at(2);
      ss << ", r: " << r;
   }

   return ss.str();
}

} // namespace Stencil
} // namespace Worland
} // namespace DenseSM
} // namespace TestSuite
} // namespace QuICC

#endif // QUICC_TESTSUITE_DENSESM_WORLAND_STENCIL_TESTER_HPP

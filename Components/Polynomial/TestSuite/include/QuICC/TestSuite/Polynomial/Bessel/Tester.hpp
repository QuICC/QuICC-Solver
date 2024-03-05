/**
 * @file Tester.hpp
 * @brief Tester for Bessel polynomials
 */

#ifndef QUICC_TESTSUITE_POLYNOMIAL_BESSEL_TESTER_HPP
#define QUICC_TESTSUITE_POLYNOMIAL_BESSEL_TESTER_HPP

// System includes
//
#include <catch2/catch.hpp>
#include <set>
#include <string>

// Project includes
//
#include "QuICC/Polynomial/Quadrature/WorlandSphEnergyRule.hpp"
#include "QuICC/TestSuite/Polynomial/TesterBase.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Typedefs.hpp"

namespace poly = ::QuICC::Polynomial;

namespace QuICC {

namespace TestSuite {

namespace Polynomial {

namespace Bessel {

template <typename TOp> class Tester : public Polynomial::TesterBase<TOp>
{
public:
   /// Typedef for parameter type
   typedef typename Polynomial::TesterBase<TOp>::ParameterType ParameterType;

   /**
    * @brief Constructor
    *
    * @param fname      file name
    * @param keepData   Dump data?
    */
   Tester(const std::string& fname, const bool keepData);

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
    * @param param   Parameters of run
    */
   virtual std::string resname(const int specN, const int physN,
      const ParameterType& param) const override;

   /**
    * @brief Test matrix operator
    *
    * @param specN   Number of spectral modes
    * @param physN   Number of physical grid points
    * @param param   Parameters of run
    * @param type    Type of test
    */
   virtual Matrix buildOperator(const int specN, const int physN,
      const ParameterType& param, const TestType type) const override;

   /**
    * @brief Format the parameters
    *
    * @param param   Parameters of run
    */
   virtual std::string formatParameter(
      const ParameterType& param) const override;

private:
   /**
    * @brief Append specific path
    */
   void appendPath();

   /**
    * @brief Test quadrature matrix operator
    *
    * @param specN      Number of spectral modes
    * @param physN      Number of physical grid points
    * @param param      Parameters of run
    * @param isWeighted Are quadrature weights included?
    */
   Matrix buildMatrix(const int specN, const int physN,
      const ParameterType& param, const bool isWeighted) const;

   /**
    * @brief Test on-the-fly inner product
    *
    * @param specN      Number of spectral modes
    * @param physN      Number of physical grid points
    * @param param      Parameters of run
    */
   Matrix buildInner(const int specN, const int physN,
      const ParameterType& param) const;

   /**
    * @brief Test on-the-fly outer product
    *
    * @param specN      Number of spectral modes
    * @param physN      Number of physical grid points
    * @param param      Parameters of run
    */
   Matrix buildOuter(const int specN, const int physN,
      const ParameterType& param) const;

   /**
    * @brief Test on-the-fly reduction
    *
    * @param specN      Number of spectral modes
    * @param physN      Number of physical grid points
    * @param param      Parameters of run
    */
   Matrix buildReduce(const int specN, const int physN,
      const ParameterType& param) const;
};

template <typename TOp>
Tester<TOp>::Tester(const std::string& fname, const bool keepData) :
    Polynomial::TesterBase<TOp>(fname, keepData)
{
   this->appendPath();
}

template <typename TOp> void Tester<TOp>::appendPath()
{
   this->mPath += "Bessel/";
}

template <typename TOp>
Matrix Tester<TOp>::buildOperator(const int specN, const int physN,
   const ParameterType& param, const TestType type) const
{
   Matrix outData;
   switch (type)
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

template <typename TOp>
Matrix Tester<TOp>::buildMatrix(const int specN, const int physN,
   const ParameterType& param, const bool isWeighted) const
{
   // Get harmonic degree from parameters
   int l = static_cast<int>(param.at(0));

   // Create quadrature
   Internal::Array igrid;
   Internal::Array iweights;
   poly::Quadrature::WorlandSphEnergyRule quad;
   quad.computeQuadrature(igrid, iweights, physN);

   TOp op;
   Matrix outData(physN, specN + 1);

   if (isWeighted)
   {
      op.template compute<MHDFloat>(outData, specN + 1, l, igrid, iweights);
   }
   else
   {
      op.template compute<MHDFloat>(outData, specN + 1, l, igrid,
         Internal::Array());
   }

   return outData;
}

template <typename TOp>
Matrix Tester<TOp>::buildOuter(const int specN, const int physN,
   const ParameterType& param) const
{
   // Get harmonic degree from parameters
   int l = static_cast<int>(param.at(0));

   // Create quadrature
   Internal::Array igrid;
   Internal::Array iweights;
   poly::Quadrature::WorlandSphEnergyRule quad;
   quad.computeQuadrature(igrid, iweights, physN);

   TOp op;
   Matrix outData(physN, specN + 1);
   Matrix inData = Matrix::Identity(specN + 1, specN + 1);

   op.template compute<MHDFloat>(outData, specN + 1, l, igrid,
      Internal::Array());

   return outData;
}

template <typename TOp>
Matrix Tester<TOp>::buildInner(const int specN, const int physN,
   const ParameterType& param) const
{
   // Get harmonic degree from parameters
   int l = static_cast<int>(param.at(0));

   // Create quadrature
   Internal::Array igrid;
   Internal::Array iweights;
   poly::Quadrature::WorlandSphEnergyRule quad;
   quad.computeQuadrature(igrid, iweights, physN);

   TOp op;
   Matrix outData(specN + 1, specN + 1);
   Matrix inData = Matrix::Identity(physN, specN + 1);
   op.template compute<MHDFloat>(inData, specN + 1, l, igrid, iweights);

   op.template compute<MHDFloat>(outData, specN + 1, l, igrid,
      Internal::Array());

   Matrix outDataT = outData.transpose();

   return outDataT;
}

template <typename TOp>
Matrix Tester<TOp>::buildReduce(const int specN, const int physN,
   const ParameterType& param) const
{
   // Get harmonic degree from parameters
   int l = static_cast<int>(param.at(0));

   // Create quadrature
   Internal::Array igrid;
   Internal::Array iweights;
   poly::Quadrature::WorlandSphEnergyRule quad;
   quad.computeQuadrature(igrid, iweights, physN);

   TOp op;
   Matrix outData(specN + 1, 1);

   op.template compute<MHDFloat>(outData, specN + 1, l, igrid,
      Internal::Array());

   return outData;
}

template <typename TOp>
std::string Tester<TOp>::resname(const int specN, const int physN,
   const ParameterType& param) const
{
   // Get harmonic degree from parameters
   int l = static_cast<int>(param.at(0));

   std::string s = "_l" + std::to_string(l) + "_n" + std::to_string(specN) +
                   "_g" + std::to_string(physN);
   return s;
}

template <typename TOp>
std::string Tester<TOp>::formatParameter(const ParameterType& param) const
{
   // Get harmonic degree from parameters
   int l = static_cast<int>(param.at(0));

   std::stringstream ss;
   ss << "l: " << l;

   return ss.str();
}

} // namespace Bessel
} // namespace Polynomial
} // namespace TestSuite
} // namespace QuICC

#endif // QUICC_TESTSUITE_POLYNOMIAL_BESSEL_TESTER_HPP

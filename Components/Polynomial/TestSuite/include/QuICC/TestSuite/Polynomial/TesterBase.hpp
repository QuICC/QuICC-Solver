/**
 * @file common.hpp
 * @brief Functions common to all Associated Legendre tests
 */

#ifndef QUICC_TESTSUITE_POLYNOMIAL_TESTERBASE_HPP
#define QUICC_TESTSUITE_POLYNOMIAL_TESTERBASE_HPP


// Configuration includes
//

// System includes
//
#include <catch2/catch.hpp>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Timers/TimerMacro.h"
#include "QuICC/TestSuite/Polynomial/TestType.hpp"
#include "TestSuite/Io.hpp"

namespace QuICC {

namespace TestSuite {

namespace Polynomial {

   template <typename TOp> class  TesterBase
   {
      public:
         /// Typedef for the parameter type
         typedef std::vector<MHDFloat> ParameterType;

         /// Typedef for the Error type
         typedef std::pair<bool,MHDFloat> ErrorType;

         /**
          * @brief Constructor
          */
         TesterBase(const std::string& fname, const bool keepData);

         /**
          * @brief Destructor
          */
         virtual ~TesterBase() = default;

         /**
          * @brief Validate implementation against reference single parameter input. Setup is read from meta data
          */
         void validate(const ParameterType& param, const TestType type) const;

         /**
          * @brief Validate implementation against reference
          */
         void validate(const int specN, const int physN, const ParameterType& param, const TestType type) const;

         /**
          * @brief Set max acceptable ulp
          */
         void setUlp(const unsigned int ulp);

      protected:
         /**
          * @brief Build filename extension with resolution information
          */
         virtual std::string resname(const int specN, const int physN, const ParameterType& param) const = 0;

         /**
          * @brief Build operator
          */
         virtual Matrix buildOperator(const int specN, const int physN, const ParameterType& param, const TestType test) const = 0;

         /**
          * @brief Format the parameters
          */
         virtual std::string formatParameter(const ParameterType& param) const = 0;

         /**
          * @brief Check if reference is special value
          */
         virtual bool isSpecial(const int i, const int j, const MHDFloat data, const MHDFloat ref) const;

         /**
          * @brief Check special values
          */
         virtual ErrorType checkSpecial(const int i, const int j, const MHDFloat data, const Matrix& ref) const;

         /**
          * @brief Check normalized values
          */
         virtual ErrorType checkNormal(const MHDFloat data, const MHDFloat ref, const MHDFloat scale = 1.0) const;

         /**
          * @brief Check denormalized values
          */
         virtual ErrorType checkSubnormal(const MHDFloat data, const MHDFloat ref, const MHDFloat scale = 1.0) const;

         /**
          * @brief Get data root
          */
         std::string dataRoot() const;

         /**
          * @brief Get reference data root
          */
         std::string refRoot() const ;

         /**
          * @brief Make filename
          */
         virtual std::string makeFilename(const int specN, const int physN, const ParameterType& param, const std::string& root, const TestType type) const;

         /**
          * @brief Test operator
          */
         Matrix testOperator(const int specN, const int physN, const ParameterType& param, const TestType type) const;

         /**
          * @brief Compute error against reference data
          */
         void computeError(const Matrix& outData, const Matrix& refData) const;

         /**
          * @brief Tolerance
          */
         MHDFloat tolerance() const;

         /**
          * @brief Save data to file
          */
         bool mKeepData;

         /**
          * @brief Max Units in the last place allowed
          */
         unsigned int mMaxUlp;

         /**
          * @brief Datatype epsilon
          */
         MHDFloat mEpsilon;

         /**
          * @brief Reference scale for the problem
          */
         MHDFloat mScale = -1.0;

         /**
          * @brief Base filename
          */
         std::string mBasename;

         /**
          * @brief Directory path
          */
         std::string mPath;
   };

   template <typename TOp> TesterBase<TOp>::TesterBase(const std::string& fname, const bool keepData)
      : mKeepData(keepData), mMaxUlp(11), mEpsilon(std::numeric_limits<MHDFloat>::epsilon()), mBasename(fname), mPath("Polynomial/")
   {
   }

   template <typename TOp> void TesterBase<TOp>::validate(const ParameterType& param, const TestType type) const
   {
      assert(param.size() == 1);

      // Read metadata
      std::string filename = this->mBasename;

      if(type == TestType::WEIGHTED_MATRIX)
      {
            filename = "weighted_" +  filename;
      }
      std::string ext = "_id" + std::to_string(static_cast<int>(param.at(0))) + "_meta";
      std::string fullname = this->refRoot() + filename.insert(filename.length()-4, ext);
      Array meta;
      readList(meta, fullname);

      assert(meta.size() == 3);
      int specN = static_cast<int>(meta(0));
      int physN = static_cast<int>(meta(1));

      ParameterType runParam;
      runParam.push_back(meta(2));

      this->validate(specN, physN, runParam, type);
   }

   template <typename TOp> void TesterBase<TOp>::validate(const int specN, const int physN, const ParameterType& param, const TestType type) const
   {
      Matrix outData;
      std::string filename = this->mBasename;

      // Compute operator
      std::string infoType = "unknown";
      switch(type)
      {
         case TestType::MATRIX:
            infoType = "matrix";
            break;
         case TestType::WEIGHTED_MATRIX:
            infoType = "weighted matrix";
            filename = "weighted_" +  filename;
            break;
         case TestType::OTF_INNER:
            infoType = "on-the-fly inner product";
            break;
         case TestType::OTF_OUTER:
            infoType = "on-the-fly outer product";
            break;
         case TestType::OTF_REDUCE:
            infoType = "on-the-fly reduction";
            break;
         case TestType::QUADRATURE:
            infoType = "quadrature rule";
            break;
      }
      outData = this->testOperator(specN, physN, param, type);

      // Read reference data
      Matrix refData(outData.rows(), outData.cols());
      std::string fullname = this->makeFilename(specN, physN, param, this->refRoot(), type);
      readData(refData, fullname);

      // Fail if reference not found
      if(refData.size() == 0)
      {
         INFO( "reference data was not found" );
         CHECK( false );
      }

      // Compare data to reference
      INFO( "type: " + infoType );
      INFO( this->formatParameter(param));
      this->computeError(outData, refData);
   }

   template <typename TOp> void TesterBase<TOp>::setUlp(const unsigned int ulp)
   {
      this->mMaxUlp = ulp;
   }

   template <typename TOp> std::string TesterBase<TOp>::makeFilename(const int specN, const int physN, const ParameterType& param, const std::string& root, const TestType type) const
   {
      std::string pre = "";
      std::string post = "";

      switch(type)
      {
         case TestType::WEIGHTED_MATRIX:
            pre = "weighted_";
            break;
         case TestType::OTF_INNER:
            post = "_Inner";
            break;
         case TestType::OTF_OUTER:
            post = "_Outer";
            break;
         case TestType::OTF_REDUCE:
            post = "_Reduce";
            break;
         default:
            break;
      }

      std::string ext = resname(specN, physN, param);
      std::string filename = pre + this->mBasename;
      std::string fullname = root + filename.insert(filename.length()-4, ext+post);

      return fullname;
   }

   template <typename TOp> Matrix TesterBase<TOp>::testOperator(const int specN, const int physN, const ParameterType& param, const TestType type) const
   {
      Matrix outData = this->buildOperator(specN, physN, param, type);

      if(this->mKeepData)
      {
         std::string fullname = this->makeFilename(specN, physN, param, this->dataRoot(), type);
         writeData(fullname, outData);
      }

      return outData;
   }

   template <typename TOp> void TesterBase<TOp>::computeError(const Matrix& outData, const Matrix& refData) const
   {
      auto refScale = std::max(std::abs(refData.maxCoeff()),
         std::abs(refData.minCoeff()));

      // Compute error
      //
      for(int j = 0; j < refData.cols(); j++)
      {
         INFO( "column: " << j );
         for(int i = 0; i < refData.rows(); i++)
         {
            auto ref = std::abs(refData(i,j));

            INFO( "position: " << i << " / " << refData.rows()-1 );
            INFO( "refData: " << std::scientific << std::setprecision(16) << refData(i,j) );

            // check special values
            if(this->isSpecial(i, j, outData(i,j), ref))
            {
               auto err = this->checkSpecial(i, j, outData(i,j), refData);

               // check error
               {
                  INFO( "checked special value" );
                  INFO( "outData: " << std::scientific << std::setprecision(16) << outData(i,j) );
                  INFO( "max ulp: " << this->mMaxUlp);
                  INFO( "measured ulp: " << err.second);
                  CHECK(err.first);
               }
            }
            // check normalized values
            else if(ref >= std::numeric_limits<MHDFloat>::min())
            {
               auto err = this->checkNormal(outData(i,j), refData(i,j), refScale);

               // check error
               {
                  INFO( "checked normal value" );
                  INFO( "outData: " << std::scientific << std::setprecision(16) << outData(i,j) );
                  INFO( "max ulp: " << this->mMaxUlp);
                  INFO( "measured ulp: " << err.second);
                  CHECK(err.first);
               }
            }
            // Check subnormal values
            else
            {
               auto err = this->checkSubnormal(outData(i,j), refData(i,j), refScale);

               // check error
               {
                  INFO( "checked subnormal value" );
                  INFO( "outData: " << std::scientific << std::setprecision(16) << outData(i,j) );
                  INFO( "max ulp: " << this->mMaxUlp);
                  INFO( "measured ulp: " << err.second);
                  CHECK(err.first);
               }
            }
         }
      }
   }

   template <typename TOp> MHDFloat TesterBase<TOp>::tolerance() const
   {
      return this->mMaxUlp*this->mEpsilon;
   }

   template <typename TOp> bool TesterBase<TOp>::isSpecial(const int i, const int j, const MHDFloat data, const MHDFloat ref) const
   {
      return false;
   }

   template <typename TOp> typename TesterBase<TOp>::ErrorType TesterBase<TOp>::checkSpecial(const int i, const int j, const MHDFloat data, const Matrix& ref) const
   {
      bool isEqual = false;

      auto ulp = 0;

      return std::make_pair(isEqual, ulp);
   }

   template <typename TOp> typename TesterBase<TOp>::ErrorType TesterBase<TOp>::checkNormal(const MHDFloat data, const MHDFloat ref, const MHDFloat scale) const
   {
      bool isEqual = false;

      auto diff = std::abs(data-ref);
      auto tol = this->tolerance();
      auto refMod = scale;
      if (refMod <= 0.0)
      {
         refMod = std::abs(ref);
      }

      if(diff < tol)
      {
         isEqual = diff < (tol * refMod);
      }
      else
      {
         isEqual = (diff / refMod ) < tol;
      }

      auto ulp = diff / (refMod * this->mEpsilon);

      return std::make_pair(isEqual, ulp);
   }

   template <typename TOp> typename TesterBase<TOp>::ErrorType TesterBase<TOp>::checkSubnormal(const MHDFloat data, const MHDFloat ref, const MHDFloat scale) const
   {
      bool isEqual = false;
      double ulp = 0;

      // check that subnormal was flushed to zero or subnormal
      if (data == 0.0 || data < std::numeric_limits<MHDFloat>::min())
      {
         isEqual = true;
         ulp = 1;
      }
      else
      {
         auto err = checkNormal(data, ref, scale);
         isEqual = err.first;
         ulp = err.second;
      }

      return std::make_pair(isEqual, ulp);
   }

   template <typename TOp> std::string TesterBase<TOp>::dataRoot() const
   {
      std::string p = "_data/" + this->mPath;
      return p;
   }

   template <typename TOp> std::string TesterBase<TOp>::refRoot() const
   {
      std::string p = "_refdata/" + this->mPath;
      return p;
   }

}
}
}

#endif //QUICC_TESTSUITE_POLYNOMIAL_TESTERBASE_HPP

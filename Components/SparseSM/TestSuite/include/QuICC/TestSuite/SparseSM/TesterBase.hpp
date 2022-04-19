/**
 * @file common.hpp
 * @brief Functions common to all SparseSM tests
 */

#ifndef QUICC_TESTSUITE_SPARSESM_TESTERBASE_HPP
#define QUICC_TESTSUITE_SPARSESM_TESTERBASE_HPP


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
#include "QuICC/Typedefs.hpp"
#include "Timers/TimerMacro.h"
#include "QuICC/TestSuite/SparseSM/TestType.hpp"
#include "QuICC/TestSuite/SparseSM/Io.hpp"

namespace QuICC {

namespace TestSuite {

namespace SparseSM {

   template <typename TOp> class TesterBase
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
          * @brief Validate implementation against reference
          */
         void validate(const ParameterType& param, const TestType type) const;

         /**
          * @brief Set max acceptable ulp
          */
         void setUlp(const unsigned int ulp);

      protected:
         /**
          * @brief Enum of data content type
          */
         enum ContentType
         {
            META = 0,
            REFERENCE
         };

         /**
          * @brief Build filename extension with resolution information
          */
         virtual std::string resname(const ParameterType& param) const = 0;

         /**
          * @brief Validate implementation against reference
          */
         template <typename TData> void validateData(const ParameterType& param, const TestType type) const;

         /**
          * @brief Build sparse operator
          */
         virtual void buildOperator(SparseMatrix& mat, const ParameterType& param, const TestType type) const = 0;

         /**
          * @brief Build banded operator
          */
         virtual void buildOperator(Matrix& mat, const ParameterType& param, const TestType type) const = 0;

         /**
          * @brief Format the parameters
          */
         virtual std::string formatParameter(const ParameterType& param) const = 0;

         /**
          * @brief Read data from file
          */
         virtual void readFile(Matrix& data, const ParameterType& param, const TestType type, const ContentType ctype) const;

         virtual void readFile(SparseMatrix& data, const ParameterType& param, const TestType type, const ContentType ctype) const;

         template <typename TData> void basicReadFile(TData& data, const ParameterType& param, const TestType type, const ContentType ctype) const;

         /**
          * @brief Check if reference is special value
          */
         virtual bool isSpecial(const int i, const int j, const MHDFloat data, const MHDFloat ref) const;

         /**
          * @brief Check special values
          */
         virtual ErrorType checkSpecial(const int i, const int j, const MHDFloat data, const Matrix& ref) const;

         /**
          * @brief Check special values
          */
         virtual ErrorType checkSpecial(const int i, const int j, const MHDFloat data, const SparseMatrix& ref) const;

         /**
          * @brief Check if reference is exact zero
          */
         virtual bool isZero(const MHDFloat ref) const;

         /**
          * @brief Check exact zero values
          */
         virtual ErrorType checkZero(const MHDFloat data, const MHDFloat ref) const;

         /**
          * @brief Check normalized values
          */
         virtual ErrorType checkNormal(const MHDFloat data, const MHDFloat ref) const;

         /**
          * @brief Check normalized values against reference scale
          */
         virtual ErrorType checkNormal(const MHDFloat data, const MHDFloat ref, const MHDFloat refScale) const;

         /**
          * @brief Check denormalized values
          */
         virtual ErrorType checkSubnormal(const MHDFloat data, const MHDFloat ref) const;

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
         virtual std::string makeFilename(const ParameterType& param, const std::string& root, const TestType type, const ContentType ctype) const;

         /**
          * @brief Test operator
          */
         template <typename TData> void testOperator(TData&, const ParameterType& param, const TestType type) const;

         /**
          * @brief Compute error against reference data
          */
         void computeError(const Matrix& outData, const Matrix& refData) const;

         /**
          * @brief Compute error against reference data
          */
         void computeError(const SparseMatrix& outData, const SparseMatrix& refData) const;

         /**
          * @brief Tolerance
          */
         MHDFloat tolerance() const;

         /**
          * @brief Save data to file
          */
         bool mKeepData;

         /**
          * @brief Datatype epsilon
          */
         unsigned int mMaxUlp;

         /**
          * @brief Datatype epsilon
          */
         MHDFloat mEpsilon;

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
      : mKeepData(keepData), mMaxUlp(11), mEpsilon(std::numeric_limits<MHDFloat>::epsilon()), mBasename(fname), mPath("SparseSM/")
   {
   }

   template <typename TOp> void TesterBase<TOp>::validate(const ParameterType& param, const TestType type) const
   {
      // Compute operator
      switch(type)
      {
         case TestType::SPARSE:
            this->validateData<SparseMatrix>(param, type);
            break;
         case TestType::BANDED:
            this->validateData<Matrix>(param, type);
            break;
      }
   }

   template <typename TOp> template <typename TData> void TesterBase<TOp>::validateData(const ParameterType& param, const TestType type) const
   {
      TData outData;
      std::string filename = this->mBasename;

      std::string infoType = "unknown";
      switch(type)
      {
         case TestType::SPARSE:
            infoType = "sparse";
            break;
         case TestType::BANDED:
            infoType = "BLAS banded";
            break;
      }

      // Compute operator
      this->testOperator(outData, param, type);

      // Read reference data
      TData refData(outData.rows(), outData.cols());
      this->readFile(refData, param, type, ContentType::REFERENCE);

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

   template <typename TOp> void TesterBase<TOp>::readFile(Matrix& data, const ParameterType& param, const TestType type, const ContentType ctype) const
   {
      this->basicReadFile(data, param, type, ctype);
   }

   template <typename TOp> void TesterBase<TOp>::readFile(SparseMatrix& data, const ParameterType& param, const TestType type, const ContentType ctype) const
   {
      this->basicReadFile(data, param, type, ctype);
   }

   template <typename TOp> template<typename TData> void TesterBase<TOp>::basicReadFile(TData& data, const ParameterType& param, const TestType type, const ContentType ctype) const
   {
      // Read reference data
      std::string fullname = this->makeFilename(param, this->refRoot(), type, ctype);
      readData(data, fullname);
   }

   template <typename TOp> void TesterBase<TOp>::setUlp(const unsigned int ulp)
   {
      this->mMaxUlp = ulp;
   }

   template <typename TOp> std::string TesterBase<TOp>::makeFilename(const ParameterType& param, const std::string& root, const TestType type, const ContentType ctype) const
   {
      std::string pre = "";
      std::string post = "";
      std::string sub = "";

      switch(type)
      {
         case TestType::SPARSE:
            sub = "Sparse/";
            break;
         case TestType::BANDED:
            sub = "Banded/";
            break;
      }

      switch(ctype)
      {
         case ContentType::META:
            post = "_meta";
            break;
         case ContentType::REFERENCE:
            post = "_ref";
            break;
      }

      std::string ext = resname(param);
      std::string filename = pre + this->mBasename;
      std::string fullname = root + sub + filename.insert(filename.length()-4, ext+post);

      return fullname;
   }

   template <typename TOp> template <typename TData> void TesterBase<TOp>::testOperator(TData& outData, const ParameterType& param, const TestType type) const
   {
      this->buildOperator(outData, param, type);

      if(this->mKeepData)
      {
         std::string fullname = this->makeFilename(param, this->dataRoot(), type, ContentType::REFERENCE);
         writeData(fullname, outData);
      }
   }

   template <typename TOp> void TesterBase<TOp>::computeError(const Matrix& outData, const Matrix& refData) const
   {
      MHDFloat globalUlp = 0;

      // Compute error
      //
      for(int j = 0; j < refData.cols(); j++)
      {
         INFO( "n: " << j );
         for(int i = 0; i < refData.rows(); i++)
         {
            ErrorType err;
            auto ref = std::abs(refData(i,j));

            INFO( "position: " << i << " / " << refData.rows()-1 );
            INFO( "refData: " << std::scientific << std::setprecision(16) << refData(i,j) );

            // check special values
            if(this->isSpecial(i, j, outData(i,j), ref))
            {
               err = this->checkSpecial(i, j, outData(i,j), refData);

               // check error
               {
                  INFO( "checked special value" );
                  INFO( "outData: " << std::scientific << std::setprecision(16) << outData(i,j) );
                  INFO( "max ulp: " << this->mMaxUlp);
                  INFO( "measured ulp: " << err.second);
                  CHECK(err.first);
               }
            }
            // check exact zero
            else if(this->isZero(ref))
            {
               err = this->checkZero(outData(i,j), refData(i,j));

               // check error
               {
                  INFO( "checked exact zero" );
                  INFO( "outData: " << std::scientific << std::setprecision(16) << outData(i,j) );
                  INFO( "max ulp: " << this->mMaxUlp);
                  INFO( "measured ulp: " << err.second);
                  CHECK(err.first);
               }
            }
            // check normalized values
            else if(ref >= std::numeric_limits<MHDFloat>::min())
            {
               err = this->checkNormal(outData(i,j), refData(i,j));

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
               err = this->checkSubnormal(outData(i,j), refData(i,j));

               // check error
               {
                  INFO( "checked subnormal value" );
                  INFO( "outData: " << std::scientific << std::setprecision(16) << outData(i,j) );
                  INFO( "max ulp: " << this->mMaxUlp);
                  INFO( "measured ulp: " << err.second);
                  CHECK(err.first);
               }
            }
            globalUlp = std::max(globalUlp, err.second);
         }
      }

      INFO( "Global max ulp: " << globalUlp );
      CHECK( globalUlp < this->mMaxUlp );

   }

   template <typename TOp> void TesterBase<TOp>::computeError(const SparseMatrix& outData, const SparseMatrix& refData) const
   {
      MHDFloat globalUlp = 0;

      // Compute error
      //
      for(int j = 0; j < refData.outerSize(); ++j)
      {
         for(SparseMatrix::InnerIterator it(refData, j); it; ++it)
         {
            int i_ = it.row();
            int j_ = it.col();
            INFO( "n: " << j_ );
            ErrorType err;
            auto ref = std::abs(it.value());

            INFO( "position: " << i_ << " / " << refData.rows()-1 );
            INFO( "refData: " << std::scientific << std::setprecision(16) << it.value() );

            auto outVal = outData.coeff(it.row(),it.col());
            // check special values
            if(this->isSpecial(i_, j_, outVal, ref))
            {
               err = this->checkSpecial(i_, j_, outVal, refData);

               // check error
               {
                  INFO( "checked special value" );
                  INFO( "outData: " << std::scientific << std::setprecision(16) << outVal );
                  INFO( "max ulp: " << this->mMaxUlp);
                  INFO( "measured ulp: " << err.second);
                  CHECK(err.first);
               }
            }
            // check exact zero
            else if(this->isZero(ref))
            {
               err = this->checkZero(outVal, it.value());

               // check error
               {
                  INFO( "checked exact zero" );
                  INFO( "outData: " << std::scientific << std::setprecision(16) << outVal );
                  INFO( "max ulp: " << this->mMaxUlp);
                  INFO( "measured ulp: " << err.second);
                  CHECK(err.first);
               }
            }
            // check normalized values
            else if(ref >= std::numeric_limits<MHDFloat>::min())
            {
               err = this->checkNormal(outVal, it.value());

               // check error
               {
                  INFO( "checked normal value" );
                  INFO( "outData: " << std::scientific << std::setprecision(16) << outVal );
                  INFO( "max ulp: " << this->mMaxUlp);
                  INFO( "measured ulp: " << err.second);
                  CHECK(err.first);
               }
            }
            // Check subnormal values
            else
            {
               err = this->checkSubnormal(outVal, it.value());

               // check error
               {
                  INFO( "checked subnormal value" );
                  INFO( "outData: " << std::scientific << std::setprecision(16) << outVal );
                  INFO( "max ulp: " << this->mMaxUlp);
                  INFO( "measured ulp: " << err.second);
                  CHECK(err.first);
               }
            }
            globalUlp = std::max(globalUlp, err.second);
         }
      }

      INFO( "Global max ulp: " << globalUlp );
      CHECK( globalUlp < this->mMaxUlp );

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

   template <typename TOp> typename TesterBase<TOp>::ErrorType TesterBase<TOp>::checkSpecial(const int i, const int j, const MHDFloat data, const SparseMatrix& ref) const
   {
      bool isEqual = false;

      auto ulp = 0;

      return std::make_pair(isEqual, ulp);
   }

   template <typename TOp> bool TesterBase<TOp>::isZero(const MHDFloat ref) const
   {
      return (ref == 0);
   }

   template <typename TOp> typename TesterBase<TOp>::ErrorType TesterBase<TOp>::checkZero(const MHDFloat data, const MHDFloat ref) const
   {
      return this->checkNormal(data, ref, 1.0);
   }

   template <typename TOp> typename TesterBase<TOp>::ErrorType TesterBase<TOp>::checkNormal(const MHDFloat data, const MHDFloat ref, const MHDFloat refMod) const
   {
      bool isEqual = false;

      auto diff = std::abs(data-ref);
      auto tol = this->tolerance();

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

   template <typename TOp> typename TesterBase<TOp>::ErrorType TesterBase<TOp>::checkNormal(const MHDFloat data, const MHDFloat ref) const
   {
      auto refMod = std::abs(ref);
      return this->checkNormal(data, ref, refMod);
   }

   template <typename TOp> typename TesterBase<TOp>::ErrorType TesterBase<TOp>::checkSubnormal(const MHDFloat data, const MHDFloat ref) const
   {
      bool isEqual = false;
      double ulp = 0;
#if 0
      auto x = ref;
      auto y = data;
      if (std::abs(x) < std::abs(y))
      {
         std::swap(x, y);
      }
      auto diff = std::abs(x-y);

      // epsilon subnormal
      auto eps = std::numeric_limits<double>::denorm_min();

      if(diff <= eps)
      {
         isEqual = true;
      }

      auto ulp = diff/eps;
#else
      // check that subnormal was flushed to zero
      if (data == 0.0)
      {
         isEqual = true;
         ulp = 1;
      }
#endif
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

#endif //QUICC_TESTSUITE_SPARSESM_TESTERBASE_HPP

/**
 * @file Tester.hpp
 * @brief Tester for Fourier transforms
 */

#ifndef QUICC_TESTSUITE_TRANSFORM_FOURIER_TESTER_HPP
#define QUICC_TESTSUITE_TRANSFORM_FOURIER_TESTER_HPP


// Configuration includes
//

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
#include "QuICC/TestSuite/Transform/TesterBase.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Setup.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Setup.hpp"

template <typename T, typename... Args>
class has_transform
{
   template <typename C,
            typename = decltype( std::declval<C>().transform(std::declval<Args>()...) )>
               static std::true_type test(int);
   template <typename C>
      static std::false_type test(...);

   public:
   static constexpr bool value = decltype(test<T>(0))::value;
};

namespace transf = ::QuICC::Transform::Fft::Fourier;

namespace QuICC {

namespace TestSuite {

namespace Transform {

namespace Fourier {

   template <typename TOp, typename TOp2 = void> class Tester: public Transform::TesterBase<TOp>
   {
      public:
         /// Typedef for parameter type
         typedef typename Transform::TesterBase<TOp>::ParameterType ParameterType;

         /**
          * @brief Constructor
          */
         Tester(const std::string& fname, const bool keepData);

         /*
          * @brief Destructor
          */
         ~Tester() = default;

      protected:
         /// Typedef ContentType from base
         typedef typename Transform::TesterBase<TOp>::ContentType ContentType;
         /// Typedef for Forward input data
         typedef std::conditional_t<std::is_same_v<typename TOp::SetupType,transf::Complex::Setup>, MatrixZ, Matrix> FwdType;
         /// Typedef for Backward input data
         typedef MatrixZ BwdType;

         /**
          * @brief Build filename extension with resolution information
          */
         std::string resname(const ParameterType& param) const override;

         /**
          * @brief Read data from file
          */
         void readFile(Matrix& data, const ParameterType& param, const TestType type, const ContentType ctype) const override;

         void readFile(MatrixZ& data, const ParameterType& param, const TestType type, const ContentType ctype) const override;

         template <typename TData> void dbReadFile(TData& data, const ParameterType& param, const TestType type, const ContentType ctype) const;

         /**
          * @brief Test operator
          */
         Matrix applyOperator(const ParameterType& param, const TestType type) const override;

         /**
          * @brief Format the parameters
          */
         std::string formatParameter(const ParameterType& param) const override;

      private:
         /**
          * @brief Append specific path
          */
         void appendPath();

         /**
          * @brief Test projector
          */
         Matrix applyProjector(const ParameterType& param) const;

         /**
          * @brief Test integrator
          */
         Matrix applyIntegrator(const ParameterType& param) const;

         /**
          * @brief Test backward-forward loop
          */
         Matrix applyBFLoop(const ParameterType& param) const;

         /**
          * @brief Initialize Poly operator
          */
         template <typename T> void initOperator(T& op, const transf::Complex::SharedSetup spSetup) const;

         /**
          * @brief Initialize FFT operator
          */
         template <typename T> void initOperator(T& op, const transf::Mixed::SharedSetup spSetup) const;

         /**
          * @brief Build transform operator setup
          */
         std::shared_ptr<typename TOp::SetupType> buildSetup(const ParameterType& param, const TestType type) const;
   };

   template <typename TOp, typename TOp2> Tester<TOp,TOp2>::Tester(const std::string& fname, const bool keepData)
      : Transform::TesterBase<TOp>(fname, keepData)
   {
      this->appendPath();
   }

   template <typename TOp, typename TOp2> void Tester<TOp,TOp2>::appendPath()
   {
      this->mPath += "Fourier/";
      if constexpr(std::is_same_v<typename TOp::SetupType,transf::Mixed::Setup>)
      {
         this->mPath += "Mixed/";
      }
      else if(std::is_same_v<typename TOp::SetupType,transf::Complex::Setup>)
      {
         this->mPath += "Complex/";
      }
   }

   template <typename TOp, typename TOp2> void Tester<TOp,TOp2>::readFile(Matrix& data, const ParameterType& param, const TestType type, const ContentType ctype) const
   {
      if(param.size() == 1)
      {
         TesterBase<TOp>::basicReadFile(data, param, type, ctype);
      }
      else
      {
         this->dbReadFile(data, param, type, ctype);
      }
   }

   template <typename TOp, typename TOp2> void Tester<TOp,TOp2>::readFile(MatrixZ& data, const ParameterType& param, const TestType type, const ContentType ctype) const
   {
      if(param.size() == 1)
      {
         TesterBase<TOp>::basicReadFile(data, param, type, ctype);
      }
      else
      {
         this->dbReadFile(data, param, type, ctype);
      }
   }

   template <typename TOp, typename TOp2> template <typename TData> void Tester<TOp,TOp2>::dbReadFile(TData& data, const ParameterType& param, const TestType type, const ContentType ctype) const
   {
      // Create setup
      auto spSetup = this->buildSetup(param, type);

      // Read database file
      TData dbData = TData::Zero(data.rows(),spSetup->specSize());
      ParameterType dbParam = {param.at(0)};
      std::string fullname = this->makeFilename(dbParam, this->refRoot(), type, ctype);
      readData(dbData, fullname);

      // Loop over indexes
      int col = 0;
      for(int j = 0; j < spSetup->slowSize(); j++)
      {
         int m = spSetup->slow(j);
         // Loop over multiplier
         for(int i = 0; i < spSetup->mult(j); i++)
         {
            data.col(col) = dbData.col(m);
            col++;
         }
      }
   }

   template <typename TOp, typename TOp2> Matrix Tester<TOp,TOp2>::applyOperator(const ParameterType& param, const TestType type) const
   {
      Matrix outData;
      switch(type)
      {
         case TestType::PROJECTOR:
            outData = this->applyProjector(param);
            break;
         case TestType::INTEGRATOR:
            outData = this->applyIntegrator(param);
            break;
         case TestType::BFLOOP:
            outData = this->applyBFLoop(param);
            break;
         default:
            throw std::logic_error("Test type not implemented");
            break;
      }

      return outData;
   }

   template <typename TOp, typename TOp2> Matrix Tester<TOp,TOp2>::applyProjector(const ParameterType& param) const
   {
      if constexpr(has_transform<TOp, FwdType&, const BwdType&>::value)
      {
         const TestType type = TestType::PROJECTOR;

         // Create setup
         auto spSetup = this->buildSetup(param, type);

         // Input data
         BwdType inData(spSetup->specSize(),spSetup->blockSize());
         this->readFile(inData, param, type, ContentType::INPUT);

         TOp op;
         this->initOperator(op, spSetup);

         FwdType outData;

         for (unsigned int i = 0; i < this->mIter; ++i)
         {
            outData = FwdType(op.outRows(), op.outCols());
            op.transform(outData, inData);
         }

         if constexpr(std::is_same_v<FwdType,MatrixZ>)
         {
            Matrix out(2*outData.rows(),outData.cols());
            out.topRows(outData.rows()) = outData.real();
            out.bottomRows(outData.rows()) = outData.imag();

            return out;
         }
         else
         {
            return outData;
         }
      } else
      {
         throw std::logic_error("This operator is not an projector");

         Matrix out;
         return out;
      }
   }

   template <typename TOp, typename TOp2> Matrix Tester<TOp,TOp2>::applyIntegrator(const ParameterType& param) const
   {
      if constexpr(has_transform<TOp, BwdType&, const FwdType&>::value)
      {
         const TestType type = TestType::INTEGRATOR;

         // Create setup
         auto spSetup = this->buildSetup(param, type);

         // Input data
         FwdType inData(spSetup->fwdSize(),spSetup->blockSize());
         this->readFile(inData, param, type, ContentType::INPUT);

         TOp op;
         this->initOperator(op, spSetup);

         BwdType outData;

         for (unsigned int i = 0; i < this->mIter; ++i)
         {
            outData = BwdType(spSetup->bwdSize(), op.outCols());
            op.transform(outData, inData);
         }

         // note, we check only "dealiased" values
         MatrixZ outDealias(op.outRows(),outData.cols());
         op.dealias(outDealias, outData);

         // extract real/im
         Matrix out(2*op.outRows(),outData.cols());
         out.topRows(op.outRows()) = outDealias.real();
         out.bottomRows(op.outRows()) = outDealias.imag();

         return out;
      } else
      {
         throw std::logic_error("This operator is not an integrator");

         Matrix out;
         return out;
      }
   }

   template <typename TOp, typename TOp2> Matrix Tester<TOp,TOp2>::applyBFLoop(const ParameterType& param) const
   {
      if constexpr(std::is_same_v<TOp2,void>)
      {
         throw std::logic_error("Bacward-forward loop can only be computed if second operator type is given");
      }
      else
      {
         const TestType type = TestType::BFLOOP;

         // Create setup
         auto spSetup = this->buildSetup(param, type);

         // Input mods data
         BwdType inData(spSetup->specSize(),spSetup->blockSize());
         this->readFile(inData, param, type, ContentType::INPUT);

         TOp opB;
         this->initOperator(opB, spSetup);

         FwdType tmpData(opB.outRows(), opB.outCols());

         opB.transform(tmpData, inData);

         TOp2 opF;
         this->initOperator(opF, spSetup);

         BwdType outData(spSetup->bwdSize(), opF.outCols());

         opF.transform(outData, tmpData);

         // note, we check only "dealiased" values
         Matrix out(2*opF.outRows(),outData.cols());
         out.topRows(opF.outRows()) = outData.real().topRows(opF.outRows());
         out.bottomRows(opF.outRows()) = outData.imag().topRows(opF.outRows());

         return out;
      }
   }

   template <typename TOp, typename TOp2> std::string Tester<TOp,TOp2>::resname(const ParameterType& param) const
   {
      auto id = param.at(0);

      std::stringstream ss;
      ss.precision(10);
      ss << "_id" << id;

      if(param.size() == 3)
      {
         int np = param.at(1);
         ss << "_np" << np;
         int r = param.at(2);
         ss << "_r" << r;
         ss << "_stage2";
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

   template <typename TOp, typename TOp2> template <typename T> void Tester<TOp,TOp2>::initOperator(T& op, const transf::Complex::SharedSetup spSetup) const
   {
      op.init(spSetup);
   }

   template <typename TOp, typename TOp2> template <typename T> void Tester<TOp,TOp2>::initOperator(T& op, const transf::Mixed::SharedSetup spSetup) const
   {
      op.init(spSetup);
   }

   template <typename TOp, typename TOp2> std::shared_ptr<typename TOp::SetupType> Tester<TOp,TOp2>::buildSetup(const ParameterType& param, const TestType type) const
   {
      // Read metadata
      Array meta(0);
      std::string fullname = this->makeFilename(param, this->refRoot(), type, ContentType::META);
      readList(meta, fullname);

      // Create setup
      int nMeta = 3;
      int specN = meta(0);
      int physN = meta(1);
      int nModes = meta(2);
      auto spSetup = std::make_shared<typename TOp::SetupType>(physN, specN, GridPurpose::SIMULATION);
      spSetup->setBoxScale(1.0);

      // Gather indices
      std::map<int,int> indices;
      int nModes2D = 0;
      int h = nMeta;
      if(param.size() == 1)
      {
         for(int i = 2; i < meta.size(); i++)
         {
            int m = static_cast<int>(meta(i));
            indices[m]++;
         }
      }
      else
      {
         assert((meta.size() - nMeta - 2*nModes) % 2 == 0);
         for(int i = 0; i < nModes; i++)
         {
            int m = static_cast<int>(meta(h));
            int mult = static_cast<int>(meta(h+1));
            indices.insert(std::pair(m,mult));
            h += 2;
            nModes2D += mult;
         }

         // Check meta data size
         if(meta.size() - nMeta - 2*nModes - 2*nModes2D != 0)
         {
            throw std::logic_error("Meta data format is not supported (file: " + fullname + ")");
         }
      }

      // Add indices with multiplier
      for(const auto& [m, mult]: indices)
      {
         spSetup->addIndex(m, mult);
      }
      spSetup->lock();

      return spSetup;
   }

} // ALegendre
} // Transform
} // TestSuite
} // QuICC

#endif //QUICC_TESTSUITE_TRANSFORM_FOURIER_TESTER_HPP

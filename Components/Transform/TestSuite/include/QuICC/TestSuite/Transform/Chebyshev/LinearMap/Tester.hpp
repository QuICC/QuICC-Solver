/**
 * @file Tester.hpp
 * @brief Tester for Chebyshev linear map transforms
 */

#ifndef QUICC_TESTSUITE_TRANSFORM_CHEBYSHEV_LINEARMAP_TESTER_HPP
#define QUICC_TESTSUITE_TRANSFORM_CHEBYSHEV_LINEARMAP_TESTER_HPP


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
#include "QuICC/Transform/Fft/Chebyshev/Setup.hpp"

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

namespace transf = ::QuICC::Transform::Fft::Chebyshev;

namespace QuICC {

namespace TestSuite {

namespace Transform {

namespace Chebyshev {

namespace LinearMap {

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
         virtual ~Tester() = default;

      protected:
         /// Typedef ContentType from base
         typedef typename Transform::TesterBase<TOp>::ContentType ContentType;

         /// Typedef for Forward input data
         typedef MatrixZ FwdType;

         /// Typedef for Backward input data
         typedef MatrixZ BwdType;

         /**
          * @brief Build filename extension with resolution information
          */
         virtual std::string resname(const ParameterType& param) const override;

         /**
          * @brief Read data from file
          */
         virtual void readFile(Matrix& data, const ParameterType& param, const TestType type, const ContentType ctype) const override;

         virtual void readFile(MatrixZ& data, const ParameterType& param, const TestType type, const ContentType ctype) const override;

         template <typename TData> void dbReadFile(TData& data, const ParameterType& param, const TestType type, const ContentType ctype) const;

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

         /**
          * @brief Test projector
          */
         virtual Matrix applyProjector(const ParameterType& param) const;

         /**
          * @brief Test integrator
          */
         virtual Matrix applyIntegrator(const ParameterType& param) const;

         /**
          * @brief Test reductor
          */
         virtual Matrix applyReductor(const ParameterType& param) const;

         /**
          * @brief Test backward-forward loop
          */
         virtual Matrix applyBFLoop(const ParameterType& param) const;

         /**
          * @brief Initialize Poly operator
          */
         template <typename T> void initOperator(T& op, const transf::SharedSetup spSetup) const;

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
      this->mPath += "Chebyshev/LinearMap/";
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
         case TestType::REDUCTOR:
            outData = this->applyReductor(param);
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

         FwdType outData(op.outRows(), op.outCols());

         op.transform(outData, inData);

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

         BwdType outData(op.outRows(), op.outCols());

         op.transform(outData, inData);

         Matrix out(2*outData.rows(),outData.cols());
         out.topRows(outData.rows()) = outData.real();
         out.bottomRows(outData.rows()) = outData.imag();

         return out;
      } else
      {
         throw std::logic_error("This operator is not an integrator");

         Matrix out;
         return out;
      }
   }

   template <typename TOp, typename TOp2> Matrix Tester<TOp,TOp2>::applyReductor(const ParameterType& param) const
   {
      if constexpr(has_transform<TOp, Matrix&, const BwdType&>::value)
      {
         const TestType type = TestType::REDUCTOR;

         // Create setup
         auto spSetup = this->buildSetup(param, type);

         // Input data
         BwdType inData(spSetup->specSize(), spSetup->blockSize());
         this->readFile(inData, param, type, ContentType::INPUT);

         TOp op;
         this->initOperator(op, spSetup);

         Matrix outData(op.outRows(), op.outCols());

         op.transform(outData, inData);

         return outData;
      }
      else
      {
         throw std::logic_error("This operator is not an reductor");

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

         // Input data
         BwdType inData(spSetup->specSize(),spSetup->blockSize());
         this->readFile(inData, param, type, ContentType::INPUT);

         TOp opB;
         this->initOperator(opB, spSetup);

         FwdType tmpData(opB.outRows(), opB.outCols());

         opB.transform(tmpData, inData);

         TOp2 opF;
         this->initOperator(opF, spSetup);

         BwdType outData(opF.outRows(), opF.outCols());

         opF.transform(outData, tmpData);

         Matrix out(2*outData.rows(),outData.cols());
         out.topRows(outData.rows()) = outData.real();
         out.bottomRows(outData.rows()) = outData.imag();

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

   template <typename TOp, typename TOp2> template <typename T> void Tester<TOp,TOp2>::initOperator(T& op, const transf::SharedSetup spSetup) const
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
      int nMeta = 4;
      int specN = meta(0);
      int physN = meta(1);
      double lb = meta(2);
      double ub = meta(3);
      auto spSetup = std::make_shared<typename TOp::SetupType>(physN, meta.size() - nMeta, specN, GridPurpose::SIMULATION);
      spSetup->setBoxScale(1.0);
      spSetup->setBounds(lb, ub);

      // Gather indices
      std::map<int,int> indices;
      if(param.size() == 1)
      {
         for(int i = nMeta; i < meta.size(); i++)
         {
            int m = static_cast<int>(meta(i));
            indices[m]++;
         }
      }
      else
      {
         assert((meta.size() - nMeta) % 2 == 0);
         for(int i = nMeta; i < meta.size(); i += 2)
         {
            int m = static_cast<int>(meta(i));
            int mult = static_cast<int>(meta(i+1));
            indices.insert(std::pair(m,mult));
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

} // LinearMap
} // Chebyshev
} // Transform
} // TestSuite
} // QuICC

#endif //QUICC_TESTSUITE_TRANSFORM_CHEBYSHEV_LINEARMAP_TESTER_HPP

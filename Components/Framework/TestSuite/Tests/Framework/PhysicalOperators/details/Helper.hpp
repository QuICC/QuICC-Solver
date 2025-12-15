#include <cstdio>
#include <memory>
#include <catch2/catch.hpp>

#include "Memory/Cpu/NewDelete.hpp"
#include "Memory/MemoryResource.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/ScalarFields/ScalarFieldSetup.hpp"
#include "QuICC/ScalarFields/FlatScalarField.hpp"
#include "QuICC/ScalarFields/ViewScalarField.hpp"
#include "QuICC/VectorFields/VectorField.hpp"
#include "QuICC/TensorFields/TensorField.hpp"

namespace details {

   enum class SetupType {
      UniformUp,
   };

   std::shared_ptr<QuICC::Datatypes::ScalarFieldSetup> createSetup(std::size_t& dim1D, std::size_t& dim2D, std::size_t& dim3D, std::vector<std::vector<std::size_t>>& idx2D, std::vector<std::size_t>& idx3D, const SetupType type, const std::uint32_t variant);

   template <typename T, template <typename> class TField> std::shared_ptr<TField<T>> createScalarField(const std::size_t dim1D, const std::size_t dim3D, const std::vector<std::vector<std::size_t>>& idx2D, const std::vector<std::size_t>& idx3D, std::shared_ptr<QuICC::Datatypes::ScalarFieldSetup> spSetup, const std::uint32_t variant);

   template <typename T, template <typename> class TField> std::shared_ptr<QuICC::Datatypes::VectorField<TField<T>, QuICC::FieldComponents::Physical::Id>> createVectorField(const std::size_t dim1D, const std::size_t dim3D, const std::vector<std::vector<std::size_t>>& idx2D, const std::vector<std::size_t>& idx3D, std::shared_ptr<QuICC::Datatypes::ScalarFieldSetup> spSetup, const std::uint32_t variant);

   template <typename T, template <typename> class TField> std::shared_ptr<QuICC::Datatypes::TensorField<TField<T>, QuICC::FieldComponents::Physical::Id>> createTensorField(const std::size_t dim1D, const std::size_t dim3D, const std::vector<std::vector<std::size_t>>& idx2D, const std::vector<std::size_t>& idx3D, std::shared_ptr<QuICC::Datatypes::ScalarFieldSetup> spSetup, const std::uint32_t variant);

   std::vector<std::uint32_t> validationVariants();
   std::vector<std::uint32_t> performanceVariants();

   template <typename T> T fieldValueGeneric(const int i, const int j, const int k, const T offset, const T base, const T c);
   template <typename T> T fieldValueA(const int i, const int j, const int k);
   template <typename T> T fieldValueB(const int i, const int j, const int k);
   template <typename T> T fieldValueC(const int i, const int j, const int k);
   template <typename T> T fieldValueD(const int i, const int j, const int k);
   template <typename T> T fieldValueE(const int i, const int j, const int k);
   template <typename T> T fieldValueF(const int i, const int j, const int k);
   template <typename T> T fieldValueG(const int i, const int j, const int k);
   template <typename T> T fieldValueH(const int i, const int j, const int k);
   template <typename T> T fieldValueI(const int i, const int j, const int k);
   template <typename T> T fieldValueK(const int i, const int j, const int k);
   template <typename T> T fieldValueBad(const int i, const int j, const int k);

   template <typename T> T fieldValueGeneric(const int i, const int j, const int k, const T offset, const T base, const T c)
   {
      return c*(offset + i + j*base + k*base*base);
   }

   template <typename T> T fieldValueA(const int i, const int j, const int k)
   {
      return fieldValueGeneric<T>(i, j, k, 0.0, 100.0, 1.0);
   }

   template <typename T> T fieldValueB(const int i, const int j, const int k)
   {
      return fieldValueGeneric<T>(i, j, k, 42.0, 101.0, -1.0);
   }

   template <typename T> T fieldValueC(const int i, const int j, const int k)
   {
      return fieldValueGeneric<T>(i, j, k, -21.0, 102.0, -1.0);
   }

   template <typename T> T fieldValueD(const int i, const int j, const int k)
   {
      return fieldValueGeneric<T>(i, j, k, 37.0, 103.0, -1.0);
   }

   template <typename T> T fieldValueE(const int i, const int j, const int k)
   {
      return fieldValueGeneric<T>(i, j, k, -56.0, 104.0, -1.0);
   }

   template <typename T> T fieldValueF(const int i, const int j, const int k)
   {
      return fieldValueGeneric<T>(i, j, k, 17.0, 33.0, 1.0);
   }

   template <typename T> T fieldValueG(const int i, const int j, const int k)
   {
      return fieldValueGeneric<T>(i, j, k, -3.0, 27.0, -1.0);
   }

   template <typename T> T fieldValueH(const int i, const int j, const int k)
   {
      return fieldValueGeneric<T>(i, j, k, -11.0, 46.0, 1.0);
   }

   template <typename T> T fieldValueI(const int i, const int j, const int k)
   {
      return fieldValueGeneric<T>(i, j, k, 9.0, 15.0, -1.0);
   }

   template <typename T> T fieldValueJ(const int i, const int j, const int k)
   {
      return fieldValueGeneric<T>(i, j, k, -5.0, 39.0, 1.0);
   }

   template <typename T> T fieldValueBad(const int i, const int j, const int k)
   {
      return fieldValueGeneric<T>(i, j, k, 42.42, 0.01, -1.0);
   }

   struct IdxFunctor
   {
      const std::vector<std::vector<std::size_t>>& _idx2D;
      const std::vector<std::size_t>& _idx3D;

      IdxFunctor(const std::vector<std::vector<std::size_t>>& idx2D, const std::vector<std::size_t>& idx3D) : _idx2D(idx2D), _idx3D(idx3D) {};

      IdxFunctor() = delete;

      ~IdxFunctor() = default;

      std::size_t dim2D(const int k) const
      {
         assert(static_cast<std::size_t>(k) < _idx2D.size());
         return _idx2D.at(k).size();
      };

      std::size_t idx2D(const int j, const int k) const
      {
         assert(static_cast<std::size_t>(k) < _idx2D.size());
         assert(static_cast<std::size_t>(j) < _idx2D.at(k).size());
         return _idx2D.at(k).at(j);
      };

      std::size_t dim3D() const
      {
         return _idx3D.size();
      };

      std::size_t idx3D(const int k) const
      {
         assert(static_cast<std::size_t>(k) < _idx3D.size());
         return _idx3D.at(k);
      };

   };

   template <typename T, template <typename> class TField> void setFieldValue(TField<T>& field, const std::size_t dim1D, const std::size_t dim3D, const std::vector<std::vector<std::size_t>>& idx2D, const std::vector<std::size_t>& idx3D, std::shared_ptr<QuICC::Datatypes::ScalarFieldSetup> spSetup, const std::uint32_t variant)
   {
      using FieldValueFct = T (*)(const int i, const int j, const int k);
      FieldValueFct valueFct = nullptr;

      if(variant < 100)
      {
         switch(variant % 10)
         {
            case 0:
               valueFct = fieldValueA;
               break;
            case 1:
               valueFct = fieldValueB;
               break;
            case 2:
               valueFct = fieldValueC;
               break;
            case 3:
               valueFct = fieldValueD;
               break;
            case 4:
               valueFct = fieldValueE;
               break;
            case 5:
               valueFct = fieldValueF;
               break;
            case 6:
               valueFct = fieldValueG;
               break;
            case 7:
               valueFct = fieldValueH;
               break;
            case 8:
               valueFct = fieldValueI;
               break;
            case 9:
               valueFct = fieldValueJ;
               break;
            default:
               throw std::logic_error("Unknown variant");
               break;
         }
      }
      else
      {
         valueFct = fieldValueBad;
      }

      assert(static_cast<std::size_t>(spSetup->nBlock()) == idx3D.size());
      assert(static_cast<std::size_t>(spSetup->nBlock()) == idx2D.size());
      for(int k = 0; k < spSetup->nBlock(); k++)
      {
         int k_ = idx3D.at(k);
         for(int j = 0; j < spSetup->blockCols(k); j++)
         {
            assert(static_cast<std::size_t>(spSetup->blockCols(k)) == idx2D.at(k).size());
            assert(static_cast<std::size_t>(spSetup->blockRows(k)) == dim1D);

            int j_ = idx2D.at(k).at(j);
            for(int i = 0; i < spSetup->blockRows(k); i++)
            {
               field.setPoint(valueFct(i,j_,k_), i, j, k);
            }
         }
      }
   }

   template <typename T, template <typename> class TField> std::shared_ptr<TField<T>> createScalarField(const std::size_t dim1D, const std::size_t dim3D, const std::vector<std::vector<std::size_t>>& idx2D, const std::vector<std::size_t>& idx3D, std::shared_ptr<QuICC::Datatypes::ScalarFieldSetup> spSetup, const std::uint32_t variant)
   {
      auto spField = std::make_shared<TField<T>>(spSetup);
      auto&& field = *spField;
      setFieldValue(field, dim1D, dim3D, idx2D, idx3D, spSetup, variant);

      return spField;
   }

   template <typename T, template <typename> class TField> std::shared_ptr<QuICC::Datatypes::VectorField<TField<T>,QuICC::FieldComponents::Physical::Id>> createVectorField(const std::size_t dim1D, const std::size_t dim3D, const std::vector<std::vector<std::size_t>>& idx2D, const std::vector<std::size_t>& idx3D, std::shared_ptr<QuICC::Datatypes::ScalarFieldSetup> spSetup, const std::uint32_t variant)
   {
      constexpr auto R = QuICC::FieldComponents::Physical::R;
      constexpr auto THETA = QuICC::FieldComponents::Physical::THETA;
      constexpr auto PHI = QuICC::FieldComponents::Physical::PHI;
      std::map<QuICC::FieldComponents::Physical::Id, bool> comps = {
         {R, true},
         {THETA, true},
         {PHI, true}
      };
      auto spField = std::make_shared<QuICC::Datatypes::VectorField<TField<T>,QuICC::FieldComponents::Physical::Id>>(spSetup, comps);
      auto&& vfield = *spField;

      std::uint32_t var = variant;
      for(auto&& [c, v]: comps)
      {
         setFieldValue(vfield.rComp(c), dim1D, dim3D, idx2D, idx3D, spSetup, var);
         var++;
      }

      return spField;
   }

   template <typename T, template <typename> class TField> std::shared_ptr<QuICC::Datatypes::TensorField<TField<T>,QuICC::FieldComponents::Physical::Id>> createTensorField(const std::size_t dim1D, const std::size_t dim3D, const std::vector<std::vector<std::size_t>>& idx2D, const std::vector<std::size_t>& idx3D, std::shared_ptr<QuICC::Datatypes::ScalarFieldSetup> spSetup, const std::uint32_t variant)
   {
      constexpr auto R = QuICC::FieldComponents::Physical::R;
      constexpr auto THETA = QuICC::FieldComponents::Physical::THETA;
      constexpr auto PHI = QuICC::FieldComponents::Physical::PHI;
      std::map<std::pair<QuICC::FieldComponents::Physical::Id,QuICC::FieldComponents::Physical::Id>, bool> comps = {
         {{R,R}, true},
         {{R,THETA}, true},
         {{R,PHI}, true},
         {{THETA,R}, true},
         {{THETA,THETA}, true},
         {{THETA,PHI}, true},
         {{PHI,R}, true},
         {{PHI,THETA}, true},
         {{PHI,PHI}, true}
      };
      auto spField = std::make_shared<QuICC::Datatypes::TensorField<TField<T>,QuICC::FieldComponents::Physical::Id>>(spSetup, comps);
      auto&& vfield = *spField;

      std::uint32_t var = variant;
      for(auto&& [c, v]: comps)
      {
         setFieldValue(vfield.rComp(c.first, c.second), dim1D, dim3D, idx2D, idx3D, spSetup, var);
         var++;
      }

      return spField;
   }

template <typename TRef, typename TOther>
void checkComputation(std::shared_ptr<QuICC::Datatypes::ScalarFieldSetup> spSetup, const TRef& ref, const TOther& other)
{
   for(int j = 0; j < spSetup->dataCols(); j++)
   {
      INFO( " j = " << j );
      for(int i = 0; i < spSetup->dataRows(); i++)
      {
         INFO( " i = " << i );
         double r = ref.data()(i,j);
         double o = other.data()(i,j);
         REQUIRE_THAT( o, Catch::Matchers::WithinULP(r, 10000) );
      }
   }
}
}

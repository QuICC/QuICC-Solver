/** 
 * @file IPhysicalPyModel.hpp
 * @brief Interface for implementation of a physical models calling into Python model
 */

#ifndef QUICC_MODEL_IPHYSICALPYMODEL_HPP
#define QUICC_MODEL_IPHYSICALPYMODEL_HPP

// System includes
//
#include "QuICC/PyQuICC/SystemHeader.hpp"
#include <string>
#include <vector>
#include <set>
#include <memory>

// Project includes
//
#include "QuICC/Model/IPhysicalModel.hpp"
#include "QuICC/Model/PyModelBackend.hpp"

namespace QuICC {

namespace Model {

   /**
    * @brief Interface for the implementation of a physical models calling into Python model
    */
   template <typename TSim, typename TState, typename TVis> class IPhysicalPyModel: public IPhysicalModel<TSim, TState, TVis>
   {
      public:
         /**
          * @brief Constructor
          */
         IPhysicalPyModel() = default;

         /**
          * @brief Destructor
          */
         virtual ~IPhysicalPyModel();

         /// Python script/module name
         virtual std::string PYMODULE() = 0;

         /// Python model class name, uses PhysicalModel by default
         virtual std::string PYCLASS();

         /**
          * @brief Initialize model
          */
         virtual void init();

      protected:

      private:
   };

   template <typename TSim, typename TState, typename TVis> IPhysicalPyModel<TSim,TState,TVis>::~IPhysicalPyModel()
   {
   }

   template <typename TSim, typename TState, typename TVis> std::string IPhysicalPyModel<TSim,TState,TVis>::PYCLASS()
   {
      return "PhysicalModel";
   }

   template <typename TSim, typename TState, typename TVis> void IPhysicalPyModel<TSim,TState,TVis>::init()
   {
      IPhysicalModel<TSim,TState,TVis>::init();

      this->mpBackend = std::make_shared<PyModelBackend>(this->PYMODULE(), this->PYCLASS());
   }

}
}

#endif // QUICC_MODEL_IPHYSICALPYMODEL_HPP

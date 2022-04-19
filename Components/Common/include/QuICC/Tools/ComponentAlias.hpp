/** 
 * @file ComponentAlias.hpp
 * @brief Implementation of simple alias for field components
 */

#ifndef QUICC_TOOLS_COMPONENTALIAS_HPP
#define QUICC_TOOLS_COMPONENTALIAS_HPP

// Configuration includes
//

// System includes
//
#include <vector>

// External includes
//

// Project includes
//

namespace QuICC {

namespace Tools {

   /**
    * @brief Implementation of a few useful output formating static functions
    */
   template <typename T> class ComponentAlias
   {
      public:
         /**
          * @brief Constructor
          */
         ComponentAlias();

         /**
          * @brief Destructor
          */
         ~ComponentAlias();

         /**
          * @brief Add component to alias
          */
         void add(const T comp);

         /**
          * @brief Alias for first component
          */
         T ONE() const;

         /**
          * @brief Alias for second component
          */
         T TWO() const;

         /**
          * @brief Alias for third component
          */
         T THREE() const;

         /**
          * @brief Get component id by index
          * 
          * @param i component index
          */
         T comp(const int i) const;

         /**
          * @brief Get component id by index
          * 
          * @param i component index
          */
         typename std::vector<T>::const_iterator cbegin() const;

         /**
          * @brief Get component id by index
          * 
          * @param i component index
          */
         typename std::vector<T>::const_iterator cend() const;

      protected:

      private:
         /**
          * @brief Storage for the component Aliases
          */
         std::vector<T> mComp;
   };

   template <typename T>
      ComponentAlias<T>::ComponentAlias()
   {
   }

   template <typename T>
      ComponentAlias<T>::~ComponentAlias()
   {
   }

   template <typename T>
      T ComponentAlias<T>::ONE() const
   {
      assert(this->mComp.size() > 0);

      return this->mComp.at(0);
   }

   template <typename T>
      T ComponentAlias<T>::TWO() const
   {
      assert(this->mComp.size() > 1);

      return this->mComp.at(1);
   }

   template <typename T>
      T ComponentAlias<T>::THREE() const
   {
      assert(this->mComp.size() > 2);

      return this->mComp.at(2);
   }

   template <typename T>
      T ComponentAlias<T>::comp(const int i) const
   {
      assert(this->mComp.size() > i);

      return this->mComp.at(i);
   }

   template <typename T>
      typename std::vector<T>::const_iterator ComponentAlias<T>::cbegin() const
   {
      return this->mComp.cbegin();
   }

   template <typename T>
      typename std::vector<T>::const_iterator ComponentAlias<T>::cend() const
   {
      return this->mComp.cend();
   }

   template <typename T>
      void ComponentAlias<T>::add(const T comp)
   {
      return this->mComp.push_back(comp);
   }

}
}

#endif // QUICC_TOOLS_COMPONENTALIAS_HPP

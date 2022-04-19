/**
 * @file IteratorRange.hpp
 * @brie Interface for range based for loop with iterator pair
 */

#ifndef QUICC_ITERATORRANGE_HPP
#define QUICC_ITERATORRANGE_HPP

// System includes
//
#include <map>

// External includes
//

// Project includes
//

namespace QuICC {

   /**
    * @brie Interface for range based for loop with iterator pair
    */
   template <typename T> class IteratorRange 
   {
      public:
         /**
          * @brief Constructor from two iterators
          */
         IteratorRange(T b, T e);

         /**
          * @brief Destructor
          */
         ~IteratorRange();

         /**
          * @brief Get begin of range
          */
         T begin() const;

         /**
          * @brief Get end of range
          */
         T end() const;

      private:
         /**
          * @brief Storage for begin iterator
          */
         T mB;

         /**
          * @brief Storage for end iterator
          */
         T mE;
   };

   template <typename T> IteratorRange<T> make_range(T b, T e);

   template <typename T> IteratorRange<T> make_range(const std::pair<T,T>& p);

   //
   //
   //
   template <typename T> IteratorRange<T>::IteratorRange(T b, T e)
      : mB(b), mE(e)
   {
   }

   template <typename T> IteratorRange<T>::~IteratorRange()
   {
   }

   template <typename T> T IteratorRange<T>::begin() const
   {
      return this->mB;
   }

   template <typename T> T IteratorRange<T>::end() const
   {
      return this->mE;
   }

   template <typename T> IteratorRange<T> make_range(T b, T e)
   {
      return IteratorRange<T>(b,e);
   }

   template <typename T> IteratorRange<T> make_range(const std::pair<T,T>& p)
   {
      return IteratorRange<T>(p.first,p.second);
   }

}

#endif // QUICC_ITERATORRANGE_HPP

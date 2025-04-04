/**
 * @file CacheLayerSize.hpp
 * @brief Cache layer size
 */
#pragma once

// External includes
//

#include <map>

// Project includes
//


namespace QuICC
{
namespace Graph
{
namespace details
{

/// @brief CacheLayerSize class
/// Computing the size on the device and copying to the host in order
/// to allocate the temporary is expensive.
/// We cache the size of the layers in a map based on the metadata pointer.
/// @tparam T type of ptr to cache
template <class T>
class CacheLayerSize
{
public:
    /// @brief delete copy ctor
    /// @param
    CacheLayerSize(CacheLayerSize const&) = delete;
    /// @brief delete assignment op
    /// @param
    void operator=(CacheLayerSize const&) = delete;
    /// @brief create/return singleton
    /// @return singleton
    static CacheLayerSize& getInstance();
    /// @brief invalidate the entire cache
    void invalidate();
    /// @brief invalidate a specific entry in the cache
    /// @param ptr pointer to the metadata to remove
    void invalidate(const T* ptr);
    /// @brief check if a value is in the cache and return it, otherwise return zero
    /// @param ptr pointer to the metadata to check
    /// @return cached size or zero if not found
    std::size_t getCachedSizeOrZero(const T* ptr) const;
    /// @brief add an entry to the cache
    /// @param ptr pointer to the metadata
    /// @param size size to cache
    void addToCache(const T* ptr, const std::size_t size);
private:
    /// @brief map pointe to size
    std::map<T*, std::size_t> _ptr2size;
    /// @brief Default ctor
    CacheLayerSize() = default;
    /// @brief dtor
    ~CacheLayerSize() = default;
};

template <class T>
CacheLayerSize<T>& CacheLayerSize<T>::getInstance()
{
    static CacheLayerSize instance;
    return instance;
}

template <class T>
void CacheLayerSize<T>::invalidate()
{
    _ptr2size.clear();
}

template <class T>
void CacheLayerSize<T>::invalidate(const T* ptr)
{
    _ptr2size.erase(ptr);
}

template <class T>
std::size_t CacheLayerSize<T>::getCachedSizeOrZero(const T* ptr) const
{
    auto it = _ptr2size.find(ptr);
    if (it != _ptr2size.end()) {
        return it->second;
    }
    return 0;
}

template <class T>
void CacheLayerSize<T>::addToCache(const T* ptr, const std::size_t size)
{
    _ptr2size[ptr] = size;
}

} // namespace details
} // namespace Graph
} // namespace QuICC
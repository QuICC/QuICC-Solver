/**
 * @file Pensieve.hpp
 * @brief Implements a singleton providing access to a global memory resource
 * https://en.wikipedia.org/wiki/Magical_objects_in_Harry_Potter#Pensieve
 */

#pragma once

// External includes
//
#include <memory>
#include <type_traits>

// Project includes
//
#include "Memory/MemoryResource.hpp"

namespace QuICC {
namespace Memory {

/**
* @brief Static interface to the singleton providing access to a
* global memory resource
*/
template <class Tmem>
class Pensieve
{
public:
    /// @brief delete copy ctor
    /// @param
    Pensieve(Pensieve const&) = delete;
    /// @brief delete assignment op
    /// @param
    void operator=(Pensieve const&) = delete;
    /// @brief create/return singleton
    /// @return singleton
    static Pensieve& getInstance();
    /// @brief access the memory resource
    /// @return
    Tmem& getMem();
private:
    /// @brief memory resource holder
    Tmem _mem;
    /// @brief Default ctor
    Pensieve() = default;
    /// @brief dtor
    ~Pensieve() = default;
};

template <class Tmem>
Pensieve<Tmem>& Pensieve<Tmem>::getInstance()
{
    static_assert(std::is_base_of_v<memory_resource, Tmem>, "Tmem is not a memory resource");
    static Pensieve instance;
    return instance;
}

template <class Tmem>
Tmem& Pensieve<Tmem>::getMem()
{
    return _mem;
}

} // namespace Memory
} // namespace QuICC


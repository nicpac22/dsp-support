/* 
 * Copyright (c) 2026 Nick Xenias
 * Based on 2012 software developed by Michael Ihde
 * 
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU Lesser General Public License as   
 * published by the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * Lesser General Lesser Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

// requires -D_XOPEN_SOURCE 600 or higher

#ifndef _ALIGNED_ALLOCATOR_H_
#define _ALIGNED_ALLOCATOR_H_

#include <cstddef>   // size_t, ptrdiff_t
#include <cstdlib>
#include <new>       // std::bad_alloc
#include <stdlib.h>

// Alignment and trailing pad are both specified in BYTES.  The alignment must
// be a power of two and at least sizeof(void*) (a posix_memalign requirement).
template<size_t AlignBytes, size_t PadBytes = 0>
struct aligned_allocator_traits
{
  static_assert((AlignBytes & (AlignBytes-1)) == 0 && AlignBytes >= sizeof(void*),
      "aligned_allocator_traits: alignment must be a power of two and at least "
      "sizeof(void*)");
  static constexpr size_t align_bytes = AlignBytes;
  static constexpr size_t pad_bytes = PadBytes;
};

typedef aligned_allocator_traits<16> align_16;
typedef aligned_allocator_traits<32> align_32;
typedef aligned_allocator_traits<64> align_64;

// Alignment with extra 2xalignment worth of pad bytes appended to the end,
// allowing for over-indexing on read of up to 2 SIMD registers without
// causing segfaults
typedef aligned_allocator_traits<16, 32> evm_16;
typedef aligned_allocator_traits<32, 64> evm_32;
typedef aligned_allocator_traits<64, 128> evm_64;
typedef evm_64 evm_max;

template<typename _Tp, typename traits = align_16>
class aligned_allocator
{
  public:
    typedef size_t     size_type;
    typedef ptrdiff_t  difference_type;
    typedef _Tp*       pointer;
    typedef const _Tp* const_pointer;
    typedef _Tp&       reference;
    typedef const _Tp& const_reference;
    typedef _Tp        value_type;

    template<typename _Tp1>
    struct rebind
    {
      typedef aligned_allocator<_Tp1,traits > other;
    };

    aligned_allocator() throw()
    {}

    aligned_allocator(const aligned_allocator&) throw()
    {}

    template<typename _Tp1>
    aligned_allocator(const aligned_allocator<_Tp1,traits>&) throw()
    {}

    ~aligned_allocator() throw()
    {}

    pointer address(reference __x) const
    {
      return &__x;
    }

    const_pointer address(const_reference __x) const
    {
      return &__x;
    }

    // NB: __n is permitted to be 0.  The C++ standard says nothing
    // about what the return value is when __n == 0.
    pointer allocate(size_type __n, const void* = 0)
    {
      if(__builtin_expect(__n > this->max_size(), false))
      {
	      throw std::bad_alloc();
      }

      void* tmpvalue = 0;
      int ret = posix_memalign(&tmpvalue,
                               traits::align_bytes,
                               (__n * sizeof(_Tp))+traits::pad_bytes);
      if(ret)
      {
	      throw std::bad_alloc();
      }
      if(!tmpvalue)
      {
	      throw std::bad_alloc();
      }
      
      return reinterpret_cast<_Tp*>(tmpvalue);
    }

    // __p is not permitted to be a null pointer.
    void deallocate(pointer __p, size_type)
    {
      free(static_cast<void*>(__p));
    }

    size_type max_size() const throw() 
    {
      // leave room for the pad bytes so n*sizeof(_Tp)+pad_bytes can't wrap
      return (size_t(-1) - traits::pad_bytes) / sizeof(_Tp);
    }

    // _GLIBCXX_RESOLVE_LIB_DEFECTS
    // 402. wrong new expression in [some_] allocator::construct
    void construct(pointer __p, const _Tp& __val) 
    {
      ::new(__p) value_type(__val);
    }

    void destroy(pointer __p)
    {
      __p->~_Tp();
    }
};

template<typename _Tp1, typename _Tp2, typename traits>
inline bool operator==(const aligned_allocator<_Tp1,traits>&,
                       const aligned_allocator<_Tp2,traits>&)
{
  return true;
}

template<typename _Tp1, typename _Tp2, typename traits>
inline bool operator!=(const aligned_allocator<_Tp1,traits>&,
                       const aligned_allocator<_Tp2,traits>&)
{
  return false;
}

#endif /* _ALIGNED_ALLOCATOR_H_ */

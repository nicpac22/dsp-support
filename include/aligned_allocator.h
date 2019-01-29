/* 
 * Copyright (c) 2012 Michael Ihde
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

#include <cstdlib>
#include <bits/functexcept.h>
#include <stdlib.h>

template <typename alignT>
class aligned_allocator_traits
{
  public:
    typedef alignT align_type;
    static const size_t align_bits = (sizeof(align_type)*8);
    static const size_t align_bytes = sizeof(align_type);
};

// NB: sizes are BYTE alignment, not bits, i.e. align_16 is a 16-byte alignment
typedef aligned_allocator_traits<int16_t> align_16;
typedef aligned_allocator_traits<int32_t> align_32;
typedef aligned_allocator_traits<int64_t> align_64;

template<typename _Tp, typename traits = aligned_allocator_traits<int16_t> >
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
	      std::__throw_bad_alloc();
      }

      void* tmpvalue = 0;
      int ret = posix_memalign(&tmpvalue,
                               traits::align_bits,
                               (__n * sizeof(_Tp)));
      if(ret)
      {
	      std::__throw_bad_alloc();
      }
      if(!tmpvalue)
      {
	      std::__throw_bad_alloc();
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
      return size_t(-1) / sizeof(_Tp);
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

template<typename _Tp>
inline bool operator==(const aligned_allocator<_Tp>&,
                       const aligned_allocator<_Tp>&)
{
  return true;
}

template<typename _Tp>
inline bool operator!=(const aligned_allocator<_Tp>&,
                       const aligned_allocator<_Tp>&)
{
  return false;
}

#endif /* _ALIGNED_ALLOCATOR_H_ */

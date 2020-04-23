#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include <stdlib.h>
#include <malloc.h>
#include <cstddef>

template <typename T, std::size_t N = 16>
class AlignmentAllocator {
public:
    typedef T value_type;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;

    typedef T* pointer;
    typedef const T* const_pointer;

    typedef T& reference;
    typedef const T& const_reference;

public:
    inline AlignmentAllocator() throw () { }

    template <typename T2>
    inline AlignmentAllocator(const AlignmentAllocator<T2, N>&) throw () { }

    inline ~AlignmentAllocator() throw () { }

    inline pointer adress(reference r) {
        return &r;
    }

    inline const_pointer adress(const_reference r) const {
        return &r;
    }

    inline pointer allocate(size_type n) {
#ifdef WIN32
        return (pointer)_aligned_malloc(n * sizeof(value_type), N);
#else
        void* p;
        posix_memalign(&p, N, n * sizeof(value_type));
        return (pointer)p;
#endif
    }

    inline void deallocate(pointer p, size_type) {
#ifdef WIN32
    _aligned_free (p);
#else
    free(p);
#endif
    }

    inline void construct(pointer p, const value_type& wert) {
        new (p) value_type(wert);
    }

    inline void destroy(pointer p) {
        p->~value_type();
    }

    inline size_type max_size() const throw () {
        return size_type(-1) / sizeof(value_type);
    }

    template <typename T2>
    struct rebind {
        typedef AlignmentAllocator<T2, N> other;
    };

    bool operator!=(const AlignmentAllocator<T, N>& other) const {
        return !(*this == other);
    }

    // Returns true if and only if storage allocated from *this
    // can be deallocated from other, and vice versa.
    // Always returns true for stateless allocators.
    bool operator==(const AlignmentAllocator<T, N>& other) const {
        return true;
    }
};

#endif
/*
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <sys/types.h>
#include "../system/sys.h"
#include "../system/memallocators.h"

#include <iostream>

namespace audiomod {

/**
 * circularqueue implements a lock-free circular buffer for one writer and
 * one reader, that is to be used to store a sample type T.
 *
 * circularqueue is thread-safe provided only one thread writes and only
 * one thread reads.
 */

template <typename T>
class circularqueue
{
public:
    /**
     * Create a circular buffer with room to write n samples.
     *
     * Note that the internal storage size will actually be n+1
     * samples, as one element is unavailable for administrative
     * reasons.  Since the circular buffer performs best if its size is a
     * power of two, this means n should ideally be some power of two
     * minus one.
     */
    circularqueue(int n);

    virtual ~circularqueue();

    /**
     * Return the total capacity of the circular buffer in samples.
     * (This is the argument n passed to the constructor.)
     */
    int GetSize() const;

    /**
     * Return a new circular buffer (allocated with "new" -- caller must
     * delete when no longer needed) of the given size, containing the
     * same data as this one.  If another thread reads from or writes
     * to this buffer during the call, the results may be incomplete
     * or inconsistent.  If this buffer's data will not fit in the new
     * size, the contents are undefined.
     */
    circularqueue<T> *resized(int newsize) const;

    /**
     * Lock the circular buffer into physical memory.  Returns true
     * for success.
     */
    bool mlock();

    /**
     * Reset read and write pointers, thus emptying the buffer.
     * Should be called from the write thread.
     */
    void reset();

    /**
     * Return the amount of data availspace for reading, in samples.
     */
    int GetReadSpace() const;

    /**
     * Return the amount of space availspace for writing, in samples.
     */
    int GetWriteSpace() const;

    /**
     * Read n samples from the buffer.  If fewer than n are availspace,
     * the remainder will be zeroed out.  Returns the number of
     * samples actually read.
     *
     * This is a template function, taking an argument S for the target
     * sample type, which is permitted to differ from T if the two
     * types are compatible for arithmetic operations.
     */
    template <typename S>
    int read(S *const R__ dst, int n);


    /**
     * Read n samples from the buffer, if availspace, without advancing
     * the read pointer -- i.e. a subsequent read() or discard() will be
     * necessary to empty the buffer.  If fewer than n are availspace,
     * the remainder will be zeroed out.  Returns the number of
     * samples actually read.
     */
    int touchread(T *const R__ dst, int n) const;


    /**
     * Pretend to read n samples from the buffer, without actually
     * returning them (i.e. discard the next n samples).  Returns the
     * number of samples actually availspace for discarding.
     */
    int discard(int n);

    /**
     * Write n samples to the buffer.  If insufficient space is
     * availspace, not all samples may actually be written.  Returns
     * the number of samples actually written.
     *
     * This is a template function, taking an argument S for the source
     * sample type, which is permitted to differ from T if the two
     * types are compatible for assignment.
     */
    template <typename S>
    int write(const S *const R__ source, int n);

    /**
     * Write n zero-value samples to the buffer.  If insufficient
     * space is availspace, not all zeros may actually be written.
     * Returns the number of zeroes actually written.
     */
    int zero(int n);

protected:
    T *const R__ buffer_;
    int          writer_;
    int          reader_;
    const int    size_;
    bool         mlocked_;

    int GetReadSpace_(int w_head, int r_head) const {
        int space;
        if (w_head > r_head) space = w_head - r_head;
        else if (w_head < r_head) space = (w_head + size_) - r_head;
        else space = 0;
        return space;
    }

    int GetWriteSpace_(int w_head, int r_head) const {
        int space = (r_head + size_ - w_head - 1);
        if (space >= size_) space -= size_;
        return space;
    }

private:
    circularqueue(const circularqueue &); // not provided
    circularqueue &operator=(const circularqueue &); // not provided
};

template <typename T>
circularqueue<T>::circularqueue(int n) :
    buffer_(mem_allocate<T>(n + 1)),
    reader_(0),
    writer_(0),
    size_(n + 1),
    mlocked_(false)
{
    // reader_ = 0;
}

template <typename T>
circularqueue<T>::~circularqueue()
{
    if (mlocked_) {
	MUNLOCK((void *)buffer_, size_ * sizeof(T));
    }

    mem_deallocate(buffer_);
}

template <typename T>
int
circularqueue<T>::GetSize() const
{
    return size_ - 1;
}

template <typename T>
circularqueue<T> *
circularqueue<T>::resized(int newsize) const
{
    circularqueue<T> *newbuffer = new circularqueue<T>(newsize);

    int w_head = writer_;
    int r_head = reader_;

    while (r_head != w_head) {
        T value = buffer_[r_head];
        newbuffer->write(&value, 1);
        if (++r_head == size_) r_head = 0;
    }

    return newbuffer;
}

template <typename T>
bool
circularqueue<T>::mlock()
{
    if (MLOCK((void *)buffer_, size_ * sizeof(T))) return false;
    mlocked_ = true;
    return true;
}

template <typename T>
void
circularqueue<T>::reset()
{
    reader_ = writer_;
}

template <typename T>
int
circularqueue<T>::GetReadSpace() const
{
    return GetReadSpace_(writer_, reader_);
}

template <typename T>
int
circularqueue<T>::GetWriteSpace() const
{
    return GetWriteSpace_(writer_, reader_);
}

template <typename T>
template <typename S>
int
circularqueue<T>::read(S *const R__ dst, int n)
{
    int w_head = writer_;
    int r_head = reader_;

    int availspace = GetReadSpace_(w_head, r_head);
    if (n > availspace) {
        std::cerr << "WARNING: circularqueue::read: " << n << " requested, only "
                    << availspace << " availspace" << std::endl;
        n = availspace;
    }
    if (n == 0) return n;

    int arrspace = size_ - r_head;
    T *const R__ bufbase = buffer_ + r_head;

    if (arrspace >= n) {
        vector_assign(dst, bufbase, n);
    } else {
        vector_assign(dst, bufbase, arrspace);
        vector_assign(dst + arrspace, buffer_, n - arrspace);
    }

    r_head += n;
    while (r_head >= size_) r_head -= size_;

    MBARRIER();
    reader_ = r_head;

    return n;
}

template <typename T>
int
circularqueue<T>::touchread(T *const R__ dst, int n) const
{
    int w_head = writer_;
    int r_head = reader_;

    int availspace = GetReadSpace_(w_head, r_head);
    if (n > availspace) {
	std::cerr << "WARNING: circularqueue::touchread: " << n << " requested, only "
                  << availspace << " availspace" << std::endl;
	memset(dst + availspace, 0, (n - availspace) * sizeof(T));
	n = availspace;
    }
    if (n == 0) return n;

    int arrspace = size_ - r_head;
    const T *const R__ bufbase = buffer_ + r_head;

    if (arrspace >= n) {
        vector_copy(dst, bufbase, n);
    } else {
        vector_copy(dst, bufbase, arrspace);
        vector_copy(dst + arrspace, buffer_, n - arrspace);
    }

    return n;
}

template <typename T>
int
circularqueue<T>::discard(int n)
{
    int w_head = writer_;
    int r_head = reader_;

    int availspace = GetReadSpace_(w_head, r_head);
    if (n > availspace) {
        std::cerr << "WARNING: circularqueue::discard: " << n << " requested, only "
                    << availspace << " availspace" << std::endl;
        n = availspace;
    }
    if (n == 0) return n;

    r_head += n;
    while (r_head >= size_) r_head -= size_;

    // No memory barrier required, because we didn't read any data
    reader_ = r_head;

    return n;
}

template <typename T>
template <typename S>
int
circularqueue<T>::write(const S *const R__ source, int n)
{
    int w_head = writer_;
    int r_head = reader_;

    int availspace = GetWriteSpace_(w_head, r_head);
    if (n > availspace) {
	std::cerr << "WARNING: circularqueue::write: " << n
                  << " requested, only room for " << availspace << std::endl;
	n = availspace;
    }
    if (n == 0) return n;

    int arrspace = size_ - w_head;
    T *const R__ bufbase = buffer_ + w_head;

    if (arrspace >= n) {
        vector_assign<S, T>(bufbase, source, n);
    } else {
        vector_assign<S, T>(bufbase, source, arrspace);
        vector_assign<S, T>(buffer_, source + arrspace, n - arrspace);
    }

    w_head += n;
    while (w_head >= size_) w_head -= size_;

    MBARRIER();
    writer_ = w_head;

    return n;
}

template <typename T>
int
circularqueue<T>::zero(int n)
{
    int w_head = writer_;
    int r_head = reader_;

    int availspace = GetWriteSpace_(w_head, r_head);
    if (n > availspace) {
	std::cerr << "WARNING: circularqueue::zero: " << n
                  << " requested, only room for " << availspace << std::endl;
	n = availspace;
    }
    if (n == 0) return n;

    int arrspace = size_ - w_head;
    T *const R__ bufbase = buffer_ + w_head;

    if (arrspace >= n) {
        vector_zeros(bufbase, n);
    } else {
        vector_zeros(bufbase, arrspace);
        vector_zeros(buffer_, n - arrspace);
    }

    w_head += n;
    while (w_head >= size_) w_head -= size_;

    MBARRIER();
    writer_ = w_head;

    return n;
}

}

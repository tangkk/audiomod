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

#include <cmath>
#include <cstdlib>
#include <map>

#include "../system/sys.h"
#include "../system/vectorization.h"
#include "../system/memallocators.h"

namespace audiomod {

enum windowfuncType {
    Rectangular,
    Bartlett,
    Hamming,
    Hanning,
    Blackman,
    Gaussian,
    Nuttall,
    BlackmanHarris
};

template <typename T>
class windowfunc
{
public:
    /**
     * Construct a windower of the given type.
     */
    windowfunc(windowfuncType type, int size) : type_(type), size_(size), window_(0) {
        makewindow();
    }
    windowfunc(const windowfunc &w) : type_(w.type_), size_(w.size_), window_(0) {
        makewindow();
    }
    windowfunc &operator=(const windowfunc &w) {
        if (&w == this) return *this;
        type_ = w.type_;
        size_ = w.size_;
        window_ = 0;
        makewindow();
        return *this;
    }
    virtual ~windowfunc() {
        mem_deallocate(window_);
    }
    
    inline void apply(T *const R__ block) const {
        vector_multiply(block, window_, size_);
    }

    inline void apply(const T *const R__ src, T *const R__ dst) const {
        vector_multiply(dst, src, window_, size_);
    }

    inline void add2dst(T *const R__ dst, T scale) const {
        vector_add_with_gain(dst, window_, scale, size_);
    }

    inline T GetRMS() const {
        T total = 0;
        for (int i = 0; i < size_; ++i) {
            total += window_[i] * window_[i];
        }
        T rms = sqrt(total / size_);
        return rms;
    }

    inline T GetArea() const { return area_; }
    inline T GetValue(int i) const { return window_[i]; }

    inline windowfuncType GetType() const { return type_; }
    inline int GetSize() const { return size_; }

protected:
    windowfuncType type_;
    int size_;
    T *R__ window_;
    T area_;
    
    void makewindow();
    void wincore(T *, T, T, T, T);
};

template <typename T>
void windowfunc<T>::makewindow()
{
    if (!window_) window_ = mem_allocate<T>(size_);

    const int n = size_;
    vector_set(window_, T(1.0), n);
    int i;

    switch (type_) {
		
    case Rectangular:
	for (i = 0; i < n; ++i) {
	    window_[i] *= 0.5;
	}
	break;
	    
    case Bartlett:
	for (i = 0; i < n/2; ++i) {
	    window_[i] *= (i / T(n/2));
	    window_[i + n/2] *= (1.0 - (i / T(n/2)));
	}
	break;
	    
    case Hamming:
        wincore(window_, 0.54, 0.46, 0.0, 0.0);
	break;
	    
    case Hanning:
        wincore(window_, 0.50, 0.50, 0.0, 0.0);
	break;
	    
    case Blackman:
        wincore(window_, 0.42, 0.50, 0.08, 0.0);
	break;
	    
    case Gaussian:
	for (i = 0; i < n; ++i) {
            window_[i] *= pow(2, - pow((i - (n-1)/2.0) / ((n-1)/2.0 / 3), 2));
	}
	break;

    case Nuttall:
        wincore(window_, 0.3635819, 0.4891775, 0.1365995, 0.0106411);
	break;

    case BlackmanHarris:
        wincore(window_, 0.35875, 0.48829, 0.14128, 0.01168);
        break;
    }
	
    area_ = 0;
    for (i = 0; i < n; ++i) {
        area_ += window_[i];
    }
    area_ /= n;
}

template <typename T>
void windowfunc<T>::wincore(T *mult, T a0, T a1, T a2, T a3)
{
    int n = int(size_);
    for (int i = 0; i < n; ++i) {
        mult[i] *= (a0
                    - a1 * cos(2 * M_PI * i / n)
                    + a2 * cos(4 * M_PI * i / n)
                    - a3 * cos(6 * M_PI * i / n));
    }
}

}


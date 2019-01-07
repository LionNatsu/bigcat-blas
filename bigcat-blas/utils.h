/*
    bigcat-blas is a simple C++ linear algebra library.
    Copyright (C) 2019  Daming Yang

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef BIGCAT_BLAS_UTILS_H
#define BIGCAT_BLAS_UTILS_H

#include <cmath>

namespace blas {
    float fast_inv_sqrt(float x) {
        /*
         * This could even return a negative number.
         * Don't use this algorithm.
        float half_x = 0.5f * x;
        float y = x;
        long i = *(long *) &y;
        i = 0x5f3759df - (i >> 1);
        y = *(float *) &i;
        y = y * (1.5f - (half_x * y * y));
        y = y * (1.5f - (half_x * y * y));
        if(y <= 0) {
                std::cerr << x << " ! " << y << std::endl;
        }*/
        return 1 / sqrtf(x);
    }

    const float PI = 3.141593f;
}

#endif //BIGCAT_BLAS_UTILS_H

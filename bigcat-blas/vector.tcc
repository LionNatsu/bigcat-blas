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

#ifndef BIGCAT_BLAS_VECTOR_TCC
#define BIGCAT_BLAS_VECTOR_TCC

#include "matrix.h"
#include "vector.h"


#include <cstring>

namespace blas {

    template<size_t M>
    vector<M>::vector(const matrix<M, 1> &rhs) {
        operator=(rhs);
    }

    template<size_t M>
    vector <M> &vector<M>::operator=(const matrix<M, 1> &rhs) {
        matrix<M, 1>::operator=(rhs);
        return *this;
    }

    template<size_t M>
    vector<M>::vector(std::initializer_list<float> list) { operator=(list); }

    template<size_t M>
    vector <M> &vector<M>::operator=(std::initializer_list<float> list) {
        size_t m = 0;
        for (auto e : list) {
            operator()(m) = e;
            m++;
        }
        return *this;
    }

    template<size_t M>
    float &vector<M>::operator()(size_t m) {
        return matrix<M, 1>::operator()(m, 0);
    }

    template<size_t M>
    float vector<M>::operator()(size_t m) const {
        return matrix<M, 1>::operator()(m, 0);
    }

    template<size_t M>
    const std::string vector<M>::to_string() const {
        std::string print;
        for (size_t m = 0; m < M; m++) {
            if (m != 0) print += ' ';
            print += std::to_string(operator()(m));
        }
        return print;
    }

    template<size_t M>
    float vector<M>::norm_square() {
        float sum = 0;
        for(size_t m = 0; m < M; m++)
            sum += operator()(m) * operator()(m);
        return sum;
    }

    template<size_t M>
    const vector<M> vector<M>::normalized() {
        return matrix<M, 1>::operator*(fast_inv_sqrt(norm_square()));
    }
}

#endif //BIGCAT_BLAS_VECTOR_TCC
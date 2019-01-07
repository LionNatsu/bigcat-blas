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
#include <type_traits>

#include <cstring>

namespace blas {

    template<size_t M>
    vector<M>::vector(const matrix<M, 1> &rhs) {
        operator=(rhs);
    }

    template<size_t M>
    vector<M>::vector(const quaternion &rhs) {
        auto v = rhs.to_vector();
        if (M == 3) {
            operator()(0) = v(1);
            operator()(1) = v(2);
            operator()(2) = v(3);
        } else if (M == 4) {
            operator()(0) = v(0);
            operator()(1) = v(1);
            operator()(2) = v(2);
            operator()(3) = v(3);
        }
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
    const vector <M> vector<M>::normalized() const {
        return matrix<M, 1>::operator*(fast_inv_sqrt((*this) * (*this)));
    }

    template<size_t M>
    float vector<M>::operator*(const vector &rhs) const {
        float sum = 0;
        for (size_t m = 0; m < M; m++)
            sum += operator()(m) * rhs(m);
        return sum;
    }

    template<size_t M>
    const vector <M> vector<M>::cross(const vector &rhs) const {
        static_assert(M == 3);
        auto u1 = operator()(0), u2 = operator()(1), u3 = operator()(2);
        auto v1 = rhs(0), v2 = rhs(1), v3 = rhs(2);
        return vector<3>{
                u2 * v3 - u3 * v2,
                u3 * v1 - u1 * v3,
                u1 * v2 - u2 * v1
        };
    }

    template<size_t M>
    const vector <M> vector<M>::operator*(float k) const {
        return matrix<M, 1>::operator*(k);
    }

    template<size_t M>
    const vector <M> vector<M>::operator/(float k) const {
        return matrix<M, 1>::operator/(k);
    }

    template<size_t M>
    vector <M> &vector<M>::operator+=(const vector &rhs) {
        matrix<M, 1>::operator+=(rhs);
        return *this;
    }

    template<size_t M>
    vector <M> &vector<M>::operator-=(const vector &rhs) {
        matrix<M, 1>::operator-=(rhs);
        return *this;
    }

    template<size_t M>
    const vector <M> operator*(float k, const vector <M> &vec) {
        return vec * k;
    }
}

#endif //BIGCAT_BLAS_VECTOR_TCC
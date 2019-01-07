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

#ifndef BIGCAT_BLAS_MATRIX_TCC
#define BIGCAT_BLAS_MATRIX_TCC

#include "matrix.h"

#include <cstring>

namespace blas {

    template<size_t M, size_t N>
    matrix<M, N>::matrix(const matrix &rhs) {
        operator=(rhs);
    }

    template<size_t M, size_t N>
    matrix<M, N>::matrix(std::initializer_list<std::initializer_list<float>> list) { operator=(list); }

    template<size_t M, size_t N>
    matrix <M, N> &matrix<M, N>::operator=(const matrix &rhs) {
        copy(rhs);
        return *this;
    }

    template<size_t M, size_t N>
    matrix <M, N> &matrix<M, N>::operator=(std::initializer_list<std::initializer_list<float>> list) {
        size_t m = 0;
        for (auto row : list) {
            size_t n = 0;
            for (auto e : row) {
                operator()(m, n) = e;
                n++;
            }
            m++;
        }
        return *this;
    }

    template<size_t M, size_t N>
    float &matrix<M, N>::operator()(size_t m, size_t n) {
        return data[m][n];
    }

    template<size_t M, size_t N>
    float matrix<M, N>::operator()(size_t m, size_t n) const {
        return data[m][n];
    }

    template<size_t M, size_t N>
    const std::string matrix<M, N>::to_string() const {
        std::string print;
        for (size_t m = 0; m < M; m++) {
            if (m != 0) print += '\n';
            for (size_t n = 0; n < N; n++) {
                if (n != 0) print += ' ';
                print += std::to_string(operator()(m, n));
            }
        }
        return print;
    }

    template<size_t M, size_t N>
    float matrix<M, N>::det() const {
        static_assert(M == N);
        if (M == 1) return operator()(0, 0);
        else if (M == 2) return operator()(0, 0) * operator()(1, 1) - operator()(0, 1) * operator()(1, 0);
        else if (M == 3) {
            auto a1 = operator()(0, 0), a2 = operator()(0, 1), a3 = operator()(0, 2);
            auto b1 = operator()(1, 0), b2 = operator()(1, 1), b3 = operator()(1, 2);
            auto c1 = operator()(2, 0), c2 = operator()(2, 1), c3 = operator()(2, 2);
            return a1 * b2 * c3 + b1 * c2 * a3 + c1 * a2 * b3 - a3 * b2 * c1 - b3 * c2 * a1 - c3 * a2 * b1;
        }
        static_assert(M <= 3);
        return 0;
    }

    template<size_t M, size_t N>
    const matrix <M, N> matrix<M, N>::operator*(float k) const {
        matrix<M, N> result;
        for (size_t m = 0; m < M; m++)
            for (size_t n = 0; n < N; n++)
                result(m, n) = operator()(m, n) * k;
        return result;
    }

    template<size_t M, size_t N>
    const matrix <M, N> matrix<M, N>::operator/(float k) const {
        return operator*(1 / k);
    }

    template<size_t M, size_t N>
    const matrix <M, N> operator*(float k, const matrix <M, N> &mat) {
        return mat * k;
    }

    template<size_t M, size_t N>
    template<size_t R>
    const matrix <M, R> matrix<M, N>::operator*(const matrix <N, R> &rhs) const {
        matrix<M, R> result;
        for (size_t m = 0; m < M; m++)
            for (size_t r = 0; r < R; r++)
                for (size_t n = 0; n < N; n++)
                    result(m, r) += operator()(m, n) * rhs(n, r);
        return result;
    }

    template<size_t M, size_t N>
    const matrix <N, M> matrix<M, N>::trans() const {
        matrix<N, M> result;
        for (size_t m = 0; m < M; m++)
            for (size_t n = 0; n < N; n++)
                result(n, m) = operator()(m, n);
        return result;
    }

    template<size_t M, size_t N>
    void matrix<M, N>::copy(const matrix &rhs) {
        std::memcpy(data, rhs.data, sizeof(data));
    }

    template<size_t M, size_t N>
    const matrix <M, N> matrix<M, N>::operator-() const {
        return operator*(-1);
    }

    template<size_t M, size_t N>
    const matrix <M, N> matrix<M, N>::operator+(const matrix &rhs) const {
        matrix<M, N> result;
        for (size_t m = 0; m < M; m++)
            for (size_t n = 0; n < N; n++)
                result(m, n) = operator()(m, n) + rhs(m, n);
        return result;
    }

    template<size_t M, size_t N>
    const matrix <M, N> matrix<M, N>::operator-(const matrix &rhs) const {
        return operator+(rhs * -1);
    }

    template<size_t M, size_t N>
    matrix <M, N> &matrix<M, N>::operator+=(const matrix &rhs) {
        for (size_t m = 0; m < M; m++)
            for (size_t n = 0; n < N; n++)
                operator()(m, n) += rhs(m, n);
        return *this;
    }

    template<size_t M, size_t N>
    matrix <M, N> &matrix<M, N>::operator-=(const matrix &rhs) {
        for (size_t m = 0; m < M; m++)
            for (size_t n = 0; n < N; n++)
                operator()(m, n) -= rhs(m, n);
        return *this;
    }
}

#endif //BIGCAT_BLAS_MATRIX_TCC
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

#ifndef BIGCAT_BLAS_MATRIX_H
#define BIGCAT_BLAS_MATRIX_H

#include <initializer_list>
#include <string>

namespace blas {
    template<size_t M, size_t N = M>
    class matrix {
    public:
        matrix() = default;

        matrix(const matrix &);

        matrix(std::initializer_list<std::initializer_list<float>>);

        matrix &operator=(std::initializer_list<std::initializer_list<float>>);

        matrix &operator=(const matrix &);

        // Accessing Elements

        float &operator()(size_t, size_t);

        float operator()(size_t, size_t) const;

        // Arithmetic

        const matrix operator-() const;

        const matrix operator+(const matrix &) const;

        matrix &operator+=(const matrix &);

        matrix &operator-=(const matrix &);

        const matrix operator-(const matrix &) const;

        const matrix operator*(float) const;

        const matrix operator/(float) const;

        template<size_t R>
        const matrix<M, R> operator*(const matrix<N, R> &) const;

        const matrix<N, M> trans() const;

        // Square Matrices

        float det() const;

        const matrix adj() const;

        const matrix inv() const;

        // Debugging

        virtual const std::string to_string() const;

    private:
        void copy(const matrix &);

    private:
        float data[M][N] = {};
    };

    template<size_t M, size_t N>
    const matrix<M, N> operator*(float, const matrix<M, N> &);

}

#endif //BIGCAT_BLAS_MATRIX_H

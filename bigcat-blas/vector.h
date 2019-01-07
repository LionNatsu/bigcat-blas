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

#ifndef BIGCAT_BLAS_VECTOR_H
#define BIGCAT_BLAS_VECTOR_H

namespace blas {
    template<size_t M>
    class vector : public matrix<M, 1> {
    public:
        vector() = default;

        vector(const quaternion &);

        vector(const matrix<M, 1> &);

        vector(std::initializer_list<float>);

        vector &operator=(std::initializer_list<float>);

        vector &operator=(const matrix<M, 1> &);

        // Accessing Elements

        float &operator()(size_t);

        float operator()(size_t) const;

        // Vectors

        vector &operator+=(const vector &);

        vector &operator-=(const vector &);

        const vector operator*(float) const;

        const vector operator/(float) const;

        float operator*(const vector &) const;

        const vector cross(const vector &) const;

        const vector normalized() const;

        // Debugging

        const std::string to_string() const;
    };

    template<size_t M>
    const vector<M> operator*(float, const vector<M> &);
}

#endif //BIGCAT_BLAS_VECTOR_H

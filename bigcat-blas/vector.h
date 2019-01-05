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

        vector(const matrix<M, 1> &);

        vector(std::initializer_list<float> list);

        vector &operator=(std::initializer_list<float> list);

        vector &operator=(const matrix<M, 1> &);

        // Accessing Elements

        float &operator()(size_t m);

        float operator()(size_t m) const;

        // Vectors

        float norm_square();

        const vector normalized();

        // Debugging

        const std::string to_string() const;
    };
}

#endif //BIGCAT_BLAS_VECTOR_H

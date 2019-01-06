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

#ifndef BIGCAT_BLAS_H
#define BIGCAT_BLAS_H

namespace blas {
    template<size_t M, size_t N>
    class matrix;

    template<size_t M>
    class vector;

    class quaternion;
}

#include "bigcat-blas/utils.h"
#include "bigcat-blas/matrix.h"
#include "bigcat-blas/vector.h"
#include "bigcat-blas/quaternion.h"

#include "bigcat-blas/matrix.tcc"
#include "bigcat-blas/vector.tcc"
#include "bigcat-blas/quaternion.tcc"

#endif
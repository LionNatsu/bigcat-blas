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

#ifndef BIGCAT_BLAS_QUATERNION_TCC
#define BIGCAT_BLAS_QUATERNION_TCC

#include "matrix.h"
#include "vector.h"
#include "quaternion.h"

#include <cmath>
#include <cstring>

namespace blas {

    quaternion::quaternion(std::initializer_list<float> list) { operator=(list); }

    quaternion &quaternion::operator=(std::initializer_list<float> list) {
        if (list.size() == 4)
            data = list;
        else if (list.size() == 3) {
            auto iter = list.begin();
            data(0) = 0;
            data(1) = *iter++;
            data(2) = *iter++;
            data(3) = *iter;
        }
        return *this;
    }

    float &quaternion::a() {
        return data(0);
    }

    float quaternion::a() const {
        return data(0);
    }

    float &quaternion::b() {
        return data(1);
    }

    float quaternion::b() const {
        return data(1);
    }

    float &quaternion::c() {
        return data(2);
    }

    float quaternion::c() const {
        return data(2);
    }

    float &quaternion::d() {
        return data(3);
    }

    float quaternion::d() const {
        return data(3);
    }

    const std::string quaternion::to_string() const {
        return data.to_string();
    }

    quaternion &quaternion::operator=(const vector<4> &vec) {
        data = vec;
        return *this;
    }

    quaternion &quaternion::operator=(const vector<3> &vec) {
        data(0) = 0;
        data(1) = vec(0);
        data(2) = vec(1);
        data(3) = vec(2);
        return *this;
    }

    quaternion::quaternion(const vector<4> &vec) { operator=(vec); }

    quaternion::quaternion(const vector<3> &vec) { operator=(vec); }

    const quaternion quaternion::normalized() const {
        return data.normalized();
    }

    const quaternion quaternion::operator*(float k) const {
        return data * k;
    }

    const quaternion operator*(float k, const quaternion &rhs) {
        return rhs * k;
    }

    const quaternion quaternion::operator/(float k) const {
        return data / k;
    }

    const quaternion quaternion::operator*(const quaternion &rhs) const {
        auto a = this->a(), b = this->b(), c = this->c(), d = this->d();
        auto e = rhs.a(), f = rhs.b(), g = rhs.c(), h = rhs.d();
        return quaternion({
                                  a * e - b * f - c * g - d * h,
                                  b * e + a * f - d * g + c * h,
                                  c * e + d * f + a * g - b * h,
                                  d * e - c * f + b * g + a * h
                          });
    }

    const quaternion quaternion::conj() const {
        return quaternion({a(), -b(), -c(), -d()});
    }

    const quaternion quaternion::rotate(const quaternion &u) const {
        return (*this) * u * (conj());
    }

    const quaternion quaternion::rotate_inv(const quaternion &u) const {
        return (conj()) * u * (*this);
    }

    const vector<4> quaternion::to_vector() const {
        return {a(), b(), c(), d()};
    }

    const vector<3> quaternion::to_vector3() const {
        return {b(), c(), d()};
    }

    const quaternion quaternion::operator+(const quaternion &rhs) const {
        return vector<4>(data + rhs.data);
    }

    quaternion &quaternion::operator+=(const quaternion &rhs) {
        data += rhs;
        return *this;
    }

    const matrix<4> quaternion::rotate_inv_jacobian(const quaternion &u) const {
        auto dxq1 = u.b() * a(), dxq2 = u.b() * b(), dxq3 = u.b() * c(), dxq4 = u.b() * d();
        auto dyq1 = u.c() * a(), dyq2 = u.c() * b(), dyq3 = u.c() * c(), dyq4 = u.c() * d();
        auto dzq1 = u.d() * a(), dzq2 = u.d() * b(), dzq3 = u.d() * c(), dzq4 = u.d() * d();
        return matrix<4>{
                {0,            0,                      0,                       0},
                {dyq4 - dzq3,  dyq3 + dzq4,            -2 * dxq3 + dyq2 - dzq1, -2 * dxq4 + dyq1 + dzq2},
                {-dxq4 + dzq2, dxq3 - 2 * dyq2 + dzq1, dxq2 + dzq4,             -dxq1 - 2 * dyq4 + dzq3},
                {dxq3 - dyq2,  dxq4 - dyq1 - 2 * dzq2, dxq1 + dyq4 - 2 * dzq3,  dxq2 + dyq3}
        } * 2;
    }

    const quaternion quaternion::operator-(const quaternion &rhs) const {
        return vector<4>(data - rhs.data);
    }

    quaternion &quaternion::operator-=(const quaternion &rhs) {
        data -= rhs;
        return *this;
    }

    const vector<3> quaternion::to_euler_angles() const {
        const float ARC2DEG = 180 / PI;
        return {
                ARC2DEG * atan2f(b() * c() + d() * a(), 0.5f - c() * c() - d() * d()),
                ARC2DEG * asinf(-2.0f * (b() * d() - a() * c())),
                ARC2DEG * atan2f(a() * b() + c() * d(), 0.5f - b() * b() - c() * c())
        };
    }
}

#endif //BIGCAT_BLAS_QUATERNION_TCC
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

#ifndef BIGCAT_BLAS_QUATERNION_H
#define BIGCAT_BLAS_QUATERNION_H

namespace blas {
    class quaternion {
    public:
        quaternion() = default;

        quaternion(const vector<3> &);

        quaternion &operator=(const vector<3> &);

        quaternion(const vector<4> &);

        quaternion &operator=(const vector<4> &);

        quaternion(std::initializer_list<float>);

        quaternion &operator=(std::initializer_list<float>);

        // Accessing Elements

        float &a();

        float a() const;

        float &b();

        float b() const;

        float &c();

        float c() const;

        float &d();

        float d() const;

        // Arithmetic

        const quaternion operator*(float) const;

        const quaternion operator/(float) const;

        const quaternion operator*(const quaternion &) const;

        const quaternion operator+(const quaternion &) const;

        const quaternion operator-(const quaternion &) const;

        quaternion &operator+=(const quaternion &);

        quaternion &operator-=(const quaternion &);

        const quaternion rotate(const quaternion &) const;

        const quaternion rotate_inv(const quaternion &) const;

        const matrix<4> rotate_inv_jacobian(const quaternion &) const;

        const quaternion conj() const;

        const quaternion normalized() const;

        // Debugging

        const std::string to_string() const;

        const vector<4> to_vector() const;

        const vector<3> to_vector3() const;

        const vector<3> to_euler_angles() const;

    private:
        vector<4> data;
    };

    const quaternion operator*(float, const quaternion &);
}

#endif //BIGCAT_BLAS_QUATERNION_H

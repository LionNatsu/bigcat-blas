/*
    This file is an example and a part of bigcat-blas.
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

#include <iostream>
#include <cmath>
#include <random>
#include "bigcat-blas.h"

namespace physical_model {
    const blas::vector<3> accelerometer();

    const blas::vector<3> magnetometer();

    const blas::vector<3> gyroscope();
}

void Madgwick_update(blas::quaternion &q);

void Mahony_update(blas::quaternion &q);

const float duration = 0.01;

int main() {
    int tick = 0;
    blas::quaternion q1 = {1, 0, 0, 0};
    blas::quaternion q2 = {1, 0, 0, 0};
    for (;;) {
        if (tick++ == 100) break;
        Mahony_update(q1);
        Madgwick_update(q2);
        std::cout << q1.to_euler_angles().to_string() << std::endl;
        std::cout << q2.to_euler_angles().to_string() << '\t' << tick * duration << 's' << std::endl << std::endl;
    }
}

/*
 * Algorithms are based on Madgwick's great work.
 * See https://x-io.co.uk/open-source-imu-and-ahrs-algorithms/
 */

void Madgwick_update(blas::quaternion &q) {
    using blas::quaternion;
    using blas::vector;
    using blas::matrix;

    const float beta = 1;

    const vector<3> acc = physical_model::accelerometer().normalized();
    // Assume the real acceleration is pointing to z+,
    // i.e. there is no accelerated motion.
    const vector<3> acc_expect = q.rotate_inv({0, 0, 1});
    const quaternion acc_error = vector<3>(acc_expect - acc);
    const matrix<4> acc_jacobian = q.rotate_inv_jacobian({0, 0, 1});

    const vector<3> mag = physical_model::magnetometer().normalized();
    const vector<3> mag_0 = q.rotate(mag);
    // Decompose into x-y plane(horizontal) and z-axis(vertical).
    const float mag_h = sqrtf(mag_0(0) * mag_0(0) + mag_0(1) * mag_0(1));
    const float mag_v = mag_0(2);
    // Assume the real direction of magnetic field is in surface y=0, pointing to x+,
    // i.e. there is no electromagnetic interference.
    const vector<3> mag_expect = q.rotate_inv({mag_h, 0, mag_v});
    const quaternion mag_error = vector<3>(mag_expect - mag);
    const matrix<4> mag_jacobian = q.rotate_inv_jacobian({mag_h, 0, mag_v});

    const quaternion gradient = vector<4>(
            acc_jacobian.trans() * acc_error.to_vector() +
            mag_jacobian.trans() * mag_error.to_vector()
    ).normalized();

    vector<3> gyr = physical_model::gyroscope();

    quaternion q_new = 0.5 * q * gyr;
    q_new -= beta * gradient;
    // Approximate the integral
    q += q_new * duration;
    q = q.normalized();
}

void Mahony_update(blas::quaternion &q) {
    using blas::quaternion;
    using blas::vector;
    using blas::matrix;

    const float Kp = 120.0, Ki = 0.02;
    static vector<3> integrated_error;

    const vector<3> acc = physical_model::accelerometer().normalized();
    // Assume the real acceleration is pointing to z+,
    // i.e. there is no accelerated motion.
    const vector<3> acc_expect = q.rotate_inv({0, 0, 1});

    const vector<3> mag = physical_model::magnetometer().normalized();
    const vector<3> mag_0 = q.rotate(mag);
    // Decompose into x-y plane(horizontal) and z-axis(vertical).
    const float mag_h = sqrtf(mag_0(0) * mag_0(0) + mag_0(1) * mag_0(1));
    const float mag_v = mag_0(2);
    // Assume the real direction of magnetic field is in surface y=0, pointing to x+,
    // i.e. there is no electromagnetic interference.
    const vector<3> mag_expect = q.rotate_inv({mag_h, 0, mag_v});

    const vector<3> error = acc.cross(acc_expect) + mag.cross(mag_expect);

    vector<3> gyr = physical_model::gyroscope();
    integrated_error += Ki * error * duration;
    gyr += integrated_error + Kp * error;

    // Approximate the integral
    q += 0.5 * q * gyr * duration;
    q = q.normalized();
}

namespace physical_model {
    std::random_device rd{};
    std::mt19937 gen{rd()};
    blas::vector<3> acc_real = {1, 0.1, 1};
    blas::vector<3> mag_real = {100, 100, 0};
    blas::vector<3> gyr_real = {0, 0, 0};

    std::normal_distribution<float> acc_noise_x{0, 0.01}, acc_noise_y{0, 0.01}, acc_noise_z{0, 0.01};
    std::normal_distribution<float> mag_noise_x{0, 1}, mag_noise_y{0, 0.1}, mag_noise_z{0, 1};
    std::normal_distribution<float> gyr_noise_x{0, 0.01}, gyr_noise_y{0, 0.01}, gyr_noise_z{0, 0.01};
    const blas::vector<3> gyr_bias = {0.01, 0.03, 0.5};

    const blas::vector<3> accelerometer() {
        return blas::vector<3>{acc_noise_x(gen), acc_noise_y(gen), acc_noise_z(gen)} + acc_real;
    }

    const blas::vector<3> magnetometer() {
        return blas::vector<3>{mag_noise_x(gen), mag_noise_y(gen), mag_noise_z(gen)} + mag_real;
    }

    const blas::vector<3> gyroscope() {
        static blas::vector<3> random_walk;
        random_walk += {gyr_noise_x(gen), gyr_noise_y(gen), gyr_noise_z(gen)};
        return gyr_real + gyr_bias + random_walk;
    }
}
#include <iostream>
#include <cmath>
#include <random>
#include "bigcat-blas.h"

const blas::vector<3> euler_angles(const blas::quaternion &q);

namespace physical_model {
    const blas::vector<3> accelerometer();

    const blas::vector<3> magnetometer();

    const blas::vector<3> gyroscope();
}

int main() {
    using blas::quaternion;
    using blas::vector;
    using blas::matrix;

    const float duration = 0.01;
    int tick = 0;
    quaternion q = {1, 0, 0, 0};

    /*
     * See http://x-io.co.uk/open-source-imu-and-ahrs-algorithms/
     */
    const float beta = 0.1;
    vector<3> integrated_error;
    for (;;) {
        if (tick++ == 200) break;
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

        // Approximate the integral
        quaternion q_new = 0.5 * q * gyr;
        q_new -= beta * gradient;
        q += q_new * duration;
        q = q.normalized();
        std::cout << euler_angles(q).to_string() << '\t' << tick * duration << 's' << std::endl;
    }
}

const blas::vector<3> euler_angles(const blas::quaternion &q) {
    const float PI = 3.141593;
    const float ARC2DEG = 180 / PI;
    return {
            ARC2DEG * atan2f(q.b() * q.c() + q.d() * q.a(), 0.5f - q.c() * q.c() - q.d() * q.d()),
            ARC2DEG * asinf(-2.0f * (q.b() * q.d() - q.a() * q.c())),
            ARC2DEG * atan2f(q.a() * q.b() + q.c() * q.d(), 0.5f - q.b() * q.b() - q.c() * q.c())
    };
}

namespace physical_model {
    std::random_device rd{};
    std::mt19937 gen{rd()};
    blas::vector<3> acc_real = {1.0, 2.0, 9.8};
    blas::vector<3> mag_real = {100.0, 100.0, 0.5};
    blas::vector<3> gyr_real = {0, 0, 0};

    std::normal_distribution<float> acc_noise_x{0, 0.01}, acc_noise_y{0, 0.01}, acc_noise_z{0, 0.01};
    std::normal_distribution<float> mag_noise_x{0, 0.01}, mag_noise_y{0, 0.01}, mag_noise_z{0, 0.01};
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
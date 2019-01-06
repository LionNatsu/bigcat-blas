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

    const float
            duration = 0.01,
            Ki = 0.02,
            Kp = 120.0;

    quaternion q = {1, 0, 0, 0};
    vector<3> integrated_error;

    for (;;) {
        static int count = 0;
        if (count++ == 100) break;

        const vector<3> acc = physical_model::accelerometer().normalized();
        const vector<3> acc_expect = q.rotate_inv({0, 0, 1});

        const vector<3> mag = physical_model::magnetometer().normalized();
        const vector<3> mag_zero = q.rotate(mag);
        const float mag_h = sqrtf(mag_zero(0) * mag_zero(0) + mag_zero(1) * mag_zero(1));
        const float mag_v = mag_zero(2);
        const vector<3> mag_expect = q.rotate_inv({mag_h, 0, mag_v});

        const vector<3> error = acc.cross(acc_expect) + mag.cross(mag_expect);

        vector<3> gyr = physical_model::gyroscope();
        integrated_error += Ki * error * duration;
        gyr += integrated_error + Kp * error;

        q += 1 / 2.0 * q * gyr * duration;
        q = q.normalized();

        std::cout << euler_angles(q).to_string() << '\t' << count * duration << 's' << std::endl;
    }
    return q.a() == 0;
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
        static blas::vector<3> random_walk;  // bias
        random_walk += {gyr_noise_x(gen), gyr_noise_y(gen), gyr_noise_z(gen)};
        return gyr_real + gyr_bias + random_walk;
    }
}
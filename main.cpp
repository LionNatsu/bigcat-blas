#include <iostream>
#include "bigcat-blas.h"

int main() {
    const blas::matrix<3> A = {{1, 2, 3}, {4, 5, 6}, {7, 0, 9}};
    const blas::matrix<3> B = A;
    const blas::matrix<3, 1> C {{1}, {2}, {3}};
    blas::vector<3> v {1, 2, 3}, q;
    std::cout << "Hello, World! " << (A.det() * 1.001)  << std::endl;
    std::cout  << v.to_string() << std::endl;
    std::cout << blas::vector<3>( (A * v) ).to_string() << std::endl;
    std::cout << blas::vector<3>( q + (A * v) ).normalized().to_string() << std::endl;
    return 0;
}

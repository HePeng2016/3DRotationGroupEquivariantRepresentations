#include <cmath>
#include <complex>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <algorithm>

class WignerD {
public:
    using Complex = std::complex<double>;
    using MatrixC = std::vector<std::vector<Complex>>;

    explicit WignerD(int maxL) : maxL_(maxL) {
        if (maxL_ < 0) {
            throw std::invalid_argument("maxL must be non-negative.");
        }

        // Need factorials up to 2*maxL.
        logFact_.resize(2 * maxL_ + 1);
        logFact_[0] = 0.0;
        for (int i = 1; i <= 2 * maxL_; ++i) {
            logFact_[i] = logFact_[i - 1] + std::log(static_cast<double>(i));
        }
    }

    // Compute full Wigner D^l(alpha, beta, gamma)
    // Matrix indices are shifted:
    //   row = mp + l, col = m + l
    // where mp,m in [-l,l].
    MatrixC D(int l, double alpha, double beta, double gamma) const {
        check_l(l);

        const int dim = 2 * l + 1;
        MatrixC mat(dim, std::vector<Complex>(dim));

        for (int mp = -l; mp <= l; ++mp) {
            Complex leftPhase = std::exp(Complex(0.0, -mp * alpha));

            for (int m = -l; m <= l; ++m) {
                Complex rightPhase = std::exp(Complex(0.0, -m * gamma));
                double dval = small_d(l, mp, m, beta);

                mat[mp + l][m + l] = leftPhase * dval * rightPhase;
            }
        }

        return mat;
    }

    // Compute reduced Wigner d^l(beta)
    // Matrix indices:
    //   row = mp + l, col = m + l
    MatrixC small_d_matrix(int l, double beta) const {
        check_l(l);

        const int dim = 2 * l + 1;
        MatrixC mat(dim, std::vector<Complex>(dim));

        for (int mp = -l; mp <= l; ++mp) {
            for (int m = -l; m <= l; ++m) {
                mat[mp + l][m + l] = Complex(small_d(l, mp, m, beta), 0.0);
            }
        }

        return mat;
    }

    // Compute one reduced Wigner d element d^l_{mp,m}(beta)
    double small_d(int l, int mp, int m, double beta) const {
        check_l(l);

        if (std::abs(m) > l || std::abs(mp) > l) {
            return 0.0;
        }

        // Formula:
        //
        // d^l_{mp,m}(beta)
        // = sum_k (-1)^{k - mp + m}
        //   sqrt((l+m)!(l-m)!(l+mp)!(l-mp)!)
        //   /[(l+m-k)! k! (mp-m+k)! (l-mp-k)!]
        //   * cos(beta/2)^{2l + m - mp - 2k}
        //   * sin(beta/2)^{mp - m + 2k}
        //
        // k range chosen so all factorials are nonnegative.

        const double cb = std::cos(0.5 * beta);
        const double sb = std::sin(0.5 * beta);

        const int kMin = std::max(0, m - mp);
        const int kMax = std::min(l + m, l - mp);

        const double logPref =
            0.5 * (
                logFact_[l + m] +
                logFact_[l - m] +
                logFact_[l + mp] +
                logFact_[l - mp]
            );

        double sum = 0.0;

        for (int k = kMin; k <= kMax; ++k) {
            const int a = l + m - k;
            const int b = k;
            const int c = mp - m + k;
            const int d = l - mp - k;

            if (a < 0 || b < 0 || c < 0 || d < 0) {
                continue;
            }

            const int powC = 2 * l + m - mp - 2 * k;
            const int powS = mp - m + 2 * k;

            if (powC < 0 || powS < 0) {
                continue;
            }

            double logTerm =
                logPref
                - logFact_[a]
                - logFact_[b]
                - logFact_[c]
                - logFact_[d];

            double term = std::exp(logTerm);

            // Handle powers robustly.
            term *= int_pow(cb, powC);
            term *= int_pow(sb, powS);

            // Sign factor: (-1)^(k - mp + m)
            if (((k - mp + m) & 1) != 0) {
                term = -term;
            }

            sum += term;
        }

        return sum;
    }

    // Rotation that aligns a direction n=(x,y,z) with the z-axis.
    //
    // If n has spherical angles theta, phi, then
    //   R_n = Rz(phi) Ry(theta)
    // maps z to n.
    //
    // The alignment rotation is
    //   R_align = R_n^{-1}.
    //
    // Therefore,
    //   D_align = D(phi, theta, 0)^dagger. 
    // Active rotation that rotates the object so that direction n aligns with z-axis.
    // R_n = Rz(phi) Ry(theta) maps z -> n
    // Therefore the alignment rotation is R_align = R_n^{-1}
    // The corresponding Wigner D matrix is D(R_align) = D(phi, theta, 0)^dagger

    MatrixC alignment_D(int l, double x, double y, double z) const {
        check_l(l);

        const double norm = std::sqrt(x * x + y * y + z * z);
        if (norm == 0.0) {
            throw std::invalid_argument("Direction vector must be nonzero.");
        }

        x /= norm;
        y /= norm;
        z /= norm;

        // Clamp z to avoid numerical errors outside [-1,1].
        z = std::max(-1.0, std::min(1.0, z));

        const double theta = std::acos(z);
        const double phi = std::atan2(y, x);

        MatrixC D_forward = D(l, phi, theta, 0.0);
        return dagger(D_forward);
    }

    static MatrixC dagger(const MatrixC& A) {
        const int rows = static_cast<int>(A.size());
        const int cols = rows > 0 ? static_cast<int>(A[0].size()) : 0;

        MatrixC B(cols, std::vector<Complex>(rows));

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                B[j][i] = std::conj(A[i][j]);
            }
        }

        return B;
    }

private:
    int maxL_;
    std::vector<double> logFact_;

    void check_l(int l) const {
        if (l < 0 || l > maxL_) {
            throw std::invalid_argument("l is outside the precomputed range.");
        }
    }

    static double int_pow(double x, int n) {
        if (n == 0) {
            return 1.0;
        }

        // std::pow is fine, but repeated squaring is faster for integer powers.
        double result = 1.0;
        double base = x;
        int exp = n;

        while (exp > 0) {
            if (exp & 1) {
                result *= base;
            }
            base *= base;
            exp >>= 1;
        }

        return result;
    }
};









/*nt main() {
    const int maxL = 10;
    WignerD wd(maxL);

    int l = 2;

    // General Euler angles, z-y-z convention.
    double alpha = 0.3;
    double beta  = 1.0;
    double gamma = -0.7;

    auto D = wd.D(l, alpha, beta, gamma);

    std::cout << "D^" << l << " matrix:\n";
    for (const auto& row : D) {
        for (const auto& val : row) {
            std::cout << "(" << val.real() << "," << val.imag() << ") ";
        }
        std::cout << "\n";
    }

    // Alignment example:
    // Direction n = (1,1,1), rotate so that it aligns with z-axis.
    auto D_align = wd.alignment_D(l, 1.0, 1.0, 1.0);

    std::cout << "\nAlignment D matrix:\n";
    for (const auto& row : D_align) {
        for (const auto& val : row) {
            std::cout << "(" << val.real() << "," << val.imag() << ") ";
        }
        std::cout << "\n";
    }

    return 0;
}
*/


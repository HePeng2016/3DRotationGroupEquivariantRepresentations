#include <Eigen/Dense>
#include <Eigen/QR>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>
#include <iostream>
#include "stdio.h"
#include "tuple"
#include "string"
#include "wigner/wigner_3nj.hpp"
#define  cdouble complex<double>
struct Index_t {
    int16_t m1;
    int16_t m2;
};


class equivalentFeatures
{
    private: int N = 3;
    public:
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> CGTableC;
    std::vector<std::vector<std::vector<std::vector<std::vector<Index_t>>>>> CGTableI;
    std::vector<std::vector<std::vector<std::vector<double>>>> W3JTableC;
    std::vector<std::vector<std::vector<std::vector<Index_t>>>> W3JTableI;
    std::vector<std::vector<uint16_t>> WignerPI;
    std::vector<std::vector<uint16_t>> SelfPI;

    void Initial();
    std::vector<double> CStoRS_Encode(const std::vector<std::complex<double>>& V1, int n);
    std::vector<std::complex<double>> RStoCS_Encode(const std::vector<double>& V1, int n);
    std::vector<std::complex<double>> IncreaseDegree(const std::vector<std::complex<double>>& V, int n);
    std::vector<double> SelfProduct(const std::vector<std::complex<double>>& V, int n, int n2);
    std::vector<std::complex<double>> W3jProduct(const std::vector<std::complex<double>>& V1,
    const std::vector<std::complex<double>>& V2,const std::vector<std::complex<double>>& V3,int n1, int n2);
    std::vector<double> W3jProductCToR(const std::vector<std::complex<double>>& InvariantV, int n1,int d2, int d3, int n2);
    std::vector<std::complex<double>> W3jProductRToC(const std::vector<double>& InvariantV, int n1,int d2, int d3, int n2);
    std::vector<std::complex<double>> W3jProductCompact(const std::vector<std::complex<double>>& V1,
    const std::vector<std::complex<double>>& V2,const std::vector<std::complex<double>>& V3,int n);
    Eigen::MatrixXcd DecodeMatrixCompact(const std::vector<std::complex<double>>& V2, const std::vector<std::complex<double>>& V3, int n);
    std::vector<double> W3jProductCompactCToR(const std::vector<std::complex<double>>& InvariantV, int n);
    std::vector<std::complex<double>> W3jProductCompactRToC(const std::vector<double>& InvariantV, int n);
    std::vector<double> SelfProductPairwise(const std::vector<std::complex<double>>& V, int n, int n2);
    std::vector< std::vector<double> > ReferencesExtract(const  std::vector<double>& VR);
    std::vector<double> ProductEncode(const std::vector<std::complex<double>>& V1, const std::vector<std::complex<double>>& V2, int n);
    std::vector<double> ProductEncodePairwise(const std::vector<std::complex<double>>& V1, const std::vector<std::complex<double>>& V2, int n);
    std::vector<std::vector<std::complex<double>>> SelfProductMatrix(const std::vector<std::complex<double>>& V, int n, int n2);
    std::vector<std::vector<std::complex<double>>> DerivativeSH(const std::vector<std::complex<double>>& Y);
    std::vector<std::vector<std::complex<double>>> DerivativeSH_XYZ(const std::vector<std::complex<double>>& Y, const std::vector<std::complex<double>>& DR, double ReciprocalRadii, double ReciprocalF);
};







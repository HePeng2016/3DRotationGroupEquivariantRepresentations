#include "rotation.h"


using namespace std;
using namespace Eigen;





void equivalentFeatures::Initial() {

    int d1 = N;
    int d2 = N;
    int d3 = std::abs(d1 + d2 - 2) - std::abs(d1 - d2) + 1;
    wigner::Wigner3jSeriesJ<double, int> w3j;

    CGTableC.resize(d1, std::vector<std::vector<std::vector<std::vector<double>>>>(d2, std::vector<std::vector<std::vector<double>>>(d3)));
    CGTableI.resize(d1, std::vector<std::vector<std::vector<std::vector<Index_t>>>>(d2, std::vector<std::vector<std::vector<Index_t>>>(d3)));

    d3 = N;
    W3JTableC.resize(d1, std::vector<std::vector<std::vector<double>>>(d2, std::vector<std::vector<double>>(d3)));
    W3JTableI.resize(d1, std::vector<std::vector<std::vector<Index_t>>>(d2, std::vector<std::vector<Index_t>>(d3)));


    for (int J1 = 0; J1 < d1; ++J1) {
        for (int J2 = J1; J2 < d2; ++J2) {
            for (int j3 = (J2 - J1); j3 <= (J1 + J2); ++j3) {
                int j1 = J1;
                int j2 = J2;
                int J = j3 - (J2 - J1);


                if (!((j2==j1)&&((j3 % 2)==1)))
                {
                    CGTableC[J1][J2][J].resize(2 * j3 + 1);
                    CGTableI[J1][J2][J].resize(2 * j3 + 1);

                    for (int i = 0; i < (2 * j3 + 1); ++i) {
                       CGTableC[J1][J2][J][i].clear();
                       CGTableI[J1][J2][J][i].clear();
                    }

                }
            }
        }
    }



    for (int J1 = 1; J1 <= d1; ++J1) {
      for (int J2 = 1; J2 <= d2; ++J2) {

        int j1 = J1 - 1;
        int j2 = J2 - 1;

        for (int m1 = -j1; m1 <= j1; ++m1) {
          for (int m2 = -j2; m2 <= j2; ++m2) {



            w3j.compute(j1,j2,m1,m2);

            for (int j3 = w3j.nmin();j3 <= w3j.nmax(); ++j3) {

                double w3j_result = w3j.get(j3);

                if (j3 < d3) {
                            int m3 = -(m1 + m2);



                            if (std::abs(m3) <= j3) {
                                W3JTableI[j1][j2][j3].push_back(Index_t{m1,m2});
                                W3JTableC[j1][j2][j3].push_back(w3j_result);
                            }
                        }
                if(j2>=j1) {
                            int m3 = m2 + m1;
                            int J  = j3 - (j2 - j1);

                             if (abs(m3) <= j3)
                               {
                                  if (!((j2 == j1)&&(j3 % 2 == 1))) {
                                     int i = m3 + j3;
                                     CGTableI[j1][j2][J][i].push_back(Index_t{m1, m2});
                                     double clebschgorda = std::sqrt(2 * j3 + 1)*std::pow(-1,j1 - j2 - m3)*w3j_result;//wigner3J coefficient to Clebsh-Gordan coefficient
                                     CGTableC[j1][j2][J][i].push_back(clebschgorda); // Placeholder for clebschgordan
                                    }

                               }
                           }
                    }
                }
              }
           }
      }

    for (int J1 = 0; J1 < d1; ++J1) {
        for (int J2 = 0; J2 < d2; ++J2) {
            for (int J3 = std::abs(J2 - J1); J3 <= std::abs(J1 + J2); ++J3) {
                if (J3 < d3) {
                    if (W3JTableI[J1][J2][J3].size()!=0) {
                        WignerPI.push_back({J1 + 1, J2 + 1, J3+1}); // Adjusting for 1-based indexing
                    }
                }
            }
        }
    }
    for (int J1 = 0; J1 < d1; ++J1) {
        for (int J2 = J1; J2 < d2; ++J2) {
            for (int J3 = std::abs(J2 - J1); J3 <= std::abs(J1 + J2); ++J3) {
                if  (!((J2 == J1)&&(J3 % 2 == 1))) {
                    int J = J3 - (J2 - J1);
                    if (CGTableI[J1][J2][J].size()!=0) {
                        SelfPI.push_back({J1 + 1, J2 + 1, J+1}); // Adjusting for 1-based indexing
                    }
                }
            }
        }
    }


};



//Complex spherical harmonics tensor to real spherical harmonics tensor

std::vector<double> equivalentFeatures::CStoRS_Encode(const std::vector<std::complex<double>>& V1, int n) {
    std::vector<double> V3(n*n, 0.0);

    for (int I = 1; I <= n; ++I) {
        int Start = (I - 1) * (I - 1);
        for (int I2 = 1; I2 <= I; ++I2) {
            if (I2 == I) {
                V3[Start + I2 - 1] = std::real(V1[Start + I2 - 1]);
            } else {
                V3[Start + I2 - 1] = std::pow(-1, I - I2) * std::imag(V1[Start + 2 * I - I2 - 1]) * std::sqrt(2);
                V3[Start + 2 * I - I2 - 1] = std::pow(-1, I - I2) * std::real(V1[Start + 2 * I - I2 - 1]) * std::sqrt(2);
            }
        }
    }

    return V3;
}



//Real spherical harmonics tensor to Complex spherical harmonics tensor
std::vector<std::complex<double>> equivalentFeatures::RStoCS_Encode(const std::vector<double>& V1, int n) {
    std::vector<std::complex<double>> V3(n*n, std::complex<double>(0, 0));

    for (int I = 1; I <= n; ++I) {
        int Start = (I - 1) * (I - 1);
        for (int I2 = 1; I2 <= I; ++I2) {
            if (I2 == I) {
                V3[Start + I2 - 1] = V1[Start + I2 - 1];
            } else {
                V3[Start + I2 - 1] = (-V1[Start + I2 - 1] * std::complex<double>(0, 1) + V1[Start + 2 * I - I2 - 1]) / std::sqrt(2.0);
                V3[Start + 2 * I - I2 - 1] = std::pow(-1, I - I2) * (V1[Start + I2 - 1] * std::complex<double>(0, 1) + V1[Start + 2 * I - I2 - 1]) / std::sqrt(2.0);
            }
        }
    }

    return V3;
}


// length of V >=4
// According to the recursion of Associated Legendre Functions
// P_n+1 = ((2n+1)/(n+1)) * (sqrt((2n+3)/(2n+1))Cos*P_n - (n/n+1)*sqrt((2n+3)/(2n-1))P_n-1
// P_n+1^(m+1) = 2CosP_n^(m+1)(sqrt((2*n+3)/(2*n+1))sqrt((n-m)/(n+m+2))  ) - P_n-1^(m+1)sqrt((2*n+3)/(2*n-1))sqrt(((n-m)(n-m-1))/((n+m+2)(n+m+1)))+ (2m+1)P_n^m*Sin*(-1)^m*sqrt((2*n+3)/(2*n+1))*sqrt(1/((n+m+2)(n+m+1)))

  std::vector<std::complex<double>> equivalentFeatures::IncreaseDegree(const std::vector<std::complex<double>>& V, int n) {
    int N = static_cast<int>(std::sqrt(V.size()));
    int N_new = N + n;
    int J0 = 1;

    std::vector<std::complex<double>> V3(N_new * N_new, std::complex<double>(0, 0));
    for (int i = 0; i < N * N; ++i) {
        V3[i] = V[i];
    }

    std::complex<double> Cos = V3[J0 * J0 + J0] * 2.046653415892977; // 2*sqrt(pi/3);
    std::complex<double> Sina_h = V3[J0 * J0 + J0 + 1] * 2.8944050182330705; // 2*(sqrt(2*pi/3));
    std::complex<double> Sina_l = V3[J0 * J0 + J0 - 1] * 2.8944050182330705; // 2*(sqrt(2*pi/3));

    for (int I__ = 1; I__ <= n; ++I__) {
        int J3 = (N - 1 + I__);
        int J1 = J3 - 1;
        int J2 = J3 - 2;
        int m1 = 0;
        int m2 = 0;
        int m3 = 0;
        int I1 = (J1 * J1 + m1 + J1);
        int I2 = (J2 * J2 + m2 + J2);
        int I3 = (J3 * J3 + m3 + J3);
        V3[I3] = ((2 * J1 + 1) / (double)(J1 + 1)) * std::sqrt((2 * J1 + 3) / (double)(2 * J1 + 1)) * Cos * V3[I1] - ((J1) / (double)(J1 + 1)) * std::sqrt((2 * J1 + 3) / (double)(2 * J1 - 1)) * V3[I2];

        for (int m = 0; m <= J1; ++m) {
            m3 = m + 1;
            m1 = m + 1;
            m2 = m + 1;

            I1 = (J1 * J1 + m1 + J1);
            I2 = (J2 * J2 + m2 + J2);
            I3 = (J3 * J3 + m3 + J3);
            int Im = (J1 * J1 + m + J1);

            if (m2 <= J2) {
                V3[I3] =(std::complex<double>)2 * Cos * std::sqrt(((2 * J1 + 3) / (double)(2 * J1 + 1)) * ((J1 - m) / (double)(J1 + m + 2))) * V3[I1] - std::sqrt(((2 * J1 + 3) / (double)(2 * J1 - 1)) * (((J1 - m) * (J1 - m - 1)) / (double)((J1 + m + 2) * (J1 + m + 1)))) * V3[I2];
            } else if (m1 <= J1) {
                V3[I3] =(std::complex<double>)2 * Cos * std::sqrt(((2 * J1 + 3) / (double)(2 * J1 + 1)) * ((J1 - m) / (double)(J1 + m + 2))) * V3[I1];
            }
            V3[I3] += std::sqrt(((2 * J1 + 3) /(double)(2 * J1 + 1)) * (1 /(double)((J1 + m + 2) * (J1 + m + 1)))) * (2 * m + 1) * V3[Im] * Sina_h;

            I1 = (J1 * J1 - m1 + J1);
            I2 = (J2 * J2 - m2 + J2);
            I3 = (J3 * J3 - m3 + J3);
            Im = (J1 * J1 - m + J1);

            if (m2 <= J2) {
                V3[I3] = (std::complex<double>)2 * Cos * std::sqrt(((2 * J1 + 3) / (double)(2 * J1 + 1)) * ((J1 - m) / (double)(J1 + m + 2))) * V3[I1] - std::sqrt(((2 * J1 + 3) /(double)(2 * J1 - 1)) * (((J1 - m) * (J1 - m - 1)) /(double)((J1 + m + 2) * (J1 + m + 1)))) * V3[I2];
            } else if (m1 <= J1) {
                V3[I3] = (std::complex<double>)2 * Cos * std::sqrt(((2 * J1 + 3) / (double)(2 * J1 + 1)) * ((J1 - m) / (double)(J1 + m + 2))) * V3[I1];
            }
            V3[I3] += std::sqrt(((2 * J1 + 3) /(double)(2 * J1 + 1)) * (1 /(double)((J1 + m + 2) * (J1 + m + 1)))) * (2 * m + 1) * V3[Im] * Sina_l;
        }
    }

    return V3;
}


// Length of V >=9
// Extract reference vector from the Quadrupole moment

  std::vector< std::vector<double> > equivalentFeatures::ReferencesExtract(const  std::vector<double>& VR) {
    if (VR.size() < 9) {
        throw std::invalid_argument("Length of V must be >= 9");
    }


    double a12 = VR[5-1]/1.0925484305920792;
    double a23 = VR[6-1]/1.0925484305920792;
    double a13 = VR[8-1]/1.0925484305920792;
    double r3 = 4.1887902047863905 * (VR[2-1] * VR[2-1] + VR[3-1] * VR[3-1] + VR[4-1] * VR[4-1]);
    double r1 = VR[7-1];
    double r2 = VR[9-1];
    double a33 = (-0.31539156525252005 * r3 - r1) / (-0.31539156525252005 - 0.63078313050504);
    double a11 = (0.5462742152960396 * (r1 - 0.63078313050504 * a33) - 0.31539156525252005 * r2) / (2 * 0.5462742152960396 * (-0.31539156525252005));
    double a22 = (0.5462742152960396 * (r1 - 0.63078313050504 * a33) + 0.31539156525252005 * r2) / (2 * 0.5462742152960396 * (-0.31539156525252005));
    double alpha = a11 + a22 + a33;
    double beta = a12 * a12 + a13 * a13 + a23 * a23 - a11 * a22 - a22 * a33 - a33 * a11;
    double gamma = a11 * a22 * a33 + 2 * a12 * a23 * a13 - a11 * a23 * a23 - a12 * a12 * a33 - a13 * a13 * a22;
    double p = std::abs(-((3 * beta + alpha * alpha) / 3));
    double q = -(gamma + 2 * alpha * alpha * alpha / 27 + alpha * beta / 3);
    double x = -q / (2 * std::sqrt(std::pow(p / 3, 3)));
    double a = -3.0 / 4.0;
    double b = -x / 4.0;
    std::complex<double> im(0.0, 1.0);
    double cos_div_3 = 2.0 * std::real(std::pow(-(b / 2.0) + std::sqrt(std::complex<double>(-(a / 3.0) * (a / 3.0) * (a / 3.0) - (b / 2.0) * (b / 2.0))) * im, 1.0 / 3.0));
    double lambda1 = alpha / 3.0 + 2.0 * std::sqrt(p / 3.0) * cos_div_3;
    a11 = a11 - lambda1;
    a22 = a22 - lambda1;
    a33 = a33 - lambda1;

    std::vector<double> T1 = {a22 * a33 - a23 * a23, -(a12 * a33 - a23 * a13), a12 * a23 - a22 * a13};
    double norm_T1 = std::sqrt(T1[0] * T1[0] + T1[1] * T1[1] + T1[2] * T1[2]);
    for (auto& val : T1) {
        val /= norm_T1;
    }

    std::vector<double> T2 = {VR[4-1], VR[2-1], VR[3-1]};
    double norm_T2 = std::sqrt(T2[0] * T2[0] + T2[1] * T2[1] + T2[2] * T2[2]);
    for (auto& val : T2) {
        val /= norm_T2;
    }

    std::vector< std::vector<double> > Array_return;
    Array_return.push_back(T1);
    Array_return.push_back(T2);
    return Array_return;
}





std::vector<std::complex<double>> equivalentFeatures::W3jProduct(const std::vector<std::complex<double>>& V1,
                                              const std::vector<std::complex<double>>& V2,
                                              const std::vector<std::complex<double>>& V3,
                                              int n1, int n2) {
    int RInvariantVSize = 0;
    int32_t d2 = static_cast<int32_t>(std::floor(std::pow(V2.size(), 0.5)));
    int32_t d3 = static_cast<int32_t>(std::floor(std::pow(V3.size(), 0.5)));



    for (int I = 1; I <= d2; ++I) {
            RInvariantVSize += std::max(std::min(d3 - 1, n1 + I - 2) - std::abs(n1 - I) + 1, 0)
                          - std::max(0, std::min(I + n1 - n2 - 1, std::min(d3 - 1, n1 + I - 2) + 1)
                          - std::max(I + n2 - n1 + 1, std::abs(n1 - I) + 1) + 1);
    }



    std::vector<std::complex<double>> RInvariantV(RInvariantVSize);
    int Index_ = 0;

    int J1    =  n1-1;

    for (int J2 = 0; J2 < d2; ++J2) {
        for (int J3 = std::abs(J2 - J1); J3 <= std::min((J1 + J2),d2-1); ++J3) {

        if( J1 >= std::abs(J2-J3) + n2)
        {
            J3 = std::max(J3,J2 + n1 - n2-2);
            continue;
        }
         auto TempI = W3JTableI[J1][J2][J3];
         auto TempC = W3JTableC[J1][J2][J3];
         std::complex<double> Value = 0;

         for (int I_1 = 0; I_1 < TempI.size(); ++I_1) {
            int m1 = TempI[I_1].m1;
            int m2 = TempI[I_1].m2;
            int m3 = -(m1 + m2);
            int I1 = (J1 * J1 + m1 + J1);
            int I2 = (J2 * J2 + m2 + J2);
            int I3 = (J3 * J3 + m3 + J3);
            Value += V1[I1] * V2[I2] * V3[I3] * TempC[I_1];
         }
         RInvariantV[Index_] = Value;
         Index_++;
        }
    }

    return RInvariantV;
}


std::vector<double> equivalentFeatures::W3jProductCToR(const std::vector<std::complex<double>>& InvariantV, int n1, int d2, int d3,int n2) {

    int RInvariantVSize = 0;

    for (int I = 1; I <= d2; ++I) {
            RInvariantVSize += std::max(std::min(d3 - 1, n1 + I - 2) - std::abs(n1 - I) + 1, 0)
                          - std::max(0, std::min(I + n1 - n2 - 1, std::min(d3 - 1, n1 + I - 2) + 1)
                          - std::max(I + n2 - n1 + 1, std::abs(n1 - I) + 1) + 1);
    }

    std::vector<double> RInvariantV(RInvariantVSize);
    int Index =  0;
    int J1    =  n1-1;


    for (int J2 = 0; J2 < d2; ++J2) {
      for (int J3 = std::abs(J2 - J1); J3 <= std::min((J1 + J2),d2-1); ++J3) {

        if( J1 >= std::abs(J2-J3) + n2)
        {
            J3 = std::max(J3,J2 + n1 - n2-2);
            continue;
        }

        RInvariantV[Index] = std::real(std::pow(std::complex<double>(0, 1), J1 - J2 - J3 ) * InvariantV[Index]);
        Index++;
        }
    }

   return RInvariantV;
}


std::vector<std::complex<double>> equivalentFeatures::W3jProductRToC(const std::vector<double>& InvariantV, int n1,int d2, int d3, int n2) {
    int CInvariantVSize = 0;


    for (int I = 1; I <= d2; ++I) {
            CInvariantVSize += std::max(std::min(d3 - 1, n1 + I - 2) - std::abs(n1 - I) + 1, 0)
                          - std::max(0, std::min(I + n1 - n2 - 1, std::min(d3 - 1, n1 + I - 2) + 1)
                          - std::max(I + n2 - n1 + 1, std::abs(n1 - I) + 1) + 1);
    }


    std::vector<std::complex<double>> CInvariantV(CInvariantVSize);
    int Index =  0;
    int J1    =  n1-1;

    for (int J2 = 0; J2 < d2; ++J2) {
      for (int J3 = std::abs(J2 - J1); J3 <= std::min((J1 + J2),d2-1); ++J3) {

         if( J1 >= std::abs(J2-J3) + n2)
         {
               J3 = std::max(J3,J2 + n1 - n2-2);
               continue;
         }

         CInvariantV[Index] = std::pow(std::complex<double>(0, 1), -J1 + J2 + J3 ) * InvariantV[Index];
         Index++;

        }
    }

    return CInvariantV;
}





std::vector<double> equivalentFeatures::SelfProduct(const std::vector<std::complex<double>>& V, int n, int n2) {
    int size__ = 0;
    int d2 = static_cast<int>(std::floor(std::pow(V.size(), 0.5)));
    int d1 = n;

    for (int J1 = 1; J1 <= d1; ++J1) {
        for (int J2 = std::max(J1-d2+1,1);J2 <= std::min(d2-1+J1,d2); ++J2) {
            int Begin = std::max({(J1 - J2 + 1), J2, (J2 - n2 + J1)});
            int End = std::min(J1 + J2 - 1, d2);


            for (int I = Begin; I <= End; I += 1) {

                if (((I == J2)&&((J1%2)==0)))
                    continue;
                size__++;
            }
        }
    }

    std::vector<double> RInvariantV(size__);
    int Index_ = 0;

    for (int J1 = 1; J1 <= d1; ++J1) {
        for (int J2 = std::max(J1-d2+1,1); J2 <= std::min(d2-1+J1,d2); ++J2) {
            int Begin = std::max({J1 - J2 + 1, J2, J2 - n2 + J1});
            int End = std::min(J1 + J2 - 1, d2);



            for (int I = Begin; I <= End; I += 1) {
                int j1 = J1 - 1;
                int j2 = J2 - 1;
                int j3 = I - 1;
                int J = j1 - (I - J2) + 1;

                if (((I == J2)&&((J%2)==0)))
                    continue;

                auto TempI = CGTableI[J2-1][I-1][J-1];
                auto TempC = CGTableC[J2-1][I-1][J-1];

                for (size_t I_1 = 0; I_1 < TempI.size(); ++I_1) {
                    std::complex<double> Value(0, 0);
                    for (size_t I_2 = 0; I_2 < TempI[I_1].size(); ++I_2) {
                        int m1 = TempI[I_1][I_2].m1;
                        int m2 = TempI[I_1][I_2].m2;
                        int I1 = (j2 * j2 + m1 + j2);
                        int I2 = (j3 * j3 + m2 + j3);
                        Value += V[I1] * V[I2] * TempC[I_1][I_2];
                    }
                    RInvariantV[Index_] += std::real(Value * std::conj(Value));
                }
                Index_++;
            }
        }
    }

    return RInvariantV;
}






std::vector<std::complex<double>> equivalentFeatures::W3jProductCompact(const std::vector<std::complex<double>>& V1,
                                                    const std::vector<std::complex<double>>& V2,
                                                    const std::vector<std::complex<double>>& V3,
                                                    int n) {
    int RInvariantVSize = n * n;
    std::vector<std::complex<double>> RInvariantV(RInvariantVSize, std::complex<double>(0, 0));
    int Index_ = 0;


    for(int J1 = 0;J1<n ;J1++)
    {
       for (int J2 = 0; J2 <= J1; ++J2) {
        for (int J3 = std::abs(J2 - J1); J3 <= std::min({(J1 + J2),J1,J1-J2+1}); ++J3) {

         auto TempI = W3JTableI[J1][J2][J3];
         auto TempC = W3JTableC[J1][J2][J3];

        std::complex<double> Value(0, 0);

        for (int I_1 = 0; I_1 < TempI.size(); I_1++) {
            int m1 = TempI[I_1].m1;
            int m2 = TempI[I_1].m2;
            int m3 = -(m1 + m2);
            int I1 = (J1 * J1 + m1 + J1);
            int I2 = (J2 * J2 + m2 + J2);
            int I3 = (J3 * J3 + m3 + J3);
            Value += V1[I1] * V2[I2] * V3[I3] * TempC[I_1];
        }

        RInvariantV[Index_] = Value;
        Index_++;

        }}
   }

    return RInvariantV;
}




std::vector<double> equivalentFeatures::W3jProductCompactCToR(const std::vector<std::complex<double>>& InvariantV, int n) {
    int RInvariantVSize = n * n;
    std::vector<double> RInvariantV(RInvariantVSize,0.0);
    int Index_ = 0;


    for(int J1 = 0;J1<n ;J1++)
    {
       for (int J2 = 0; J2 <= J1; ++J2) {
        for (int J3 = std::abs(J2 - J1); J3 <= std::min({(J1 + J2),J1,J1-J2+1}); ++J3) {

        std::complex<double> Value(0, 0);

        RInvariantV[Index_] = std::real(std::pow(std::complex<double>(0.0, 1.0), J1 - J2 -J3) * InvariantV[Index_]);
        Index_++;

        }}
   }

    return RInvariantV;
}



std::vector<std::complex<double>> equivalentFeatures::W3jProductCompactRToC(const std::vector<double>& InvariantV, int n) {
    int CInvariantVSize = n * n;
    std::vector<std::complex<double>> CInvariantV(CInvariantVSize, std::complex<double>(0.0, 0.0));
    int Index_ = 0;

    for(int J1 = 0;J1<n ;J1++)
    {
       for (int J2 = 0; J2 <= J1; ++J2) {
        for (int J3 = std::abs(J2 - J1); J3 <= std::min({(J1 + J2),J1,J1-J2+1}); ++J3) {

        CInvariantV[Index_] = std::pow(std::complex<double>(0.0, 1.0), -J1 + J2 +J3) * InvariantV[Index_];
        Index_++;
    }}}
    return CInvariantV;
}









Eigen::MatrixXcd equivalentFeatures::DecodeMatrixCompact(const std::vector<std::complex<double>>& V2, const std::vector<std::complex<double>>& V3, int n) {
    int VSize = n * n;
    MatrixXcd LTm = MatrixXcd::Zero(VSize, VSize);
    int Index_ = 0;



    for(int J1 = 0;J1<n ;J1++)
    {
      for (int J2 = 0; J2 <= J1; ++J2) {
        for (int J3 = std::abs(J2 - J1); J3 <= std::min({(J1 + J2),J1,J1-J2+1}); ++J3) {

         for (int I_1 = 0; I_1 < W3JTableI[J1][J2][J3].size(); I_1++) {
            int m1 = W3JTableI[J1][J2][J3][I_1].m1;
            int m2 = W3JTableI[J1][J2][J3][I_1].m2;
            int m3 = -(m1 + m2);
            int I1 = (J1 * J1 + m1 + J1);
            int I2 = (J2 * J2 + m2 + J2);
            int I3 = (J3 * J3 + m3 + J3);
            LTm(Index_, I1) += V2[I2] * V3[I3] * W3JTableC[J1][J2][J3][I_1];
        }
        Index_++;
        }}
   }

   Eigen::MatrixXcd G_inv = LTm.completeOrthogonalDecomposition().pseudoInverse();
   return G_inv;
}




std::vector<double> equivalentFeatures::SelfProductPairwise(const std::vector<std::complex<double>>& V, int n, int n2) {
    int32_t d2 = static_cast<int32_t>(std::floor(std::pow(V.size(), 0.5)));
    int d1 = std::min(n,d2);
    int Max_size = 0;
    int size_sum = 0;

    for (int J1 = 1; J1 <= d1; ++J1) {
        int size__ = 0;
        for (int J2 = std::max(J1-d2+1,1); J2 <= std::min(d2-1+J1,d2); ++J2) {
            int Begin = std::max({J1 - J2 + 1, J2, J2-n2+J1});
            int End = std::min(J1 + J2 - 1, d2);

            for (int I = Begin; I <= End; I += 1) {

                if (((I == J2)&&((J1%2)==0)))
                    continue;
                size__++;
            }
        }
        if (size__ > Max_size) {
            Max_size = size__;
        }

        if(J1>1)
          size_sum = size__ * (size__ + 1) / 2 + size_sum;
        else
          size_sum = size__  + size_sum;
    }


    std::vector<std::vector<std::complex<double>>> LTm(Max_size, std::vector<std::complex<double>>(d2 * d2 - (d2 - 1) * (d2 - 1), 0.0));
    std::vector<std::complex<double>> LTmConj(d2 * d2 - (d2 - 1) * (d2 - 1), std::complex<double>(0.0,0.0));
    std::vector<int> DiffJ(d2 * d2 - (d2 - 1) * (d2 - 1), 0.0);
    std::vector<double> Result_(size_sum, 0);
    size_sum = 0;

    for (int J1 = 1; J1 <= d1; ++J1) {
        int Index = 0;
        for (int J2 = std::max(J1-d2+1,1); J2 <= std::min(d2-1+J1,d2); ++J2) {
            int Begin = std::max({J1 - J2 + 1, J2, J2 - n2 + J1});
            int End = std::min(J1 + J2 - 1, d2);


            for (int I = Begin; I <= End; I += 1) {

                int j1 = J1 - 1;
                int j2 = J2 - 1;
                int j3 = I - 1;
                int J  = j1 - (I - J2) + 1;



                if (((I == J2)&&((J%2)==0)))
                    continue;


                auto TempI = CGTableI[J2-1][I-1][J-1];
                auto TempC = CGTableC[J2-1][I-1][J-1];

                std::fill(LTm[Index].begin(), LTm[Index].end(), std::complex<double>(0.0, 0.0));
                for (size_t I_1 = 0; I_1 < TempI.size(); ++I_1) {
                    for (size_t I_2 = 0; I_2 < TempI[I_1].size(); ++I_2) {
                        auto m1 = TempI[I_1][I_2].m1;
                        auto m2 = TempI[I_1][I_2].m2;
                        auto m3 = (m1 + m2);
                        int I1 = (j2 * j2 + m1 + j2 );
                        int I2 = (j3 * j3 + m2 + j3 );
                        int I3 = (m3 + j1);
                        LTm[Index][I3] += V[I1] * V[I2] * TempC[I_1][I_2];
                    }
                }

                DiffJ[Index] = J;
                Index++;
            }
        }
        if (Index > 1) {

            int Len = J1 * J1 - (J1 - 1) * (J1 - 1);

            for (int I_3 = 1; I_3 <= Index; ++I_3) {

                  for (int i =0;i<Len;i++){
                       LTmConj[i] = std::conj(LTm[I_3-1][i]);
                    }

                  if (Len > 1){

                     for (int I_2 = I_3; I_2 <= Index; ++I_2) {
                      if ( abs(DiffJ[I_3-1]-DiffJ[I_2-1]) != 1 )
                      {
                          std::complex<double> Self_Product = std::inner_product(LTmConj.begin(), LTmConj.begin() + Len, LTm[I_2-1].begin(), std::complex<double>(0.0, 0.0));
                          Result_[size_sum] = std::real(Self_Product);
                          size_sum++;
                      }
                    }
                  }else
                   {
                          std::complex<double> Self_Product = std::inner_product(LTmConj.begin(), LTmConj.begin() + Len, LTm[I_3-1].begin(), std::complex<double>(0.0, 0.0));
                          Result_[size_sum] = std::real(Self_Product);
                          size_sum++;
                   }


            }
        }
    }
    return std::vector<double>(Result_.begin(), Result_.begin() + size_sum);
}





std::vector<double> equivalentFeatures::ProductEncode(const std::vector<std::complex<double>>& V1, const std::vector<std::complex<double>>& V2, int n) {
    int d1 = static_cast<int>(std::floor(std::sqrt(V1.size())));
    int d2 = static_cast<int>(std::floor(std::sqrt(V2.size())));
    int SIZE = std::min({d1 + d2 - 1, N});

    std::vector<std::complex<double>> LTm(SIZE * SIZE, 0.0);
    std::vector<std::complex<double>> LTmConj(SIZE * SIZE, 0.0);

    int LENGTH = 0;
    for (int J3 = std::abs(d1 - d2) + 1; J3 <= std::min({d1 + d2 - 1, N}); ++J3) {
        for (int J1 = std::max({J3 - d2 + 1, 1}); J1 <= std::min({J3 + d2 - 1, d1}); ++J1) {
            for (int J2 = std::max({J3 - J1 + 1, J1 - J3 + 1, 1}); J2 <= std::min({J3 + J1 - 1,d2}); ++J2) {
                   if( J3 > std::abs(J1-J2) + n)
                    {
                         J2 = std::max(J2,J1 + N - n-1);
                         continue;
                    }

                LENGTH++;
            }
        }
    }

    std::vector<double> pseudoInput(LENGTH);
    int I = 0;
    for (int J3 = std::abs(d1 - d2) + 1; J3 <= std::min({d1 + d2 - 1, N}); ++J3) {
        int Index = 1;
        int BASE = (J3 - 1) * (J3 - 1);

        for (int J1 = std::max({J3 - d2 + 1, 1}); J1 <= std::min({J3 + d2 - 1, d1}); ++J1) {
            for (int J2 = std::max({J3 - J1 + 1, J1 - J3 + 1, 1}); J2 <= std::min({J3 + J1 - 1,d2}); ++J2) {

                if( J3 > std::abs(J1-J2) + n)
                {
                    J2 = std::max(J2,J1 + N - n-1);
                    continue;
                }
                const auto& TempI = W3JTableI[J1-1][J2-1][J3-1];
                const auto& TempC = W3JTableC[J1-1][J2-1][J3-1];
                for (size_t I_1 = 0; I_1 < TempI.size(); ++I_1) {
                    int m1 = TempI[I_1].m1;
                    int m2 = TempI[I_1].m2;
                    int m3 = -(m1 + m2);
                    int I1 = (J1 - 1) * (J1 - 1) + m1 + J1;
                    int I2 = (J2 - 1) * (J2 - 1) + m2 + J2;
                    int I3 = (J3 - 1) * (J3 - 1) + m3 + J3;
                    LTm[I3 - BASE - 1] += V1[I1-1] * V2[I2-1] * TempC[I_1] * std::sqrt(2 * J3 - 1) * std::pow(-1, J1 - J2 - m3);
                }
                int Len = J3*J3 - BASE;
                for (int i =0;i<Len;i++)
                {
                    LTmConj[i] = std::conj(LTm[i]);
                }
                pseudoInput[I] = std::real(std::inner_product(LTmConj.begin(),LTmConj.begin()+Len,LTm.begin(),std::complex<double>(0.0, 0.0)));
                I++;
                std::fill(LTm.begin(), LTm.end(), 0.0);
            }
        }
    }

    return pseudoInput;
}






std::vector<double> equivalentFeatures::ProductEncodePairwise(const std::vector<std::complex<double>>& V1, const std::vector<std::complex<double>>& V2, int n) {
    int d1 = static_cast<int>(std::floor(std::sqrt(V1.size())));
    int d2 = static_cast<int>(std::floor(std::sqrt(V2.size())));
    int SIZE = std::min({d1 + d2 - 1, N});
    int size__ = 0;
    int Max_size = 0;
    int size_sum = 0;
    int LENGTH = 0;

    for (int J3 = std::abs(d1 - d2) + 1; J3 <= std::min({d1 + d2 - 1, N}); J3++) {
        size__ = 0;
        for (int J1 = std::max({J3 - d2 + 1,1}); J1 <= std::min({J3 + d2 - 1, d1}); J1++) {
            for (int J2 = std::max({J3 - J1 + 1, J1 - J3 + 1, 1}); J2 <= std::min({J3 + J1 - 1,d2}); J2++) {
                    if( J3 > std::abs(J1-J2) + n)
                    {
                        J2 = std::max(J2,J1 + N - n-1);
                        continue;
                    }
                size__++;
            }
        }
        if (Max_size < size__) {
            Max_size = size__;
        }
        if (J3 > 1) {
            LENGTH = static_cast<int>(size__ * (size__ + 1) / 2) + LENGTH;
        } else {
            LENGTH = size__ + LENGTH;
        }
        if (J3 == 2) {
            LENGTH = LENGTH - 2;
        }
    }

    std::vector<std::vector<std::complex<double>>> LTm(Max_size, std::vector<std::complex<double>>(SIZE * SIZE, 0));
    std::vector<std::complex<double>> LTmConj(SIZE * SIZE, 0.0);
    std::vector<double> DiffJ(Max_size, 0);
    std::vector<double> pseudoInput(LENGTH, 0);

    int I = 1;
    for (int J3 = std::abs(d1 - d2) + 1; J3 <= std::min({d1 + d2 - 1,N}); J3++) {
        int Index = 1;
        int BASE = (J3 - 1) * (J3 - 1);
        for (int J1 = std::max({J3 - d2 + 1, 1}); J1 <= std::min({J3 + d2 - 1, d1}); J1++) {
            for (int J2 = std::max({J3 - J1 + 1, J1 - J3 + 1, 1}); J2 <= std::min({J3 + J1 - 1, d2}); J2++) {
                if( J3 > std::abs(J1-J2) + n)
                {
                    J2 = std::max(J2,J1 + N - n-1);
                    continue;
                }
                auto TempI = W3JTableI[J1-1][J2-1][J3-1];
                auto TempC = W3JTableC[J1-1][J2-1][J3-1];
                for (size_t I_1 = 0; I_1 < TempI.size(); I_1++) {
                    int m1 = TempI[I_1].m1;
                    int m2 = TempI[I_1].m2;
                    int m3 = -(m1 + m2);
                    int I1 = (J1 - 1) * (J1 - 1) + m1 + J1;
                    int I2 = (J2 - 1) * (J2 - 1) + m2 + J2;
                    int I3 = (J3 - 1) * (J3 - 1) + m3 + J3;
                    LTm[Index-1][I3 - BASE-1] += V1[I1-1] * V2[I2-1] * TempC[I_1] * (std::sqrt(2 * J3-1) * std::pow(-1, J1 - J2 - m3)); // wigner3J coefficient to Clebsh-Gordan coefficient
                }
                DiffJ[Index-1] = std::abs(J1 - J2);
                Index++;
                I++;
            }
        }

        if (Index > 1) {
            int Len = J3 * J3 - (J3 - 1) * (J3 - 1);
            for (int I_3 = 1; I_3 < Index; I_3++) {


                for (int i =0;i<Len;i++)
                {
                    LTmConj[i] = std::conj(LTm[I_3-1][i]);
                }

                if (Len > 1) {
                    for (int I_2 = I_3; I_2 < Index; I_2++) {
                        if (!(((I_3 == 2) || (I_3 == 1)) && (I_2 == 3) && (Len == 3))) {
                            auto Self_Product = std::inner_product(LTmConj.begin(),LTmConj.begin()+Len,LTm[I_2-1].begin(),std::complex<double>(0.0, 0.0));
                            if (std::abs(DiffJ[I_3-1] - DiffJ[I_2-1]) != 1) {
                                pseudoInput[size_sum] = std::real(Self_Product);
                            } else {
                                pseudoInput[size_sum] = std::imag(Self_Product);
                            }
                            size_sum++;
                        }
                    }
                } else {
                    auto Self_Product = std::inner_product(LTmConj.begin(),LTmConj.begin()+Len,LTm[I_3-1].begin(),std::complex<double>(0.0, 0.0));
                    pseudoInput[size_sum] = std::real(Self_Product);
                    size_sum++;
                }
            }
        }
        std::fill(LTm.begin(), LTm.end(), std::vector<std::complex<double>>(SIZE * SIZE, 0));
    }
    return pseudoInput;
}






std::vector<std::vector<std::complex<double>>> equivalentFeatures::SelfProductMatrix(const std::vector<std::complex<double>>& V, int n, int n2) {
    int d2 = static_cast<int>(floor(pow(V.size(), 0.5)));
    int d1 = n;
    int size_sum = 0;
    int size_add = 0;
    int odd_sum  = 0;

    int J1 = d1;
    int size__ = 0;
    size_add = 0;

    for (int J2 = max({J1 - d2 + 1, 1}); J2 <= min({d2 - 1 + J1, d2}); J2++) {
        int Begin = max({(J1 - J2 + 1), J2, (J2 - n2 + J1)});
        int End = min(J1 + J2 - 1, d2);
        for (int I = Begin; I <= End; I++) {
            if ((I == J2) && ((J1 % 2) == 0)) {
                continue;
            }
            size__++;
        }
    }

    vector<vector<complex<double>>> LTm(size__, vector<complex<double>>(d1 * d1 - (d1 - 1) * (d1 - 1), complex<double>(0.0, 0.0)));

    J1 = d1;
    int Index = 0;
    for (int J2 = max({J1 - d2 + 1, 1}); J2 <= min({d2 - 1 + J1, d2}); J2++) {
        int Begin = max({(J1 - J2 + 1), J2, (J2 - n2 + J1)});
        int End = min(J1 + J2 - 1, d2);

        for (int I = Begin; I <= End; I++) {
            int j1 = J1 - 1;
            int j2 = J2 - 1;
            int j3 = I - 1;
            int J = j1 + 1;

            if (((I == J2) && ((J % 2) == 0))) {
                continue;
            }

            auto TempI = W3JTableI[J2-1][I-1][J-1];
            auto TempC = W3JTableC[J2-1][I-1][J-1];
            fill(LTm[Index].begin(), LTm[Index].end(), complex<double>(0.0, 0.0));
            for (size_t I_1 = 0; I_1 < TempI.size(); I_1++) {
                int m1 = TempI[I_1].m1;
                int m2 = TempI[I_1].m2;
                int m3 = (m1 + m2);
                int I1 = (j2 * j2 + m1 + j2);
                int I2 = (j3 * j3 + m2 + j3);
                int I3 = (m3 + j1);
                LTm[Index][I3] += V[I1] * V[I2] * TempC[I_1];
            }
            Index++;
        }
    }

    return LTm;
}


//std::sqrt(2 * j3 + 1)*std::pow(-1,j1 - j2 - m3)
using namespace std;

vector<vector<complex<double>>> equivalentFeatures::DerivativeSH(const vector<complex<double>>& Y) {
    double SinTheta = sqrt(pow(Y[1].real(), 2) + pow(Y[1].imag(), 2));
    complex<double> ExpVarphi = SinTheta != 0.0 ? Y[1] / SinTheta : complex<double>(0, 0);
    double CotTheta = SinTheta != 0.0 ? Y[2].real() / (SinTheta * sqrt(2.0)) : 0.0;

    vector<complex<double>> DerivationTheta(Y.size(), complex<double>(0, 0));
    vector<complex<double>> DerivationVarphi(Y.size(), complex<double>(0, 0));

    int N = static_cast<int>(sqrt(Y.size()));
    int Base_ = 0;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < i*2; ++j) {
            int m = j - i;
            DerivationTheta[Base_ + j] = m * CotTheta * Y[Base_ + j] + sqrt((i - m ) * (i+1+m)) * ExpVarphi * Y[Base_ +j+1];
            DerivationVarphi[Base_ + j] = complex<double>(0, m) * Y[Base_ + j];
        }
        int j = i * 2;
        int m = j - i;
        DerivationTheta[Base_ + j] = m * CotTheta * Y[Base_ + j];
        DerivationVarphi[Base_ + j] = complex<double>(0, m) * Y[Base_ + j];
        Base_ = (i+1)*(i+1);
    }

    return {DerivationTheta, DerivationVarphi};
}





vector<vector<cdouble>> equivalentFeatures::DerivativeSH_XYZ(const vector<cdouble>& Y, const vector<cdouble>& DR, double ReciprocalRadii, double ReciprocalF) {

    double SinTheta = sqrt(pow(Y[1].real(), 2) + pow(Y[1].imag(), 2));
    cdouble ExpVarphi = SinTheta != 0.0 ? Y[1] / SinTheta : cdouble(0, 0);
    double CotTheta = SinTheta != 0.0 ? (Y[2].real()) / (SinTheta * sqrt(2.0)) : 0.0;

    double SinVarphi = SinTheta != 0.0 ? -Y[1].imag() / SinTheta : 0.0;
    double CosVarphi = SinTheta != 0.0 ? Y[1].real() / SinTheta : 0.0;
    double CosTheta = ReciprocalF * (Y[2].real() / (0.5 * sqrt(3.0 / M_PI)));



    size_t lenY = Y.size();
    vector<cdouble> DerivationX(lenY, cdouble(0, 0));
    vector<cdouble> DerivationY(lenY, cdouble(0, 0));
    vector<cdouble> DerivationZ(lenY, cdouble(0, 0));
    vector<cdouble> DRc = DR;

    SinTheta = ReciprocalF * (SinTheta / (0.5 * sqrt(3.0 / (M_PI * 2.0))));
    double SinVarphiSlashSinTheta = SinTheta != 0.0 ? SinVarphi / SinTheta : 0.0;
    double CosVarphiSlashSinTheta = SinTheta != 0.0 ? CosVarphi / SinTheta : 0.0;

    int N = static_cast<int>(sqrt(Y.size()));
    int Base_ = 0;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < i*2; ++j) {

            int m = j - i;
            cdouble DerivationTheta  = m * CotTheta * Y[Base_ + j] + sqrt((i - m ) * (i+1+m)) * ExpVarphi * Y[Base_ +j+1];
            cdouble DerivationVarphi = complex<double>(0, m) * Y[Base_ + j];

            DerivationX[Base_ + j] = DRc[Base_ + j] * Y[Base_ + j] * SinTheta * CosVarphi + ReciprocalRadii * (CosTheta * CosVarphi * DerivationTheta - SinVarphiSlashSinTheta * DerivationVarphi);
            DerivationY[Base_ + j] = DRc[Base_ + j] * Y[Base_ + j] * SinTheta * SinVarphi + ReciprocalRadii * (CosTheta * SinVarphi * DerivationTheta + CosVarphiSlashSinTheta * DerivationVarphi);
            DerivationZ[Base_ + j] = DRc[Base_ + j] * Y[Base_ + j] * CosTheta + ReciprocalRadii * (-SinTheta * DerivationTheta);
        }
        int j = i * 2;
        int m = j - i;
        cdouble DerivationTheta     = m * CotTheta * Y[Base_ + j];
        cdouble DerivationVarphi    = complex<double>(0, m) * Y[Base_ + j];
        DerivationX[Base_ + j] = DRc[Base_ + j] * Y[Base_ + j] * SinTheta * CosVarphi + ReciprocalRadii * (CosTheta * CosVarphi * DerivationTheta - SinVarphiSlashSinTheta * DerivationVarphi);
        DerivationY[Base_ + j] = DRc[Base_ + j] * Y[Base_ + j] * SinTheta * SinVarphi + ReciprocalRadii * (CosTheta * SinVarphi * DerivationTheta + CosVarphiSlashSinTheta * DerivationVarphi);
        DerivationZ[Base_ + j] = DRc[Base_ + j] * Y[Base_ + j] * CosTheta + ReciprocalRadii * (-SinTheta * DerivationTheta);

        Base_ = (i+1)*(i+1);
    }

    // Return concatenated matrix: each vector is a column
    // So result is vector of vectors: rows = 3, cols = lenY
    vector<vector<cdouble>> result(3, vector<cdouble>(lenY));
    for (size_t i = 0; i < lenY; ++i) {
        result[0][i] = DerivationX[i];
        result[1][i] = DerivationY[i];
        result[2][i] = DerivationZ[i];
    }

    return result;
}

















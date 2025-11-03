#include "include/rotation.h"
#include <fstream>
#include <iostream>

int main()
{
    equivalentFeatures Test;
    Test.Initial();
    std::vector<double> V2 = {0.2820947917738782,0.12324306613285871,-0.41829957442835725,0.220383818192421,0.12430010752649816,-0.2359278576176394,0.3780897816045623,-0.42188728105548623,0.07638145146836267};
    std::vector<double> V1 = {0.2820947917738782,-0.29294906671981563,0.34398535815919695,-0.18597669831976563,0.24933311274303446,-0.4611703555484876,0.1535728097664055,-0.2927708254824036,-0.11723009694482389};
    std::vector<double> V3 = {0.2820947917738782,0.40547817509140544,0.2458934345224644,-0.1177127139315595,-0.21843396179773633,0.45629291253978804,-0.07575460360437422,-0.13246453294480323,-0.3445070791340872};
    std::vector<std::complex<double>> V_2 = Test.RStoCS_Encode(V2,3);
    std::vector<std::complex<double>> V_1 = Test.RStoCS_Encode(V1,3);
    std::vector<std::complex<double>> V_3 = Test.RStoCS_Encode(V3,3);
    std::vector<std::complex<double>> V__ = Test.W3jProductCompact(V_2,V_1,V_3,3);
    std::vector<double> V___ = Test.W3jProductCompactCToR(V__,3);
    std::vector<std::complex<double>> V____ = Test.W3jProductCompactRToC(V___,3);
                                          V__   = Test.W3jProduct(V_2,V_1,V_3,3,3);
                                          V___  = Test.W3jProductCToR(V__,3,3,3,3);
                                          V____  = Test.W3jProductRToC(V___,3,3,3,3);
    Eigen::MatrixXcd M = Test.DecodeMatrixCompact(V_2,V_3,3);
    for (size_t i = 0; i < V1.size(); ++i) {
        V3[i] = V1[i] + V2[i] + V3[i];
    }
    std::vector<double> V_sum = {0.8462843753216346,0.2357721745044485,0.1715792182533041,-0.08330559405890413,0.15519925847179628,-0.24080530062633895,0.4559079877665936,-0.847122639482693,-0.3853557246105484};
//    std::vector<double> TT = Test.ProductEncode(Test.RStoCS_Encode(V1,3),Test.RStoCS_Encode(V2,3),5);
//    std::vector<double> TT = Test.SelfProductPairwise(Test.RStoCS_Encode(V1,3),6,2);
      std::vector<double> TT = Test.ProductEncode(Test.RStoCS_Encode(V1,3),Test.RStoCS_Encode(V2,3),3);
      std::vector<std::vector<std::complex<double>>> TTT =  Test.SelfProductMatrix(V_2,3,2);
      std::vector<std::complex<double>> S2 = {
	          {0.28209479177387814, 0.0},
		      {0.10713552128272251, -0.048368181166014684},
		          {0.5748596526258442, 0.0},
			      {-0.10713552128272251, -0.048368181166014684},
			          {0.029572706782343073, -0.033538029951111326},
				      {0.28185427219504855, -0.12724797841762917},
				          {0.8366481389933829, 0.0},
					      {-0.28185427219504855, -0.12724797841762917},
					          {0.029572706782343073, 0.033538029951111326}
      };

      std::vector< std::complex<double> > D  = {0.0, 1.0, 1.0, 1.0, 2.449489742783178, 2.449489742783178, 2.449489742783178, 2.449489742783178, 2.449489742783178};
      double ReciprocalRadii = 0.8164797044666512;
      double ReciprocalF     = 0.8164965809277261;
      std::vector<std::vector<std::complex<double>>> Result  = Test.DerivativeSH(S2);
      std::vector<std::vector<std::complex<double>>> Result2 = Test.DerivativeSH_XYZ(S2,D,ReciprocalRadii,ReciprocalF);

//    TT = Test.SelfProduct(V_3,2,2);
//    std::cout << Test.RStoCS_Encode(V3,3) << std::endl;
    printf("Test");
    //std::vector<std::complex<double>> V3(4);
    //for (size_t i = 0; i < 4; ++i) {
    //    V3[i] = V_[i];
    //};
    //std::vector< std::vector<double> > V__ = Test.ReferencesExtract(V_);
    /*std::vector<double> V1 =  Test.CStoRS_Encode(V_,3);
    std::vector<double> V4 = Test.IncreaseDegree(V3,1);
    double norm_V3 = 0.0;
    for (double val : V3) {
        norm_V3 += val * val;
    };
    norm_V3 = std::sqrt(norm_V3);*/
    return 0;
}

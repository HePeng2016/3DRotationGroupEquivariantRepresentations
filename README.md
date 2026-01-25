# Diagram


![Image](https://github.com/user-attachments/assets/1766cb40-c146-4f60-94d7-701053fce36a)


# Introduction

This software package is used to represent the spherical harmonic tensors in a rotation equivariance or invariance way. The key functions are implemented in two different versions. One version uses Julia,the other version uses C++, they yield the same result. 
# Installation (julia)
     include("./equivalentFeatures.jl")
     using  .equivalentFeatures 
# Initial
     equivalentFeatures.setN(N_new);
*N_new* indicates the maximum orbital quantum number *l* of the input spherical harmonic representation *Y(lm)*, 

Because in julia the index begins from 1, so we define *N_new* as *l+1*. In the following, the definition of degree n is the quantum number *l* plus one. If this command is missing, then the default maximum orbital quantum number *l* is 2, so *N_new* is 3. 
     
     equivalentFeatures.Initial();
   Once this function had been carried out, the tables for the Clebsch Gordan coefficients and the Wigner 3j coefficients were generated. And the data structures for storing the relations between coefficients and variables were generated.
# Usage
    equivalentFeatures.CtoS_Encode(V1);

  This function changes the cartesian coordination terms to the spherical harmonic tensor terms. 
  
  *V1* is the vector that records the multipole moment terms of the multipole expansion in the Cartesian system of coordinates.
  The order of the terms in this vector is determined by the length of the term, if the terms are equally long, they are ordered according to the lexicographical order. 
  
  e.g.
  
   *x<y<z*,
   *xx<x* ,
   *xx<xy<xz<yy<yz<zz*
   
*V1* is formatted as: 

N=1: [monopole]

*[1]*
    
N=2: [monopole,dipole]

*[1,x,y,z]*
    
N=3: [monopole,dipole,quadruple]

*[1,x,y,z,xx,xy,xz,yy,yz,zz]*

The output is a spherical harmonic tensor stored as a complex vector.

e.g.

N=3:

*[Y(l=0,m=0),Y(l=1,m=-1),Y(l=1,m=0),Y(l=1,m=1),Y(l=2,m=-2),Y(l=2,m=-1),Y(l=2,m=0),Y(l=2,m=1),Y(l=2,m=2)]*

*Y(l,m)* is the term for spherical multipole moment in the spherical coordinate system. *l* is the orbital quantum number and *m* is the azimuthal quantum number.

     equivalentFeatures.IncreaseDegree(S2,n);
     
This function increases the degree of the spherical harmonic tensor by n. S2 should be a spherical harmonic tensor [monopole,dipole] derived from a normalized Cartesian coordinate vector.

     equivalentFeatures.SelfProduct(V1)

This function will change a spherical harmonic tensor to a rotation invariance embedding vector.

     equivalentFeatures.SelfProduct(V1,n,n2)

For the Clebsch–Gordan (CG) product of any two shells in the V1 spherical harmonic tensor, only top n2  with smaller degree are selected and the maximum degree is n.  
     
     equivalentFeatures.SelfProductPairwise(V1,n,n2)
     
The definition of n,n2 is the same as SelfProduct function. For the rotation-invariance embedding, both norm and pairwise product are considered. 

e.g. 

      include("./equivalentFeatures.jl")
      using  .equivalentFeatures
      using LinearAlgebra
	  
      equivalentFeatures.Initial();  
      C1 =  [0.0043477849927746155,0.0,0.9999905483381614];
      R1 =  norm(C1,2);
      C1 =  C1/R1; # Cartesian coordinate normalization.
      V1 =  vcat(1, C1);
      
      C2 =  [0.772027518982468,0.33454525822573616,0.5404192632877276];
      R2 =  norm(C2,2);
      C2 =  C2/R2;# Cartesian coordinate normalization.
      V2 =  vcat(1, C2);
      #c2 is the vector obtained by rotating c1.
      
      S1 = equivalentFeatures.CtoS_Encode(V1,2);
      S2 = equivalentFeatures.CtoS_Encode(V2,2);    # Convert the vector to a dipole spherical harmonic tensor(*l=1*), because n indicates *l+1*, the input *n = 2*.  
      S1 = equivalentFeatures.IncreaseDegree(S1,1);  
      S2 = equivalentFeatures.IncreaseDegree(S2,1); # Increase the dipole spherical harmonic tensor into quadrupole spherical harmonic tensor. 
	  S1[2:4] = S1[2:4]*R1; 
	  S1[5:9] = S1[5:9]*R1*R1; 
	  S2[2:4] = S2[2:4]*R2; 
	  S2[5:9] = S2[5:9]*R2*R2; 
      E1 = equivalentFeatures.SelfProduct(S1); 
      E2 = equivalentFeatures.SelfProduct(S2);
      Loss = sum(abs.(E1 - E2))
      E1 = equivalentFeatures.SelfProduct(S1,2,2); 
      E2 = equivalentFeatures.SelfProduct(S2,2,2);
      Loss = sum(abs.(E1 - E2))
      E1 = equivalentFeatures.SelfProductPairwise(S1,2,2);
      E2 = equivalentFeatures.SelfProductPairwise(S2,2,2);
      Loss = sum(abs.(E1 - E2)) 

     
    
V1 and V2 are two identical spherical harmonic tensors with different rotations. 

     equivalentFeatures.W3jProduct(V1,V2,V3) 

This function transforms three spherical harmonic tensors into a vector that is invariant to rotation. These three spherical harmonic tensors could be the tensors that represent the density of two nodes with different centers and the tensor that is derived from the subtraction between the coordinates of two nodes.

      equivalentFeatures.DecodeMatrix(V1,V2)
      
This function will return a matrix that can be used to convert the invariant coding calculated by the W3jProduct function into one of the original spherical harmonic tensors (V3). And the V1, V2 are two other original spherical harmonic tensors.

e.g.
      using LinearAlgebra


      C1 =  [0.0043477849927746155,0.0,0.9999905483381614];
      R1 =  norm(C1,2);
      C1 =  C1/R1; # Cartesian coordinate normalization.
      V1 =  vcat(1, C1);
      
      C2 =  [0.772027518982468,0.33454525822573616,0.5404192632877276];
      R2 =  norm(C2,2);
      C2 =  C2/R2;# Cartesian coordinate normalization.
      V2 =  vcat(1, C2);

      C3 =  [0.654676578,0.7654334543987,0.765432216890];
      R3 =  norm(C3,2);
      C3 =  C3/R3;# Cartesian coordinate normalization.
      V3 =  vcat(1, C3);
      

      
      V1 = equivalentFeatures.CtoS_Encode(V1,2);
      V2 = equivalentFeatures.CtoS_Encode(V2,2);
      V3 = equivalentFeatures.CtoS_Encode(V3,2);
      
      V1 = equivalentFeatures.IncreaseDegree(V1,1);
      V2 = equivalentFeatures.IncreaseDegree(V2,1);
      V3 = equivalentFeatures.IncreaseDegree(V3,1);  
       
       W3 = equivalentFeatures.W3jProduct(V1,V2,V3)
       M  = equivalentFeatures.DecodeMatrix(V2,V3) 
       norm(M*W3-V1,2)

 V1, V2, V3 are three original spherical harmonic tensors.
 W3 is the invariant encoding derived from the W3j product of these three tensors. 
 M is the matrix that converts the invariant coding (W3) into the original spherical harmonic tensor (V3).
 
 
      equivalentFeatures.W3jProduct(V1,V2,V3,n1,n2); 
       
This function transforms three spherical harmonic tensors into a vector that is invariant to rotation, only shell n1 for V1 is selected, and for three shells wigner production from V1,V2,V3 with degrees s1,s2,s3 satisfied that s1 = n1, and s1 < abs(s2-s3)+n2.

     equivalentFeatures.W3jProductCompact(V1,V2,V3,n1)； 
     
This function transforms three spherical harmonic tensors into a vector that is invariant to rotation，n1 is the maximum degree of V1. this function will give invariance vector with the same length of V3. 
   
     equivalentFeatures.DecodeMatrixCompact(V1,V2,V3,n1); 

This function will return a matrix that can be used to convert the invariant coding calculated by the W3jProductCompact function into one of the original spherical harmonic tensors (V1).  Because the length of V ,therefore this matrix is a square matrix.


 e.g.

 
      W3 = equivalentFeatures.W3jProductCompact(V1,V2,V3)
      M  = equivalentFeatures.DecodeMatrixCompact(V2,V3) 
      norm(M*W3-V1_,2);

The length of V1 is equal to the length of W3, therefore M is a square matrix. V2,V3 are two reference vectors. 

      equivalentFeatures.ReferencesExtract(V_input)
      
This function will yield two reference vectors from the diople and quadruple shells of V_input. 

    
      equivalentFeatures.SelfProductMatrix( V_input, int n, int n2);  
      
      
For the Clebsch–Gordan (CG) product of any two shells in the V_input spherical harmonic tensor, only top n2  with smaller degree are selected. The spherial harmonic tensors with degree n are returned.

e.g 

     V_input = [ 0.8462843753216345 + 0.0im, 0.49442005281239043 - 0.38003620975475416im,1.1266402071841877 + 0.0im,-0.49442005281239043 - 0.38003620975475416im,0.12624905371373815 - 0.586665146512387im,0.712812513540789 - 0.5922990542960181im,0.8306662145480092 + 0.0im,-0.712812513540789 - 0.5922990542960181im,0.12624905371373815 + 0.586665146512387im ]
    V_output = [ 0.5641895835477564 + 0.0im,0.046180362258809796-0.5151788499369871im,0.563454914288809,-0.046180362258809796-0.5151788499369871im,-0.35831018662719905 -0.036901839488305416im,0.09441424373213234-0.6546737533363572im,0.008536735434823653+0.0im,-0.09441424373213234-0.6546737533363572im,-0.35831018662719905+0.036901839488305416im]
      
      
      MM =  equivalentFeatures.ReferencesExtract(V_input)
      v1 = vcat(1,MM[:,1]);
      v2 = vcat(1,MM[:,2]);# MM[:,1] and MM[:,2] have been normalize.
      v2 = equivalentFeatures.RStoCS_Encode(v2,2);
      v1 = equivalentFeatures.RStoCS_Encode(v1,2);
      v1 = equivalentFeatures.IncreaseDegree(v1,1);
      v2 = equivalentFeatures.IncreaseDegree(v2,1);
      V_Encode = equivalentFeatures.W3jProductCompact(V_output,v1,v2,3);
      M = equivalentFeatures.DecodeMatrixCompact(v1,v2,3);
      norm(M*V_Encode - V_output,2);
      
      
      MM=equivalentFeatures.SelfProductMatrix(V_input,2,2)
      ReciprocalRadii   = (MM[1,2]*MM[1,2]-MM[1,1]*MM[1,3]*2)/(0.25*(3/pi));
      ReciprocalRadii   = (ReciprocalRadii)^(-(1/2));
      v1 = vcat(1,MM[1,:].*ReciprocalRadii);
      ReciprocalRadii   = (MM[2,2]*MM[2,2]-MM[2,1]*MM[2,3]*2)/(0.25*(3/pi));
      ReciprocalRadii   = (ReciprocalRadii)^(-(1/2));
      v2 = vcat(1,MM[2,:].*ReciprocalRadii);   # normalize MM[1,:] and MM[2,:]. 
      v2 = equivalentFeatures.RStoCS_Encode(v2,2);
      v1 = equivalentFeatures.RStoCS_Encode(v1,2);
      v1 = equivalentFeatures.IncreaseDegree(v1,1);
      v2 = equivalentFeatures.IncreaseDegree(v2,1);
      V_Encode = equivalentFeatures.W3jProductCompact(V_output,v1,v2,3);
      M = equivalentFeatures.DecodeMatrixCompact(v1,v2,3);
      norm(M*V_Encode - V_output,2);

ReferencesExtract,SelfProductMatrix can be used to generate the reference vectors for W3jProductCompact and DecodeMatrixCompact functions. 

      equivalentFeatures.ProductEncode(V1,V2)


This function will return the invariant coding of the tensor product of two spherical harmonic tensors (V1, V2).

e.g.

      V1 =[1.0,0.0043477849927746155,0.0,0.9999905483381614,0.11]; 
      V2 =[1.0,0.772027518982468,0.33454525822573616,0.5404192632877276];
      include("./equivalentFeatures.jl")
      using  .equivalentFeatures
      equivalentFeatures.setN(2); 
      equivalentFeatures.Initial();
      S1 = equivalentFeatures.CtoS_Encode(V1);
      S2 = equivalentFeatures.CtoS_Encode(V2);
      equivalentFeatures.ProductEncode(S1,S2);
S1 and S2 are two spherical harmonic tensors share the same centre.

       equivalentFeatures.DerivativeSH(Y)

This function will return the derivative of a Spherical Harmonic tensor with respect to $\theta$ and $\phi$.

      equivalentFeatures.DerivativeSH_XYZ(Y,DR,ReciprocalRadii,ReciprocalF)
      
This function returns the derivative of a spherical harmonic tensor with respect to Cartesian coordinates. DR is a vector whose elements are fractions of the derivative of radius function and radius function((df(r)/dr)/f(r)). ReciprocalRadii is the reciprocal of the radius in Cartesian coordinates (1/r), ReciprocalF is the reciprocal of the radius function. 

e.g.

      x=0.31
      y=0.14
      z=0.41
      r = (x^2+y^2+z^2)^0.5;
      z  = (z/r);
      x  = (x/r);
      y  = (y/r);
      V3 =[1,x,y,z]
      S = equivalentFeatures.CtoS_Encode(V3,2);
      S = equivalentFeatures.IncreaseDegree(S,1);
      S = equivalentFeatures.IncreaseDegree(S,1);
      f0 = 1;
      f1 = r;
      f2 = r*r;
      f3 = r*r*r; 
      S[1]     = S[1]*f0;
      S[2:4]   = S[2:4]*f1; 
      S[5:9]   = S[5:9]*f2; 
      S[10:16] = S[10:16]*f3; 
      equivalentFeatures.DerivativeSH(S);
      D = zeros(16);
      D[2:4] .= 1.0;  #f1'
      D[5:9] .= 2*r;  #f2'
      D[10:16] .= 3*r*r; #f3'
      ReciprocalRadii = 1/r; 
      ReciprocalF     = 1/f1;
      equivalentFeatures.DerivativeSH(Y,DR,ReciprocalRadii,ReciprocalF)
The DerivativeSH function returns a matrix with two columns. One column indicates the θ and the other indicates the φ.  DerivativeSH_XYZ returns a three-column matrix indicating X, Y, and Z, respectively. 

# Installation (c++)


   Eigen package is needed for install. 

       cd  include
       make
   
   rotation.o will be generated. 
   
   In your project file: 
   
       #include "include/rotation.h"
       equivalentFeatures Test; 
       Test.Initial();

   In your Makefile file: 

        objects = include/rotation.o 
        CC = g++  -std=c++11  -g  -I include 
        $(CC)  -o yourproject yourproject.cpp $(objects)
  
# Usage   
    std::vector<std::complex<double>> RStoCS_Encode(const std::vector<double>& V1, int n); 
    std::vector<double> CStoRS_Encode(const std::vector<std::complex<double>>& V1, int n);

RStoCS_Encode will convert a real spherical harmonic tensor into a complex spherical harmonic tensor, CStoRS_Encode will convert a complex spherical harmonic tensor back to a real spherical harmonic tensor. 
    std::vector<std::complex<double>> IncreaseDegree(const std::vector<std::complex<double>>& V, int n);
IncreaseDegree function can increase the degree of spherical harmonic tensor by n. 
e.g.

    std::vector<double> V1 = { 1.0,0.2820947917738782,0.1558348923076441,-0.41829957442835725};
    std::vector< std::complex<double> > V2 = Test.RStoCS_Encode(V1,2);
    std::vector< std::complex<double> > V3 = Test.IncreaseDegree(V2,1);
    std::vector<double> V4 = Test.CStoRS_Encode(V3,3);
V1 is a real spherical harmonic tensor to be converted as a complex spherical harmonic V2 via RStoCS_Encode function, and the complex spherical harmonic V2 is increased by one degree using IncreaseDegree function, increased spherical harmonic tensor V3 is converted back as a  complex spherical harmonic tensor V4 using CStoRS_Encode. 

    std::vector<double> SelfProduct(const std::vector<std::complex<double>>& V, int n, int n2); 
 V is a spherical harmonic tensor to be encoded. 
 Two shells of V will be Clebsch–Gordan (CG) producted to generate a set of spherical harmonic tensors,the maximum degree of spherical harmonic tensors in this set is n, and the maximum size of the set is n2. The norm of each spherical harmonic tensor from this set derived from the shells (CG) production of V is calculated as return.  
 
    std::vector<double> SelfProductPairwise(const std::vector<std::complex<double>>& V, int n, int n2);
   
 V is a spherical harmonic tensor to be encoded. 
 Two shells of V will be  Clebsch–Gordan (CG)  producted to generate a set of spherical harmonic tensors, the maximum degree of spherical harmonic tensors in this set is n, and the maximum size of the set is n2. For this set derived from the shells (CG) production of V, the norm of each spherical harmonic tensor and the product of two spherical harmonic tensors are caculated as return. 
e.g. 


    std::vector<std::complex<double>> V = {
        {0.8462843753216346, 0.0},
        {-0.05890595046982487, -0.16671610340719356},
        {0.1715792182533041, 0.0},
        {0.05890595046982487, -0.16671610340719356},
        {-0.2724876460411745, -0.10974244810053088},
        {-0.5990061628748591, 0.17027506101854945},
        {0.4559079877665936, 0.0},
        {0.5990061628748591, 0.17027506101854945},
        {-0.2724876460411745, 0.10974244810053088}
    };
    std::vector<std::complex<double>> V_R = {
        {0.8462843753216346, 0.0},
        {-0.035869334374274975, 0.19458805470434914},
        {0.11689964946194337, 0.0},
        {0.035869334374274975, 0.19458805470434914},
        {-0.40579691028618065, -0.3759821000777476},
        {0.36937292043771863, -0.25410952237718903},
        {0.37677386504280463, 0.0},
        {-0.36937292043771863, -0.25410952237718903},
        {-0.40579691028618065, 0.3759821000777476}
    };
    std::vector<double> N1 = Test.SelfProduct(V,2,2);
    std::vector<double> N2 = Test.SelfProduct(V_R,2,2);
    double norm = 0.0;
    for (size_t i = 0; i < N1.size(); ++i) {
        norm += std::abs(N1[i] - N2[i]);
    }
   
    N1 = Test.SelfProductPairwise(V,2,2);
    N2 = Test.SelfProductPairwise(V_R,2,2);
    norm = 0.0;
    for (size_t i = 0; i < N1.size(); ++i) {
        norm += std::abs(N1[i] - N2[i]);
    }
    
V is a spherical harmonic vector,  V_R is a rotated spherical harmonic vector for the V spherical harmonic vector. For SelfProduct, SelfProductPairwise functions, V and V_R will yield the same result. 
 
    std::vector<std::complex<double>> W3jProduct(const std::vector<std::complex<double>>& V1, const std::vector<std::complex<double>>& V2,const std::vector<std::complex<double>>& V3,int n1, int n2); 
   
V1, V2, V3 are spherical harmonic tensors with different orientation to be encoded with wigner 3J. Only shell n1 for V1 is selected, and for three shells wigner production from V1,V2,V3 with degrees s1,s2,s3 satisfied that s1 = n1, and s1 < abs(s2-s3)+n2. 

    std::vector<double> W3jProductCToR(const std::vector<std::complex<double>>& InvariantV, int n1,int d2, int d3, int n2);

The definition of n1 and n2 is the same of W3jProduct, d2,d3 are the dimensions of V2 and V3 that are as inputs of W3jProduct.This funcation will convert the spherical harmonic tensor that is the result of W3jProduct into a real vector. 

     std::vector<std::complex<double>> W3jProductRToC(const std::vector<double>& InvariantV, int n1,int d2, int d3, int n2);

The definition of n1,n2,d2,d3 is the same of W3jProductCToR, the input is a real vector converted by W3jProductCToR, the output is spherical harmonic tensor that is the same as input of W3jProductCToR. This function is the inverse process of the W3jProductCToR function.

      std::vector<std::complex<double>> W3jProductCompact(const std::vector<std::complex<double>>& V1,const std::vector<std::complex<double>>& V2,const std::vector<std::complex<double>>& V3,int n);

 V1, V2, V3 are spherical harmonic tensors with different orientation to be encoded with wigner 3J.  The length of result of wigner 3j product for these three spherical harmonic should be the same as V1, and this result is a rotation invariance vector. 
     
     Eigen::MatrixXcd DecodeMatrixCompact(const std::vector<std::complex<double>>& V2, const std::vector<std::complex<double>>& V3, int n);

The spherical harmonic tensors V2 and V3 have the same reference coordinate as V1 the input of W3jProductCompact. This function will return a matrix M. This matrix will recover an invariance vector yielded by W3jProductCompact into an spherical harmonic tensor that share the same reference coordinate as V2, V3, the orientation of this spherical harmonic tensor is provide by V2, V3.  
e.g.

       equivalentFeatures Test;
       Test.Initial();

       std::vector<std::complex<double>> V1 = {
            {0.2820947917738782, 0.0}, {0.1558348923076441, -0.08714600779676653},
            {-0.41829957442835725, 0.0}, {-0.1558348923076441, -0.08714600779676653},
            {0.05400984229015042, -0.08789344893420385}, {-0.29831935733068915, 0.16682618799224708},
            {0.3780897816045623, 0.0}, {0.29831935733068915, 0.16682618799224708},
            {0.05400984229015042, 0.08789344893420385}
       };
      std::vector<std::complex<double>> V2 = {
        {0.2820947917738782, 0.0}, {-0.13150538452459107, 0.20714627161985197},
        {0.34398535815919695, 0.0}, {0.13150538452459107, 0.20714627161985197},
        {-0.08289419650884133, -0.17630513479494964}, {-0.20702023603219083, 0.32609668569054673},
        {0.1535728097664055, 0.0}, {0.20702023603219083, 0.32609668569054673},
        {-0.08289419650884133, 0.17630513479494964}
    };
     std::vector<std::complex<double>> V3 = {
        {0.2820947917738782, 0.0}, {-0.08323545825287791, -0.286716367230279},
        {0.2458934345224644, 0.0}, {0.08323545825287791, -0.286716367230279},
        {-0.24360329182248358, 0.15445613562862262}, {-0.09366656951197919, -0.32264781266424436},
        {-0.07575460360437422, 0.0}, {0.09366656951197919, -0.32264781266424436},
        {-0.24360329182248358, -0.15445613562862262}
    };
     std::vector<std::complex<double>> V1_ = {
        {0.2820947917738782, 0.0}, {0.1294158205116877, -0.22846248270670816},
        {0.3175614797663446, 0.0}, {-0.1294158205116877, -0.22846248270670816},
        {-0.11470689480471548, -0.19135797511692804}, {0.18808081324411152, -0.3320259406721129},
        {0.08429133903919787, 0.0}, {-0.18808081324411152, -0.3320259406721129},
        {-0.11470689480471548, 0.19135797511692804}
    };
    std::vector<std::complex<double>> V2_ = {
        {0.2820947917738782, 0.0}, {-0.11794378575971845, 0.12037014272293214},
        {-0.42653604644964793, 0.0}, {0.11794378575971845, 0.12037014272293214},
        {-0.0018711923316910401, -0.09188362988533172}, {0.23022904373214756, -0.23496534959002102},
        {0.4056684557059928, 0.0}, {-0.23022904373214756, -0.23496534959002102},
        {-0.0018711923316910401, 0.09188362988533172}
    };
     std::vector<std::complex<double>> V3_ = {
        {0.2820947917738782, 0.0}, {-0.04734136912624448, 0.3026803946881253},
        {0.22587421614524694, 0.0}, {0.04734136912624448, 0.3026803946881253},
        {-0.28921882314977415, -0.09274049507548826}, {-0.04893693653854026, 0.3128817678849446},
        {-0.11318592970238676, 0.0}, {0.04893693653854026, 0.3128817678849446},
        {-0.28921882314977415, 0.09274049507548826}
    };

    std::vector<std::complex<double>> V__ = Test.W3jProductCompact(V1,V2,V3,3); 
    Eigen::MatrixXcd M = Test.DecodeMatrixCompact(V2_,V3_,3); 
    std::vector<std::complex<double>> result(M.rows(), 0);

    for (size_t i = 0; i < M.rows(); ++i) {
        for (size_t j = 0; j < M.cols(); ++j) {
            result[i] += M(i,j) * V__[j];
          }
    }


    double norm = 0.0;
    for (size_t i = 0; i < result.size(); ++i) {
        norm += (std::pow(result[i] - V1_[i], 2)).real();
    }
    norm = std::sqrt(norm);
  V1, V2, V3  are same as  V1_, V2_, V3_  in rotation equivalention.  V1 is embedded as rotation invariant vector V__ using V2, V3.  V2_, V3_ generated a auxiliary matrix that can convert invariant vector V__ into a spherical harmonic tensor that is the same as V1_ and rotational equivalent to V1. 
  
       std::vector<double> W3jProductCompactCToR(const std::vector<std::complex<double>>& InvariantV, int n);
       std::vector<std::complex<double>> W3jProductCompactRToC(const std::vector<double>& InvariantV, int n);
    
W3jProductCompactRToC will change the result of W3jProductCompact into real vector, and W3jProductCompactRToC will recover the this real vector into original spherical harmonic tensor. 

       std::vector< std::vector<double> > ReferencesExtract(const  std::vector<double>& VR);
VR is a real spherical harmonic tensor, two reference vectors are extracted from dipole and quadrupole shells of this real spherical harmonic tensor using ReferencesExtract.  

       std::vector<std::vector<std::complex<double>>> SelfProductMatrix(const std::vector<std::complex<double>>& V, int n, int n2);
V is a complex spherical harmonic tensor, the pair of shells from this spherical harmonic tensor are producted as a set of spherical harmonic tensors, the maximum size of each set is n2, in each set, only spherical harmonic tensors with degree n are returned.   

The reference vectors for the W3jProductCompact and DecodeMatrixCompact functions can be generated by either the ReferencesExtract or the SelfProductMatrix functions.
e.g.

       std::vector<std::complex<double>> V_input = {
        {0.5641895835477564, 0.0},
        {0.02432950778305304, 0.12000026382308544},
        {-0.0743142162691603, 0.0},
        {-0.02432950778305304, 0.12000026382308544},
        {-0.028884354218690914, -0.2641985837291535},
        {-0.5053395933628799, 0.4929228736827938},
        {0.5316625913709678, 0.0},
        {0.5053395933628799, 0.4929228736827938},
        {-0.028884354218690914, 0.2641985837291535}
    };
    std::vector<std::complex<double>> V_output = {
        {0.5641895835477564, 0.0},
        {0.046180362258809796, -0.5151788499369871},
        {0.563454914288809, 0.0},
        {-0.046180362258809796, -0.5151788499369871},
        {-0.35831018662719905, -0.036901839488305416},
        {0.09441424373213234, -0.6546737533363572},
        {0.008536735434823653, 0.0},
        {-0.09441424373213234, -0.6546737533363572},
        {-0.35831018662719905, 0.036901839488305416}
    };
    /* reference vectors extraction via ReferencesExtract function*/
    std::vector<std::complex<double>>  V_input_ =  V_input;
    std::vector< std::vector<double> > V__ = Test.ReferencesExtract(Test.CStoRS_Encode(V_input_,3));
    std::vector<double> v1_r = {1};
    v1_r.insert(v1_r.end(), V__[0].begin(), V__[0].end()); 
    std::vector<double> v2_r = {1};
    v2_r.insert(v2_r.end(), V__[1].begin(), V__[1].end()); 
    std::vector<std::complex<double>>  v1 = Test.RStoCS_Encode(v1_r,2); 
    std::vector<std::complex<double>>  v2 = Test.RStoCS_Encode(v2_r,2); 
    v1 = Test.IncreaseDegree(v1,1);
    v2 = Test.IncreaseDegree(v2,1);
    std::vector<std::complex<double>> V_Encode = Test.W3jProductCompact(V_output,v1,v2,3);//Embedding as rotation-invariant coding 
    Eigen::MatrixXcd M = Test.DecodeMatrixCompact(v1,v2,3);//Generate an auxiliary matrix for recovering embedded. 
    std::vector<std::complex<double>> V_Re(M.rows(), 0);
    for (size_t i = 0; i < M.rows(); ++i) {
        for (size_t j = 0; j < M.cols(); ++j) {
            V_Re[i] += M(i,j) * V_Encode[j];
          }
    }// V_Re=M*V_Encode;
    double norm = 0.0;
    for (size_t i = 0; i < V_Re.size(); ++i) {
        norm += (std::pow(V_Re[i] - V_output[i], 2)).real();
    }
    norm = std::abs(norm);
    /* reference vectors extraction via SelfProductMatrix function*/
    std::vector<std::vector<std::complex<double>>>  V_ = Test.SelfProductMatrix(V_input,2,2);  
    v1.resize(1); 
    v2.resize(1);
    v1[0] = {1.0, 0.0}; 
    v2[0] = {1.0, 0.0};
    v1.insert(v1.end(), V_[0].begin(), V_[0].end());
    v2.insert(v2.end(), V_[1].begin(), V_[1].end());
    v1 = Test.IncreaseDegree(v1,1);
    v2 = Test.IncreaseDegree(v2,1);
    V_Encode = Test.W3jProductCompact(V_output,v1,v2,3);//Embedding as rotation-invariant coding 
    M = Test.DecodeMatrixCompact(v1,v2,3);//Generate an auxiliary matrix for recovering embedded. 
    equivalentFeatures.DecodeMatrixCompact(v1,v2,3);
    std::fill(V_Re.begin(), V_Re.end(),0); 
    for (size_t i = 0; i < M.rows(); ++i) {
        for (size_t j = 0; j < M.cols(); ++j) {
            V_Re[i] += M(i,j) * V_Encode[j];
          }
    }// V_Re=M*V_Encode; 
    norm = 0.0;
    for (size_t i = 0; i < V_Re.size(); ++i) {
        norm += (std::pow(V_Re[i] - V_output[i], 2)).real();
    }
    norm = std::abs(norm);

Reference vectors are extracted from input spherical harmonic tensors (V_input) via the ReferencesExtract or SelfProductMatrix function. The spherical harmonic tensor (V_output) is converted into a rotation-invariant vector, which is insensitive to the rotation of (V_input) and (V_output), using the W3jProductCompact function and the reference vectors. This rotation-invariant vector is then recovered as a spherical harmonic tensor with coordinate information using an auxiliary matrix generated by the DecodeMatrixCompact function and the reference vectors.

     std::vector<double> ProductEncode(const std::vector<std::complex<double>>& V1, const std::vector<std::complex<double>>& V2, int n);
 V1 and V2 are spherical harmonic tensors to be encoded, one shell of V1 and one shell of V2 will be Clebsch–Gordan (CG) products to generate a set of spherical harmonic tensors, the maximum size of each set is n, the output is a list of the norms of each spherical harmonic tensor from these sets.    
 
    std::vector<double> ProductEncodePairwise(const std::vector<std::complex<double>>& V1, const std::vector<std::complex<double>>& V2, int n);

 V1 and V2 are spherical harmonic tensors to be encoded,  one shell of V1 and one shell of V2 will be Clebsch–Gordan (CG) producted to generate a set of spherical harmonic tensor, the maximum size of each set is n, the return is a list for the norm of each spherical harmonic tensors and pairwise products among spherical harmonic tensors from these sets.
e.g.


    std::vector<std::complex<double>> Va = {
        {0.5641895835477564, 0.0},
        {0.02432950778305304, 0.12000026382308544},
        {-0.0743142162691603, 0.0},
        {-0.02432950778305304, 0.12000026382308544},
        {-0.028884354218690914, -0.2641985837291535},
        {-0.5053395933628799, 0.4929228736827938},
        {0.5316625913709678, 0.0},
        {0.5053395933628799, 0.4929228736827938},
        {-0.028884354218690914, 0.2641985837291535}
    };
    std::vector<std::complex<double>> Vb = {
        {0.5641895835477564, 0.0},
        {-0.21474084277746897, -0.07957009561042702},
        {0.5898787926816613, 0.0},
        {0.21474084277746897, -0.07957009561042702},
        {-0.3264974883313249, -0.021848999166327016},
        {-0.30068680554417004, 0.003448873026302368},
        {0.07781820616203128, 0.0},
        {0.30068680554417004, 0.003448873026302368},
        {-0.3264974883313249, 0.021848999166327016}
    };
    std::vector<std::complex<double>> Va_ = {
        {0.5641895835477564, 0.0},
        {0.011472034751969437, -0.1080923399837762},
        {-0.10897456668330358, 0.0},
        {-0.011472034751969437, -0.1080923399837762},
        {-0.11657808713640651, -0.2832416050022595},
        {0.41830985697625883, -0.5669912902621337},
        {0.4899597947451914, 0.0},
        {-0.41830985697625883, -0.5669912902621337},
        {-0.11657808713640651, 0.2832416050022595}
    };
    std::vector<std::complex<double>> Vb_ = {
        {0.5641895835477564, 0.0},
        {-0.16528515488596268, 0.4230505374110573},
        {-0.2006618303044012, 0.0},
        {0.16528515488596268, 0.4230505374110573},
        {-0.2910900154814652, -0.1846241249608196},
        {0.1812921071936071, 0.07791641829492388},
        {0.29248252600360675, 0.0},
        {-0.1812921071936071, 0.07791641829492388},
        {-0.2910900154814652, 0.1846241249608196}
    };
    std::vector<double> N1 = Test.ProductEncode(Va,Vb,2);
    std::vector<double> N2 = Test.ProductEncode(Va_,Vb_,2);
    double norm = 0.0;
    for (size_t i = 0; i < N1.size(); ++i) {
        norm += std::abs(N1[i] - N2[i]);
    }

    N1 = Test.ProductEncodePairwise(Va,Vb,2);
    N2 = Test.ProductEncodePairwise(Va_,Vb_,2);
    norm = 0.0;
    for (size_t i = 0; i < N1.size(); ++i) {
        norm += std::abs(N1[i] - N2[i]);
    }

  Va,Vb are same as Va_,Vb_ in rotational equivariant way,  ProductEncode and ProductEncodePairwise functions are applied in Va,Vb and Va_,Vb_, the same results are yielded. 

     vector< vector< complex<double> > > DerivativeSH(const vector< complex<double>>& Y) 
     
Y is the spherical harmonic tensors for obtaining the derivatives, the result is the derivatives of this spherical harmonic tensors respect to $\theta$ and $\phi$.   

    vector< vector< complex<double > > > DerivativeSH_XYZ(const vector< complex<double>>& Y, const vector< complex<double>>& DR, double ReciprocalRadii, double ReciprocalF)
    
  The S2 spherical harmonic tensors is the input for the the derivatives extraction. D is a list in which each element is the ratio of the derivative of the radii function to the radii function itself. 'ReciprocalRadii' is the multiplicative inverse of 'Radii', and 'ReciprocalF' is the multiplicative inverse of the 'Radii Function'.  The returns are the derivatives of spherical harmonics with respect to the Cartesian coordinates X, Y and Z.
  e.g.

  
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



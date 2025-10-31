module  equivalentFeatures

          using SphericalHarmonicExpansions
          using WignerSymbols
          using DataStructures
          using MultivariatePolynomials
          using LinearAlgebra
          using SphericalHarmonics

          struct Index
            m1::Int16
            m2::Int16
          end

          N=3;
          CTS_Closed = false;# Whether Cartesian to spherical convertion setting is closed

          #CGTableV = fill(Vector{Float64}(undef,0),(d1,d2,d3));
          export CGTableC
          export CGTableI
          export W3JTableC
          export W3JTableI
          export SelfPI
          export SelfPI_test
          export WignerPI
          export CtoS_C
          export  S_cache

          function setN(N_new)
            global N = N_new;
          end

          function setCTS_Closed( Adjust )
             global CTS_Closed  = Adjust;
          end


          function Initial()

             d1 = N
             d2 = N
             d3 = abs(d1+d2-2)-abs(d1-d2)+1
             global S_cache = SphericalHarmonics.cache(N-1);
             global CGTableC = fill(Vector{ Vector{Float64} }(undef,0),(d1,d2,d3));
             global CGTableI = fill(Vector{ Vector{Index} }(undef,0),(d1,d2,d3));



             for J1 = 1:d1
               for J2 = J1:d2
                  for j3 = (J2-J1):(J1+J2-2)
                    j1 = J1-1;
                    j2 = J2-1;
                    J  = j3-(J2-J1)+1;

                    # CGTableV[J1,J2,J+1]=Vector{Float64}(undef,2*J+1);
                    CGTableC[J1,J2,J]=Vector{ Vector{Float64} }(undef,2*j3+1);
                    CGTableI[J1,J2,J]=Vector{ Vector{Index} }(undef,2*j3+1);

                    for i = 1:(2*j3+1)
                      CGTableC[J1,J2,J][i]=Vector{Float64}(undef,0);
                      CGTableI[J1,J2,J][i]=Vector{Index}(undef,0);
                    end

                    for m1 = -j1:j1
                       for m2 = -j2:j2
                          m3=m2+m1;
                          if abs(m3) <= j3
                             i  = m3+j3+1;
                             push!(CGTableI[J1,J2,J][i],Index(m1,m2));
                             push!(CGTableC[J1,J2,J][i],clebschgordan(Float64,j1,m1,j2,m2,j3,m1+m2));
                          end
                       end
                    end
                  end
               end
             end



             d3 = N
             global W3JTableC = fill(Vector{Float64}(undef,0),(d1,d2,d3));
             global W3JTableI = fill(Vector{Index}(undef,0),(d1,d2,d3));

             for J1 = 1:d1
                for J2 = 1:d2 #for J2 = J1:d2
                   for J3 = abs(J2-J1):abs(J1+J2-2)
                      if J3 < d3
                         j1 = J1-1;
                         j2 = J2-1;
                         j3 = J3;
                         J3 = J3+1;

                         W3JTableC[J1,J2,J3]=Vector{Float64}(undef,0);
                         W3JTableI[J1,J2,J3]=Vector{Index}(undef,0);


                         for m1 = -j1:j1
                            for m2 = -j2:j2
                               m3 =-(m1+m2);
                               if abs(m3) <= j3
                                  push!(W3JTableI[J1,J2,J3],Index(m1,m2));
                                  push!(W3JTableC[J1,J2,J3],wigner3j(Float64,j1,j2,j3,m1,m2,m3));
                               end
                            end
                         end

                      end
                   end
                end
             end


             dictT = Dict{Index,Float64}();
             empty!(dictT);
             List  = Any[];
             AlwaysZero=true;



             for J1 = 1:d1
                for J2 = J1:d2
                   for j3 = (J2-J1):(J1+J2-2)
                      J3 = j3-(J2-J1)+1;

                     if J1 == J2
                         AlwaysZero=true;

                         for I = 1:length(CGTableC[J1,J2,J3])

                                  empty!(dictT);

                                  for I2 = 1:length(CGTableC[J1,J2,J3][I])

                                        m1 = max(CGTableI[J1,J2,J3][I][I2].m1,CGTableI[J1,J2,J3][I][I2].m2);
                                        m2 = min(CGTableI[J1,J2,J3][I][I2].m1,CGTableI[J1,J2,J3][I][I2].m2);
                                        M_Index = Index(m1,m2);
                                        if haskey(dictT,M_Index) == false
                                          dictT[M_Index] = 0.0;
                                        end
                                          dictT[M_Index] = dictT[M_Index]+CGTableC[J1,J2,J3][I][I2];
                                  end

                                  for (k, v) in dictT
                                     if v != 0
                                        AlwaysZero=false;
                                     end
                                  end
                         end
                     end

                    if (AlwaysZero != true) && (length(CGTableC[J1,J2,J3])!=0)
                         push!(List,[J1,J2,J3]);
                    end

                   end
                end
             end

             global SelfPI   = Array{Int}(undef,length(List),3);


             for I = 1:length(List)
               SelfPI[I,:]=List[I];
             end

             List = Any[];
                for J1 = 1:d1
                  for J2 = 1:d2 # for J2 = J1:d2
                      for J3 = abs(J2-J1):abs(J1+J2-2)
                         if J3 < d3

                            J3 = J3+1;
                            if (length(W3JTableI[J1,J2,J3]) != 0)
                               push!(List,[J1,J2,J3]);
                            end
                         end
                      end
                   end
                end

             global WignerPI = Array{Int}(undef,length(List),3);

             for I = 1:length(List)
               WignerPI[I,:]=List[I];
             end

             if CTS_Closed == true
                return;
             end

             @polyvar x y z
             global CtoS_C = Vector{ Vector{Float64} }(undef,N*N);# Cartesian to spherical
             I_= 0;
             for J = 1:N
               j=J-1;
               dict = DictInitial(j);
               for m =-j:j
                 I_ = I_+1;
                 p=rlylm(j,m,x,y,z);
                 CtoS_C[I_] =  zeros(length(dict));
                 List = collect(terms(p));
                   for I = 1:length(List)
                     i1 = degree(monomial(List[I]),x);
                     i2 = degree(monomial(List[I]),y);
                     i3 = degree(monomial(List[I]),z);
                     CtoS_C[I_][dict[[i1,i2,i3]]]=coefficient(List[I]);
                   end
               end
             end
          end

         function DictInitial(N)

             s = Stack{Int}();
             push!(s,2);
             Index = 0;
             dict = Dict{Vector{UInt32},Int64}();

             if N<=0
               dict[[0,0,0]]=1.0;
               return dict;
             end

             while length(s)!= 0
                n = first(s);
                if length(s) < N
                  if n >= 0
                     push!(s,n);
                  else
                     pop!(s);
                     if(length(s)== 0)
                       break;
                     end
                     n=pop!(s);
                     push!(s,n-1);
                  end
                elseif length(s) == N
                  Index = Index+1;
                  x=0;y=0;z=0;
                  for i in Iterators.reverse(s)
                     if i == 2
                       x=x+1;
                     elseif i == 1
                       y=y+1;
                     elseif i == 0
                       z=z+1;
                     end
                  end
                     dict[[x,y,z]]=Index;
                  n=pop!(s);
                  n=n-1;
                    if n>=0
                      push!(s,n);
                    else
                      if(length(s)== 0)
                        break;
                      end
                      n=pop!(s);
                      push!(s,n-1);
                    end
                end
             end
             return dict;
         end

         function SelfProduct(V)
             RInvariantV = Complex.(zeros(size(SelfPI)[1]));

             for I =1:size(SelfPI)[1]

                 TempI = CGTableI[SelfPI[I,1],SelfPI[I,2],SelfPI[I,3]];
                 TempC = CGTableC[SelfPI[I,1],SelfPI[I,2],SelfPI[I,3]];
                 J1 = SelfPI[I,1]-1;
                 J2 = SelfPI[I,2]-1;

                 for I_1 = 1:length(TempI)
                   Value = Complex.(0);
                   for I_2 = 1:length(TempI[I_1])
                      m1 = TempI[I_1][I_2].m1;
                      m2 = TempI[I_1][I_2].m2;
                      I1 = (J1^2 + m1+J1+1);
                      I2 = (J2^2 + m2+J2+1);
                      Value=Value+V[I1]*V[I2]*TempC[I_1][I_2];
                   end
                   RInvariantV[I] = RInvariantV[I] + Value*conj(Value);
                 end
             end
             return RInvariantV;
         end

         function SelfProduct(V,n)
             size__ = 0;

             for I in 1:size(SelfPI)[1]
                if (equivalentFeatures.SelfPI[I,1]<= n && equivalentFeatures.SelfPI[I,2]<= n)
                    size__=size__ +1;
                end
             end

             RInvariantV = Complex.(zeros(size__));
             Index = 1;
             for I =1:size(SelfPI)[1]

                if (equivalentFeatures.SelfPI[I,1]<= n && equivalentFeatures.SelfPI[I,2]<= n)

                 TempI = CGTableI[SelfPI[I,1],SelfPI[I,2],SelfPI[I,3]];
                 TempC = CGTableC[SelfPI[I,1],SelfPI[I,2],SelfPI[I,3]];
                 J1 = SelfPI[I,1]-1;
                 J2 = SelfPI[I,2]-1;

                 for I_1 = 1:length(TempI)
                   Value = Complex.(0);
                   for I_2 = 1:length(TempI[I_1])
                      m1 = TempI[I_1][I_2].m1;
                      m2 = TempI[I_1][I_2].m2;
                      I1 = (J1^2 + m1+J1+1);
                      I2 = (J2^2 + m2+J2+1);
                      Value=Value+V[I1]*V[I2]*TempC[I_1][I_2];
                   end
                   RInvariantV[Index] = RInvariantV[Index] + Value*conj(Value);
                 end
                  Index = Index +1;
               end
            end
             return RInvariantV;
         end

              function SelfProduct(V,n,n2)

                 size__ = 0;
                 d2 = Int32(floor((length(V))^0.5));
                 d1 = n;


                 for J1 = 1:d1
                      for J2 = max(J1-d2+1,1):min(d2-1+J1,d2)
                                 Begin = max((J1-J2+1),J2,(J2-n2+J1));
                                 End   = min(J1+J2-1,d2);
                   
                                 for I in Begin:End
                                      if ((I == J2)&&((J1%2)==0))
                                          continue;
                                      end 
                                      size__ = size__ + 1;
                                 end
                       end
                 end

                 RInvariantV = zeros(size__);
                 Index = 1;

                 for J1 = 1:d1
                      for J2 = max(J1-d2+1,1):min(d2-1+J1,d2)
                                 Begin = max((J1-J2+1),J2,(J2-n2+J1));
                                 End   = min(J1+J2-1,d2);

                                 for I in Begin:End

                                   j1 = J1-1;
                                   j2 = J2-1;
                                   j3 = I -1;
               			               J  = j1 -(I-J2)+1;
                                  #print( [J1,J2,I]);
                                  #print("\n");
                                   TempI = CGTableI[J2,I,J];
                                   TempC = CGTableC[J2,I,J];
                                   if (((I == J2)&&((J%2)==0)))
                                        continue;
                                   end 

                                   for I_1 = 1:length(TempI)
                                        Value = Complex.(0);
                                      for I_2 = 1:length(TempI[I_1])
                                        m1 = TempI[I_1][I_2].m1;
                                        m2 = TempI[I_1][I_2].m2;
                                        I1 = (j2^2 + m1+j2+1);
                                        I2 = (j3^2 + m2+j3+1);
                                        Value=Value+V[I1]*V[I2]*TempC[I_1][I_2];
                                      end
                                         RInvariantV[Index] = RInvariantV[Index] + real(Value*conj(Value));
                                    end
                                         Index = Index+1;
                                end
                        end
                    end

              return RInvariantV;

              end




                   function  W3jProduct(V1,V2,V3)
                           RInvariantV = Complex.(zeros(size(WignerPI)[1]));

                           for I =1:size(WignerPI)[1]
                              TempI = W3JTableI[WignerPI[I,1],WignerPI[I,2],WignerPI[I,3]];
                              TempC = W3JTableC[WignerPI[I,1],WignerPI[I,2],WignerPI[I,3]];
                              J1 = WignerPI[I,1]-1;
                              J2 = WignerPI[I,2]-1;
                              J3 = WignerPI[I,3]-1;
                              Value = Complex.(0);
                              for I_1 = 1:length(TempI)
                                 m1 = TempI[I_1].m1;
                                 m2 = TempI[I_1].m2;
                                 m3 = -(m1+m2);
                                 I1 = (J1^2 + m1+J1+1);
                                 I2 = (J2^2 + m2+J2+1);
                                 I3 = (J3^2 + m3+J3+1);
                                 Value = Value + V1[I1]*V2[I2]*V3[I3]*TempC[I_1];
                              end
                              RInvariantV[I]=Value;
                           end
                           return RInvariantV;
                    end


                    function  W3jProductCToR(InvariantV)
                           RInvariantV = zeros(size(WignerPI)[1]);

                           for I =1:size(WignerPI)[1]
                              RInvariantV[I]= real(im^(WignerPI[I,1]-WignerPI[I,2]-WignerPI[I,3]+1.0)*InvariantV[I]);
                           end
                           return RInvariantV;
                    end


            function  W3jProductRToC(InvariantV)
                   CInvariantV = Complex.(zeros(size(WignerPI)[1]));

                   for I =1:size(WignerPI)[1]
                      CInvariantV[I]= im^(-WignerPI[I,1]+WignerPI[I,2]+WignerPI[I,3]-1.0)*InvariantV[I];
                   end
                   return CInvariantV;
            end







             function  W3jProduct(V1,V2,V3,n)

                RInvariantVSize = 0;

                for I in 1:n
                   RInvariantVSize = RInvariantVSize +(min((n-1),(n+I-2))-abs(n-I))+1;
                end

                RInvariantV = Complex.(zeros(RInvariantVSize));
                Index =1;
                 for I =1:size(WignerPI)[1]
                   if WignerPI[I,1] != n
                      continue;
                   end
                   if WignerPI[I,2] > n || WignerPI[I,3] > n
                      continue;
                   end
                    TempI = W3JTableI[WignerPI[I,1],WignerPI[I,2],WignerPI[I,3]];
                    TempC = W3JTableC[WignerPI[I,1],WignerPI[I,2],WignerPI[I,3]];
                    J1 = WignerPI[I,1]-1;
                    J2 = WignerPI[I,2]-1;
                    J3 = WignerPI[I,3]-1;
                    Value = Complex.(0);
                    for I_1 = 1:length(TempI)
                       m1 = TempI[I_1].m1;
                       m2 = TempI[I_1].m2;
                       m3 = -(m1+m2);
                       I1 = (J1^2 + m1+J1+1);
                       I2 = (J2^2 + m2+J2+1);
                       I3 = (J3^2 + m3+J3+1);
                       Value = Value + V1[I1]*V2[I2]*V3[I3]*TempC[I_1];
                    end
                    RInvariantV[Index]=Value;
                    Index = Index+1;
                 end
                 return RInvariantV;
             end


        function  W3jProductCToR(InvariantV,n)

                    RInvariantVSize = 0;

                    for I in 1:n
                       RInvariantVSize = RInvariantVSize +(min((n-1),(n+I-2))-abs(n-I))+1;
                    end

                    RInvariantV = zeros(RInvariantVSize);
                    Index =1;
                     for I =1:size(WignerPI)[1]
                       if WignerPI[I,1] != n
                          continue;
                       end
                       if WignerPI[I,2] > n || WignerPI[I,3] > n
                          continue;
                       end
                        RInvariantV[Index]=real(im^(WignerPI[I,1]-WignerPI[I,2]-WignerPI[I,3]+1.0)*InvariantV[Index]);
                        Index = Index+1;
                     end
                     return RInvariantV;
                 end

        function  W3jProductRToC(InvariantV,n)

                    CInvariantVSize = 0;

                    for I in 1:n
                       CInvariantVSize = CInvariantVSize +(min((n-1),(n+I-2))-abs(n-I))+1;
                    end

                    CInvariantV = Complex.(zeros(CInvariantVSize));
                    Index =1;
                     for I =1:size(WignerPI)[1]
                       if WignerPI[I,1] != n
                          continue;
                       end
                       if WignerPI[I,2] > n || WignerPI[I,3] > n
                          continue;
                       end
                        CInvariantV[Index]=im^(-WignerPI[I,1]+WignerPI[I,2]+WignerPI[I,3]-1.0)*InvariantV[Index];
                        Index = Index+1;
                     end
                     return CInvariantV;
                 end


        function  W3jProductCToR(InvariantV,n)

                    RInvariantVSize = 0;

                    for I in 1:n
                       RInvariantVSize = RInvariantVSize +(min((n-1),(n+I-2))-abs(n-I))+1;
                    end

                    RInvariantV = zeros(RInvariantVSize);
                    Index =1;
                     for I =1:size(WignerPI)[1]
                       if WignerPI[I,1] != n
                          continue;
                       end
                       if WignerPI[I,2] > n || WignerPI[I,3] > n
                          continue;
                       end
                        RInvariantV[Index]= real(im^(WignerPI[I,1]-WignerPI[I,2]-WignerPI[I,3]+1.0)*InvariantV[Index]);
                        Index = Index+1;
                     end
                     return RInvariantV;
                 end




        function  W3jProduct(V1,V2,V3,n1,n2)

                    RInvariantVSize = 0;
                    d2 = Int32(floor((length(V2))^0.5));
                    d3 = Int32(floor((length(V3))^0.5));
                


                    for I in 1:d2
                       RInvariantVSize = RInvariantVSize + max(min((d3-1),(n1+I-2))-abs(n1-I)+1,0) -  max(0, min(I+n1-n2-1,min((d3-1),(n1+I-2))+1) - max(I+n2-n1+1,abs(n1-I)+1)+1);
                    end

                    RInvariantV = Complex.(zeros(RInvariantVSize));
                    Index =1;
                    for I =1:size(WignerPI)[1]
                       if WignerPI[I,1] != n1
                          if WignerPI[I,1] > n1
                             break;
                          else
                             continue;
                          end
                       end
                       if  WignerPI[I,1] > abs(WignerPI[I,2]-WignerPI[I,3])+n2 
                          continue; 
                       end 
                       if WignerPI[I,2] > d2 || WignerPI[I,3] > d3
                          continue;
                       end
                       print([WignerPI[I,1],WignerPI[I,2],WignerPI[I,3]]);
                       print("\n");

                       TempI = W3JTableI[WignerPI[I,1],WignerPI[I,2],WignerPI[I,3]];
                       TempC = W3JTableC[WignerPI[I,1],WignerPI[I,2],WignerPI[I,3]];
                       J1 = WignerPI[I,1]-1;
                       J2 = WignerPI[I,2]-1;
                       J3 = WignerPI[I,3]-1;
                       Value = Complex.(0);
                        for I_1 = 1:length(TempI)
                           m1 = TempI[I_1].m1;
                           m2 = TempI[I_1].m2;
                           m3 = -(m1+m2);
                           I1 = (J1^2 + m1+J1+1);
                           I2 = (J2^2 + m2+J2+1);
                           I3 = (J3^2 + m3+J3+1);
                           Value = Value + V1[I1]*V2[I2]*V3[I3]*TempC[I_1];
                        end
                        RInvariantV[Index]=Value;
                        Index = Index+1;
                     end
                     return RInvariantV;
                 end

        function  W3jProductCToR(InvariantV,n1,n2,d2,d3)

                    RInvariantVSize = 0;

                    for I in 1:d2
                       RInvariantVSize = RInvariantVSize + max(min((d3-1),(n1+I-2))-abs(n1-I)+1,0) -  max(0, min(I+n1-n2-1,min((d3-1),(n1+I-2))+1) - max(I+n2-n1+1,abs(n1-I)+1)+1);
                    end

                    RInvariantV = zeros(RInvariantVSize);
                    Index =1;
                     for I =1:size(WignerPI)[1]
                       if WignerPI[I,1] != n1
                          if WignerPI[I,1] > n1
                             break;
                          else
                             continue;
                          end
                       end
                       if  WignerPI[I,1] > abs(WignerPI[I,2]-WignerPI[I,3])+n2
                          continue;
                       end 

                       if WignerPI[I,2] > d2 || WignerPI[I,3] > d3
                          continue;
                       end
                       RInvariantV[Index]= real(im^(WignerPI[I,1]-WignerPI[I,2]-WignerPI[I,3]+1.0)*InvariantV[Index]);
                       Index = Index+1;
                     end
                     return RInvariantV;
                 end



        function  W3jProductRToC(InvariantV,n1,n2,d2,d3)

                    CInvariantVSize = 0;

                    for I in 1:d2
                       CInvariantVSize = CInvariantVSize + max(min((d3-1),(n1+I-2))-abs(n1-I)+1,0) -  max(0, min(I+n1-n2-1,min((d3-1),(n1+I-2))+1) - max(I+n2-n1+1,abs(n1-I)+1)+1);
                    end

                    CInvariantV = Complex.(zeros(CInvariantVSize));
                    Index =1;
                     for I =1:size(WignerPI)[1]
                       if WignerPI[I,1] != n1
                          if WignerPI[I,1] > n1
                             break;
                          else
                             continue;
                          end
                       end
                       if  WignerPI[I,1] > abs(WignerPI[I,2]-WignerPI[I,3])+n2
                          continue;
                       end 
                       if WignerPI[I,2] > d2 || WignerPI[I,3] > d3
                          continue;
                       end
                        CInvariantV[Index]= im^(-WignerPI[I,1]+WignerPI[I,2]+WignerPI[I,3]-1.0)*InvariantV[Index];
                        Index = Index+1;
                     end
                     return CInvariantV;
                 end


        function  W3jProductCompact(V1,V2,V3,n)

             RInvariantVSize = n*n;


                    RInvariantV = Complex.(zeros(RInvariantVSize));
                    Index =1;
                    for I =1:size(WignerPI)[1]
                       if  WignerPI[I,1] > n   || WignerPI[I,2] > WignerPI[I,1] || WignerPI[I,3] > WignerPI[I,1]
                           continue;
                       end
                       if (WignerPI[I,2] + WignerPI[I,3]) >= (WignerPI[I,1] +3)
                          continue;
                       end

                        TempI = W3JTableI[WignerPI[I,1],WignerPI[I,2],WignerPI[I,3]];
                        TempC = W3JTableC[WignerPI[I,1],WignerPI[I,2],WignerPI[I,3]];
                        J1 = WignerPI[I,1]-1;
                        J2 = WignerPI[I,2]-1;
                        J3 = WignerPI[I,3]-1;
                        Value = Complex.(0);
                        for I_1 = 1:length(TempI)
                           m1 = TempI[I_1].m1;
                           m2 = TempI[I_1].m2;
                           m3 = -(m1+m2);
                           I1 = (J1^2 + m1+J1+1);
                           I2 = (J2^2 + m2+J2+1);
                           I3 = (J3^2 + m3+J3+1);
                           Value = Value + V1[I1]*V2[I2]*V3[I3]*TempC[I_1];
                        end
                        RInvariantV[Index]=Value;
                        Index = Index+1;
                     end
                     return RInvariantV;
        end

        function  W3jProductCompactCToR(InvariantV,n)
           RInvariantVSize = n^2;
           RInvariantV = Complex.(zeros(RInvariantVSize));
           Index =1;
           for I =1:size(WignerPI)[1]
               if  WignerPI[I,1] > n   || WignerPI[I,2] > WignerPI[I,1] || WignerPI[I,3] > WignerPI[I,1]
                      continue;
                end
               if (WignerPI[I,2] + WignerPI[I,3]) >= (WignerPI[I,1] +3)
                      continue;
                end

                RInvariantV[Index]= real(im^(WignerPI[I,1]-WignerPI[I,2]-WignerPI[I,3]+1.0)*InvariantV[Index]);

                Index = Index+1;
           end
           return  RInvariantV;
        end


        function  W3jProductCompactRToC(InvariantV,n)

              CInvariantVSize = n^2;

              CInvariantV = Complex.(zeros(CInvariantVSize));
                    Index =1;
                     for I =1:size(WignerPI)[1]
                       if  WignerPI[I,1] > n || WignerPI[I,2] > WignerPI[I,1] || WignerPI[I,3] > WignerPI[I,1]
                          continue;
                       end
                       if (WignerPI[I,2] + WignerPI[I,3]) >= (WignerPI[I,1] +3)
                          continue;
                       end

                        CInvariantV[Index]= im^(-WignerPI[I,1]+WignerPI[I,2]+WignerPI[I,3]-1.0)*InvariantV[Index];
                        Index = Index+1;
                     end
                     return CInvariantV;
        end





         function  CtoS_Encode( V1 )
            V2 = zeros(length(CtoS_C));
            Base_ = 0;
            PreviousLength = 0;
            for I = 1:length(CtoS_C)

               if length(CtoS_C[I]) != PreviousLength
                  Base_ = Base_+length(CtoS_C[I]);
                  PreviousLength = length(CtoS_C[I]);
               end
               Start = Base_ - length(CtoS_C[I]);
               for I2 = 1:length(CtoS_C[I])
                  V2[I] =V2[I]+V1[Start+I2]*CtoS_C[I][I2];
               end
            end
            V3 = Complex.(zeros(length(CtoS_C)));

            for I = 1:N
               Start = (I-1)*(I-1);
               for I2 = 1:I
                 if I2 == I
                   V3[Start+I2] = V2[Start+I2];
                 else
                   V3[Start+I2]   = (-V2[Start+I2]*im + V2[Start+2*I-I2])/(2^0.5);
                   V3[Start+2*I-I2] = ((-1)^(I-I2))*(V2[Start+I2]*im + V2[Start+2*I-I2])/(2^0.5);
                 end
               end
            end

            return V3;
         end

        function  CtoS_Encode( V1,n)

            V2 = zeros(n*n);
            Base_ = 0;
            PreviousLength = 0;
            Size_ = 0;
            Times = 0;

            for I = 1:length(V2)
               if length(CtoS_C[I]) != PreviousLength
                  Base_ = Base_+length(CtoS_C[I]);
                  PreviousLength = length(CtoS_C[I]);
                  Times = Times +1;
               end

               Start = Base_ - length(CtoS_C[I]);
               for I2 = 1:length(CtoS_C[I])
                  V2[I] =V2[I]+V1[Start+I2]*CtoS_C[I][I2];
               end
            end

            V3 = Complex.(zeros(length(V2)));

            for I = 1:n
               Start = (I-1)*(I-1);
               for I2 = 1:I
                 if I2 == I
                   V3[Start+I2] = V2[Start+I2];
                 else
                   V3[Start+I2] = (-V2[Start+I2]*im + V2[Start+2*I-I2])/(2^0.5);
                   V3[Start+2*I-I2] = ((-1)^(I-I2))*(V2[Start+I2]*im + V2[Start+2*I-I2])/(2^0.5);
                 end
               end
            end

            return V3;
         end

         function  XYZtoYlm_Encode( Cartesian_ ,n )

           n = n-1;
           Vector_Norm1 = (Cartesian_[1].^2 + Cartesian_[2].^2 + Cartesian_[3].^2)^0.5;
           Angle1 = acos(Cartesian_[3]/Vector_Norm1);
           Vector_Norm2 = (Cartesian_[1].^2 + Cartesian_[2].^2)^0.5;
           Angle2 = (Vector_Norm2 == 0) ? 0 : acos(Cartesian_[1]/Vector_Norm2)*(1-2*(Cartesian_[2]<0.0));
           return computeYlm(Angle1,Angle2,lmax = n)[:,1];

         end

           function  XYZtoYlm_Encode( Cartesian_ )

           n = N-1;
           Vector_Norm1 = (Cartesian_[1].^2 + Cartesian_[2].^2 + Cartesian_[3].^2)^0.5;
           Angle1 = acos(Cartesian_[3]/Vector_Norm1);
           Vector_Norm2 = (Cartesian_[1].^2 + Cartesian_[2].^2)^0.5;
           Angle2 = (Vector_Norm2 == 0) ? 0 : acos(Cartesian_[1]/Vector_Norm2)*(1-2*(Cartesian_[2]<0.0));
           computePlmcostheta!(S_cache, Angle1,n);
           return  computeYlm!(S_cache,Angle1,Angle2,n)[:,1];

         end



         function  CStoRS_Encode( V1 )

            V3 = zeros(length(V1));

            for I = 1:N
               Start = (I-1)*(I-1);
               for I2 = 1:I
                 if I2 == I
                   V3[Start+I2] = real(V1[Start+I2]);
                 else
                   V3[Start+I2]   = ((-1)^(I-I2))*(imag(V1[Start+2*I-I2]))*(2^0.5);
                   V3[Start+2*I-I2] = ((-1)^(I-I2))*(real(V1[Start+2*I-I2]))*(2^0.5);
                 end
               end
           end

           return V3;
         end


         function  CStoRS_Encode( V1,n )

            V3 = zeros(length(V1));

            for I = 1:n
               Start = (I-1)*(I-1);
               for I2 = 1:I
                 if I2 == I
                   V3[Start+I2] = real(V1[Start+I2]);
                 else
                   V3[Start+I2]   = ((-1)^(I-I2))*(imag(V1[Start+2*I-I2]))*(2^0.5);
                   V3[Start+2*I-I2] = ((-1)^(I-I2))*(real(V1[Start+2*I-I2]))*(2^0.5);
                 end
               end
           end

           return V3;
         end



         function  CStoRS_Encode_sg(V1,n)
            V3 = zeros(length(V1));
            for I2 = 1:n
                 if I2 == n
                   V3[I2] = real(V1[I2]);
                 else
                   V3[I2]   = ((-1)^(n-I2))*(imag(V1[2*n-I2]))*(2^0.5);
                   V3[2*n-I2] = ((-1)^(n-I2))*(real(V1[2*n-I2]))*(2^0.5);
                 end
            end
            return V3;
         end

         function  RStoCS_Encode_sg(V1,n)
            V3 = Complex.(zeros(length(V1)));
            for I2 = 1:n
                 if I2 == n
                   V3[I2] = V1[I2];
                 else
                   V3[I2]     = (-V1[I2]*im + V1[2*n-I2])/(2^0.5);
                   V3[2*n-I2] = ((-1)^(n-I2))*(V1[I2]*im + V1[2*n-I2])/(2^0.5);
                 end
            end
            return V3;
         end

         function  RStoCS_Encode( V1 )

           V3 = Complex.(zeros(length(V1)));

           for I = 1:N
               Start = (I-1)*(I-1);
               for I2 = 1:I
                 if I2 == I
                   V3[Start+I2] = V1[Start+I2];
                 else
                   V3[Start+I2]   = (-V1[Start+I2]*im + V1[Start+2*I-I2])/(2^0.5);
                   V3[Start+2*I-I2] = ((-1)^(I-I2))*(V1[Start+I2]*im + V1[Start+2*I-I2])/(2^0.5);
                 end
               end
           end

           return V3;
         end


         function  RStoCS_Encode( V1,n )

           V3 = Complex.(zeros(length(V1)));

           for I = 1:n
               Start = (I-1)*(I-1);
               for I2 = 1:I
                 if I2 == I
                   V3[Start+I2] = V1[Start+I2];
                 else
                   V3[Start+I2]   = (-V1[Start+I2]*im + V1[Start+2*I-I2])/(2^0.5);
                   V3[Start+2*I-I2] = ((-1)^(I-I2))*(V1[Start+I2]*im + V1[Start+2*I-I2])/(2^0.5);
                 end
               end
           end

           return V3;
         end




         function DecodeMatrix(V2,V3)

               LTm =  Complex.(zeros(size(WignerPI)[1],N*N));
               pseudoInput = zeros(size(WignerPI)[1]);

               for I =1:size(WignerPI)[1]
                  TempI = W3JTableI[WignerPI[I,1],WignerPI[I,2],WignerPI[I,3]];
                  TempC = W3JTableC[WignerPI[I,1],WignerPI[I,2],WignerPI[I,3]];
                  J1 = WignerPI[I,1]-1;
                  J2 = WignerPI[I,2]-1;
                  J3 = WignerPI[I,3]-1;
                     for I_1 = 1:length(TempI)
                        m1 = TempI[I_1].m1;
                        m2 = TempI[I_1].m2;
                        m3 = -(m1+m2);
                        I1 = (J1^2 + m1+J1+1);
                        I2 = (J2^2 + m2+J2+1);
                        I3 = (J3^2 + m3+J3+1);
                        LTm[I,I1] = LTm[I,I1]+V2[I2]*V3[I3]*TempC[I_1];
                     end
               end

               #for I =1:size(WignerPI)[1]
               #   pseudoInput[I]=(LTm[I,:]'LTm[I,:]).re;
               #end

               G_inv=LinearAlgebra.pinv(LTm);
               return G_inv;
         end


          function ProductEncode(V1,V2)


            d1 = Int32(floor((length(V1))^0.5));
            d2 = Int32(floor((length(V2))^0.5));
            SIZE = min((d1+d2-1),N);
            Base = (abs(d1-d2)*abs(d1-d2));
            LTm  = Complex.(zeros(1,SIZE*SIZE));

            LENGTH = 0;
            for J3 = (abs(d1-d2)+1):min((d1+d2-1),N)  
              for J1 = max(J3-d2+1,2-J3,1):min(J3+d2-1,d1) 
                for J2 = max(J3-J1+1,J1-J3+1,1):min(J3+J1-1,d2)
                  LENGTH = LENGTH+1; 
                end
              end 
            end  

            pseudoInput = zeros(LENGTH);
            I = 1;
            for J3 = (abs(d1-d2)+1):min((d1+d2-1),N)  
              for J1 = max(J3-d2+1,2-J3,1):min(J3+d2-1,d1) 
                for J2 = max(J3-J1+1,J1-J3+1,1):min(J3+J1-1,d2)
                        TempI = W3JTableI[J1,J2,J3];
                        TempC = W3JTableC[J1,J2,J3];
                        for I_1 = 1:length(TempI)
                           m1 = TempI[I_1].m1;
                           m2 = TempI[I_1].m2;
                           m3 = -(m1+m2);
                           I1 = ((J1-1)^2 + m1+J1);
                           I2 = ((J2-1)^2 + m2+J2);
                           I3 = ((J3-1)^2 + m3+J3);
                           LTm[1,I3-Base] = LTm[1,I3-Base]+V1[I1]*V2[I2]*TempC[I_1]*((2*J3-1)^0.5*((-1)^(J1-J2-m3)));#wigner3J coefficient to Clebsh-Gordan coefficient;
                        end  
                        pseudoInput[I]=(LTm[1,:]'LTm[1,:]).re;
                        I = I+1;                        
                        LTm.=0.0;
                  end 
              end 
            end 

            return pseudoInput;
          end





          function SelfProductPairwise(V,n,n2)

             d2 = Int32(floor((length(V))^0.5));
             d1 = n;
             Max_size = 0;
             size_sum = 0;
             size_add = 0; 
             odd_sum  = 0;


             for J1 = 1:d1
                size__   = 0;
                size_add = 0;
                for J2 = max(J1-d2+1,1):min(d2-1+J1,d2)
                  Begin = max((J1-J2+1),J2,(J2-n2+J1));
                  End   = min(J1+J2-1,d2);
                  for I in Begin:End
                   #   print( [J1,J2,I]);
                   #   print("\n");
                      if ((I == J2)&&((J1%2)==0))
                          continue;
                      end 
                      size__ = size__ + 1;
                  end    
               end

               if ( Max_size <  size__ ) 
                    Max_size = size__; 
               end 
               if J1 > 1
                 size_sum = size__*(size__+1)/2  + size_sum; 
               else 
                 size_sum = size__  + size_sum;
               end
            end 

            LTm =  Complex.(zeros(Max_size,d2*d2-(d2-1)*(d2-1)));
            DiffJ = zeros(d2*d2-(d2-1)*(d2-1),1);
            Result_ = zeros(Int16(size_sum),1);
            size_sum = 1; 


            for J1 = 1:d1
                Index = 1;
                for J2 = max(J1-d2+1,1):min(d2-1+J1,d2)
                   Begin = max((J1-J2+1),J2,(J2-n2+J1));
                   End   = min(J1+J2-1,d2);

                   for I in Begin:End

                       j1  = J1-1;
                       j2  = J2-1;
                       j3  = I -1;
                       J   = j1 -(I-J2)+1;

                      if (((I == J2)&&((J%2)==0)))
                         continue;
                      end 

                      #print( [J2,I,J]);  
                      #print("\n");

                      TempI = CGTableI[J2,I,J];
                      TempC = CGTableC[J2,I,J];
                      LTm[Index,:] .= 0.0im;
                      for I_1 = 1:length(TempI)
                          for I_2 = 1:length(TempI[I_1])
                            m1 = TempI[I_1][I_2].m1;
                            m2 = TempI[I_1][I_2].m2;
                            m3 = (m1+m2);
                            I1 = (j2^2 + m1+j2+1);
                            I2 = (j3^2 + m2+j3+1);
                            I3 = (m3+j1+1);
                            LTm[Index,I3] = LTm[Index,I3] + V[I1]*V[I2]*TempC[I_1][I_2];
                         end
                      end
                          DiffJ[Index] = J;
                          Index = Index+1;
                   end
                end

                 
                if Index >1 
                      Len = J1*J1-(J1-1)*(J1-1);

                      for I_3 in 1:(Index-1)

                         if Len > 1 
                           for I_2 in I_3:(Index-1) #for I_2 in I_3:min((Index-1),J1+I_3-1)
                           
                              if ( abs(DiffJ[I_3]-DiffJ[I_2]) != 1 )
                                Self_Product = LTm[I_3,1:Len]'*LTm[I_2,1:Len]; 
                                Result_[size_sum] = real(Self_Product);
                                size_sum = size_sum+1;
                               end 
                                #print(abs(DiffJ[I_3]-DiffJ[I_2]));
                                #print("  "); 
                                #print(Self_Product);
                                #print("  ");  
                                #print("\n");    
                           end 
                         else 
                           Self_Product = LTm[I_3,1:Len]'*LTm[I_3,1:Len];  
                           Result_[size_sum] = real(Self_Product);
                           size_sum = size_sum+1;
                         end 
                      end 
                end 


             end
             return Result_[1:size_sum-1];
          end



          function SelfProductMatrix(V,n,n2)

             d2 = Int32(floor((length(V))^0.5));
             d1 = n;
             size_sum = 0;
             size_add = 0; 
             odd_sum  = 0;


             J1 = d1
             size__   = 0;
             size_add = 0;
             for J2 = max(J1-d2+1,1):min(d2-1+J1,d2)
                Begin = max((J1-J2+1),J2,(J2-n2+J1));
                End   = min(J1+J2-1,d2);
                  for I in Begin:End
                   #   print( [J1,J2,I]);
                   #   print("\n");
                      if ((I == J2)&&((J1%2)==0))
                          continue;
                      end 
                      size__ = size__ + 1;
                end    
            end

            
            LTm =  Complex.(zeros(size__,d1*d1-(d1-1)*(d1-1)));

                J1 = d1;
                Index = 1;
                for J2 = max(J1-d2+1,1):min(d2-1+J1,d2)
                   Begin = max((J1-J2+1),J2,(J2-n2+J1));
                   End   = min(J1+J2-1,d2);

                   for I in Begin:End

                       j1  = J1-1;
                       j2  = J2-1;
                       j3  = I -1;
                       J   = j1+1;

                      if (((I == J2)&&((J%2)==0)))
                         continue;
                      end 

                      #print( [J2,I,J]);  
                      #print("\n");

                      TempI = W3JTableI[J2,I,J];
                      TempC = W3JTableC[J2,I,J];
                      LTm[Index,:] .= 0.0im;
                      for I_1 = 1:length(TempI)
                            m1 = TempI[I_1].m1;
                            m2 = TempI[I_1].m2;
                            m3 = (m1+m2);
                            I1 = (j2^2 + m1+j2+1);
                            I2 = (j3^2 + m2+j3+1);
                            I3 = (m3+j1+1);
                            LTm[Index,I3] = LTm[Index,I3] + V[I1]*V[I2]*TempC[I_1];
                      end
                          Index = Index+1;
                   end
                end
             
             return LTm;
          end










         # length of V >=4
         # According to the recursion of Associated Legendre Functions
         # P_n+1 = ((2n+1)/(n+1)) * (sqrt((2n+3)/(2n+1))Cos*P_n - (n/n+1)*sqrt((2n+3)/(2n-1))P_n-1
         # P_n+1^(m+1) = 2CosP_n^(m+1)(sqrt((2*n+3)/(2*n+1))sqrt((n-m)/(n+m+2))  ) - P_n-1^(m+1)sqrt((2*n+3)/(2*n-1))sqrt(((n-m)(n-m-1))/((n+m+2)(n+m+1)))+ (2m+1)P_n^m*Sin*(-1)^m*sqrt((2*n+3)/(2*n+1))*sqrt(1/((n+m+2)(n+m+1)))
         function IncreaseDegree(V,n)
	
	     N  = length(V);
             N  = Int16(N^0.5);
             N_new = (N+n);
             J0 = 1;

             V3 = Complex.(zeros(N_new*N_new));
		     V3[1:N*N] = V[1:N*N];

             Cos    = V3[J0^2+J0+1]*2.046653415892977;#2*sqrt(pi/3);
             Sina_h = V3[J0^2+J0+1+1]*2.8944050182330705; #2*(sqrt(2*pi/3));
		     Sina_l = V3[J0^2+J0+1-1]*2.8944050182330705; #2*(sqrt(2*pi/3));

            for I__ in 1:n
	          J3 = (N-1+I__);
              J1 = J3-1;
              J2 = J3-2;
		      m1 = 0;
              m2 = 0;
              m3 = 0;
		      I1 = (J1^2 + m1+J1+1);
              I2 = (J2^2 + m2+J2+1);
              I3 = (J3^2 + m3+J3+1);
              V3[I3] = ((2*J1+1)/(J1+1))*sqrt((2*J1+3)/(2*J1+1))*Cos*V3[I1]-((J1)/(J1+1))*sqrt((2*J1+3)/(2*J1-1))*V3[I2];

            for m in 0:J1
		          m3 = m+1;
		          m1 = m+1;
                  m2 = m+1;

		          I1 = (J1^2 + m1+J1+1);
                  I2 = (J2^2 + m2+J2+1);
                  I3 = (J3^2 + m3+J3+1);
                  Im = (J1^2 + m +J1+1);

             if m2 <= J2
 		            V3[I3] = 2*Cos*sqrt(((2*J1+3)/(2*J1+1))*((J1-m)/(J1+m+2)))*V3[I1] - sqrt(((2*J1+3)/(2*J1-1))*(((J1-m)*(J1-m-1))/((J1+m+2)*(J1+m+1))))*V3[I2];
 		     elseif m1 <= J1
 		            V3[I3] = 2*Cos*sqrt(((2*J1+3)/(2*J1+1))*((J1-m)/(J1+m+2)))*V3[I1];
             end
                V3[I3] = V3[I3] + sqrt(((2*J1+3)/(2*J1+1))*(1/((J1+m+2)*(J1+m+1))))*(2*m+1)*V3[Im]*Sina_h;


   		     I1 = (J1^2 - m1+J1+1);
             I2 = (J2^2 - m2+J2+1);
             I3 = (J3^2 - m3+J3+1);
             Im = (J1^2 - m +J1+1);

             if m2 <= J2
 		            V3[I3] = 2*Cos*sqrt(((2*J1+3)/(2*J1+1))*((J1-m)/(J1+m+2)))*V3[I1] - sqrt(((2*J1+3)/(2*J1-1))*(((J1-m)*(J1-m-1))/((J1+m+2)*(J1+m+1))))*V3[I2];
 		     elseif m1 <= J1
 		            V3[I3] = 2*Cos*sqrt(((2*J1+3)/(2*J1+1))*((J1-m)/(J1+m+2)))*V3[I1];
             end
                V3[I3] = V3[I3] + sqrt(((2*J1+3)/(2*J1+1))*(1/((J1+m+2)*(J1+m+1))))*(2*m+1)*V3[Im]*Sina_l;
             end
           end

              return  V3;
           end

        # Length of V >=9
        # Extract reference vector from the Quadrupole moment
        function ReferencesExtract(V)

          VR  =  equivalentFeatures.CStoRS_Encode(V,3)
          a12 = VR[5]/1.0925484305920792;
          a23 = VR[6]/1.0925484305920792;
          a13 = VR[8]/1.0925484305920792;
          r3    = 4.1887902047863905*(VR[2]^2 + VR[3]^2 + VR[4]^2);
          r1    =  VR[7];
          r2    =  VR[9];
          a33  = (-0.31539156525252005*r3 -r1)/(-0.31539156525252005-0.63078313050504);
          a11  = (0.5462742152960396*(r1-0.63078313050504*a33)-0.31539156525252005*r2)/(2*0.5462742152960396*(-0.31539156525252005));
          a22  = (0.5462742152960396*(r1-0.63078313050504*a33)+0.31539156525252005*r2)/(2*0.5462742152960396*(-0.31539156525252005));
          alpha  =  a11 + a22+a33;
          beta   =  a12*a12+a13*a13+a23*a23 -a11*a22-a22*a33-a33*a11;
gamma    =  a11*a22*a33 +2*a12*a23*a13-a11*a23*a23-a12*a12*a33-a13*a13*a22;
          p = abs(-((3*beta+alpha*alpha)/3));
          q = -(gamma+2*alpha*alpha*alpha/27+alpha*beta/3);
          x = -q/( 2*((p/3)^3)^0.5 )
          a = -3/4;
          b = -x/4;
          cos_div_3 = 2*real((-(b/2) +  (Complex( -(a/3)^3 -(b/2)^2)^0.5)*im)^(1/3))
          lambda1 =  alpha/3  +2*(p/3)^0.5*cos_div_3;
          a11 = a11 - lambda1;
          a22 = a22 - lambda1;
          a33 = a33 - lambda1;
          T1 = [a22*a33-a23*a23,-(a12*a33-a23*a13),a12*a23-a22*a13 ];
          T1 = T1/(T1'*T1)^0.5;
          T2 = [VR[4], VR[2],VR[3]];
          T2 = T2/(T2'*T2)^0.5;
	    hcat(T1,T2);
        end

        # Derivative of Spherical Harmonic with Respect to $\theta$ and $\phi$

        function  DerivativeSH(Y)

	          SinTheta  = ((real(Y[2])^2 + imag(Y[2])^2))^0.5;
	          ExpVarphi = SinTheta != 0.0 ? Y[2]/SinTheta : 0;
	          CotTheta  = SinTheta != 0.0 ? (Y[3])/(SinTheta*(2^0.5)) : 0; 
	          DerivationTheta   = Complex.(zeros(length(Y))); 
	          DerivationVarphi  = Complex.(zeros(length(Y)));
	
	          N = Int64(sqrt(length(Y)));
	          Base_ = 0; 
	          for I in 1:N  
	            for J in 1:((I-1)*2)
	                m = J - I;
	                DerivationTheta[Base_+J] = m*CotTheta*Y[Base_+J] + ((I-m-1)*(I+m))^0.5*ExpVarphi*Y[Base_+J+1];
	                DerivationVarphi[Base_+J]= m*im*Y[Base_+J];
	            end
	            J = ((I-1)*2)+1; 
	            m = J - I;
	            DerivationTheta[Base_+J]=m*CotTheta*Y[Base_+J];
	            DerivationVarphi[Base_+J]= m*im*Y[Base_+J];
	            Base_ = I*I; 
	          end 
	          
	          return  hcat(DerivationTheta,DerivationVarphi); 
       end

        # Derivative of Spherical Harmonic with Respect to Cartesian coordinates, DR is (df(r)/dr)/f(r)

         function  DerivativeSH_XYZ(Y,DR)

		          ReciprocalRadii   = (Y[3]*Y[3]-Y[2]*Y[4]*2)/(0.25*(3/pi));
		          ReciprocalRadii   = ReciprocalRadii != 0.0 ? (ReciprocalRadii)^(-(1/2)) : 0;
		          SinTheta  = ((real(Y[2])^2 + imag(Y[2])^2))^0.5;
		          ExpVarphi = SinTheta != 0.0 ? Y[2]/SinTheta : 0;
		          CotTheta  = SinTheta != 0.0 ? (Y[3])/(SinTheta*(2^0.5)) : 0; 
		
		          SinVarphi         = SinTheta != 0.0 ? -imag(Y[2])/SinTheta : 0;
		          CosVarphi         = SinTheta != 0.0 ? real(Y[2])/SinTheta : 0;
		          CosTheta          = ReciprocalRadii*(Y[3]/((1/2)*((3/pi)^(1/2)))); 
		          DerivationX       = Complex.(zeros(length(Y))); 
		          DerivationY       = Complex.(zeros(length(Y)));
		          DerivationZ       = Complex.(zeros(length(Y)));
		          DR                = Complex.(DR);
		          SinTheta  = ReciprocalRadii*(SinTheta/(0.5*((3/(pi*2))^0.5)));
		          SinVarphiSlashSinTheta  = SinTheta  != 0.0 ? SinVarphi/SinTheta : 0; 
		          CosVarphiSlashSinTheta  = SinTheta != 0.0 ?  CosVarphi/SinTheta : 0; 
		
		          N = Int64(sqrt(length(Y)));
		          Base_ = 0; 
		          for I in 1:N  
		            for J in 1:((I-1)*2)
		                m = J - I;
		                DerivationTheta  = m*CotTheta*Y[Base_+J] + ((I-m-1)*(I+m))^0.5*ExpVarphi*Y[Base_+J+1]
		                DerivationVarphi = m*im*Y[Base_+J];
		                DerivationX[Base_+J] = DR[Base_+J]*Y[Base_+J]*SinTheta*CosVarphi + ReciprocalRadii*(CosTheta*CosVarphi*DerivationTheta - SinVarphiSlashSinTheta*DerivationVarphi);
		                DerivationY[Base_+J] = DR[Base_+J]*Y[Base_+J]*SinTheta*SinVarphi + ReciprocalRadii*(CosTheta*SinVarphi*DerivationTheta + CosVarphiSlashSinTheta*DerivationVarphi);
		                DerivationZ[Base_+J] = DR[Base_+J]*Y[Base_+J]*CosTheta + ReciprocalRadii*(-SinTheta*DerivationTheta);
		            end
		            J = ((I-1)*2)+1; 
		            m = J - I;
		            DerivationTheta  = m*CotTheta*Y[Base_+J];
		            DerivationVarphi = m*im*Y[Base_+J];
		            DerivationX[Base_+J] = DR[Base_+J]*Y[Base_+J]*SinTheta*CosVarphi + ReciprocalRadii*(CosTheta*CosVarphi*DerivationTheta - SinVarphiSlashSinTheta*DerivationVarphi);
		            DerivationY[Base_+J] = DR[Base_+J]*Y[Base_+J]*SinTheta*SinVarphi + ReciprocalRadii*(CosTheta*SinVarphi*DerivationTheta + CosVarphiSlashSinTheta*DerivationVarphi);
		            DerivationZ[Base_+J] = DR[Base_+J]*Y[Base_+J]*CosTheta + ReciprocalRadii*(-SinTheta*DerivationTheta);
		
		            Base_ = I*I; 
		          end 
		          
		          return  hcat(DerivationX,DerivationY,DerivationZ); 
          end





         function DecodeMatrix(V2,V3,n)

            VSize = 0;
            for I in 1:n
               VSize = VSize +(min((n-1),(n+I-2))-abs(n-I))+1;
            end

            LTm =  Complex.(zeros(VSize,2*(n-1)+1));
            # pseudoInput = zeros(VSize);
               Index = 1;
               for I =1:size(WignerPI)[1]
                  if WignerPI[I,1] != n
                     continue;
                  end
                  if WignerPI[I,2] > n || WignerPI[I,3] > n
                     continue;
                  end
                  TempI = W3JTableI[WignerPI[I,1],WignerPI[I,2],WignerPI[I,3]];
                  TempC = W3JTableC[WignerPI[I,1],WignerPI[I,2],WignerPI[I,3]];
                  J1 = WignerPI[I,1]-1;
                  J2 = WignerPI[I,2]-1;
                  J3 = WignerPI[I,3]-1;
                     for I_1 = 1:length(TempI)
                        m1 = TempI[I_1].m1;
                        m2 = TempI[I_1].m2;
                        m3 = -(m1+m2);
                        I1 = (m1+J1+1);
                        I2 = (J2^2 + m2+J2+1);
                        I3 = (J3^2 + m3+J3+1);
                        LTm[Index,I1] = LTm[Index,I1]+V2[I2]*V3[I3]*TempC[I_1];
                     end
                  Index = Index+1;
               end

              #for I =1:VSize
              #    pseudoInput[I]=(LTm[I,:]'LTm[I,:]).re;
              #end

              G_inv=LinearAlgebra.pinv(LTm);
              return G_inv;
         end


function DecodeMatrixCompact(V2,V3,n)

            VSize = n^2;


            LTm =  Complex.(zeros(VSize,VSize));
            # pseudoInput = zeros(VSize);
               Index = 1;
               for I =1:size(WignerPI)[1]
                  if  WignerPI[I,1] > n   || WignerPI[I,2] > WignerPI[I,1] || WignerPI[I,3] > WignerPI[I,1]
                     continue;
                  end
                  if (WignerPI[I,2] + WignerPI[I,3]) >= (WignerPI[I,1] +3)
                     continue;
                  end
                  TempI = W3JTableI[WignerPI[I,1],WignerPI[I,2],WignerPI[I,3]];
                  TempC = W3JTableC[WignerPI[I,1],WignerPI[I,2],WignerPI[I,3]];
                  J1 = WignerPI[I,1]-1;
                  J2 = WignerPI[I,2]-1;
                  J3 = WignerPI[I,3]-1;
                     for I_1 = 1:length(TempI)
                        m1 = TempI[I_1].m1;
                        m2 = TempI[I_1].m2;
                        m3 = -(m1+m2);
                        I1 = (J1^2 + m1+J1+1);
                        I2 = (J2^2 + m2+J2+1);
                        I3 = (J3^2 + m3+J3+1);
                        LTm[Index,I1] = LTm[Index,I1]+V2[I2]*V3[I3]*TempC[I_1];
                     end
                  Index = Index+1;
               end

              #for I =1:VSize
              #    pseudoInput[I]=(LTm[I,:]'LTm[I,:]).re;
              #end

              G_inv=LinearAlgebra.pinv(LTm);
              return G_inv;
         end





          function ProductEncode(V1,V2,n)


            d1 = Int32(floor((length(V1))^0.5));
            d2 = Int32(floor((length(V2))^0.5));
            SIZE = min((d1+d2-1),N);
            BASE = ((abs(d1-d2))*(abs(d1-d2)));
            LTm  = Complex.(zeros(1,SIZE*SIZE));

            LENGTH = 0;
            for J3 = (abs(d1-d2)+1):min((d1+d2-1),N)  
              for J1 = max(J3-d2+1,1):min(J3+d2-1,d1) 
                for J2 = max(J3-J1+1,J1-J3+1,1):min(J3+J1-1,d2)
                 
                  if  J3 > abs(J1-J2) + n
                      J2 = max(J2,J1 + N - n -1);
                      if J2 > min(J3+J1-1,d2)
                         break; 
                      end 
                      continue; 
                  end 
                  LENGTH = LENGTH+1;
                  print( [J1,J2,J3]); 
                  print("\n");
                end
              end 
            end  

            pseudoInput = zeros(LENGTH);
            I = 1;
            for J3 = (abs(d1-d2)+1):min((d1+d2-1),N)  
              BASE  = (J3-1)*(J3-1);
              for J1 = max(J3-d2+1,1):min(J3+d2-1,d1) 
                for J2 = max(J3-J1+1,J1-J3+1,1):min(J3+J1-1,d2)
                 
                        if  J3 > abs(J1-J2) + n
                            J2 = max(J2,J1 + N - n-1);
                            if J2 > min(J3+J1-1,d2)
                               break; 
                            end 
                            continue; 
                        end 

                        TempI = W3JTableI[J1,J2,J3];
                        TempC = W3JTableC[J1,J2,J3];
                        for I_1 = 1:length(TempI)
                           m1 = TempI[I_1].m1;
                           m2 = TempI[I_1].m2;
                           m3 = -(m1+m2);
                           I1 = ((J1-1)^2 + m1+J1);
                           I2 = ((J2-1)^2 + m2+J2);
                           I3 = ((J3-1)^2 + m3+J3);
                           LTm[1,I3-BASE] = LTm[1,I3-BASE]+V1[I1]*V2[I2]*TempC[I_1]*((2*J3-1)^0.5*((-1)^(J1-J2-m3)));#wigner3J coefficient to Clebsh-Gordan coefficient;
                        end  
                        pseudoInput[I]=(LTm[1,:]'LTm[1,:]).re;
                        I = I+1;                        
                        LTm.=0.0;
                  end 
              end 
            end 

            return pseudoInput;
          end




          function ProductEncodePairwise(V1,V2,n)


            d1 = Int32(floor((length(V1))^0.5));
            d2 = Int32(floor((length(V2))^0.5));
            SIZE = min((d1+d2-1),N);
            size__   = 0;
            Max_size = 0; 
            size_sum = 1; 
            LENGTH =   0;



            for J3 = (abs(d1-d2)+1):min((d1+d2-1),N)  
              size__   = 0;
              for J1 = max(J3-d2+1,1):min(J3+d2-1,d1) 
                for J2 = max(J3-J1+1,J1-J3+1,1):min(J3+J1-1,d2)

                  if  J3 > abs(J1-J2) + n
                      J2 = max(J2,J1 + N - n -1);
                      if J2 > min(J3+J1-1,d2)
                        break; 
                      end  
                      continue; 
                  end 

                  size__ = size__ + 1;
                  print( [J1,J2,J3]); 
                  print("\n");
                end
              end 
              if ( Max_size <  size__ ) 
                  Max_size = size__; 
              end 
              if J3 > 1
                 LENGTH = Int32(size__*(size__+1)/2)  + LENGTH; 
              else 
                 LENGTH = size__  + LENGTH;
              end
              if (J3 == 2)&&(size__>2)  
                 LENGTH = LENGTH - 2; 
              end 
            end  

            LTm =  Complex.(zeros(Max_size,SIZE*SIZE));
            DiffJ = zeros(1,Max_size);
            pseudoInput = zeros(LENGTH);
            

            I = 1;
            for J3 = (abs(d1-d2)+1):min((d1+d2-1),N)
              Index = 1;  
              BASE  = (J3-1)*(J3-1);
              for J1 = max(J3-d2+1,1):min(J3+d2-1,d1) 
                for J2 = max(J3-J1+1,J1-J3+1,1):min(J3+J1-1,d2)

         
                  if  J3 > abs(J1-J2) + n
                      J2 = max(J2,J1 + N - n -1);
                      if J2 > min(J3+J1-1,d2)
                        break; 
                      end  
                      continue; 
                  end 


                        TempI = W3JTableI[J1,J2,J3];
                        TempC = W3JTableC[J1,J2,J3];
                        for I_1 = 1:length(TempI)
                           m1 = TempI[I_1].m1;
                           m2 = TempI[I_1].m2;
                           m3 = -(m1+m2);
                           I1 = ((J1-1)^2 + m1+J1);
                           I2 = ((J2-1)^2 + m2+J2);
                           I3 = ((J3-1)^2 + m3+J3);
                           LTm[Index,I3-BASE] = LTm[Index,I3-BASE]+V1[I1]*V2[I2]*TempC[I_1]*((2*J3-1)^0.5*((-1)^(J1-J2-m3)));#wigner3J coefficient to Clebsh-Gordan coefficient;   
                        end  
                        DiffJ[Index]= abs(J1-J2); 
                        Index = Index+1;
                        I = I+1;                        
                  end 
              end 
                
                if Index >1 
                      Len = J3*J3-(J3-1)*(J3-1);

                      for I_3 in 1:(Index-1)

                         if Len > 1 
                          for I_2 in I_3:(Index-1) 


                           if !(((I_3 == 2)||(I_3 == 1))&&(I_2 == 3)&&(Len ==3))
                               Self_Product = LTm[I_3,1:Len]'*LTm[I_2,1:Len]; 
                               if ( abs(DiffJ[I_3]-DiffJ[I_2]) != 1 ) 
                                 pseudoInput[size_sum] = real(Self_Product);
                               else 
                                 pseudoInput[size_sum] = imag(Self_Product);
                               end 
                               size_sum = size_sum+1; 
                            end 
                               print(Self_Product); 
                               print("  ");
                               print(I_2);
                               print("  ");
                               print(I_3);
                               print("\n");
                                   
                          end 
                         else 
                               Self_Product = LTm[I_3,1:Len]'*LTm[I_3,1:Len];  
                               pseudoInput[size_sum] = real(Self_Product);
                               size_sum = size_sum+1;
                               print("\n");
                         end 
                      end 
                end 
                LTm.=0.0;
            end 
            return pseudoInput;
          end




   end






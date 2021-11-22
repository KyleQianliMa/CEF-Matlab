function[calcscattering,CalcEnergy] = sample(C20,C40,C44,C60,C64,h_mf,v_mf)
%KyleNd-LSCO
%calling the function


%total angular momentum found by Hund's rule for Er3+
J=9/2;

%iteration = length(C20).*length(C40).*length(C44).*length(C60).*length(C64);

%set down the temperature (in K)
%T = 5;

%boltzman constant in meV*K-1
%k=8.6173324*10^(-2);

[O20,O40,O44,O60,O64,Jx,Jy,Jz,Jplus,Jminus] = OperatorCuprate(J);


bestchi = 1e99; %make sure that the first one will be lower than that.


counter=0;
step=1;
k20=-284/0.912;
k40=-344/(1.25e-2);
k60=-88/(2.09e-4);
k44=93/(-2.82e-2);
k64=104/(-2.77e-3);

                            D20=C20/k20;
                            D40=C40/k40;
                            D60=C60/k60;
                            D44=C44/k44;
                            D64=C64/k64;
                            H=D20*O20+D40*O40+D60*O60+D44*O44+D64*O64+h_mf*(Jx+Jy)/(sqrt(2))+v_mf*Jz;

                            [eigenvector,SolveEnergy] = eig(H,'vector');
                            Energytemp = sort(SolveEnergy);
                            Energy = Energytemp + abs(min(SolveEnergy(:,1)));
                            CalcEnergy = [Energy(3,1),Energy(5,1),Energy(7,1),Energy(9,1)];
                            
                            scattering1=scattering_CEF(eigenvector(:,1),eigenvector(:,3),Jx,Jy,Jz)...
                                        +scattering_CEF(eigenvector(:,1),eigenvector(:,4),Jx,Jy,Jz);
                            
                            scattering2=scattering_CEF(eigenvector(:,1),eigenvector(:,5),Jx,Jy,Jz)...
                                        +scattering_CEF(eigenvector(:,1),eigenvector(:,6),Jx,Jy,Jz);
                            
                            scattering3=scattering_CEF(eigenvector(:,1),eigenvector(:,7),Jx,Jy,Jz)...
                                        +scattering_CEF(eigenvector(:,1),eigenvector(:,8),Jx,Jy,Jz); 
        
                            scattering4=scattering_CEF(eigenvector(:,1),eigenvector(:,9),Jx,Jy,Jz)...
                                       +scattering_CEF(eigenvector(:,1),eigenvector(:,10),Jx,Jy,Jz);
                            
                            %scattering6=scattering_CEF(eigenvector(:,1),eigenvector(:,9),Jx,Jy,Jz)...
                            %           +scattering_CEF(eigenvector(:,1),eigenvector(:,10),Jx,Jy,Jz);
                            
                            %scattering6=scattering_CEF(eigenvector(:,1),eigenvector(:,9),Jx,Jy,Jz)...
                            %            +scattering_CEF(eigenvector(:,1),eigenvector(:,10),Jx,Jy,Jz);%24meV
                            
                            %scattering7=scattering_CEF(eigenvector(:,1),eigenvector(:,9),Jx,Jy,Jz)...
                            %           +scattering_CEF(eigenvector(:,1),eigenvector(:,10),Jx,Jy,Jz);%26.76meV
                            
                            %scattering7=scattering_CEF(eigenvector(:,1),eigenvector(:,9),Jx,Jy,Jz)...
                            %           +scattering_CEF(eigenvector(:,1),eigenvector(:,10),Jx,Jy,Jz);  %93meV
                            
                             N=scattering2;        
                             s1=100*scattering1/N;
                             s2=100*scattering2/N;
                             s3=100*scattering3/N;
                             s4=100*scattering4/N;
                            
                             
                             calcscattering=[s1,s2,s3,s4];
                             
                           

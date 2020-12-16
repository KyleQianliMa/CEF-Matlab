

%Boothroyad
%B20=x;%-269;
%B40=-144;%;-360;
%B60=81;
%B44=298;
%B64=94;
%h=0;%0.43;
%v=0;

%Infrared paper
B20=212;%188;
B40=211;
B60=-206;
B44=-111;
B64=79;
h=0;
v=0.0;

k20=-281/0.912;
k40=-344/(1.25e-2);
k60=-88/(2.09e-4);
k44=93/(-2.82e-2);
k64=104/(-2.77e-3);
kh=0.43/0.31;
v=0;

%Nd2CuO4 10meV
%B20=0.873;
%B40=1.31e-2;
%B60=1.59e-4;
%B44=-2.59e-2;
%B64=-3.30e-3;
%h=0.31;

[O20,O40,O44,O60,O64,Jx,Jy,Jz,Jplus,Jminus] = OperatorCuprate(9/2);
H=B20*O20/k20+B40*O40/k40+B60*O60/k60+B44*O44/k44+B64*O64/k64+h*(Jx+Jy)/(sqrt(2))+v*Jz;
[eigenvector,SolveEnergy] = eig(H,'vector');
[Energy,index] = sort(SolveEnergy);
k=8.6173324*10^(-2);
Ei=0;
Energy = Energy + abs(min(SolveEnergy(:,1)));
                            scattering1=scattering_CEF(eigenvector(:,1),eigenvector(:,2),Jx,Jy,Jz);
        
                            scattering2=scattering_CEF(eigenvector(:,1),eigenvector(:,3),Jx,Jy,Jz);
                            
                            scattering3=scattering_CEF(eigenvector(:,1),eigenvector(:,4),Jx,Jy,Jz);        

                            scattering4=scattering_CEF(eigenvector(:,1),eigenvector(:,5),Jx,Jy,Jz);
                            
                            scattering5=scattering_CEF(eigenvector(:,1),eigenvector(:,6),Jx,Jy,Jz);
                   
                            scattering6=scattering_CEF(eigenvector(:,1),eigenvector(:,7),Jx,Jy,Jz);
                            
                            scattering7=scattering_CEF(eigenvector(:,1),eigenvector(:,8),Jx,Jy,Jz);
                    
                            scattering8=scattering_CEF(eigenvector(:,1),eigenvector(:,9),Jx,Jy,Jz);
                            
                            scattering9=scattering_CEF(eigenvector(:,1),eigenvector(:,10),Jx,Jy,Jz);
                            
                            
                             N=scattering4+scattering5;           
                             s1=(scattering2+scattering3)/N;
                             s2=(scattering4+scattering5)/N;
                             s3=(scattering6+scattering7)/N;
                             s4=(scattering8+scattering9)/N;
                             
                             %g-tensor
                             S=3/2;L=6;J=9/2;
                             gj=3/2+(S*(S+1)-L*(L+1))/(2*J*(J+1));
                             i=eigenvector(:,1);
                             j=eigenvector(:,2);
                             gll=2*gj*abs(transpose(j)*Jz*j);
                             gperp=gj*(transpose(i)*Jplus*j);
                             
                             
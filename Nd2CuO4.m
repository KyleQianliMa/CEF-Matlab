%B20=-28;%188;
%B40=-263;
%B60=34;
%B44=199;
%B64=183;
%h=0.18;
%v=0.0;

%k20=-284/0.912;
%k40=-344/(1.25e-2);
%k60=-88/(2.09e-4);
%k44=93/(-2.82e-2);
%k64=104/(-2.77e-3);
%kh=0.43/0.31;
%v=0;

%Nd2CuO4 10meV
B20=0.873;
B40=1.31e-2;
B60=1.59e-4;
B44=-2.59e-2;
B64=-3.30e-3;
h=0.31;
v=0;

[O20,O40,O44,O60,O64,Jx,Jy,Jz,Jplus,Jminus] = OperatorCuprate(9/2);
H=B20*O20+B40*O40+B60*O60+B44*O44+B64*O64+h*(Jx+Jy)/(sqrt(2))+v*Jz;
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
                             s1=scattering1/N;
                             s2=scattering2/N;
                             s3=scattering3/N;
                             s4=scattering4/N;
                             %s5=scattering5/N;
                             %s6=scattering6/N;
                             %s7=scattering7/(scattering4+scattering3);
                             
                             
                             calcscattering=[s1,s2,s3,s4];
                             
                             
                             %g-tensor
                             S=3/2;L=6;J=9/2;
                             gj=3/2+(S*(S+1)-L*(L+1))/(2*J*(J+1));
                             i=eigenvector(:,1);
                             j=eigenvector(:,2);
                             gll=2*gj*abs(transpose(i)*Jz*i);
                             gperp=gj*abs(transpose(i)*Jplus*j);
                             
                            
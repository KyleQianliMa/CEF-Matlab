function [testEnergy,testscattering]=Testsolutions(solution)
[O20,O40,O44,O60,O64,Jx,Jy,Jz,Jplus,Jminus] = OperatorCuprate(9/2);
k20=-281/0.912;
k40=-344/(1.25e-2);
k60=-88/(2.09e-4);
k44=93/(-2.82e-2);
k64=104/(-2.77e-3);
D20=solution(1,1)/k20;
D40=solution(1,2)/k40;
D60=solution(1,3)/k60;
D44=solution(1,4)/k44;
D64=solution(1,5)/k64;
h=solution(1,6);
%v=0.1;
v=solution(1,7);
%h=-1.3; %-0.01,0.6 
H=D20*O20+D40*O40+D60*O60+D44*O44+D64*O64+(h*(Jx+Jy)/(sqrt(2)))+v*Jz;

[eigenvector,SolveEnergy] = eig(H,'vector');
                                                        
                            scattering1=scattering_CEF(eigenvector(:,1),eigenvector(:,2),Jx,Jy,Jz);
        
                            scattering2=scattering_CEF(eigenvector(:,1),eigenvector(:,3),Jx,Jy,Jz);
                            
                            scattering3=scattering_CEF(eigenvector(:,1),eigenvector(:,4),Jx,Jy,Jz);        

                            scattering4=scattering_CEF(eigenvector(:,1),eigenvector(:,5),Jx,Jy,Jz);
                            
                            scattering5=scattering_CEF(eigenvector(:,1),eigenvector(:,6),Jx,Jy,Jz);
                   
                            scattering6=scattering_CEF(eigenvector(:,1),eigenvector(:,7),Jx,Jy,Jz);
                            
                            scattering7=scattering_CEF(eigenvector(:,1),eigenvector(:,8),Jx,Jy,Jz);
                    
                            scattering8=scattering_CEF(eigenvector(:,1),eigenvector(:,9),Jx,Jy,Jz);
                            
                            scattering9=scattering_CEF(eigenvector(:,1),eigenvector(:,10),Jx,Jy,Jz);
                            
                            
                             N=scattering6+scattering7;           
                             s1=scattering1/N;
                             s2=scattering2/N;
                             s3=scattering3/N;
                             s4=scattering4/N;
                             s5=scattering5/N;
                             s6=scattering6/N;
                             s7=scattering7/N;
                             s8=scattering8/N;
                             s9=scattering9/N;
                            testscattering=transpose([s1,s2,s3,s4,s5,s6,s7,s8,s9]);
                            Energytemp = sort(SolveEnergy);
                            testEnergy = Energytemp + abs(min(SolveEnergy(:,1)));
%L+2S=gJ
load('/SNS/users/kyleqianlima/Desktop/CEF/Susceptibilitymeasured.mat')
C=(6.02214086e23)*(9.274e-21)^2/(3*1.38e-16);
Xd=-2.5*98.64e-6;

Xion=zeros(1,69);
%Xion=zeros(1,300);

%T=linspace(1,300,300);
T=Temp';

%Z=zeros(1,300);
Z=zeros(1,69);
Energycgs=Energysol/0.0862;

for n = 1:10
    Z=Z+exp(-Energycgs(n)./T);
end


for n = 1:10
    i = eigenvector(:,n);
    Xion=Xion+((abs((gj*ctranspose(i)*Jz*i))^2).*(exp(-Energycgs(n)./T))./T);
    for m = 1:10
           j = eigenvector(:,m);
        if abs(Energycgs(m)- Energycgs(n)) == 0
            continue;
        else
            Xion=Xion+(((abs(gj*ctranspose(j)*Jz*i))^2).*(exp(-Energycgs(n)./T)-exp(-Energycgs(m)./T))./(Energycgs(m)-Energycgs(n)));
        end
    end
    
end
for n = 1:10
    i = eigenvector(:,n);
    Xion=Xion+((abs((gj*ctranspose(i)*Jx*i))^2).*(exp(-Energycgs(n)./T))./T);
    for m = 1:10
           j = eigenvector(:,m);
        if abs(Energycgs(m)- Energycgs(n)) == 0
            continue;
        else
            Xion=Xion+(((abs(gj*ctranspose(j)*Jx*i))^2).*(exp(-Energycgs(n)./T)-exp(-Energycgs(m)./T))./(Energycgs(m)-Energycgs(n)));
        end
    end
    
end
for n = 1:10
    i = eigenvector(:,n);
    Xion=Xion+((abs((gj*ctranspose(i)*Jy*i))^2).*(exp(-Energycgs(n)./T))./T);
    for m = 1:10
           j = eigenvector(:,m);
        if abs(Energycgs(m)- Energycgs(n)) == 0
            continue;
        else
            Xion=Xion+(((abs(gj*ctranspose(j)*Jy*i))^2).*(exp(-Energycgs(n)./T)-exp(-Energycgs(m)./T))./(Energycgs(m)-Energycgs(n)));
        end
    end
    
end


Xion=(C./Z).*Xion;
X=Xion+Xd;
plot(T,X,'o',Temp,Xmeasure,'g');
legend('calculated','measured');

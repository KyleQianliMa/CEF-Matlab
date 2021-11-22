%L+2S=gJ
C=(6.02214086e23)*(9.274e-21)^2/(3*60*1.38e-16);
Xd=-2.5*98.64e-6;
Xion=zeros(1,300);
T=linspace(1,300,300);

for n = 1:10
    i = eigenvector(:,n);
    Xion=Xion+C.*((abs(gj*ctranspose(i)*Jz*i))^2).*exp(-Energysol(n)./T)./T+...
              C.*((abs(gj*ctranspose(i)*Jx*i))^2).*exp(-Energysol(n)./T)./T+...
              C.*((abs(gj*ctranspose(i)*Jy*i))^2).*exp(-Energysol(n)./T)./T;
end

for n = 1:10
    for m = 1:10
        if abs(Energysol(m)- Energysol(n)) < 2e-14
            continue;
        else
            i=eigenvector(:,n);
            j=eigenvector(:,m);
            Xion=Xion+...
                C.*((abs(gj*ctranspose(j)*Jz*i))^2).*(exp(-Energysol(n)./T)-exp(-Energysol(m)./T))./(Energysol(m)-Energysol(n))+...
                C.*((abs(gj*ctranspose(j)*Jx*i))^2).*(exp(-Energysol(n)./T)-exp(-Energysol(m)./T))./(Energysol(m)-Energysol(n))+...
                C.*((abs(gj*ctranspose(j)*Jy*i))^2).*(exp(-Energysol(n)./T)-exp(-Energysol(m)./T))./(Energysol(m)-Energysol(n));
        end       
    end
end
X=Xion+Xd;
iXion=1./(X);
X=2*3.14159265358*X;
subplot(2,1,1)
plot(T,X,Temp,Xmeasure,'g');
legend('calculated','measured');
subplot(2,1,2)
plot(T,X,Temp,Xmeasure,'g');
legend('calculated','measured');

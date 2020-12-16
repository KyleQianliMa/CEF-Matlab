%L+2S=gJ
C=(6.02214076e23)*(9.274e-24)^2/(3*60*1.38e-23);
Xion=zeros(1,300);
T=linspace(1,300,300);
for n = 1:9
    i = eigenvector(:,n);
    Xion=Xion+C.*((abs(gj*ctranspose(i)*Jz*i))^2).*exp(-Energysol(n)./T)./T+...
              C.*((abs(gj*ctranspose(i)*Jx*i))^2).*exp(-Energysol(n)./T)./T+...
              C.*((abs(gj*ctranspose(i)*Jy*i))^2).*exp(-Energysol(n)./T)./T;
end

for n = 1:9
    for m = 1:9
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
iXion=1./(Xion);
plot(T,iXion);
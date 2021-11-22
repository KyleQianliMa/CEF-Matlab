ct=296; %allsolCopy(ct,1);%


C20=[-284];%allsolCopy(ct,1);%[70 95];%[174];%[212];
C40=[-344];%allsolCopy(ct,2);%[290 330];%[322];%[211];
C60=[-88];%allsolCopy(ct,3);%[-130 -100];%[-193];%[-206];
C44=[93];%allsolCopy(ct,4);%[-295 -270];%[-108];%[-111];
C64=[104];%allsolCopy(ct,5);%[20 50];%[56];%[79];
%green is for 0.26
h_mf=[0.31];
v_mf=[0];
[calcscatteringout,Energysol,solution,bestchi,eigenvector,Jplus,Jminus,Jx,Jy,Jz,allsol,chi,bestsol] = main(C20,C40,C44,C60,C64,h_mf,v_mf);

%g-tensor
      S=3/2;L=6;J=9/2;
      gj=3/2+(S*(S+1)-L*(L+1))/(2*J*(J+1));
      i=eigenvector(:,1);
      j=eigenvector(:,2);
      gz=2*gj*abs(ctranspose(i)*Jz*i);
      gxy=-2*gj*abs(ctranspose(i)*Jplus*j);
    
Susceptibility3
C20=[174];%[174];%[212];
C40=[322];%[322];%[211];
C60=[-193];%[-193];%[-206];
C44=[-108];%[-108];%[-111];
C64=[56];%[56];%[79];
%green is for 0.26
h_mf=[0];
v_mf=[0];
[calcscatteringout,Energysol,solution,bestchi,eigenvector,Jplus,Jminus,Jx,Jy,Jz,allsol,chi,bestsol] = main(C20,C40,C44,C60,C64,h_mf,v_mf);

%g-tensor
      S=3/2;L=6;J=9/2;
      gj=3/2+(S*(S+1)-L*(L+1))/(2*J*(J+1));
      i=eigenvector(:,1);
      j=eigenvector(:,2);
      gz=2*gj*abs(ctranspose(i)*Jz*i);
      gxy=-2*gj*(ctranspose(i)*Jx*j);


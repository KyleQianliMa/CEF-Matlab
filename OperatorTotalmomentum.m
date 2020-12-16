function[O20,O40,O43,O60,O63,O66,Jx,Jy,Jz,Jplus,Jminus] = OperatorTotalmomentum(J)
% basis states are the L,S,Lz,Sz states written in such a way that
% the first 1st line : (-3,-S) 2nd line : (-3,-S+1)...(-3,S),(-2,-S),....etc

%set the size of the operators that will be use for the Hamiltonian
Jplus  = zeros((2.*J+1),(2.*J+1));
Jminus = zeros((2.*J+1),(2.*J+1));
Jz     = zeros((2.*J+1),(2.*J+1));

%set down the matrix Jz
Jz = zeros((2.*J)+1,(2.*J)+1);
for i = 1:((2.*J)+1)
    Jz(i,i) = ((-J) + (i-1));
end

%set down the matrix J+ and J-
for i = 1:(2.*J)
    Jplus(i+1,i) = sqrt( ((J) - (-(J  )+(i-1))).*((J)+((-J    )+(i-1))+1) );
    Jminus(i,i+1)  = sqrt( ((J) + (-(J-1)+(i-1))).*((J)-((-(J-1))+(i-1))+1) );
end

%definition of the Jx and Jy operator
imag = sqrt(-1);
Jx = (0.5).*(Jplus + Jminus);
%Jy = (0.5).*imag.*(Jminus - Jplus);
Jy = (0.5)*imag.*(Jminus - Jplus);

%set the L square matrix
Jsquare = (Jx*ctranspose(Jx)) + (Jy*ctranspose(Jy)) + (Jz*ctranspose(Jz));
Jplussquare = Jplus*Jplus;
Jminussquare = Jminus*Jminus;

%set the CEF hamiltonian
O20 = 3.*(Jz^2) - Jsquare;
O40 = 35.*(Jz^4) - 30.*(Jsquare*(Jz^2)) + 25.*(Jz^2) -6.*(Jsquare) + 3.*(Jsquare^2);
O60 = 231.*(Jz^6) - 315.*(Jsquare*(Jz^4)) + 735.*(Jz^4) + 105.*((Jsquare^2)*(Jz^2)) - 525.*(Jsquare*(Jz^2)) + 294.*(Jz^2) - 5.*(Jsquare^3) + 40.*(Jsquare^2) - 60.*(Jsquare);
O43 = 0.25.*(Jz*(Jplus^3+Jminus^3) + (Jplus^3+Jminus^3)*Jz);
O63 = 0.25.*( ((11.*Jz^3 -3.*(Jsquare*Jz) - 59.*Jz)*(Jplus^3+Jminus^3)) +  ((Jplus^3+Jminus^3))*(11.*Jz^3 -3.*(Jsquare*Jz) - 59.*Jz));
O66 = 0.5.*(Jplussquare^3 + Jminussquare^3);

end
R=1000000000;
index=0;
for ct=1:887
    Initiate
    theory=1./X;
    measure=1./Xmeasure;
    D=0;
    for counter= 1:69
        D=D+abs(theory(counter)-measure(counter));
    end
    if R > D
        R=D;
        index=ct;
    end
end

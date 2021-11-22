newchi=zeros(1,2);
for i = 1:1707
    if chi(i)<3
        newchi=[newchi;i,chi(i)];
    else
        continue
    end
end
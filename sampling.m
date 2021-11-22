B20=[-100 1000];
B40=[-1000 1000];
B60=[-1000 1000];
B44=[-1000 1000];
B64=[-1000 1000];
h_mf=[0];
v_mf=[0];
step=133;
%data=zeros(1,13);
counter = 0;
for C20 = min(B20):step:max(B20)
     for C40 = min(B40):step:max(B40)
         for C60 = min(B60):step:max(B60)
             for C44 = min(B44):step:max(B44)
                 for C64 = min(B64):step:max(B64)
                     for h_mf=min(h_mf):0.1:(h_mf)
                         for v_mf=min(v_mf):0.01:max(v_mf)
                             counter=counter+1;
                             if counter == 100000000000
                                 return
                             else                                 
                                 [calcscattering,CalcEnergy] = sample(C20,C40,C44,C60,C64,h_mf,v_mf);
                                 dt=[C20,C40,C60,C44,C64,CalcEnergy,calcscattering];
                                 data=[data;dt];
                             end
                         end
                     end
                 end
             end
         end
     end
end
%FileData=load('/SNS/users/kyleqianlima/Desktop/CEF/data.mat')
%csvwrite('/SNS/users/kyleqianlima/Desktop/CEF/data1.csv',FileData.data)


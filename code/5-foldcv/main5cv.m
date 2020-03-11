neighbor=3;T=20;
[KL,KM]=Data_fusion_SNF(neighbor,T); 
Random_roder_5fold;
for cv=1:2
    cv
    load('Random_order.mat')
    lambda=0.1;
    res=HGILMA(Random_order(cv,:),lambda);
    RES(cv,:)=res;
end
a=mean(RES,1);
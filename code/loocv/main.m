neighbor=3; T=1;
[KL,KM]=Data_fusion_SNF(neighbor,T); 
for cv=1:1
    cv
    lambda=0.3;
    p=SNFHGILMA_LOOCV(lambda);
    POSITION(cv,:)=p;
    Overallauc(cv)=Plot_roc_curve(p);
end
save Overallauc Overallauc
save POSITION POSITION
a=mean(Overallauc,2);
b=std(Overallauc,0,2);
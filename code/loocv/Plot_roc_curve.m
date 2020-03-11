function AUC=Plot_roc_curve(globalposition)
load('C:\Users\Luo\Desktop\one work\SNFHGILMI11\Data\Known_lncRNA_miRNA_association.mat');
A=unique(Known_lncRNA_miRNA_association,'rows');
nl=max(A(:,1));
nm=max(A(:,2));
pp=size(A,1);

for i=1:pp
    if globalposition(i)>nl*nm-pp+1
        globalposition(i)=nl*nm-pp+1;
    end
end
for k=1:nl*nm-pp+1
    tp=0;
    for t=1:pp
        if globalposition(1,t)<=k
            tp=tp+1;
        end
    end
    tpr(1,k)=tp/pp;
    fp=k*pp-tp;
    fpr(1,k)=fp/(pp*(nl*nm-pp));
end
plot(fpr,tpr)
clear area;
area(1,1)=tpr(1,1)*fpr(1,1)/2;
for k=2:nl*nm-pp+1
    area(1,k)=[tpr(1,k-1)+tpr(1,k)]*[fpr(1,k)-fpr(1,k-1)]/2;
end
AUC=sum(area);

plot(fpr,tpr,'linewidth',1)
xlabel('1-Specificity');ylabel('Sensitivity');
end
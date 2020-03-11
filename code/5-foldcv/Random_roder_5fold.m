load('C:\Users\CJ\Desktop\SNFHGILMI-master\SNFHGILMI-master\Data\Known_lncRNA_miRNA_association.mat');

A=unique(Known_lncRNA_miRNA_association,'rows');
num=size(A,1);

for cv=1:2
    x=randperm(num)';
    Random_order(cv,:)=x;
end

save Random_order.mat Random_order
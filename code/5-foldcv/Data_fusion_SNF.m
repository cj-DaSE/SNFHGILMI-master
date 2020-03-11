function [KL,KM]=Data_fusion_SNF(neighbor,T)

% Import data
load('C:\Users\CJ\Desktop\SNFHGILMI-master\SNFHGILMI-master\Data\Known_lncRNA_miRNA_association.mat');
load('C:\Users\CJ\Desktop\SNFHGILMI-master\SNFHGILMI-master\Data\Sequence_data\invalid_lnc_seq.mat');
load('C:\Users\CJ\Desktop\SNFHGILMI-master\SNFHGILMI-master\Data\Sequence_data\lnc_seq_similarity_matrix.mat');
load('C:\Users\CJ\Desktop\SNFHGILMI-master\SNFHGILMI-master\Data\Sequence_data\mi_seq_similarity_matrix.mat');
load('C:\Users\CJ\Desktop\SNFHGILMI-master\SNFHGILMI-master\Data\Gaussian_kernel_data\lnc_gaussian_similarity_matrix.mat');
load('C:\Users\CJ\Desktop\SNFHGILMI-master\SNFHGILMI-master\Data\Gaussian_kernel_data\mi_gaussian_similarity_matrix.mat');

A=unique(Known_lncRNA_miRNA_association,'rows');
nl=max(A(:,1));
nm=max(A(:,2));
pp=size(A,1);
interaction=zeros(nl,nm);
for i=1:pp
    interaction(A(i,1),A(i,2))=1;
end

lnc_gaussian_similarity=lnc_gaussian_similarity_matrix;
mi_gaussian_similarity=mi_gaussian_similarity_matrix;
mi_seq_similarity=mi_seq_similarity;
lnc_seq_similarity=lnc_seq_similarity;
lnc_seq_similarity(invalid_lnc_seq,:)=[];
lnc_seq_similarity(:,invalid_lnc_seq)=[];

% Normalization of lncRNA/miRNA sequence similarity matrix
renorm = @(M) bsxfun(@rdivide, M, sum(M,2)); 
lnc_seq_similarity1= renorm(lnc_seq_similarity); 
mi_seq_similarity1=renorm(mi_seq_similarity); 

% SNF fusion process
I1=eye(nl,nl);
I2=eye(nm,nm);
Sseq_lnc=(lnc_seq_similarity1+lnc_seq_similarity1')/2+0.5*I1;
Sseq_mi=(mi_seq_similarity1+mi_seq_similarity1')/2+0.5*I2;
Sgauss_lnc=(lnc_gaussian_similarity+lnc_gaussian_similarity')/2;
Sgauss_mi=(mi_gaussian_similarity+mi_gaussian_similarity')/2;

kl1=renorm(Sseq_lnc);
km1=renorm(Sseq_mi);
kgl1=renorm(Sgauss_lnc);
kgm1=renorm(Sgauss_mi);

p1=(kl1+kl1')/2;
p2=(km1+km1')/2;
p3=(kgl1+kgl1')/2;
p4=(kgm1+kgm1')/2;

% obtain k nearest neighbors
k=neighbor;
for n=1:4
    l=[];
    eval(['matrice','=','p', num2str(n)]);
    [row,col]=size(matrice);
    k=neighbor;
    near_array=zeros(row,col);
    near_array_index=zeros(row,col);
    index_array=zeros(row,col);
    for i=1:row
        index_array(i,:)=1:1:col;
    end
    for i=1:row
        [sortRow,index]=sort(matrice(i,:),'descend');
        index_array_row=index_array(i,:);
        index_array_row=index_array_row(index);
        flag=1;
        curmax=sortRow(1);
        for j=1:length(sortRow)
            if flag<=k||curmax==sortRow(j)
                flag=flag+1;
                curmax=sortRow(j);
                near_array(i,j)=sortRow(j);
                near_array_index(i,j)=index_array_row(j);
            end
        end
    end
    for i=1:row
        if sum(near_array(i,:))==0
            sum_near_row=1;
        else
            sum_near_row=sum(near_array(i,:));
        end
        for j=1:col
            if near_array_index(i,j)~=0
                l(i,near_array_index(i,j))=near_array(i,j)/sum_near_row;
            end
        end     
    end
    eval(['l',num2str(n),'=','l']);
    eval(['save ' 'l',num2str(n) '.mat']); 
end

load('l1.mat');%lncRNA 
load('l2.mat');%miRNA
load('l3.mat');%lncRNA
load('l4.mat');%miRNA
for t=1:T
    pl1=l1*p3*l1';
    pgl3=l3*p1*l3';
    pm2=l2*p4*l2';
    pgm4=l4*p2*l4';
    p1=pl1;
    p3=pgl3;
    p2=pm2;
    p4=pgm4;
end

KL=(p1+p3)/2;
KM=(p2+p4)/2;

norm_KL= renorm(KL);
norm_KM= renorm(KM);

Il=eye(nl,nl);
Im=eye(nm,nm);
KL=(norm_KL+norm_KL'+Il)/2;
KM=(norm_KM+norm_KM'+Im)/2;

save KL KL;
save KM KM;
end






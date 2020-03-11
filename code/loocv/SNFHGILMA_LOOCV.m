function globalposition=SNFHGILMA_LOOCV(lambda)
%%predict lncRNA-related miRNA based on SNFHGILMI in the term of LOOCV 
%A: Binary relations between lncRNA and miRNA, 1st column:lncRNA, 2nd column:miRNA

load('C:\Users\Luo\Desktop\one work\SNFHGILMI11\Data\Known_lncRNA_miRNA_association.mat');
load('C:\Users\Luo\Desktop\one work\SNFHGILMI11\code\loocv\KL.mat');
load('C:\Users\Luo\Desktop\one work\SNFHGILMI11\code\loocv\KM.mat');
A=unique(Known_lncRNA_miRNA_association,'rows');
nl=max(A(:,1));
nm=max(A(:,2));
pp=size(A,1);
interaction=zeros(nl,nm);
for i=1:pp
    interaction(A(i,1),A(i,2))=1;
end
save interaction interaction;
original_interaction=interaction;


%implement leave-one-out cross validation
for ccv=1:pp
    ccv
    % obtain training sample
    load interaction;
    interaction(A(ccv,1),A(ccv,2))=0;

    %SNFHGILMI
    for i1=1:nl
        for j1=1:nl
            SL(i1,j1)=KL(i1,j1)/(sqrt(sum(KL(i1,:)))*sqrt(sum(KL(:,j1))));
        end
    end
    for i2=1:nm
        for j2=1:nm
            SM(i2,j2)=KM(i2,j2)/(sqrt(sum(KM(i2,:)))*sqrt(sum(KM(:,j2))));
        end
    end
    Pt=interaction;
    for t=1:30
        Ptt=lambda*SL*Pt*SM+(1-lambda)*interaction;
        delta=norm((Pt-Ptt),1);
        Pt=Ptt;
        if  delta < (1e-6) 
            %t
            break;
        end
    end
    F_ensemble=Ptt;
    
    % obtain the score of tested  lncRNA-miRNA interaction
    finalscore2=F_ensemble(A(ccv,1),A(ccv,2));
%end
% make the score of training lncRNA-miRNA interactions as zero
for i=1:nl
    for j=1:nm
        if original_interaction(i,j)==1
            F_ensemble(i,j)=-999999999999999999;
        end
    end
end
% obtain the position of tested lncRNA-miRNA interaction as variable globalposition(1,cv),

ll1=size(find(F_ensemble>=finalscore2),1);
ll2=size(find(F_ensemble>finalscore2),1);
globalposition(1,ccv)=ll2+1+(ll1-ll2-1)/2;
end
save globalposition.mat globalposition
end
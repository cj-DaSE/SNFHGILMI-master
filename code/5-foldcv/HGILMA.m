function result=HGILMA(x,lambda)
%predict disease-related microbe based on SNFHGILMI  in the term of 5-fold cross validation
%A: Binary relations between disease and microbe, 1st column:lncRNA, 2nd
%column:miRNA

load('C:\Users\CJ\Desktop\SNFHGILMI-master\SNFHGILMI-master\Data\Known_lncRNA_miRNA_association.mat');
load('C:\Users\CJ\Desktop\SNFHGILMI-master\SNFHGILMI-master\code\5-foldcv\KL.mat');
load('C:\Users\CJ\Desktop\SNFHGILMI-master\SNFHGILMI-master\code\5-foldcv\KM.mat');
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

%implement 5-fold cross validation

result=zeros(1,4);
for ccv=1:5
    %ccv
    
    load interaction interaction;
    if ccv<5
        AA=A(x((ccv-1)*floor(pp/5)+1:floor(pp/5)*ccv),:);
        % obtain training sample
        for i=1:floor(pp/5)
            interaction(AA(i,1),AA(i,2))=0;
        end
    else
        AA=A(x((ccv-1)*floor(pp/5)+1:pp),:);
        % obtain training sample
        for i=1:pp-floor(pp/5)*4
            interaction(AA(i,1),AA(i,2))=0;
        end
    end
    
    % SNFHGILMI
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
    F=Ptt; 
    save F F;
    result=result+model_evaluate(original_interaction,F,interaction);
    
end
result=result/ccv;
end



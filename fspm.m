%Subspace-based Preprocessing Module for Fast Hyperspectral Endmember Subset Selection

function [R, time,bIdxRecIner,oIdxRecIner,oR]=fspm(HIM,q,ratio,Kn,delta,showResults,type)

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%       Input:
%        HIM: [nr,nc,nb] dataset
%        q:   endmember numbers
%        ratio: percentatge of pixels retained from subsquent EE
%       showResults: visually display convex hull points
%        Kn: number of neighborhood
%        delta: local density score threshold
%--------------------------------------------------------------------------
% Copyright (May, 2020):    Wenxing Bao (bwx71@163.com)
%                           Xiangfei Shen (xfshen95@163.com)
%
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ----------------------------------------------------------------------
%
% More details in:
%
% [1] Xiangfei Shen, Wenxing Bao, Kewen Qu. Spatial-spectral hyperspectral 
% endmember extraction using a spatial energy prior constrained maximum simplex 
% volume approach[J]. IEEE Journal of Selected Topics in Applied Earth Observations 
% and Remote Sensing, 2020, 13: 1347-1361.
%
% [2] Xiangfei Shen, Wenxing Bao, Kewen Qu. Subspace-based preprocessing 
% module for fast hyperspectral endmember selection[J]. IEEE Journal of Selected 
% Topics in Applied Earth Observations and Remote Sensing, 2021, 14:3386-3402.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%
[nr,nc,nb]=size(HIM);
Y=reshape(HIM,nr*nc,nb)';
[~,N]=size(Y);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% subspace transforming
eta=10;type3=0;
tic
[Ud, ~, ~] = svds((Y*Y.')/N, q);
Xd = Ud.'*Y;
u = mean(Xd, 2);
Y =  Xd ./ repmat( sum( Xd .* repmat(u,[1 N]) ) ,[q 1]);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% display results
Yshow=Y;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stopNum=ceil(N*ratio);  
stopNumAux=0;  
x=linspace(1,q,q);
r=nchoosek(x,2);
[rr,~]=size(r);
iter=1;
%% 
Ycov=Y;
bIdxRecIner=[];
oIdxRecIner=[];
idxRecIner=linspace(1,N,N);

while stopNumAux<stopNum

    for i=1:rr
 
        yAux=Ycov(r(i,:),:);
        
        K=convhull(yAux');
        K=unique(K);
        
        
        %% standardization
        yAux=[(yAux(1,:)-mean(yAux(1,:)))/std(yAux(1,:));(yAux(2,:)-mean(yAux(2,:)))/std(yAux(2,:))];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if type3==1
            figure
            plot(yAux(1,:),yAux(2,:),'o');
            hold on
            plot(yAux(1,K),yAux(2,K),'b>')
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [minVal1,~]=min(yAux(1,:));
        [minVal2,~]=min(yAux(2,:));
        [maxVal1,~]=max(yAux(1,:));
        [maxVal2,~]=max(yAux(2,:));
        
        sigmaA=(maxVal1-minVal1)/eta;
        sigmaB=(maxVal2-minVal2)/eta;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        kY=yAux(:,K);
        
        [~,minIdx1]=min(kY(1,:));
        [~,minIdx2]=min(kY(2,:));
        [~,maxIdx1]=max(kY(1,:));
        [~,maxIdx2]=max(kY(2,:));
        
        ex=K([minIdx1,minIdx2,maxIdx1,maxIdx2]);
        ex=unique(ex);
        exk=length(ex);
        x=yAux(:,ex);
        
        isOutlier=zeros(1,exk);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if type3==1
            plot(yAux(1,ex),yAux(2,ex),'ro')
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for jj=1:exk
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if type3==1
                plot(x(1,jj),x(2,jj),'r*')
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            idxAux=yAux(1,:)<x(1,jj)+sigmaA & yAux(1,:) > x(1,jj)-sigmaA;
            idxAux2  = yAux(2,idxAux==1)<x(2,jj)+sigmaB & yAux(2,idxAux==1) > x(2,jj)-sigmaB;
            idx3=find(idxAux==1);
            idx4=idx3(idxAux2==1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if type3==1
                plot(yAux(1,idx4),yAux(2,idx4),'g.')
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            leng=length(idx4);
            
            if leng==1
                
                density=0;
                densityN=inf;
            else
                aux=sqrt(sum((repmat(x(:,jj),1,leng)-yAux(:,idx4)).^2));
                [aux2,aux2idx]=sort(aux);
                aux2idx=idx4(aux2idx);
                if leng<=Kn
                    density=(leng-1)/(Kn*sum(aux2));
                    maxidx=aux2idx(leng);
                else
                    density=1/sum(aux2(1:Kn+1));
                    maxidx=aux2idx(Kn+1);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if type3==1
                    if length(aux2idx)<=Kn
                        plot(yAux(1,aux2idx(1:end)),yAux(2,aux2idx(1:end)),'b^')
                    else
                        plot(yAux(1,aux2idx(1:Kn)),yAux(2,aux2idx(1:Kn)),'b^')
                    end
                    plot(yAux(1,maxidx),yAux(2,maxidx),'r.')
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                idxAux=yAux(1,:)<yAux(1,maxidx)+sigmaA & yAux(1,:) > yAux(1,maxidx)-sigmaA;
                idxAux2  = yAux(2,idxAux==1)<yAux(2,maxidx)+sigmaB & yAux(2,idxAux==1) > yAux(2,maxidx)-sigmaB;
                idx3=find(idxAux==1);
                idx4=idx3(idxAux2==1);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if type3==1
                    plot(yAux(1,idx4),yAux(2,idx4),'r.')
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                leng=length(idx4); 
                aux=sqrt(sum((repmat(yAux(:,maxidx),1,leng)-yAux(:,idx4)).^2));
                [aux2,aux2idx]=sort(aux);
                aux2idx=idx4(aux2idx);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if type3==1
                    if length(aux2idx)<=Kn
                        plot(yAux(1,aux2idx(1:end)),yAux(2,aux2idx(1:end)),'b^')
                    else
                        plot(yAux(1,aux2idx(1:Kn)),yAux(2,aux2idx(1:Kn)),'b^')
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if leng<=Kn
                    densityN=(leng-1)/(Kn*sum(aux2));
                else
                    densityN=1/sum(aux2(1:Kn+1));
                end

            end
            
            if density/densityN < delta
                isOutlier(jj)=1;
            end
        end

        Ko=ex(isOutlier==1)';
        if isempty(Ko)
            Ycov(:,K)=[];
            bIdxRecIner=[bIdxRecIner,idxRecIner(K)];
            idxRecIner(K)=[];
        else
            K(K==Ko)=[];
            Ycov(:,[K;Ko])=[];
            bIdxRecIner=[bIdxRecIner,idxRecIner(K)];
            oIdxRecIner=[oIdxRecIner,idxRecIner(Ko)];
            idxRecIner([K;Ko])=[];
        end

        %%
        fprintf(' *** %i pixels (%i are removed) in %i th pair subspace (%i pairs) in %i iter ***\n',length(K),length(Ko),i,rr, iter);

    end

    
    
    %% 
    
    [~,stopNumAux]=size(bIdxRecIner);

    
    
    %fprintf(' >>> %i pixels are retained in %i iter: targeted %i pixels <<<\n',length(idx),iter,stopNum);
    %fprintf('-------  -------  -------  -------  -------  -------  -------\n')
    
    iter=iter+1;
    %     idx=[];
    %     outIdx=[];
end

R = Ud*Xd(:,bIdxRecIner);
oR =Ud*Xd(:,oIdxRecIner);
time=toc;
fprintf('-------  -------  -------  -------  -------  -------  -------\n')
fprintf('------- %i pixels retained (%.2f%)     -------\n',length(R(1,:)), length(R(1,:))/N*100)
fprintf('------- %.2f total execution time (s)  -------\n',time)
fprintf('-------  -------  -------  -------  -------  -------  -------\n')
%fprintf('**** execution   time: %.2f seconds    ****\n',Execution_Time_in_Seconds)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if showResults==1
    if type==1
        for i=1:rr
            K=convhull(Yshow(r(i,:),:)');
            
            YY=Yshow(r(i,:),:)';
            figure(i)
            plot(YY(:,1),YY(:,2),'o')
            set(gcf,'unit','normalized','position',[0.4,0.4,0.22,0.3]);
            %xlabel('1st Component')
            %ylabel('2nd Component')
            
            hold on
            
            plot(YY(K,1),YY(K,2),'o','MarkerFaceColor','r','MarkerSize',5)
            %plot(YY(outIdx,1),YY(outIdx,2),'x','MarkerSize',5)
            %plot(YY(1,1),YY(1,2),'^','MarkerFaceColor','g','MarkerSize',3)
            %plot(YY(2,1),YY(2,2),'d','MarkerFaceColor','m','MarkerSize',3)
            %plot(YY(3,1),YY(3,2),'s','MarkerFaceColor','k','MarkerSize',3)
            hold off
        end
    else
        YY=Yshow([3,2],:)';
        figure
        plot(YY(:,1),YY(:,2),'o')
        set(gcf,'unit','normalized','position',[0.4,0.4,0.22,0.3]);
        hold on
        plot(YY(bIdxRecIner,1),YY(bIdxRecIner,2),'ro')
        plot(YY(oIdxRecIner,1),YY(oIdxRecIner,2),'go')
    end
    %
    % end
end
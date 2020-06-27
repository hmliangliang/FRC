function [ J,R,FM,CD,K,RT,NMI] =Evaluation(target,cluster,U,Point,signal)
%此函数主要功能是评价聚类的结果  计算Jaccard coefficient、Rand statistic、Fowlkes and Mallows
%index、Purity、Mean square error:
%输入 cluster:表示聚类的结果,类型为cell型,每一个cell单元即为一个类簇,保存的是数据的序号    target:数据的真实类标签
%U:为所有样本的数据   Point:表示是聚类类簇的中心点   signal:表示是采用分类型计算距离,还是使用数值型距离,signal=1表示采用数值型方法计算距离,signal=0表示使用分类型方法计算距离
SS=0;
SD=0;
DS=0;
DD=0;
C=zeros(size(target,1),1);%保存与target顺序相对应的类标签
for i=1:size(cluster,2)
    C(cluster{i},:)=i;
end
for i=1:size(target,1)-1
    for j=i+1:size(target,1)
        if C(i,1)==C(j,1)%两个数据在聚类中处于同一类簇
            if target(i,1)==target(j,1)%在真实情况中也是处于同一类
                SS=SS+1;
            else%在真实情况中是处于不同的类
                SD=SD+1;
            end
        else%两个数据在聚类中处于不同的类簇
            if target(i,1)==target(j,1)%在真实情况中也是处于同一类
                DS=DS+1;
            else%在真实情况中是处于不同的类
                DD=DD+1;
            end
        end
    end
end
J=SS/(SS+SD+DS);%计算Jaccard coefficient
R=(SS+DD)/(SS+DD+SD+DS);%计算Rand statistic
FM=sqrt(SS*SS/((SS+SD)*(SS+DS)));%计算Fowlkes and Mallows index
%以下主要MSE
% mse=0;
% if signal==0%采用分类型方法计算距离
%     for i=1:size(cluster,2)
%         for j=1:length(cluster{i})
%            mse=mse+size(U,2)-sum(Point(i,:)==U(cluster{i}(j),:));
%         end
%     end
% elseif signal==1%采用数值型方法计算距离
%      for i=1:size(cluster,2)
%         for j=1:length(cluster{i})
%            mse=mse+norm(Point(i,:)-U(cluster{i}(j),:),2);
%         end
%     end
% end
% mse=mse/(size(target,1)*(size(U,2)));
% disp(['本算法的指标结果：J=',num2str(J)]);
% disp(['本算法的指标结果：R=',num2str(R)]);
% disp(['本算法的指标结果：FM=',num2str(FM)]);
% disp(['本算法的指标结果：P=',num2str(P)]);
% disp(['本算法的指标结果：MSE=',num2str(mse)]);
CD=2*SS/(2*SS+DS+SD);
K=0.5*(SS/(SS+SD)+SS/(SS+DS));
RT=(DD+SS)/(DD+SS+2*(DS+SD));
NMI=nmi( C,target);
end

function MIhat = nmi( A, B )
%NMI Normalized mutual information
% http://en.wikipedia.org/wiki/Mutual_information
% http://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html
% Author: http://www.cnblogs.com/ziqiao/   [2011/12/15]
 
if length( A ) ~= length( B)
    error('length( A ) must == length( B)');
end
if iscolumn(A)
    A=A';
end
if iscolumn(B)
    B=B';
end
total = length(A);
A_ids = unique(A);
A_class = length(A_ids);
B_ids = unique(B);
B_class = length(B_ids);
% Mutual information
idAOccur = double (repmat( A, A_class, 1) == repmat( A_ids', 1, total )); %得到节点社区矩阵N*C
idBOccur = double (repmat( B, B_class, 1) == repmat( B_ids', 1, total ));
idABOccur = idAOccur * idBOccur';
Px = sum(idAOccur') / total;
Py = sum(idBOccur') / total;
Pxy = idABOccur / total;
MImatrix = Pxy .* log2(Pxy ./(Px' * Py)+eps);
MI = sum(MImatrix(:));
% Entropies
Hx = -sum(Px .* log2(Px + eps),2);
Hy = -sum(Py .* log2(Py + eps),2);
%Normalized Mutual information
MIhat = 2 * MI / (Hx+Hy);
end

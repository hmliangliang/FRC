function [ J,R,FM,CD,K,RT,NMI] =Evaluation(target,cluster,U,Point,signal)
%�˺�����Ҫ���������۾���Ľ��  ����Jaccard coefficient��Rand statistic��Fowlkes and Mallows
%index��Purity��Mean square error:
%���� cluster:��ʾ����Ľ��,����Ϊcell��,ÿһ��cell��Ԫ��Ϊһ�����,����������ݵ����    target:���ݵ���ʵ���ǩ
%U:Ϊ��������������   Point:��ʾ�Ǿ�����ص����ĵ�   signal:��ʾ�ǲ��÷����ͼ������,����ʹ����ֵ�;���,signal=1��ʾ������ֵ�ͷ����������,signal=0��ʾʹ�÷����ͷ����������
SS=0;
SD=0;
DS=0;
DD=0;
C=zeros(size(target,1),1);%������target˳�����Ӧ�����ǩ
for i=1:size(cluster,2)
    C(cluster{i},:)=i;
end
for i=1:size(target,1)-1
    for j=i+1:size(target,1)
        if C(i,1)==C(j,1)%���������ھ����д���ͬһ���
            if target(i,1)==target(j,1)%����ʵ�����Ҳ�Ǵ���ͬһ��
                SS=SS+1;
            else%����ʵ������Ǵ��ڲ�ͬ����
                SD=SD+1;
            end
        else%���������ھ����д��ڲ�ͬ�����
            if target(i,1)==target(j,1)%����ʵ�����Ҳ�Ǵ���ͬһ��
                DS=DS+1;
            else%����ʵ������Ǵ��ڲ�ͬ����
                DD=DD+1;
            end
        end
    end
end
J=SS/(SS+SD+DS);%����Jaccard coefficient
R=(SS+DD)/(SS+DD+SD+DS);%����Rand statistic
FM=sqrt(SS*SS/((SS+SD)*(SS+DS)));%����Fowlkes and Mallows index
%������ҪMSE
% mse=0;
% if signal==0%���÷����ͷ����������
%     for i=1:size(cluster,2)
%         for j=1:length(cluster{i})
%            mse=mse+size(U,2)-sum(Point(i,:)==U(cluster{i}(j),:));
%         end
%     end
% elseif signal==1%������ֵ�ͷ����������
%      for i=1:size(cluster,2)
%         for j=1:length(cluster{i})
%            mse=mse+norm(Point(i,:)-U(cluster{i}(j),:),2);
%         end
%     end
% end
% mse=mse/(size(target,1)*(size(U,2)));
% disp(['���㷨��ָ������J=',num2str(J)]);
% disp(['���㷨��ָ������R=',num2str(R)]);
% disp(['���㷨��ָ������FM=',num2str(FM)]);
% disp(['���㷨��ָ������P=',num2str(P)]);
% disp(['���㷨��ָ������MSE=',num2str(mse)]);
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
idAOccur = double (repmat( A, A_class, 1) == repmat( A_ids', 1, total )); %�õ��ڵ���������N*C
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

function [Z,result] = FRC(data,k,gamma)
%此算法执行的FRC算法  data:数据的样本,每一行代表一个样本
%result是一个k*1的结构体
col=size(data,2);%样本的列数
w=zeros(1,col);%保存属性的权值
N_MAX = 100;%最大的迭代次数
for i=1:col
    value=unique(data(:,i));%获取每列的数据，即每个属性的不同值
    for j=1:size(value,1)
        w(i)=w(i)+(length(find(data(:,i)==value(j,1))))/(size(data,1))^2;
    end
    w(i)=w(i)/size(data,1);%计算每个属性划分的粒度
end
value=sum(w);
for i=1:col%计算每个属性的粒度比值
    w(i)=value/w(i);
end
value=sum(w);
for i=1:col%计算每个属性的权值
    w(i)=w(i)/value;
end
red=[];%记录冗余的属性
for i=1:size(w,2)
    if w(i)<mean(w)-3*std(w)
        red=[red,i];
    end
end
data(red,:)=[];%消除冗余属性
Z=data(randperm(size(data,1),k),:);%随机选取k个mode点进行聚类
f=inf;%初始化目标函数为无穷大
M=M_center(data,Z,w);%初始化样本与聚类中心的相似度矩阵
M_before=M;
Z_before=Z;
f_before=f_obj(data,M,Z,w,gamma);
for i=1:N_MAX%进行迭代
    if f<=f_obj(data,M,Z,w,gamma)
       M=M_center(data,Z,w);%更新M
       [~,re]=max(M,[],2);%确定各个样本的归属
       %更新modes
       for j=1:k
           Z(j,:)=mode(data(find(re(:,1)==j),:));
       end
       f=f_obj(data,M,Z,w,gamma);%计算当前目标函数值
       if f<f_before%目标函数值减小
          M_before=M;
          Z_before=Z;
          f_before=f;
       else
           M=M_before;
           Z=Z_before;
           break;
       end
    else
        M=M_before;
        Z=Z_before;
        break;
    end
end
[~,cresult]=max(M,[],2);%确定各个样本的归属
% for i=1:size(M,2)
%     maxvalue=max(M(:,i));
%     minvalue=min(M(:,i));
%     for j=1:size(M,1)
%         if maxvalue>minvalue
%             M(j,i)=(M(j,i)-minvalue)/(maxvalue-minvalue);
%         end
%     end
% end
%cresult=kmeans(M,k);
value=unique(cresult);
result=cell(1,k);
for i=1:k
    result{i}=find(cresult(:,1)==value(i,1))';
end
end

function f=f_obj(data,M,Z,w,gamma)%计算目标函数值
f=0;
for i=1:size(data,1)
    for j=1:size(M,2)
        f=f+(1-M(i,j)*simlarity(data,data(i,:),Z(j,:),w))^2+gamma*M(i,j)*log(M(i,j));
    end
end
end

function sim=simlarity(data,data1,data2,w)%按照Eq.(12)计算样本之间的相似度
sim=0;
for i=1:size(data,2)%每个维度
    sim=sim+w(i)*length(intersect(find(data(:,i)==data1(i)),find(data(:,i)==data2(i))))/length(union(find(data(:,i)==data1(i)),find(data(:,i)==data2(i))));
end
sim=sim/size(data,2);%确保sim范围在[0,1]之间
end

function M=M_center(data,Z,w)%按照Eq.(13)计算样本与中心点之间的相似度
distance=zeros(size(data,1),size(Z,1));%初始化距离矩阵
M=zeros(size(data,1),size(Z,1));
for i=1:size(M,1)%每一个样本
    for j=1:size(M,2)%每一个中心点
        distance(i,j)=simlarity(data,data(i,:),Z(j,:),w);
    end
end
for i=1:size(M,1)%每一个样本
    value=sum(distance(i,:));
    for j=1:size(M,2)%每一个中心点
        M(i,j)=distance(i,j)/value;
    end
end
end
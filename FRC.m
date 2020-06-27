function [Z,result] = FRC(data,k,gamma)
%���㷨ִ�е�FRC�㷨  data:���ݵ�����,ÿһ�д���һ������
%result��һ��k*1�Ľṹ��
col=size(data,2);%����������
w=zeros(1,col);%�������Ե�Ȩֵ
N_MAX = 100;%���ĵ�������
for i=1:col
    value=unique(data(:,i));%��ȡÿ�е����ݣ���ÿ�����ԵĲ�ֵͬ
    for j=1:size(value,1)
        w(i)=w(i)+(length(find(data(:,i)==value(j,1))))/(size(data,1))^2;
    end
    w(i)=w(i)/size(data,1);%����ÿ�����Ի��ֵ�����
end
value=sum(w);
for i=1:col%����ÿ�����Ե����ȱ�ֵ
    w(i)=value/w(i);
end
value=sum(w);
for i=1:col%����ÿ�����Ե�Ȩֵ
    w(i)=w(i)/value;
end
red=[];%��¼���������
for i=1:size(w,2)
    if w(i)<mean(w)-3*std(w)
        red=[red,i];
    end
end
data(red,:)=[];%������������
Z=data(randperm(size(data,1),k),:);%���ѡȡk��mode����о���
f=inf;%��ʼ��Ŀ�꺯��Ϊ�����
M=M_center(data,Z,w);%��ʼ��������������ĵ����ƶȾ���
M_before=M;
Z_before=Z;
f_before=f_obj(data,M,Z,w,gamma);
for i=1:N_MAX%���е���
    if f<=f_obj(data,M,Z,w,gamma)
       M=M_center(data,Z,w);%����M
       [~,re]=max(M,[],2);%ȷ�����������Ĺ���
       %����modes
       for j=1:k
           Z(j,:)=mode(data(find(re(:,1)==j),:));
       end
       f=f_obj(data,M,Z,w,gamma);%���㵱ǰĿ�꺯��ֵ
       if f<f_before%Ŀ�꺯��ֵ��С
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
[~,cresult]=max(M,[],2);%ȷ�����������Ĺ���
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

function f=f_obj(data,M,Z,w,gamma)%����Ŀ�꺯��ֵ
f=0;
for i=1:size(data,1)
    for j=1:size(M,2)
        f=f+(1-M(i,j)*simlarity(data,data(i,:),Z(j,:),w))^2+gamma*M(i,j)*log(M(i,j));
    end
end
end

function sim=simlarity(data,data1,data2,w)%����Eq.(12)��������֮������ƶ�
sim=0;
for i=1:size(data,2)%ÿ��ά��
    sim=sim+w(i)*length(intersect(find(data(:,i)==data1(i)),find(data(:,i)==data2(i))))/length(union(find(data(:,i)==data1(i)),find(data(:,i)==data2(i))));
end
sim=sim/size(data,2);%ȷ��sim��Χ��[0,1]֮��
end

function M=M_center(data,Z,w)%����Eq.(13)�������������ĵ�֮������ƶ�
distance=zeros(size(data,1),size(Z,1));%��ʼ���������
M=zeros(size(data,1),size(Z,1));
for i=1:size(M,1)%ÿһ������
    for j=1:size(M,2)%ÿһ�����ĵ�
        distance(i,j)=simlarity(data,data(i,:),Z(j,:),w);
    end
end
for i=1:size(M,1)%ÿһ������
    value=sum(distance(i,:));
    for j=1:size(M,2)%ÿһ�����ĵ�
        M(i,j)=distance(i,j)/value;
    end
end
end
%������
tic;
data=Horse;%��������
k=2;
JJ=[];
RR=[];
FFM=[];
PP=[];
MSE=[];
RRT=[];
NNMI=[];
PPP=[];
N_MAX=50;
label=data(:,size(data,2));
data=data(:,1:(size(data,2)-1));
for i=1:N_MAX
    %[ C,p ] = DKmodes(data,k);P=p;
    %[ C,w,pnew ] = WKModes(data,k);P=pnew;
    %[C,Z]=kmodes(data,k);P=Z;
    gamma=2;
    %gamma=0.2;
    [Z,C] = FRC(data,k,gamma);P=Z;
    [ J,R,FM,CD,K,RT,NMI] =Evaluation(label,C,data,P,0);%0:��ʾ���������÷�����  1:��ʾ������������ֵ��
    JJ=[JJ,J];
    RR=[RR,R];
    FFM=[FFM,FM];
    PP=[PP,CD];
    MSE=[MSE,K];
    RRT=[RRT,RT];
    NNMI=[NNMI,NMI];
end
disp(['���㷨����J��ƽ��ֵΪ��',num2str(mean(JJ)),'$\pm$',num2str(std(JJ))]);
disp(['���㷨����R��ƽ��ֵΪ��',num2str(mean(RR)),'$\pm$',num2str(std(RR))]);
disp(['���㷨����FM��ƽ��ֵΪ��',num2str(mean(FFM)),'$\pm$',num2str(std(FFM))]);
disp(['���㷨����CD��ƽ��ֵΪ��',num2str(mean(PP)),'$\pm$',num2str(std(PP))]);
disp(['���㷨����K��ƽ��ֵΪ��',num2str(mean(MSE)),'$\pm$',num2str(std(MSE))]);
disp(['���㷨����RT��ƽ��ֵΪ��',num2str(mean(RRT)),'$\pm$',num2str(std(RRT))]);
disp(['���㷨����NMI��ƽ��ֵΪ��',num2str(mean(NNMI)),'$\pm$',num2str(std(NNMI))]);
toc;
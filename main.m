%主函数
tic;
data=Horse;%导入数据
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
    [ J,R,FM,CD,K,RT,NMI] =Evaluation(label,C,data,P,0);%0:表示计算距离采用分类型  1:表示计算距离采用数值型
    JJ=[JJ,J];
    RR=[RR,R];
    FFM=[FFM,FM];
    PP=[PP,CD];
    MSE=[MSE,K];
    RRT=[RRT,RT];
    NNMI=[NNMI,NMI];
end
disp(['本算法性能J的平均值为：',num2str(mean(JJ)),'$\pm$',num2str(std(JJ))]);
disp(['本算法性能R的平均值为：',num2str(mean(RR)),'$\pm$',num2str(std(RR))]);
disp(['本算法性能FM的平均值为：',num2str(mean(FFM)),'$\pm$',num2str(std(FFM))]);
disp(['本算法性能CD的平均值为：',num2str(mean(PP)),'$\pm$',num2str(std(PP))]);
disp(['本算法性能K的平均值为：',num2str(mean(MSE)),'$\pm$',num2str(std(MSE))]);
disp(['本算法性能RT的平均值为：',num2str(mean(RRT)),'$\pm$',num2str(std(RRT))]);
disp(['本算法性能NMI的平均值为：',num2str(mean(NNMI)),'$\pm$',num2str(std(NNMI))]);
toc;
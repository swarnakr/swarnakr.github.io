clc
close all
clear all
NoClasses=3
naive=0
Data=load('./1D/group4.txt');

plot(Data(:,1),'*','Color','g');
for i=1:NoClasses
    Class(:,:,i) = Data( (i-1)*500 + 1: i*500 ,:);
end

% TAKING 75% DATA FOR TRAINING
t=500*0.75; train=[];
for i=1:NoClasses
    train=[train; Class(1:t,:,i)];
end

%TAKING 25% DATA FOR TESTING
t1=0.25*500; test=[];
for i=1:NoClasses
    test=[test; Class(t+1:end,:,i)];
end

%CALCULATING MEAN (ML)
for i=1:NoClasses
    M(i,:)=mean(Class(1:t,:,i));
end
%CALCULATING COVARIANCE MATRIX (ML)
for i=1:NoClasses
    C(:,:,i)=cov(Class(1:t,:,i));
end

%% NAIVE BAYES
if naive

    for i=1:NoClasses
        tempdiag=diag(Class(:,:,i));
        Class(:,:,i) = diag(tempdiag);
    end
end

%% TRAINING
[rows,cols]=size(train);
for i=1:NoClasses
    x = train - repmat(M(i,:),rows,1);
    g(:,i) = -0.5 * diag((x*inv(C(:,:,i))*x'))  + repmat(log(1/3),rows,1) - repmat(0.5*log(det(C(:,:,i))),rows,1); %log-natural logarithm
end

FinalClass={};
[maxval,index]=max(g,[],2);
for i=1:NoClasses
    u=find(index==i);   
    FinalClass{i}=train(u,:);
end
[m,n]=size(g)

%PLOTTING Classes
figure; xlabel('X'),ylabel('Y'); hold on;
for i=1:NoClasses
    u=find(index==i);
    plot(u,FinalClass{i}(:,1),'*','Color',rand(1,3)), hold on;
 
end
% xL = get(gca,'XLim');
% yL = get(gca,'YLim');
% pts=[];
% %pts=pts'
% for i=1:10
%     for j=1:10
%                pts=[pts;i j]
%     end
% end
%  [rows,cols]=size(pts);
% for i=1:NoClasses
%     x = pts(1) -repmat(M(i,:),rows,1);
%      %x(2) = pts(2) -repmat(M(i,:),rows,1);
%         g(:,i) = -0.5 * diag((x*inv(C(:,:,i))*x'))  + repmat(log(1/3),rows,1) - repmat(0.5*log(det(C(:,:,i))),rows,1); %log-natural logarithm
% end
% 
% FinalClass={};
% [maxval,index]=max(g,[],2);
% for i=1:NoClasses
%     u=find(index==i);   
%     FinalClass{i}=train(u,:);
% end
% 
% 
% %PLOTTING Classes
% figure; xlabel('X'),ylabel('Y'); hold on;
% for i=1:NoClasses
%     u=find(index==i); 
%     plot(pts,FinalClass{i}(:,1),'o','Color',rand(1,3)), hold on;
% end
%    
%   
% legend('Class1','Class2','Class3');

%% TESTING

[rows,cols]=size(test);

for i=1:NoClasses
    x = test - repmat(M(i,:),rows,1);
    gTest(:,i) = -0.5 * diag((x*inv(C(:,:,i))*x'))  + repmat(log(1/3),rows,1) - repmat(0.5*log(det(C(:,:,i))),rows,1); %log-natural logarithm
end

FinalTestClass={};
[maxval,index]=max(gTest,[],2);
classBelongs = index;
for i=1:NoClasses
    u=find(index==i);   
    FinalTestClass{i}=test(u,:);
end


%PLOTTING Classes
figure; xlabel('X'),ylabel('Y'); hold on;
for i=1:NoClasses
    u=find(index==i);   
    plot(u,FinalTestClass{i}(:,1),'*','Color',rand(1,3)), hold on;
end
 xlimits = get(gca,'XLim');
ylimits = get(gca,'YLim');
a=xlimits
b=ylimits
legend('Class1','Class2','Class3','Class4');
x=a(1,1):.1:a(1,2)
y=b(1,1):.1:b(1,2)
i=size(x)
j=size(y)
for i=1:size(x)
    for j=1:size(y)
plot(x(i),y(j),'*','color','g')
    end
end

%PLOTTING DECISION REGIONS


%%

%Must be changed
for i=1:NoClasses
trgt((i-1)*t1 + 1: i*t1 )=i;
end

%CALCULATING CONFUSION MATRIX
conf_matrix=zeros(NoClasses);
for i=1:size(trgt,2)
    conf_matrix(trgt(i),classBelongs(i)) = conf_matrix(trgt(i),classBelongs(i))+ 1;
end
% Do Similarly for impostor


%CALCULATING CLASSIFICATION ACCURACY FOR ALL CLASSES
accuracy=sum(diag(conf_matrix))/(size(trgt,2))*100;

% CALCULATING CLASSIFICATION ACCURACY FOR EACH CLASS
for i=1:NoClasses
    class_acc(i)=(sum(conf_matrix(i,i))/(t1))*100;
end



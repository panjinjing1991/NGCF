clearvars;
%% Construct series
% 构建序列

% jd = [1:365,1:365,1:365];

% a=1.5;
% b=0.4;
% c=0.5;

a = 1.8;

%len = numel(jd);

len = 1000;

% normal distribution
k1 = normrnd(0,0.02,[len,1]);
k2 = normrnd(0,0.01,[len,1]);

% possion distribution
%lamda1 = 0.05; lamda2 = 0.04;
%k1 = poissrnd(lamda1, [len,1]);
%k2 = poissrnd(lamda2, [len,1]);

% uniform distribution
%k1 = unifrnd(-0.01,0.03,[len,1]);
%k2 = unifrnd(0,0.01,[len,1]);

% gamma distribution
%k1 = gamrnd(1,2,[len,1]);
%k2 = gamrnd(2,2,[len,1]);

phi1 = zeros(len,51);
phi2 = zeros(len,51);
x1 = zeros(len,51);
x2 = zeros(len,51);
x1(1,:) = repmat(rand(1),[51,1]);
x2(1,:) = repmat(rand(1),[51,1]);

for j = 1:51
    
    e = (j-1)*0.02;
    
    for i = 1:len-1
        
        %         x1(i+1,j) = (1-e)*a*x1(i,j)*(1-x1(i,j))+e*a*x2(i,j)*(1-x2(i,j));
        %         x2(i+1,j) = a*x2(i,j)*(1-x2(i,j));
        x1(i+1,j) = (1-e)*(1-a*x1(i,j)^2)+e*((1-a*x2(i,j)^2).^2);%theta1(i,j);
        x2(i+1,j) = (1-a*x2(i,j)^2).^2;%theta2(i,j);
        %         x1(i+1,j) = (1-e)*(1-a*x1(i,j)^2)+e*((1-a*x2(i,j)^2).^2+...
        %            b*(((sin(jd(i)'*2*pi./365)).^2)*((cos(jd(i)'*2*pi./365)).^2)));
        %         x2(i+1,j) = (1-a*x2(i,j)^2).^2+c*((sin(jd(i)'*2*pi./365)).^2)*((cos(jd(i)'*2*pi./365)).^2);
        %
    end
    
    x1 = k1+x1;
    x2 = k2+x2;
    
    % x1跟x1,x2均有关系，x2只与x2有关系
    phi1(:,j) = exp(-(x1(:,j)-mean(x1(:,j))).^2./(2*(var(x1(:,j)).^2)));
    phi2(:,j) = exp(-(x2(:,j)-mean(x2(:,j))).^2./(2*(var(x2(:,j)).^2)));
    
end

% figure
% % 
% plot(x1(:,20),'LineWidth',1)
% hold on
% plot(x2(:,20),'LineWidth',1)
% legend('y1','y2')
% grid on
% axis([1 1000 -2 2])

%% deperiod



%period = sin(jd(1:end-5)'*2*pi./365);





%% Linear
% 
% p = 6;
% for i = 1:numel(x1(1,:))
% 
% x10 = x1(:,i);x11 = [x10(1:end-p+1),x10(2:end-p+2),x10(3:end-p+3),x10(4:end-p+4),x10(5:end-p+5)];x12 = x10(p:end);
% x20 = x2(:,i);x21 = [x20(1:end-p+1),x20(2:end-p+2),x20(3:end-p+3),x20(4:end-p+4),x20(5:end-p+5)];x22 = x20(p:end);
% 
% x12 = stdd(x12);
% x22 = stdd(x22);
% 
% for j = 1:5
%     x11(:,j) = stdd(x11(:,j));
%     x21(:,j) = stdd(x21(:,j));
% end
% 
% b = regress(x12,x11);eps0 = var(x12-x11*b);b = regress(x12,[x11,x21]);eps1 = var(x12-[x11,x21]*b);
% G_21(i) = log(eps0/eps1);
% 
% b = regress(x22,x21);eps0 = var(x22-x21*b);b = regress(x22,[x21,x11]);eps1 = var(x22-[x21,x11]*b);
% G_12(i) = log(eps0/eps1);
% 
% end
% 
% figure
% subplot(3,2,1)
% plot(G_21','LineWidth',3);hold on;plot(G_12','LineWidth',3);grid on;legend('y2->y1','y1->y2');
% axis([1 51 0 1]); 
% set(gca,'XTick',1:5:51)
% set(gca,'XTickLabel',0:0.1:1)
% ylabel('Granger Cause','FontSize',13); xlabel('e','FontSize',13) ; 
% title('Linear Model','FontSize',15,'Color','b')

%% Kernel fit

p = 6;

for i = 1:numel(x1(1,:))
    
phi10 = phi1(:,i); phi20 = phi2(:,i);
    
x10 = x1(:,i);x11 = [phi10(1:end-p+1),phi10(2:end-p+2),phi10(3:end-p+3),phi10(4:end-p+4),phi10(5:end-p+5)];x12 = x10(p:end);
x20 = x2(:,i);x21 = [phi20(1:end-p+1),phi20(2:end-p+2),phi20(3:end-p+3),phi20(4:end-p+4),phi20(5:end-p+5)];x22 = x20(p:end);

x12 = stdd(x12);
x22 = stdd(x22);

for j = 1:5
    x11(:,j) = stdd(x11(:,j));
    x21(:,j) = stdd(x21(:,j));
end

b = regress(x12,x11);eps0 = var(x12-x11*b);b = regress(x12,[x11,x21]);eps1 = var(x12-[x11,x21]*b);
G_21(i) = log(eps0/eps1);

b = regress(x22,x21);eps0 = var(x22-x21*b);b = regress(x22,[x21,x11]);eps1 = var(x22-[x21,x11]*b);
G_12(i) = log(eps0/eps1);

end

Kernel_1 = G_21;
% subplot(2,2,2)
% plot(G_21','LineWidth',3);hold on;plot(G_12','LineWidth',3);grid on;
% legend({'y2->y1','y1->y2'},'Location','northwest');
% legend('boxoff')
% axis([1 51 0 1]); 
% set(gca,'XTick',1:5:51)
% set(gca,'XTickLabel',0:0.1:1)
% ylabel('Granger Cause','FontSize',13); xlabel('e','FontSize',13) ;
% title('Gaussian Kernal Model','FontSize',15,'Color','b')

%% BP

p = 6;
for i = 1:length(x1(1,:))
    
x10 = x1(:,i);x11 = [x10(1:end-p+1),x10(2:end-p+2),x10(3:end-p+3),x10(4:end-p+4),x10(5:end-p+5)];x12 = x10(p:end);
x20 = x2(:,i);x21 = [x20(1:end-p+1),x20(2:end-p+2),x20(3:end-p+3),x20(4:end-p+4),x20(5:end-p+5)];x22 = x20(p:end);

x12 = stdd(x12);
x22 = stdd(x22);

for j = 1:5
    x11(:,j) = stdd(x11(:,j));
    x21(:,j) = stdd(x21(:,j));
end

net = feedforwardnet(10);
net.trainParam.showWindow = 0;
net = train(net,x11',x12');eps0 = var(x12'-sim(net,x11'));
net = feedforwardnet(10);
net.trainParam.showWindow = 0;
net = train(net,[x11,x21]',x12');eps1 = var(x12'-sim(net,[x11,x21]'));
G_21(i) = log(eps0/eps1);

net = feedforwardnet(10);
net.trainParam.showWindow = 0;
net = train(net,x21',x22');eps0 = var(x22'-sim(net,x21'));
net = feedforwardnet(10);
net.trainParam.showWindow = 0;
net = train(net,[x21,x11]',x22');eps1 = var(x22'-sim(net,[x21,x11]'));
G_12(i) = log(eps0/eps1);

end

BP_1 = G_21;

% subplot(2,2,3)
% plot(G_21','LineWidth',3);hold on;plot(G_12','LineWidth',3);grid on;
% legend({'y2->y1','y1->y2'},'Location','northwest');
% legend('boxoff')
% axis([1 51 0 1]); 
% set(gca,'XTick',1:5:51)
% set(gca,'XTickLabel',0:0.1:1)
% ylabel('Granger Cause','FontSize',13); xlabel('e','FontSize',13) ;
% title('Neural Network Model','FontSize',15,'Color','b')

%% RF

p = 6;

for i = 1:length(x1(1,:))
    
x10 = x1(:,i);x11 = [x10(1:end-p+1),x10(2:end-p+2),x10(3:end-p+3),x10(4:end-p+4),x10(5:end-p+5)];x12 = x10(p:end);
x20 = x2(:,i);x21 = [x20(1:end-p+1),x20(2:end-p+2),x20(3:end-p+3),x20(4:end-p+4),x20(5:end-p+5)];x22 = x20(p:end);

x12 = stdd(x12);
x22 = stdd(x22);

for j = 1:5
    x11(:,j) = stdd(x11(:,j));
    x21(:,j) = stdd(x21(:,j));
end

% B = TreeBagger(60,[x11,period],x12,'method','regression','Surrogate','on',...
%         'OOBPredictorImportance','on','PredictorSelection','curvature');  
% eps0 = var(x12-predict(B,[x11,period]));
% 
% B = TreeBagger(60,[x11,x21,period],x12,'method','regression','Surrogate','on',...
%         'OOBPredictorImportance','on','PredictorSelection','curvature');  
% eps1 = var(x12-predict(B,[x11,x21,period]));
% G_21(i) = log(eps0/eps1);
% 
% B = TreeBagger(60,[x21,period],x22,'method','regression','Surrogate','on',...
%         'OOBPredictorImportance','on','PredictorSelection','curvature');  
% eps0 = var(x22-predict(B,[x21,period]));
% 
% B = TreeBagger(60,[x21,x11,period],x22,'method','regression','Surrogate','on',...
%         'OOBPredictorImportance','on','PredictorSelection','curvature');  
% eps1 = var(x22-predict(B,[x21,x11,period]));
% G_12(i) = log(eps0/eps1);

B = TreeBagger(60,[x11],x12,'method','regression','Surrogate','on',...
        'OOBPredictorImportance','on','PredictorSelection','curvature');  
eps0 = var(x12-predict(B,[x11]));

B = TreeBagger(60,[x11,x21],x12,'method','regression','Surrogate','on',...
        'OOBPredictorImportance','on','PredictorSelection','curvature');  
eps1 = var(x12-predict(B,[x11,x21]));
G_21(i) = log(eps0/eps1);

B = TreeBagger(60,[x21],x22,'method','regression','Surrogate','on',...
        'OOBPredictorImportance','on','PredictorSelection','curvature');  
eps0 = var(x22-predict(B,[x21]));

B = TreeBagger(60,[x21,x11],x22,'method','regression','Surrogate','on',...
        'OOBPredictorImportance','on','PredictorSelection','curvature');  
eps1 = var(x22-predict(B,[x21,x11]));
G_12(i) = log(eps0/eps1);

end
RF_1 = G_21;
% subplot(2,2,4)
% plot(G_21','LineWidth',3);hold on;plot(G_12','LineWidth',3);grid on;
% legend({'y2->y1','y1->y2'},'Location','northwest');
% legend('boxoff')
% axis([1 51 0 1]); 
% set(gca,'XTick',1:5:51)
% set(gca,'XTickLabel',0:0.1:1)
% ylabel('Granger Cause','FontSize',13); xlabel('e','FontSize',13) ;
% title('Random Forest Model','FontSize',15,'Color','b')

%% GLM

p = 6;

for i = 1:numel(x1(1,:))
    
x10 = x1(:,i);x11 = [x10(1:end-p+1),x10(2:end-p+2),x10(3:end-p+3),x10(4:end-p+4),x10(5:end-p+5)];x12 = x10(p:end);
x20 = x2(:,i);x21 = [x20(1:end-p+1),x20(2:end-p+2),x20(3:end-p+3),x20(4:end-p+4),x20(5:end-p+5)];x22 = x20(p:end);

x12 = stdd(x12);
x22 = stdd(x22);

for j = 1:5
    x11(:,j) = stdd(x11(:,j));
    x21(:,j) = stdd(x21(:,j));
end

one = ones(numel(x12),1);
% 
% [b,dev2,stats,iter,iter_flag,rank_flag,scale_flag]=...
%             glmfit_itersave([one,x11,period],x12,'binomial','link',...
%             'probit','estdisp','on','constant','off');         
% eps0 = var(x12-glmval(b,[one,x11,period],'probit','constant','off'));
% 
% [b,dev2,stats,iter,iter_flag,rank_flag,scale_flag]=...
%             glmfit_itersave([one,x11,x21,period],x12,'binomial','link',...
%             'probit','estdisp','on','constant','off');        
% eps1 = var(x12-glmval(b,[one,x11,x21,period],'probit','constant','off'));
% G_21(i) = log(eps0/eps1);
% 
% [b,dev2,stats,iter,iter_flag,rank_flag,scale_flag]=...
%             glmfit_itersave([one,x21,period],x22,'binomial','link',...
%             'probit','estdisp','on','constant','off');        
% eps0 = var(x22-glmval(b,[one,x21,period],'probit','constant','off'));
% 
% [b,dev2,stats,iter,iter_flag,rank_flag,scale_flag]=...
%             glmfit_itersave([one,x21,x11,period],x22,'binomial','link',...
%             'probit','estdisp','on','constant','off');        
% eps1 = var(x22-glmval(b,[one,x21,x11,period],'probit','constant','off'));
% G_12(i) = log(eps0/eps1);

[b,dev2,stats,iter,iter_flag,rank_flag,scale_flag]=...
            glmfit_itersave([one,x11],x12,'binomial','link',...
            'probit','estdisp','on','constant','off');         
eps0 = var(x12-glmval(b,[one,x11],'probit','constant','off'));

[b,dev2,stats,iter,iter_flag,rank_flag,scale_flag]=...
            glmfit_itersave([one,x11,x21],x12,'binomial','link',...
            'probit','estdisp','on','constant','off');        
eps1 = var(x12-glmval(b,[one,x11,x21],'probit','constant','off'));
G_21(i) = log(eps0/eps1);

[b,dev2,stats,iter,iter_flag,rank_flag,scale_flag]=...
            glmfit_itersave([one,x21],x22,'binomial','link',...
            'probit','estdisp','on','constant','off');        
eps0 = var(x22-glmval(b,[one,x21],'probit','constant','off'));

[b,dev2,stats,iter,iter_flag,rank_flag,scale_flag]=...
            glmfit_itersave([one,x21,x11],x22,'binomial','link',...
            'probit','estdisp','on','constant','off');        
eps1 = var(x22-glmval(b,[one,x21,x11],'probit','constant','off'));
G_12(i) = log(eps0/eps1);

end
GLM_1 = G_21;
% subplot(2,2,1)
% plot(G_21','LineWidth',3);hold on;plot(G_12','LineWidth',3);grid on;
% legend({'y2->y1','y1->y2'},'Location','northwest');
% legend('boxoff')
% axis([1 51 0 1]); 
% set(gca,'XTick',1:5:51)
% set(gca,'XTickLabel',0:0.1:1)
% ylabel('Granger Cause','FontSize',13); xlabel('e','FontSize',13) ;
% title('Generalized Linear Model','FontSize',15,'Color','b')

figure
plot(Kernel_1','LineWidth',2,'LineStyle','--');
hold on
plot(GLM_1','LineWidth',2,'LineStyle','--');
plot(BP_1','LineWidth',3);
plot(RF_1','LineWidth',3);
grid on

legend({'GK','GLM','BP','RF'},'Location','northwest');
legend('boxoff')

axis([1 51 0 1]); 
set(gca,'XTick',1:5:51)
set(gca,'XTickLabel',0:0.1:1)
ylabel('Nonlinear Granger Cause Value','FontSize',15);
xlabel('Nonlinear Impact Index','FontSize',15)
%xlabel('$${\alpha}$$','Interpreter','latex','FontSize',15,'Color','b')
%title('y2->y1','FontSize',15,'Color','b')

function x_standard = stdd(x)
    x_standard = (x-min(x))./(max(x)-min(x));
end

function plot_dynamic_bias_signal_hanks

% Code to fit the history-dependent drift diffusion models as described in
% Urai AE, de Gee JW, Tsetsos K, Donner TH (2019) Choice history biases subsequent evidence accumulation. eLife, in press.
%
% MIT License
% Copyright (c) Anne Urai, 2019
% anne.urai@gmail.com

% colors = [77,175,74; 55,126,184] ./ 256; % green blue
cols2 = cbrewer('qual', 'Dark2', 8);
close all;
global mypath colors
kostisPath = sprintf('%s/Anke_MEG_transition/KostisFits', mypath);

% nicer colors
subplot(3,2,1);

 % =================== %
%% start plot
 % =================== %

durs=240;

I1=.09*26; % assuming 9% coherence, this is relevant for OU only in which the state

P1=0;
P2=0;
dt=0.01/2;
% lk=4.22;
% % incr=0.1602;
% % a=[];
% % load('OUDallmodels.mat')
% % for t=1:durs;
% %
% %     P1=P1+dt.*(lk+incr)*P1+I1*dt; % accumulator with bias
% %
% %     P2=P2+dt.*lk*P2+I1*dt; % accumulator without bias
% %     a=[a;P1 P2];
% % end
% % %%
% %
% % hold on;
% % plot([1:240]*dt,(min(a(:,1)-a(:,2)-11.5,0))./11.5,'LineWidth',3);
ylabel({'Effective bias' '(fraction of bound)'});
xlabel('Time (s)');
hold on;
%set(gca,'YTicklabel',[0:0.1:1]); hold on;

 % =================== %
%%DDM Collapsing
 % =================== %

load(sprintf('%s/DynDDMCol_allmodels.mat', kostisPath));
thr=mean(params4(:,1));
bias=mean(abs(params4(:,6)));
dv=mean(params4(:,4));
timest=240;
dt=0.01/2;

tv=1:timest;
y=thr-thr.*(tv./(tv+dv));
TH1=max(y,thr/2);
TH2=min(thr/2,-y+thr);

P1=0;
P1 = mean(abs(params4(:,7)));% should be replaced by P1=stp, where stp=mean(abs(PARAMETER))
P2=0;
I1=.01*33;
I2=.01*33;
a=[];
for t=1:durs;
    
    P1=P1+dt*bias+I1*dt;
    
    P2=P2+I2*dt;
    a=[a;P1 P2];
end
plot([1:240]*dt,min((a(:,1)-a(:,2))./(TH1'-thr/2),1),'color', colors(3, :), 'linewidth', 2);

 % =================== %
%%OU Collapsing
 % =================== %

for m=1:2;
    load(sprintf('%s/OUcollapse_allmodels.mat', kostisPath));
    clear a;
    thr=mean(params2(:,1));
    bias=mean(abs(params2(:,end-1)));
    dv=mean(params2(:,4));
    timest=240;
    dt=0.01/2;
    
    tv=1:timest;
    
    
    P1=0;
    P2=0;
    I1=.09*50;
    I2=.09*50;
    a=[];
    if m==1
        indx=find(params2(:,end)<0);
        
    else
        indx=find(params2(:,end)>0);
    end
    lk=mean(params2(indx,end));
    thr=mean(params2(indx,1));
    bias=mean(abs(params2(indx,end-1)));
    dv=mean(params2(indx,4));
    y=thr-thr.*(tv./(tv+dv));
    TH1=max(y,thr/2);
    TH2=min(thr/2,-y+thr);
    for t=1:durs;
        
        P1=P1+dt*lk*P1+I1*dt+bias*dt;
        
        P2=P2+dt*lk*P2+I2*dt;
        a=[a;P1 P2];
    end
    %plot([1:240]*dt,min((a(:,1)-a(:,2))./(TH1'-thr/2+eps),1),'LineWidth',3);
    
    biasTV(:,m)=min((a(:,1)-a(:,2))./(TH1'-thr/2),1); %here I will average the signal between the 2 OU domains since they look very similar
    % it makes less sense to just average the lambda parameter without taking
    % into account the sign, i.e. averaging can lead to cancelation
end
plot([1:240]*dt,mean(biasTV'),'LineWidth',2, 'color', cols2(3, :));

% legend({'Collapsing DDM: v_{bias}','Collapsing leaky accumulator: input bias'}, ...
%     'Location', 'NorthWest');
% legend boxoff;
axis tight;
set(gca, 'ytick', [0:1]);
offsetAxes;

%subplot(423); axis off;
tightfig;
print(gcf, '-dpdf', '~/Data/serialHDDM/effective_bias_timecourse.pdf');

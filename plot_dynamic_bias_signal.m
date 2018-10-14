function plot_dynamic_bias_signal
% from a range of fitted dynamic decision models, get the average
% bias across the group and plot the effective bias signal (DV to bound)
% over the timecourse of the trial
% by Konstantinos Tsetsos
% edits by Anne Urai

clear all;
close all;
colors = [77,175,74; 55,126,184] ./ 256; % green blue
cols2 = cbrewer('qual', 'Dark2', 8);

kostisPath = '~/Data/HDDM/Anke_MEG_transition/KostisFits';

% nicer colors
% set(groot, 'defaultaxescolororder', linspecer(5));
subplot(4,3,[1 2]);

%% OU

durs=450;
I1=.09*26; % assuming 9% coherence, this is relevant for OU only in which the state

P1=0;
P2=0;
dt=0.01/2;
lk=4.22;
incr=0.1602;
a=[];
load(sprintf('%s/OUDallmodels.mat', kostisPath));
for t=1:durs;
    
    P1=P1+dt.*(lk+incr)*P1+I1*dt; % accumulator with bias
    P2=P2+dt.*lk*P2+I1*dt; % accumulator without bias
    a=[a;P1 P2];
end
%

hold on;
plot([1:450]*dt,(min(a(:,1)-a(:,2)-11.5,0))./11.5,'LineWidth',1, 'color', cols2(4, :));

%% DDM Collapsing
load(sprintf('%s/DynDDMCol_allmodels.mat', kostisPath));
thr=mean(params2(:,1));
bias=mean(abs(params2(:,end)));
dv=mean(params2(:,4));
timest=450;
dt=0.01/2;

tv=1:timest;
y=thr-thr.*(tv./(tv+dv));
TH1=max(y,thr/2);
TH2=min(thr/2,-y+thr);

P1=0;
P2=0;
I1=.01*33;
I2=.01*33;
a=[];
for t=1:durs;
    
    P1=P1+dt*bias+I1*dt;
    
    P2=P2+I2*dt;
    a=[a;P1 P2];
end
plot([1:450]*dt,min(a(:,1)-a(:,2)-TH1'+thr/2,0)./(thr/2),'LineWidth',1, 'color', colors(2, :), 'linestyle', '--');


%% OU Collapsing
for m=1:2;
    load(sprintf('%s/OUcollapse_allmodels.mat', kostisPath));
    clear a;
    thr=mean(params2(:,1));
    bias=mean(abs(params2(:,end-1)));
    dv=mean(params2(:,4));
    timest=450;
    dt=0.01/2;
    
    tv=1:timest;
    y=thr-thr.*(tv./(tv+dv));
    TH1=max(y,thr/2);
    TH2=min(thr/2,-y+thr);
    
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
    for t=1:durs;
        
        P1=P1+dt*lk*P1+I1*dt+bias*dt;
        
        P2=P2+dt*lk*P2+I2*dt;
        a=[a;P1 P2];
    end
    if m == 1,
        plot([1:450]*dt,min(a(:,1)-a(:,2)-TH1'+thr/2,0)./(thr/2),'LineWidth',2, 'color', cols2(3, :), 'linestyle', ':');
    elseif m == 2,
        plot([1:450]*dt,min(a(:,1)-a(:,2)-TH1'+thr/2,0)./(thr/2),'LineWidth',2, 'color', cols2(3, :), 'linestyle', '-.');
    end
end


%% DDM Drift-criterion
load(sprintf('%s/allmodels.mat', kostisPath));
clear a;
P1=0;
P2=0;
a=[];
thr=mean(params(:,1));
scale=mean(params(:,2));
I1=.09*scale;
I2=.09*scale;
bias=mean(abs(params(:,end)));
for t=1:durs;
    
    P1=P1+I1*dt+bias*dt;
    
    P2=P2+I2*dt;
    a=[a;P1 P2];
end
plot([1:450]*dt,min(a(:,1)-a(:,2)-thr+thr/2,0)./(thr/2),'LineWidth',1, 'color', colors(2, :));

%% LAYOUT
l = legend({'O-U: \lambda bias','Collapsing DDM: v_{bias}', ...
    'Collapsing O-U (\lambda<0): input bias','Collapsing O-U (\lambda>0): input bias','Standard DDM: v_{bias}'}, ...
    'Location', 'EastOutside', 'Box', 'off');
l.Box = 'off';

axis tight; set(gca, 'xtick', 0:0.5:2);
ylabel({'Effective bias' '(fraction of bound)'});
xlabel('Time (s)');
offsetAxes;
set(gca, 'ytick', [-1:0.5:0], 'YTicklabel',[0:0.5:1]);
xlim([-0.05 2]); 
box off;

subplot(433); axis off;
tightfig;
print(gcf, '-dpdf', '~/Data/serialHDDM/effective_bias_timecourse.pdf');

end
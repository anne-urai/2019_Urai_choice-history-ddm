function [RT,ch]= diffusion_custom2(trls,cor,pars)

  thr=pars(1)*10; %bound
  scale=pars(2)*60;%scaling parameter
  T0=pars(3)*0.5;%non-decision time
  dv=pars(4)*3;%drift-rate variability
  noise_std=1;
  bsp=pars(5)*thr*0.75;%starting point variability
  
  inp=scale*trls; %input to the diffusion model across trials
  iters=size(trls,1);
  timest=300;%maximum number of simulation time-steps
  dt=0.01/2;

  dvariance=randn(iters,1)*dv;
  I=randn(iters,timest)*noise_std*sqrt(dt)+(inp)*dt+dvariance*sqrt(dt); %create model input time-series
  P=ones(iters,1)*thr/2+(-1+2*rand(iters,1))*bsp;%initialise decision variable
  RT=[];
  ch=[];
  corct=[];
  t=1;
  while (length(P)>0 & t<timest) % run the model up until all trials have been terminated, either diffusion crossed bound or reached time-out
      P=P+I(:,t);% update decision variable
      indx=find(P>thr);% identify trials where upper bound was crossed at time t
      indx1=find(P<0); % identify trials where the lower bound was crossed at time t
      P([indx;indx1])=[];% clear the dv vector from trials that terminated
      I([indx;indx1],:)=[]; % same for the input vector
      corct=[corct cor(indx)' cor(indx1)']; % store the stimulus identities associated with choices

      cor([indx;indx1])=[]; % clear stimulus identity vector
      ch=[ch ones(1,length(indx))]; %store the choices for 'right'
      ch=[ch zeros(1,length(indx1))];%store the chocies for 'left'
      RT=[RT t*ones(1,length(indx)+length(indx1))]; % store the RT's 
      t=t+1;
  end
  
  ch=[ch double(P'>thr/2)];% for those trials that didn't cross a threshold (reached time-out) assign choice based on whether the dv was closer to the upper or lower bound.
  RT=[RT ones(1,length(P))*timest];% register maximum simulation time as the RT's of the "time-out" trials;
  corct=[corct cor']; %update stimulus identity structure
  
  %%here transform choices into correctness, with 1--> correct and 0-->incorrect
ch2=ch;
ch2((ch==1 & corct==1) | (ch==0 & corct==-1))=1;
ch2((ch==0 & corct==1) | (ch==1 & corct==-1))=0;

ch=ch2;
%%
RT=((RT)*dt)+T0; %add non-decision time to RT's
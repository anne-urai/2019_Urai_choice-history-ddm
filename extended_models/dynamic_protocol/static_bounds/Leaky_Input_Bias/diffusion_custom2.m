function [RT,ch]= diffusion_custom(trls,cor,stim,pars)



  thr=pars(1)*60;
  scale=pars(2)*100;
  T0=pars(3)*1.5;
  dv=pars(4)*15; % leak parameter
  noise_std=1;
  bsp=pars(5)*thr*0.75;
  spbias=pars(6)*5;
  
  inp=scale*trls;
  iters=size(trls,1);
  timest=450;
  dt=0.01/2;

  I=[(inp)*dt+randn(iters,150)*noise_std*sqrt(dt)+(spbias.*stim).*dt randn(iters,300)*noise_std*sqrt(dt)+spbias.*stim.*dt]; %randn(iters,timest)*noise_std*sqrt(dt)];

  
P=ones(iters,1)*thr/2+(-1+2*rand(iters,1))*bsp;
 swtch=[zeros(iters,450)] ;
 swtch(:,150)=1;
C=0*P; 
  RT=[];
  ch=[];
  corct=[];
  t=1;
  while (length(P)>0 & t<timest)
      P=P+P*dv*dt+I(:,t)+C;
      C=C+swtch(:,t).*(P/t);
      indx=find(P>thr);
      indx1=find(P<0);
      P([indx;indx1])=[];
      C([indx;indx1])=[];
      I([indx;indx1],:)=[];
      corct=[corct cor(indx)' cor(indx1)'];
      swtch([indx;indx1],:)=[];
      cor([indx;indx1])=[];
      ch=[ch ones(1,length(indx))];
      ch=[ch zeros(1,length(indx1))];
      RT=[RT t*ones(1,length(indx)+length(indx1))];
      t=t+1;
  end
  
  ch=[ch double(P'>thr/2)];
  RT=[RT ones(1,length(P))*timest];
  corct=[corct cor'];
ch2=ch;
ch2((ch==1 & corct==1) | (ch==0 & corct==-1))=1;
ch2((ch==0 & corct==1) | (ch==1 & corct==-1))=0;

ch=ch2;
RT=((RT)*dt)+T0;

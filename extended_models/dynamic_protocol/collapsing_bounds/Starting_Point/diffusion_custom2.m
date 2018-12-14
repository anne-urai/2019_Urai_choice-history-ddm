function [RT,ch]= diffusion_custom(trls,cor,stim,pars)


  thr=pars(1)*60;
  scale=pars(2)*100;
  T0=pars(3)*1.5;
  dv=pars(4)*2000;
  noise_std=1;
  spbias=pars(6)*1;
  bsp=pars(5)*thr*0.75;

  inp=scale*trls;
  iters=size(trls,1);
  timest=450;
  dt=0.01/2;
 tv=1:timest;
  y=thr-thr.*(tv./(tv+dv));
  TH1=max(y,thr/2);
  TH2=min(thr/2,-y+thr);
  I=[(inp)*dt+randn(iters,150)*noise_std*sqrt(dt) randn(iters,300)*noise_std*sqrt(dt)];  
P=ones(iters,1)*thr/2+(spbias.*stim)+(-1+2*rand(iters,1))*bsp;
  
swtch=[zeros(iters,150) ones(iters,300)];
  RT=[];
  ch=[];
  corct=[];
  t=1;
  while (length(P)>0 & t<timest)
      P=P+I(:,t)+swtch(:,t).*(P/t);
      indx=find(P>TH1(t));
      indx1=find(P<TH2(t));
      P([indx;indx1])=[];
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
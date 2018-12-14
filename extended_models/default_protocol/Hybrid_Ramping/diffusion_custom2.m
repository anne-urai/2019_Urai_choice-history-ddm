function [RT,ch]= diffusion_custom(trls,cor,stim,pars)



  thr=pars(1)*10;
  scale=pars(2)*60;
  T0=pars(3)*.5;
  dv=pars(4)*3;
  noise_std=1;
  spbias2=pars(7)*2.5;
  bsp=pars(5)*thr*0.75;
  spbias=spbias2+sign(spbias2)*pars(6)*5;
  
timest=300;
  tt=[1:timest]./timest;
  inp=repmat(scale*(trls),1,timest)+(spbias2.*stim)+repmat((spbias-spbias2).*stim,1,timest).*repmat(tt,length(trls),1);
  
  iters=length(trls);
  dt=0.01/2;

  dvariance=randn(iters,1)*dv;
  I=randn(iters,timest)*noise_std*sqrt(dt)+(inp)*dt+repmat(dvariance,1,timest)*sqrt(dt);
  P=ones(iters,1)*thr/2+(-1+2*rand(iters,1))*bsp;%+(spbias.*stim);
  RT=[];
  ch=[];
  corct=[];
  t=1;
  while (length(P)>0 & t<timest)
      P=P+I(:,t);
      indx=find(P>thr);
      indx1=find(P<0);
      P([indx;indx1])=[];
      I([indx;indx1],:)=[];
      corct=[corct cor(indx)' cor(indx1)'];

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
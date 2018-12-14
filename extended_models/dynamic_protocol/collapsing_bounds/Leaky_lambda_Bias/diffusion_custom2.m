function [RT,ch]= diffusion_custom(trls,cor,stim,pars)



  thr=pars(1)*60;
  scale=pars(2)*100;
  T0=pars(3)*1.5;
  dv=pars(4)*2000;
  noise_std=1;
  spbias=pars(6)*5;
  bsp=pars(5)*thr*0.75;  
 lk=pars(7)*15; 
 inp=scale*(trls);%+(spbias.*stim);
  iters=size(trls,1);
  timest=450;
  dt=0.01/2;
  incr1=double(stim>0)*spbias;
  incr2=double(stim<0)*spbias;

 
 tv=1:timest;
  y=thr-thr.*(tv./(tv+dv));
  TH1=max(y,thr/2);
  TH2=min(thr/2,-y+thr);

 I=[(inp)*dt]; 
  temp=I;  
  I1=[double(temp>0).*temp+randn(iters,150)*noise_std*sqrt(dt)/sqrt(2) randn(iters,300)*noise_std*sqrt(dt)/sqrt(2)];
  I2=[double(temp<0).*(-temp)+randn(iters,150)*noise_std*sqrt(dt)/sqrt(2) randn(iters,300)*noise_std*sqrt(dt)/sqrt(2)];

P1=zeros(iters,1)+(-1+2*rand(iters,1))*bsp;
 P=0*P1;
  

 
RT=[];
  ch=[];
  corct=[];
  t=1;
swtch1=[zeros(iters,150) ones(iters,300)];
swtch2=[zeros(iters,150) ones(iters,300)];

  while (length(P)>0 & t<timest)
      P1=P1+P1.*(lk+incr1)*dt+I1(:,t)+swtch1(:,t).*(P1/t);
      P2=P2+P2.*(lk+incr2)*dt+I2(:,t)+swtch2(:,t).*(P2/t);
      P=P1-P2;

      indx=find(P>TH1(t));
      indx1=find(P<TH2(t));
      P([indx;indx1])=[];
      P1([indx;indx1])=[];
      P2([indx;indx1])=[];
      incr1([indx;indx1])=[];
      incr2([indx;indx1])=[];
      swtch1([indx;indx1],:)=[];
      swtch2([indx;indx1],:)=[];
      I1([indx;indx1],:)=[];
      I2([indx;indx1],:)=[];


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
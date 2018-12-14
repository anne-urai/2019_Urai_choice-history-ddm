function [RT,ch]= diffusion_custom(trls,cor,stim,pars)

  thr=pars(1)*60;
  scale=pars(2)*100;
  T0=pars(3)*1.5;
  dv=pars(4)*15;
  noise_std=1;
  bsp=pars(5)*thr*0.75;
  sbias=pars(6)*5;

  
  incr1=double(stim>0)*sbias;
  incr2=double(stim<0)*sbias;

  inp=scale*trls;
  iters=size(trls,1);
  timest=450;
  dt=0.01/2;
  
  I=[(inp)*dt]; 
  temp=I;  
  I1=[double(temp>0).*temp+randn(iters,150)*noise_std*sqrt(dt)/sqrt(2) randn(iters,300)*noise_std*sqrt(dt)/sqrt(2)];
  I2=[double(temp<0).*(-temp)+randn(iters,150)*noise_std*sqrt(dt)/sqrt(2) randn(iters,300)*noise_std*sqrt(dt)/sqrt(2)];

P1=zeros(iters,1)+(-1+2*rand(iters,1))*bsp;
P2=zeros(iters,1);
 P=0*P1;
  RT=[];
  ch=[];
  corct=[];
  t=1;
swtch1=[zeros(iters,450) ];
swtch2=[zeros(iters,450) ];
swtch1(:,150)=1;
swtch2(:,150)=1;
C1=0*P1;
C2=0*P2;
  while (length(P)>0 & t<timest)
      P1=P1+P1.*(dv+incr1)*dt+I1(:,t)+C1;
      P2=P2+P2.*(dv+incr2)*dt+I2(:,t)+C2;
      
      P=P1-P2;
      C1=C1+swtch1(:,t).*(P1/t);
      C2=C2+swtch2(:,t).*(P2/t);
      indx=find(P>thr/2);
      indx1=find(P<-thr/2);
      P([indx;indx1])=[];
      P1([indx;indx1])=[];
      P2([indx;indx1])=[];
       C1([indx;indx1])=[];
      C2([indx;indx1])=[];
      
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
  
  ch=[ch double(P'>0)];
  RT=[RT ones(1,length(P))*timest];
  corct=[corct cor'];
ch2=ch;
ch2((ch==1 & corct==1) | (ch==0 & corct==-1))=1;
ch2((ch==0 & corct==1) | (ch==1 & corct==-1))=0;

ch=ch2;
RT=((RT)*dt)+T0;
v=600e3;r=1e5;L=1e9;
L=[1e5 1e6 1e7 1e8 1e9 1e10 1e11];
f=logspace(-4,1,20); for L=[1e5 1e6 1e7 1e8 1e9 1e10 1e11]; P(:,:,i)=cbspec(f,[0 1 3 5 10 20 30 40 50 60 70 80 90],2000,1,L,0,v,0,1,[-12 0]); end



load ruhe

[f2,Pl4e3]=psdtf(LFP_ruhe(4,:),1/sf,0,EMG_ruhe(3,:));

for i=1:up-3-low;
    x=LFP_ruhe(4,(i-1)*sf+1:(i+3)*sf);
    y=EMG_ruhe(3,(i-1)*sf+1:(i+3)*sf); 
    [f,pxy(i,:)]=psdtf(x,1/sf,0,y);
    [f,py(i,:)]=psdtf(y,1/sf);
    [f,px(i,:)]=psdtf(x,1/sf);
end

ca=abs(mean(pxy)).^2./mean(px)./mean(py);
cr=abs(mean(real(pxy)).^2)./mean(px)./mean(py);
ci=abs(mean(imag(pxy)).^2)./mean(px)./mean(py);



for i=1:up-3-low;
    x=LFP_ruhe(4,(i-1)*sf+1:(i+3)*sf);
    y=EMG_ruhe(3,(i-1)*sf+1:(i+3)*sf); 
    [f,pxy(i,:)]=psdtf(x,1/sf,0,y,1);
    [f,py(i,:)]=psdtf(y,1/sf,0,0,1);
    [f,px(i,:)]=psdtf(x,1/sf,0,0,1);
end

cah=abs(mean(pxy)).^2./mean(px)./mean(py);
crh=abs(mean(real(pxy)).^2)./mean(px)./mean(py);
cih=abs(mean(imag(pxy)).^2)./mean(px)./mean(py);



for i=1:7;
    x=LFP_ruhe(4,(i-1)*3*sf+1:i*3*sf);
    y=EMG_ruhe(3,(i-1)*3*sf+1:i*3*sf); 
    [f,pxy(i,:)]=psdtf(x,1/sf,0,y);
    [f,py(i,:)]=psdtf(y,1/sf,0,0);
    [f,px(i,:)]=psdtf(x,1/sf,0,0);
end

cah=abs(mean(pxy)).^2./mean(px)./mean(py);
crh=abs(mean(real(pxy)).^2)./mean(px)./mean(py);
cih=abs(mean(imag(pxy)).^2)./mean(px)./mean(py);
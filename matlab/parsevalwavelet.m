%% Import Data:
A=importdata('60m_04301_010.dat',' ',13);
t=A.data(:,1);
Bx=A.data(:,2);
Et=var(Bx);

%% Parameters
dj=0;
dt=0.14;
pad=0;
s0=2*dt*(param+sqrt(2+6^2))/(4*pi);
n=length(t);

for u=1:50
    dj(u)=0.02*u;
    J1(u)=round(log2(n*dt/s0)/dj(u));
    [wave,period,scale,coi]=wavelet(Bx,dt,pad,dj(u),s0,J1(u));
    k=size(wave);

    C=repmat(coi,k(1),1);
    P=repmat(period',1,k(2));
    coi=C./P>1;
    
    Ew=0;
    for i=1:k(1)
        N=length(find(coi(i,:))); %number of ones
        if N~=0
            Ew=Ew+sum(abs(wave(i,:).*coi(i,:)).^2,2)'./scale(i)*dj(u)*dt/(0.776*N);
        end
    end
    a(u)=Ew/Et;
    
end

clear C P k Ew Et wave period scale

%% 
plot(dj,a,dj,b,dj,paddy);
xlabel('\delta_j');
ylabel('a');

print -depsc parseval_v2
function C = mycoh ( x, y )

[wavex,period,scale,coix] = wavelet(x,dt);
wavey = wavelet(y,dt);

pxy=wavex.*conj(wavey);
pyy=abs(wavey).^2;
pxx=abs(wavex).^2;
spxy=smooth_mat(pxy,dt,scale);
spxx=smooth_mat(pxx,dt,scale);
spyy=smooth_mat(pyy,dt,scale);
C2=abs(spxy2).^2./(spxx.*spyy);
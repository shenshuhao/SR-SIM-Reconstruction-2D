clc;
clear all;
close all;

w1 = 1000;
w2 = 1000;
wo1 = w1/2;
wo2 = w2/2;
x = linspace(0,w1-1,w1);
y = linspace(0,w2-1,w2);
[X,Y] = meshgrid(x,y);

object = imread('Beads-obj.tif');
object = double(object);
figure,imshow(object,[]);title('Object');
freq_obj=abs(fftshift(fft2(object)));
figure,imshow((freq_obj.^0.1),[]);

h = fspecial('gaussian',[w1,w2],5); 
I = conv2(object,h,'same');
figure,imshow(I,[]);title('Filter');
otf=(fftshift(fft2(h)));
figure,imshow(abs(otf.^0.1),[]);title('OTF');

m = 0.8;
f = 0.3;
thita0=  0;
thita1 = 2*pi/3;
thita2 = 1*pi/3;
phase1 = 0*pi/3;
phase2 = 2*pi/3;
phase3 = 4*pi/3;

mask1 = 1+m*sin(f*(X*sin(thita0)+Y*cos(thita0))+phase1);
% figure,imagesc(mask1);
mask2 = 1+m*sin(f*(X*sin(thita0)+Y*cos(thita0))+phase2);
% figure,imagesc(mask2);
mask3 = 1+m*sin(f*(X*sin(thita0)+Y*cos(thita0))+phase3);
% figure,imagesc(mask3);
mask4 = 1+m*sin(f*(X*sin(thita1)+Y*cos(thita1))+phase1);
% figure,imagesc(mask4);
mask5 = 1+m*sin(f*(X*sin(thita1)+Y*cos(thita1))+phase2);
% figure,imagesc(mask5);
mask6 = 1+m*sin(f*(X*sin(thita1)+Y*cos(thita1))+phase3);
% figure,imagesc(mask6);
mask7 = 1+m*sin(f*(X*sin(thita2)+Y*cos(thita2))+phase1);
% figure,imagesc(mask7);
mask8 = 1+m*sin(f*(X*sin(thita2)+Y*cos(thita2))+phase2);
% figure,imagesc(mask8);
mask9 = 1+m*sin(f*(X*sin(thita2)+Y*cos(thita2))+phase3);
% figure,imagesc(mask9);

I1 = object.*mask1;
I1 = conv2(I1,h,'same');
I2 = object.*mask2;
I2 = conv2(I2,h,'same');
I3 = object.*mask3;
I3 = conv2(I3,h,'same');
I4 = object.*mask4;
I4 = conv2(I4,h,'same');
I5 = object.*mask5;
I5 = conv2(I5,h,'same');
I6 = object.*mask6;
I6 = conv2(I6,h,'same');
I7 = object.*mask7;
I7 = conv2(I7,h,'same');
I8 = object.*mask8;
I8 = conv2(I8,h,'same');
I9 = object.*mask9;
I9 = conv2(I9,h,'same');

FS1aT=fftshift(fft2(I1));
FS2aT=fftshift(fft2(I2));
FS3aT=fftshift(fft2(I3));
FS4aT=fftshift(fft2(I4));
FS5aT=fftshift(fft2(I5));
FS6aT=fftshift(fft2(I6));
FS7aT=fftshift(fft2(I7));
FS8aT=fftshift(fft2(I8));
FS9aT=fftshift(fft2(I9));

ModFac=m/2;
M = [1 -ModFac*exp(-1i*phase1) -ModFac*exp(1i*phase1);
     1 -ModFac*exp(-1i*phase2) -ModFac*exp(1i*phase2);
     1 -ModFac*exp(-1i*phase3) -ModFac*exp(1i*phase3);];

Minv = inv(M);
FS1 = Minv(1,1)*FS1aT + Minv(1,2)*FS2aT + Minv(1,3)*FS3aT;
FS2 = Minv(2,1)*FS1aT + Minv(2,2)*FS2aT + Minv(2,3)*FS3aT;
FS3 = Minv(3,1)*FS1aT + Minv(3,2)*FS2aT + Minv(3,3)*FS3aT;
FS4 = Minv(1,1)*FS4aT + Minv(1,2)*FS5aT + Minv(1,3)*FS6aT;
FS5 = Minv(2,1)*FS4aT + Minv(2,2)*FS5aT + Minv(2,3)*FS6aT;
FS6 = Minv(3,1)*FS4aT + Minv(3,2)*FS5aT + Minv(3,3)*FS6aT;
FS7 = Minv(1,1)*FS7aT + Minv(1,2)*FS8aT + Minv(1,3)*FS9aT;
FS8 = Minv(2,1)*FS7aT + Minv(2,2)*FS8aT + Minv(2,3)*FS9aT;
FS9 = Minv(3,1)*FS7aT + Minv(3,2)*FS8aT + Minv(3,3)*FS9aT;

% figure,imagesc(abs(FS1.^(1/10)));
% figure,imagesc(abs(FS2.^(1/10)));
% figure,imagesc(abs(FS3.^(1/10)));
% figure,imagesc(abs(FS4.^(1/10)));
% figure,imagesc(abs(FS5.^(1/10)));
% figure,imagesc(abs(FS6.^(1/10)));
% figure,imagesc(abs(FS7.^(1/10)));
% figure,imagesc(abs(FS8.^(1/10)));
% figure,imagesc(abs(FS9.^(1/10)));

dist=f/pi/2/(x(2)-x(1))*(w1-1);
%  dist=64.5;
kA=[dist,0,dist*cos(pi/3),dist*sin(pi/3)];
Shift_F2 = fft2( (ifft2(FS2)).*exp( +1i.*2*pi*(kA(1)/w1.*(X-wo1) + kA(2)/w2.*(Y-wo2)) ));
Shift_F3 = fft2( (ifft2(FS3)).*exp( -1i.*2*pi*(kA(1)/w1.*(X-wo1) + kA(2)/w2.*(Y-wo2)) ));
Shift_F5 = fft2( (ifft2(FS5)).*exp( +1i.*2*pi*(kA(3)/w1.*(X-wo1) - kA(4)/w2.*(Y-wo2)) ));
Shift_F6 = fft2( (ifft2(FS6)).*exp( -1i.*2*pi*(kA(3)/w1.*(X-wo1) - kA(4)/w2.*(Y-wo2)) ));
Shift_F8 = fft2( (ifft2(FS8)).*exp( +1i.*2*pi*(kA(3)/w1.*(X-wo1) + kA(4)/w2.*(Y-wo2)) ));
Shift_F9 = fft2( (ifft2(FS9)).*exp( -1i.*2*pi*(kA(3)/w1.*(X-wo1) + kA(4)/w2.*(Y-wo2)) ));
% 
% figure,imagesc(abs(FS1.^(1/10)));
% figure,imagesc(abs(Shift_F2.^(1/10)));
% figure,imagesc(abs(Shift_F3.^(1/10)));
% figure,imagesc(abs(FS4.^(1/10)));
% figure,imagesc(abs(Shift_F5.^(1/10)));
% figure,imagesc(abs(Shift_F6.^(1/10)));
% figure,imagesc(abs(FS7.^(1/10)));
% figure,imagesc(abs(Shift_F8.^(1/10)));
% figure,imagesc(abs(Shift_F9.^(1/10)));

FSum1=FS1+Shift_F2+Shift_F3;
% figure,imagesc(abs(FSum1.^(1/10)));
FSum2=FS4+Shift_F5+Shift_F6;
% figure,imagesc(abs(FSum2.^(1/10)));
FSum3=FS7+Shift_F8+Shift_F9;
% figure,imagesc(abs(FSum3.^(1/10)));

otf2 = fft2( (ifft2(otf)).*exp( +1i.*2*pi*(kA(1)/w1.*(X-wo1) + kA(2)/w2.*(Y-wo2)) ));
otf3 = fft2( (ifft2(otf)).*exp( -1i.*2*pi*(kA(1)/w1.*(X-wo1) + kA(2)/w2.*(Y-wo2)) ));
otf5 = fft2( (ifft2(otf)).*exp( +1i.*2*pi*(kA(3)/w1.*(X-wo1) - kA(4)/w2.*(Y-wo2)) ));
otf6 = fft2( (ifft2(otf)).*exp( -1i.*2*pi*(kA(3)/w1.*(X-wo1) - kA(4)/w2.*(Y-wo2)) ));
otf8 = fft2( (ifft2(otf)).*exp( +1i.*2*pi*(kA(3)/w1.*(X-wo1) + kA(4)/w2.*(Y-wo2)) ));
otf9 = fft2( (ifft2(otf)).*exp( -1i.*2*pi*(kA(3)/w1.*(X-wo1) + kA(4)/w2.*(Y-wo2)) ));


sigma=0.5;
FS10=FS1./(abs(otf).^2+sigma);
Shift_F20=Shift_F2./(abs(otf2).^2+sigma);
Shift_F30=Shift_F3./(abs(otf3).^2+sigma);
FS40=FS4./(abs(otf).^2+sigma);
Shift_F50=Shift_F5./(abs(otf5).^2+sigma);
Shift_F60=Shift_F6./(abs(otf6).^2+sigma);
FS70=FS7./(abs(otf).^2+sigma);
Shift_F80=Shift_F8./(abs(otf8).^2+sigma);
Shift_F90=Shift_F9./(abs(otf9).^2+sigma);

FSum=FS10+Shift_F20+Shift_F30+FS40+Shift_F50+Shift_F60+FS70+Shift_F80+Shift_F90;
% figure,imagesc(abs(FSum.^(1/10)));



figure,imshow(abs(FSum.^(1/10)),[]);title('FT of SIM image');
cm=I1+I2+I3+I4+I5+I6+I7+I8+I9;
freqcm=fftshift(fft2(cm));
figure,imshow(abs(freqcm.^(1/10)),[]);title('FT of CM image');


DSim=abs(ifft2(fftshift(FSum)));
figure,imshow(DSim,[]);title('SIM image');
figure,imshow(cm,[]);title('CM image');
DSim=DSim/max(max(DSim))*65536;
DSim=uint16(DSim);
cm=cm/max(max(cm))*65536;
cm=uint16(cm);
imwrite(DSim,'Super_resolution_beads.tif');
imwrite(cm,'Normal_resolution_beads.tif');

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

filename='star-chart-bars';
object = imread([filename '.tif']);
object = double(object);
figure,imshow(object,[]);title('Object');
freq_obj=abs(fftshift(fft2(object)));
figure,imshow((freq_obj.^0.1),[]);

h = fspecial('gaussian',[w1,w2],5); 
I = conv2(object,h,'same');
figure,imshow(I,[]);title('Filter');
otf=(fftshift(fft2(h)));
figure,imshow(abs(otf.^0.1),[]);title('OTF');

m = 0.8;%Modulation depth
f = 0.7;%Pattern frequency
thita0=  0;
thita1 = 5*pi/6;%Orientation
thita2 = 1*pi/6;
phase1 = 0*pi/3;%Phase
phase2 = 2*pi/3;
phase3 = 4*pi/3;

grating1 = 1+m*cos(f*x+phase1);
grating2 = 1+m*cos(f*x+phase2);
grating3 = 1+m*cos(f*x+phase3);

for i = 1:w1
    mask1(i,:) = grating1;
    mask2(i,:) = grating2;
    mask3(i,:) = grating3;
end
mask4 = 1+m*cos(f*(X*sin(thita1)+Y*cos(thita1))+phase1);
% figure,imagesc(mask4);
mask5 = 1+m*cos(f*(X*sin(thita1)+Y*cos(thita1))+phase2);
% figure,imagesc(mask5);
mask6 = 1+m*cos(f*(X*sin(thita1)+Y*cos(thita1))+phase3);
% figure,imagesc(mask6);
mask7 = 1+m*cos(f*(X*sin(thita2)+Y*cos(thita2))+phase1);
% figure,imagesc(mask7);
mask8 = 1+m*cos(f*(X*sin(thita2)+Y*cos(thita2))+phase2);
% figure,imagesc(mask8);
mask9 = 1+m*cos(f*(X*sin(thita2)+Y*cos(thita2))+phase3);
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

c1=1+2/m*cos(f*x+phase1);
c2=1+2/m*cos(f*x+phase2);
c3=1+2/m*cos(f*x+phase3);

for i = 1:w1
    cimage1(i,:) = c1;
    cimage2(i,:) = c2;
    cimage3(i,:) = c3;
end
cimage4 = 1+2/m*cos(f*(X*sin(thita1)+Y*cos(thita1))+phase1);
% figure,imagesc(cimage4);
cimage5 = 1+2/m*cos(f*(X*sin(thita1)+Y*cos(thita1))+phase2);
% figure,imagesc(cimage5);
cimage6 = 1+2/m*cos(f*(X*sin(thita1)+Y*cos(thita1))+phase3);
% figure,imagesc(cimage6);
cimage7 = 1+2/m*cos(f*(X*sin(thita2)+Y*cos(thita2))+phase1);
% figure,imagesc(cimage7);
cimage8 = 1+2/m*cos(f*(X*sin(thita2)+Y*cos(thita2))+phase2);
% figure,imagesc(cimage8);
cimage9 = 1+2/m*cos(f*(X*sin(thita2)+Y*cos(thita2))+phase3);
% figure,imagesc(cimage9);


sim=I1.*cimage1+I2.*cimage2+I3.*cimage3+I4.*cimage4+I5.*cimage5+I6.*cimage6+I7.*cimage7+I8.*cimage8+I9.*cimage9;

figure,imshow(I,[]);title('WF image');
figure,imshow(sim,[]);title('SIM image');

Fwf=fftshift(fft2(I));
Fsim=fftshift(fft2(sim));
figure,imagesc(abs(Fwf.^(1/10)));title('FT of WF image');
figure,imagesc(abs(Fsim.^(1/10)));title('FT of SIM image');

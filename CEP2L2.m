% This Demo repeats the results of the follwing papers 
% Please cite them if you find them useful for your research 

% Reference:
% 1: W Lu, J Duan, Z Qiu, Z Pan, W Ryan Liu, L Bai
%    Implementation of high order variational models made easy for image processing

% 2: J Duan, Z Qiu, W Lu, G Wang, Z Pan, L Bai
%    An edge-weighted second order variational model for image decomposition

% code Writen by 
% Wenqi Lu and Jinming Duan
% contact email: j.duan@imperial.ac.uk
% March 2018

function CEP2L2()
clc
clear
close all

f0=imread('castle.bmp');
f0=f0(:,:,1);
figure; imagesc(f0); colormap(gray); axis off; axis equal;

padNum=10;
f0=padarray(f0,[padNum,padNum],'symmetric'); 

[m,n]=size(f0);
f0=double(f0);

u2 =zeros(m,n);
w1 =zeros(m,n);
w2 =zeros(m,n);
b11=zeros(m,n);
b12=zeros(m,n);
v  =zeros(m,n);
b2 =zeros(m,n);

lamda =30; 
alfa  =4*lamda; 
theta1=5;
theta2=5;

[Y,X]=meshgrid(0:n-1,0:m-1);
G=cos(2*pi*X/m)+cos(2*pi*Y/n)-2;

tic
for step=1:100
    step
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update u1 using FFT
    div_w_b=Bx(b11-w1)+By(b12-w2);
    g1=theta1*div_w_b+(f0-u2);
    u1=real(ifft2(fft2(g1)./(1-2*theta1*G)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update w using solt thresholding
    c1=Fx(u1)+b11;
    c2=Fy(u1)+b12;
    abs_c=sqrt(c1.^2+c2.^2+eps);
    w1=max(abs_c-lamda/theta1,0).*c1./abs_c;
    w2=max(abs_c-lamda/theta1,0).*c2./abs_c;   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % update u2 using FFT
    div_v_b=Bx(Fx(v-b2))+By(Fy(v-b2));
    g2=f0-u1+theta2*div_v_b;
    u2=real(ifft2(fft2(g2)./(1+4*theta2*G.^2)));
 
    % update v using solt thresholding
    s=Bx(Fx(u2))+By(Fy(u2))+b2;
    v=max(abs(s)-alfa/theta2,0).*sign(s);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update Bregman iterative parameters
    b11=c1-w1;
    b12=c2-w2;
    b2 =s-v;

end

u=u1+u2;
figure; imagesc(u(padNum+1:m-padNum,padNum+1:n-padNum)); colormap(gray); axis off; axis equal;

% Forward derivative operator on x with periodic boundary condition
function Fxu = Fx(u)
Fxu = circshift(u,[0 -1])-u;

% Forward derivative operator on y with periodic boundary condition
function Fyu = Fy(u)
Fyu = circshift(u,[-1 0])-u;

% Backward derivative operator on x with periodic boundary condition
function Bxu = Bx(u)
Bxu = u - circshift(u,[0 1]);

% Backward derivative operator on y with periodic boundary condition
function Byu = By(u)
Byu = u - circshift(u,[1 0]);
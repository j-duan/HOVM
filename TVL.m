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

function CTVL()
clc
close all

f0=imread('castle.bmp');
f0=f0(:,:,1);
figure; imagesc(f0); colormap(gray); axis off; axis equal;
padNum=10;
f0=padarray(f0,[padNum,padNum],'symmetric'); % extend the boundaries 

[m,n]=size(f0);
f0=double(f0);

w1=zeros(m,n);
w2=zeros(m,n);
b11=zeros(m,n);
b12=zeros(m,n);
v=zeros(m,n);
b2=zeros(m,n);

alfa=20; 
beta=20; 
theta1=5;
theta2=5;

[Y,X]=meshgrid(0:n-1,0:m-1);
G=cos(2*pi*X/m)+cos(2*pi*Y/n)-2;

tic
for step=1:100
    step
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update u using FFT
    div_w_b=Bx(w1-b11)+By(w2-b12);
    lap_v_b=Bx(Fx(v-b2))+By(Fy(v-b2));
    g=f0-theta1*div_w_b+theta2*lap_v_b;
    u=real(ifftn(fftn(g)./(1-2*theta1*G+4*theta2*G.^2)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update w using solt thresholding
    c1=Fx(u)+b11;
    c2=Fy(u)+b12;
    
    abs_c=sqrt(c1.^2+c2.^2+eps);
    w1=max(abs_c-alfa/theta1,0).*c1./abs_c;
    w2=max(abs_c-alfa/theta1,0).*c2./abs_c;
    
    % update v using solt thresholding
    s=Bx(Fx(u))+By(Fy(u))+b2;
    v=max(abs(s)-beta/theta2,0).*sign(s);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update Bregman iterative parameters
    b11=c1-w1;
    b12=c2-w2;
    b2=s-v;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
toc
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
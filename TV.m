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

function splitBregmanFFTROF()
clc
close all

f0=imread('castle.bmp');
 
f0=f0(:,:,1);
figure; imagesc(f0); colormap(gray); axis off; axis equal;

% to ease the artefact cause by periodic boundary condition 
padNum=10; 
f0=padarray(f0,[padNum,padNum],'symmetric');

[m,n]=size(f0);
f0=double(f0);
u=f0;

% initialisation of auxiliary variable w=(w1 w2) and bregman iterative
% parameter b=(b1 b2)
b1=zeros(m,n);
b2=zeros(m,n);
w1=zeros(m,n);
w2=zeros(m,n);

lamda=30; % smooth parameter the larger the smoother
theta=5; % penalty parameter, change of which will affect convergence rate

% FFT cofficients 
[Y,X]=meshgrid(0:n-1,0:m-1);
G=cos(2*pi*X/m)+cos(2*pi*Y/n)-2;

tic
for step=1:100
    
    temp_u = u;
    
    % update u using FFT
    div_w_b=Bx(w1-b1)+By(w2-b2);
    g=f0-theta*div_w_b;
    u=real(ifft2(fft2(g)./(1-2*theta*G)));
    
    % update w using solt thresholding
    ux=Fx(u);
    uy=Fy(u);
    c1=ux+b1;
    c2=uy+b2;
    abs_c=sqrt(c1.^2+c2.^2+eps);
    w1=max(abs_c-lamda/theta,0).*c1./abs_c;
    w2=max(abs_c-lamda/theta,0).*c2./abs_c;
    
    % update Bregman iterative parameters
    b1=c1-w1;
    b2=c2-w2;
    
    % stopping condition (one can also use energy to break the loop)
    if sum(abs(u-temp_u))/sum(abs(temp_u)) < 0.000001
        break;
    end
    
end
toc
figure; imagesc(u(padNum+1:m-padNum,padNum+1:n-padNum)); colormap(gray); axis off; axis equal;
figure; imagesc(f0(padNum+1:m-padNum,padNum+1:n-padNum)-u(padNum+1:m-padNum,padNum+1:n-padNum)); colormap(gray); axis off; axis equal;

% Forward derivative operator on x with periodic boundary condition
function Fxu = Fx(u)
Fxu = circshift(u,[0 -1])-u;

% Forward derivative operator on y with periodic boundary condition
function Fyu = Fy(u)
Fyu = circshift(u,[-1 0])-u;

% Backward derivative operator on x with periodic boundary condition
function Bxu = Bx(u)
Bxu = u-circshift(u,[0 1]);

% Backward derivative operator on y with periodic boundary condition
function Byu = By(u)
Byu = u-circshift(u,[1 0]);


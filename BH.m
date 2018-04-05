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

function splitBregmanSecondOrderTV()
clc
close all

f0=imread('castle.bmp');
f0=f0(:,:,1);
figure; imagesc(f0); colormap(gray); axis off; axis equal;

padNum=5;
f0=padarray(f0,[padNum,padNum],'symmetric');

[m,n]=size(f0);
f0=double(f0);
u=f0;

w1=zeros(m,n);
w2=zeros(m,n);
w3=zeros(m,n);

b1=zeros(m,n);
b2=zeros(m,n);
b3=zeros(m,n);

alfa=40; theta=2;

[Y,X]=meshgrid(0:n-1,0:m-1);
G=cos(2*pi*X/m)+cos(2*pi*Y/n)-2;

tic
for step=1:100
    step
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update u using FFT
    div_w_b=Bx(Fx(w1-b1))+By(Fy(w2-b2))+2*By(Bx(w3-b3));
    g=f0+theta*div_w_b;
    u=real(ifftn(fftn(g)./(1+4*theta*G.^2)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update w using solt thresholding
    c1=Bx(Fx(u))+b1;
    c2=By(Fy(u))+b2;
    c3=Fy(Fx(u))+b3;
    
    abs_c=sqrt(c1.^2+c2.^2+2*c3.^2+eps);
    w1=max(abs_c-alfa/theta,0).*c1./abs_c;
    w2=max(abs_c-alfa/theta,0).*c2./abs_c;
    w3=max(abs_c-alfa/theta,0).*c3./abs_c;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update Bregman iterative parameters
    b1=c1-w1;
    b2=c2-w2;
    b3=c3-w3;
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
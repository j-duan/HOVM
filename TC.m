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

function TotalCurvature()
clc
clear
close all

f=imread('castle.bmp');
figure; imshow(f); colormap(gray); axis off; axis equal;
f=double(f(:,:,1));

padNum=5;
f=padarray(f,[padNum,padNum],'symmetric');

[m,n]=size(f);

mu1=10;
mu2=10;
mu3=1e4;
mu4=5e2;
lamda=2e2;

u =f;
p1=zeros(m,n);
p2=zeros(m,n);
n1=zeros(m,n);
n2=zeros(m,n);
m1=zeros(m,n);
m2=zeros(m,n);

lamda1 =zeros(m,n);
lamda21=zeros(m,n);
lamda22=zeros(m,n);
lamda3 =zeros(m,n);
lamda41=zeros(m,n);
lamda42=zeros(m,n);

[Y,X]=meshgrid(0:n-1,0:m-1);
G=cos(2*pi*Y/n)+cos(2*pi*X/m)-2;
theta2=mu4/mu3;
a11=theta2-2*(cos(2*pi*Y/n)-1);
a22=theta2-2*(cos(2*pi*X/m)-1);
i=sqrt(-1);
a12=-(-1+cos(2*pi*Y/n)+i*sin(2*pi*Y/n)).*(1-cos(2*pi*X/m)+i*sin(2*pi*X/m));
a21=-(-1+cos(2*pi*X/m)+i*sin(2*pi*X/m)).*(1-cos(2*pi*Y/n)+i*sin(2*pi*Y/n));
D=theta2^2-2*theta2*G;

tic
for step=1:100
    step
    
    div_p=Bx(p1)+By(p2);
    div_lamda2=Bx(lamda21)+By(lamda22);
    g=f-mu2*div_p-div_lamda2;
    u=real(ifft2(fft2(g)./(1-2*mu2*G)));
    u=max(min(u,255),0);
    
    div_n=Bx(n1)+By(n2);
    temp=div_n-lamda3/mu3;
    q=max(abs(temp)-lamda/mu3,0).*sign(temp);
    
    temp1=Fx(u)+((mu1+lamda1).*m1-lamda21)/mu2;
    temp2=Fy(u)+((mu1+lamda1).*m2-lamda22)/mu2;
    abs_temp=sqrt(temp1.^2+temp2.^2+eps);
    p1=max(abs_temp-(lamda1+mu1)/mu2,0).*temp1./abs_temp;
    p2=max(abs_temp-(lamda1+mu1)/mu2,0).*temp2./abs_temp;
    
    m1=n1+((mu1+lamda1).*p1+lamda41)/mu4;
    m2=n2+((mu1+lamda1).*p2+lamda42)/mu4;
    abs_m=sqrt(m1.^2+m2.^2);
    m1=m1./max(abs_m,1);
    m2=m2./max(abs_m,1);
    
    h1=-Fx(q+lamda3/mu3)+mu4*m1/mu3-lamda41/mu3;
    h2=-Fy(q+lamda3/mu3)+mu4*m2/mu3-lamda42/mu3;
    Fh1=fft2(h1);
    Fh2=fft2(h2);
    n1=real(ifft2((a22.*Fh1-a12.*Fh2)./D));
    n2=real(ifft2((a11.*Fh2-a21.*Fh1)./D));
    
    lamda1 =lamda1 +mu1*(sqrt(p1.^2+p2.^2)-p1.*m1-p2.*m2);
    lamda21=lamda21+mu2*(p1-Fx(u));
    lamda22=lamda22+mu2*(p2-Fy(u));
    lamda3 =lamda3 +mu3*(q-Bx(n1)-By(n2));
    lamda41=lamda41+mu4*(n1-m1);
    lamda42=lamda42+mu4*(n2-m2);
    
end
toc
u=u(padNum+1:m-padNum,padNum+1:n-padNum);
figure; imshow(u, []); colormap(gray); axis off; axis equal;

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

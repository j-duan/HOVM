% This Demo shows a second order Vese-Osher image decomposition model 
% which is solved by the fast FFT-based split Bregman

% Please cite the following papers if you find them useful for your research 
% Reference:
% 1: W Lu, J Duan, Z Qiu, Z Pan, W Ryan Liu, L Bai
%    Implementation of high order variational models made easy for image processing

% 2: J Duan, Z Qiu, W Lu, G Wang, Z Pan, L Bai
%    An edge-weighted second order variational model for image decomposition

% code Writen by 
% Jinming Duan
% contact email: j.duan@imperial.ac.uk
% March 2018

function splitBregmanSecondOrderVO()
clc
close all

f0=imread('barbara.png');
f0=f0(:,:,1);
figure; imagesc(f0); colormap(gray); axis off; axis equal;
padNum=10;
f0=padarray(f0,[padNum,padNum],'symmetric');
[m,n]=size(f0);
f0=double(f0);

b11=zeros(m,n);
b12=zeros(m,n);
% b13=zeros(m,n);
b14=zeros(m,n);

w1=zeros(m,n);
w2=zeros(m,n);
% w3=zeros(m,n);
w4=zeros(m,n);

m1=zeros(m,n);
m2=zeros(m,n);
b21=zeros(m,n);
b22=zeros(m,n);

g1=zeros(m,n);
g2=zeros(m,n);
u=f0;

% see the above paper to understand how to tune 
% the following paper
lamda=30; % smooth parameter
mu=10; % for edge weights
ganm=3; % texture parameter
theta1=5; % convergence parameter 1
theta2=3; % convergence parameter 2

[Y,X]=meshgrid(0:n-1,0:m-1);
G=cos(2*pi*Y/n)+cos(2*pi*X/m)-2;
a11=theta2-2*(cos(2*pi*Y/n)-1);
a22=theta2-2*(cos(2*pi*X/m)-1);
i=sqrt(-1);
a12=-(-1+cos(2*pi*Y/n)+i*sin(2*pi*Y/n)).*(1-cos(2*pi*X/m)+i*sin(2*pi*X/m));
a21=-(-1+cos(2*pi*X/m)+i*sin(2*pi*X/m)).*(1-cos(2*pi*Y/n)+i*sin(2*pi*Y/n));
D=theta2^2-2*theta2*G;

% set parameters for edge detector
sigma=1; % Gaussian standard deviation
scalarP=100; % scalar parameter controlling thickness of edge
gb=edgeFinder(f0,sigma,scalarP,m,n);% edge detector
figure; imagesc(gb); colormap(gray); axis off; axis equal;

rr=10;
Cm = Cmcalc(rr);
tic
for step=1:100
    step
    
    temp_u=u;
    temp_g1=g1;
    temp_g2=g2;
    
    temp_b11=b11;
    temp_b12=b12;
    temp_b14=b14;
    temp_b21=b21;
    temp_b22=b22;
    
    % update u using FFT
    % div_w_b=Fx(Bx(w1-b11))+By(Bx(w2-b12))+Bx(By(w3-b13))+Fy(By(w4-b14));
    div_w_b=Fx(Bx(w1-b11))+2*By(Bx(w2-b12))+Fy(By(w4-b14));
    div_g=Bx(g1)+By(g2);
    P=f0-div_g+theta1*div_w_b;
    u=real(ifft2(fft2(P)./(1+4*theta1*G.^2)));
    
    % update g1 and g2 using FFT
    h1=Fx(u-f0)+theta2*(m1-b21);
    h2=Fy(u-f0)+theta2*(m2-b22);
    Fh1=fft2(h1);
    Fh2=fft2(h2);
    g1=real(ifft2((a22.*Fh1-a12.*Fh2)./D));
    g2=real(ifft2((a11.*Fh2-a21.*Fh1)./D));
    
     % update weights
    [gradx, grady] = gsderiv(u, sigma, 1);
    grad = sqrt(gradx.^2 + grady.^2 +eps);
    gb = 1 - exp(-Cm./( grad/mu).^rr);
    
    % update w using solt thresholding
    uxx=Bx(Fx(u));
    uxy=Fy(Fx(u));
    % uyx=Fx(Fy(u));
    uyy=By(Fy(u));
    s1=uxx+b11;
    s2=uxy+b12;
    % s3=uyx+b13;
    s4=uyy+b14;
    % abs_s=sqrt(s1.^2+s2.^2+s3.^2+s4.^2+eps);
    abs_s=sqrt(s1.^2+2*s2.^2+s4.^2+eps);
    w1=max(abs_s-gb*lamda/theta1,0).*s1./abs_s;
    w2=max(abs_s-gb*lamda/theta1,0).*s2./abs_s;
    % w3=max(abs_s-gb*lamda/theta1,0).*s3./abs_s;
    w4=max(abs_s-gb*lamda/theta1,0).*s4./abs_s;
    
    % update v using solt thresholding
    c1=g1+b21;
    c2=g2+b22;
    abs_c=sqrt(c1.^2+c2.^2+eps);
    m1=max(abs_c-ganm/theta2,0).*c1./abs_c;
    m2=max(abs_c-ganm/theta2,0).*c2./abs_c;
    
    % update Bregman iterative parameters
    b11=s1-w1;
    b12=s2-w2;
    % b13=s3-w3;
    b14=s4-w4;
    b21=c1-m1;
    b22=c2-m2;
    
    temp1=abs(u-temp_u);
    temp2=abs(temp_u);
    Err_u(step)=log(sum(temp1(:))/sum(temp2(:)));
    
    temp1=(temp_g1-g1).^2+(temp_g2-g2).^2;
    temp1=sqrt(temp1);
    temp2=temp_g1.^2+temp_g2.^2;
    temp2=sqrt(temp2);
    Err_g(step)=log(sum(temp1(:))/sum(temp2(:)));
    
    A=m*n;
    abs_w=sqrt((w1-uxx).^2+2*(w2-uxy).^2+(w4-uyy).^2);
    abs_v=sqrt((m1-g1).^2+(m2-g2).^2);
    Err_w(step)=log(sum(abs_w(:))/A);
    Err_v(step)=log(sum(abs_v(:))/A);
    
    temp1=(temp_b11-b11).^2+2*(temp_b12-b12).^2+2*(temp_b14-b14).^2;
    temp1=sqrt(temp1);
    temp2=temp_b11.^2+2*temp_b12.^2+temp_b14.^2;
    temp2=sqrt(temp2);
    Err_b1(step)=log(sum(temp1(:))/sum(temp2(:)));
    
    temp1=(temp_b21-b21).^2+(temp_b22-b22).^2;
    temp1=sqrt(temp1);
    temp2=temp_b21.^2+temp_b22.^2;
    temp2=sqrt(temp2);
    Err_b2(step)=log(sum(temp1(:))/sum(temp2(:)));
    
    % update energy
    abs_u=sqrt(uxx.^2+2*uxy.^2+uyy.^2);
    abs_g=sqrt(g1.^2+g2.^2);
    v=Bx(g1)+By(g2);
    E=lamda*gb.*abs_u+ganm*abs_g+0.5*(f0-u-v).^2;
    Energy(step)=(sum(E(:)));
    
end
toc

outPutU=u(padNum+1:m-padNum,padNum+1:n-padNum);
outPutV=v(padNum+1:m-padNum,padNum+1:n-padNum);
outPutN=f0(padNum+1:m-padNum,padNum+1:n-padNum)-outPutU-outPutV;

imwrite(uint8(outPutU),['..\' num2str(1) '.bmp']);
imwrite(uint8(outPutV+150),['..\' num2str(2) '.bmp']);

figure; imagesc(outPutU); colormap(gray); axis off; axis equal;
figure; imagesc(outPutV); colormap(gray); axis off; axis equal;
figure; imshow(uint8(outPutN+150)); colormap(gray); axis off; axis equal;

% plots quantities showing convergence of the algorithm
stp = 5;
idxs = 1:stp:numel(Energy);
figure;

figure;
plot(idxs,Err_v(idxs),'bx','LineWidth',1)
hold on;
plot(idxs,Err_w(idxs),'rs','LineWidth',1)
hold on;
plot(Err_w,'r','LineWidth',1);
hold on;
plot(Err_v,'b','LineWidth',1);
hold on;
xlabel('Iterations');
ylabel('Residuals (log)');
legend('R1','R2');
grid on;

figure;
plot(idxs,Err_b1(idxs),'rx','LineWidth',1)
hold on;
plot(idxs,Err_b2(idxs),'bs','LineWidth',1)
hold on;
plot(Err_b1,'r','LineWidth',1);
hold on;
plot(Err_b2,'b','LineWidth',1);
xlabel('Iterations');
ylabel('Relative error in Bregman iterative parameters (log)');
legend('B1','B2');
grid on;

figure;
plot(idxs,Err_u(idxs),'rx','LineWidth',1)
hold on;
plot(idxs,Err_g(idxs),'bs','LineWidth',1)
hold on;
plot(Err_u,'r','LineWidth',1);
hold on;
plot(Err_g,'b','LineWidth',1);
xlabel('Iterations');
ylabel('Relative error in u and g (log)');
legend('R_u','g');
grid on;

figure;
plot(idxs,Energy(idxs),'rs','LineWidth',1)
hold on;
plot(Energy,'r-','LineWidth',1);
xlabel('Iterations');
ylabel('Energy (log)');
legend('split Bregman algorithm');
grid on

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

% Compute edge detector
function gb = edgeFinder(Img,sigma,scalarP,rows,cols)
Imin  = min(Img(:));
Imax  = max(Img(:));
normImg = (Img-Imin)/(Imax-Imin);  % Normalize Img to the range [0,1]
Gauss = fspecial('gaussian',[rows,cols],sigma);
Im0s = real(ifftshift(ifft2(fft2(normImg).* fft2(Gauss))));
[Gx,Gy] = gradient(Im0s);
NormGrad = sqrt(Gx.^2 + Gy.^2);
gb = 1./ (1 + scalarP*NormGrad.^2);% edge detector;

function Cm = Cmcalc(m)
if m <= 1
    error('Use m > 1')
    return
else
    Cm = fzero(strcat('1-exp(-x)-x*exp(-x)*',num2str(m)),[1e-10 1e100]);
end
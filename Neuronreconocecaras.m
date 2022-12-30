% Entrenamiento patron
clc
clear all
close all
%Imagen1
img_1=imread('roy.jpg');
img_1_bin=imbinarize(img_1(:,:,1));
img_1_canny=edge(img_1(:,:,1),'Canny');
figure(1)
image(img_1)

%Imagen1
img_2=imread('alex.jpg');
img_2_bin=imbinarize(img_2(:,:,1));
img_2_canny=edge(img_2(:,:,1),'Canny');
figure(2)
image(img_2)


%%
% [M1,N1]=size(img_1_bin);
% [M2,N2]=size(img_2_bin);
% MM_NN=[M1,M2,N1,N2];
% min_MM_NN=min(MM_NN);
%min_MM_NN=20;
min_MM_NN=300;
%convertir la imagen en un cuadrado
img_1_resize=imresize(img_1_canny,[min_MM_NN min_MM_NN]);
img_2_resize=imresize(img_2_canny,[min_MM_NN min_MM_NN]);

figure(3)
image(img_1_resize);
imshow(img_1_resize);
figure(4)
image(img_2_resize);
imshow(img_2_resize);
%%
% Pixel blanco = 0 
% Pixel negro  = 1

cara1=img_1_resize;
cara2=img_2_resize;

% Pixel blanco = 0 
% Pixel negro  = 1


[nf nc] = size(cara1);
x(1,:) = cara1(1,:);
x(2,:) = cara2(1,:);

for k = 2:nf
   caras = [ cara1(k,:)
             cara2(k,:)
             ];
   x = [ x  caras ];
end
[ nxf nxc ] = size(x);
nx = nxf;

yb(1,:) = [ 1 0 ]; 
yb(2,:) = [ 0 1 ]; 

[ nyf nyc ] = size(yb);
ny = nyf;

ne = nxc;
nm = 400;
nm1= 400;
ns = nyc;
bias = input('Bias:  SI = 1 : ');
if(bias == 1)
   ne = ne + 1;
   x = [ x ones(nx,1) ];   
end
v = 0.2*(rand(ne,nm) - 0.5);
v1= 0.2*(rand(nm,nm1) - 0.5);
w = 0.2*(rand(nm,ns) - 0.5);


% load pesoscaras;

eta = input('eta pesos : ');

for iter = 1:1000
count(iter,1) = iter;
dJdw = 0;
dJdv = 0;
dJdv1=0;
for k = 1:nx   
in = (x(k,:))';
m = v'*in;
% n = 1.0./(1+exp(-m));    % Sigmoidea 1
n = 2.0./(1+exp(-m)) - 1; % sigmoidea 2
% n = exp(-m.^2);         % Gaussiana
%SEGUNDA CAPA
o=v1'*n;
% n1 = 1.0./(1+exp(-o));    % Sigmoidea 1
n1 = 2.0./(1+exp(-o)) - 1; % sigmoidea 2
% n1 = exp(-o.^2);         % Gaussiana

out = w'*n1;

y(k,:) = out';
er = out - (yb(k,:))';
error(k,:) = er';
% dn1dm = n1.*(1 - o);       % Sigmoidea 1
dn1dm = (1 - n1.*n1)/2;    % Sigmoidea 2
%  dn1dm = -2.0*(n1.*o);     % Gaussiana

% dndm = n.*(1 - n);       % Sigmoidea 1
dndm = (1 - n.*n)/2;    % Sigmoidea 2
%  dndm = -2.0*(n.*m);     % Gaussiana
dydw = n1; 
dJdw = 1*dJdw + dydw*er';    
dJdv1 = 1*dJdv1 + n*((w*er).*dn1dm)'; 
dJdv = 1*dJdv + in*((w*er).*dndm)'; 
%w = w - eta*dJdw/nx;   
%v = v - eta*dJdv/nx;      
end
w = w - eta*dJdw/nx;
v1 = v1 - eta*dJdv1/nx; 
v  = v - eta*dJdv/nx; 

JJ = 0.5*sum(sum(error.*error))
J(iter,1) = JJ;
end

save pesoscaras1 ne ns nm nm1 v v1 w  bias; 

figure();
plot(y(2,:),'or');   % Se grafica sin redondeo
hold on;
plot(yb(2,:),'*b');

figure();
plot(count,J);
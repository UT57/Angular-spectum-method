%% Задаём поле Z=0
%Парамеры источника:
global Lx Ly D f c
Lx =10;Ly=10; f =10^6; c =1500000;

%Расчёт шага сетки,k,lymda:
global k lymda step
lymda =c/f; k = 2*pi/lymda; step = lymda/6;

%Задаём поле источника
global N M
N = Lx/step;
M = Ly/step;
R = N/2;
IST = zeros(2*N,2*M);
I1 = 1:2*N; 
I2 = 1:2*M;
x = I1-N;                 % mask x-coordinates 
y = I1-M;                 % mask y-coordinates
[X,Y] = meshgrid(x,y);    % create 2-D mask grid
A = (X.^2 + Y.^2 <= R^2); % circular aperture of radius R
IST(A) = 1; 
%pcolor(IST);


%% УГЛОВОЙ СПЕКТР
global z 
z = 2000;
hh = (pi/Lx)^2;


%Расчёт углового спектра(БПФ):
%F =(fft2(IST));
F = fftshift(fft2(IST));

%Фильтрация:
A = (X.^2 + Y.^2 <= (k*Lx/pi)^2); 
IST(A) = 1;


%Домножение на член распространения:
for n = 1:2*N
    for m = 1:2*M
         kx = (n-N)*pi/Lx;
         ky = (m-M)*pi/Ly;
        S(n,m) = F(n,m)*exp(1i*z*sqrt(k^2-kx^2-ky^2)-i*k*z)*A(n,m);     
    end
end
%Фильтрация:
A = (X.^2 + Y.^2 <= (k*Lx/pi)^2); 
IST(A) = 1;


%Обратное БПФ:     
 Res =ifft2(S);
% pcolor(abs(Res)); %Результат прямой задачи

%% Обратная задача
Res_0 = conj(Res);
Fobr = fftshift(fft2(Res_0));
%Домножение на член распространения:

for n = 1:2*N
    for m = 1:2*M
        kx = (n-N)*pi/Lx;
        ky = (m-M)*pi/Ly;
        Sobr(n,m) = Fobr(n,m)*exp(1i*z*sqrt(k^2-kx^2-ky^2)-i*k*z)*A(n,m);     
    end
end

%Обратное БПФ:     
Res_obr =(ifft2(Sobr));
pcolor(abs(Res_obr));   %Результат обратной задачи


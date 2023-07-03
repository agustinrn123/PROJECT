clc,clear

%Data x1=Susceptible (S) x2=Infected (I), x3=Recovered (R)

%pendefinisian parameter
% A=0.0256;
m    = 0.12;
b    = 0.731;
ro   = 0.0713;
myu  = 0.0031/365;
gamma= 0.002;
dt   = 1;
T    = 244;
R    = 0.1;
Q    = 0.01;

%pendefinisian matriks parameter
H  = [0,0,0;0,1,0;0,0,1];
Rk = eye(3)*R;
Qk = eye(3)*Q;
G  = eye(3);

%pengalokasian variabel
f = zeros(3,T);
x = zeros(3,T);
z = zeros(3,T);
h = zeros(3,T);
hhat = zeros(3,T);
xhat = zeros(3,T);
xhat0= zeros(3,T);

%Pendefinisian nilai awal
p0    = [1,0,0;0,1,0;0,0,1];
x(:,1)= [667173,14,0];

%Model Sistem & Model Pengukuran
for i=1:T-1
N = sum(x(:,i));
A = 0.0261*N/365;
%Model Sistem
%x(k+1)=f(xk,uk)+wk
    f(:,i) = [(A-m*b/N*x(1,i)*x(2,i)-myu*x(1,i))*dt+x(1,i);...
       (m*b/N*x(1,i)*x(2,i)-ro*x(2,i)-gamma*x(2,i))*dt+x(2,i);...
        (ro*x(2,i)-myu*x(3,i))*dt+x(3,i)];
    x(:,i+1) = f(:,i);
%Model Pengukuran
%z(k)=Mxk+vk
       h(:,i+1) = [(A-m*b/N*x(1,i+1)*x(2,i+1)-myu*x(1,i+1))*dt+x(1,i+1);...
       (m*b/N*x(1,i+1)*x(2,i+1)-ro*x(2,i+1)-gamma*x(2,i+1))*dt+x(2,i+1);...
        (ro*x(2,i+1)-myu*x(3,i+1))*dt+x(3,i+1)];
       z(:,i) = h(:,i);

%Tahap Prediksi (Time Update)
%Estimasi
    A2 = [(-m*b/N*x(2,i)-myu)*dt+1,-m*b/N*x(1,i)*dt,0;...
        m*b/N*x(2,i)*dt,(m*b/N*x(1,i)-ro-gamma)*dt+1,0;...
        0,ro*dt,-myu*dt+1];
    
%x_hat(k+1)=Akxk + Bkuk
     xhat(:,i) = f(:,i);
%Kovariansi Error
%pMin(k+1)=A2kPkA2k'+PkA2k+GkQkGk'
    pMin = A2*p0+p0*A2'+G*Qk*G';

%Tahap Koreksi (Measurement Update)
%Kalman Gain
%Kmin(k+1)=pMin(k+1)H'(HpMin+R
    K = pMin*H'*(inv(H*pMin*H'+Rk));
%Estimasi
%xhat(k+1)=x_hat(k+1)+Kmin(k+1)(z(k+1)-Hx_hat(k+1))
    hhat(:,i+1) = [(A-m*b/N*xhat(1,i+1)*xhat(2,i+1)-myu*xhat(1,i+1))*dt+xhat(1,i+1);...
        (m*b/N*xhat(1,i+1)*xhat(2,i+1)-ro*xhat(2,i+1)-gamma*xhat(2,i+1))*dt+xhat(2,i+1);...
        (ro*xhat(2,i+1)-myu*xhat(3,i+1))*dt+xhat(3,i+1)];
    
    xhat0(:,i+1)=(xhat(:,i+1)+K*(z(:,i+1)-hhat(:,i+1)));
%Kovarian Error
    KE = (eye(3)-K*H)*pMin;
end
xhat0;

xreal = xlsread("DataMatlab");

for j = 1:T
E1(j) = (xreal(j,1)')-(xhat0(1,j));
SE1(j)= E1(j).^2;
RMSE_x1(j) = sqrt(sum(SE1(j))/T);
E2(j) = (xreal(j,2)')-(xhat0(2,j));
SE2(j)= E2(j).^2;
RMSE_x2(j) = sqrt(sum(SE2(j))/T);
E3(j) = (xreal(j,3)')-(xhat0(3,j));
SE3(j)= E3(j).^2;
RMSE_x3(j) = sqrt(sum(SE3(j))/T);
end
rataRMSE_x1 = mean(RMSE_x1)
rataRMSE_x2 = mean(RMSE_x2)
rataRMSE_x3 = mean(RMSE_x3)

%kesalahan Error Root Mean Square Error (RMSE)
% RMSE_x1 = sqrt(sum((xhat0(1,:)-xreal(:,1)').^2))/T
% RMSE_x2 = sqrt(sum((xhat0(2,:)-xreal(:,2)').^2))/T
% RMSE_x3 = sqrt(sum((xhat0(3,:)-xreal(:,3)').^2))/T

%MAPE
% Mape_x1 = (1/T*(sum(abs(xreal1(:,1)-xhat0(:,1))/xreal1(:,1))))*100
% Mape_x2 = (1/T*(sum(abs(xreal1(:,2)-xhat0(:,2))/xreal1(:,2))))*100
% Mape_x3 = (1/T*(sum(abs(xreal1(:,3)-xhat0(:,3))/xreal1(:,3))))*100

figure(1)
plot(xhat(1,:),'m','LineWidth',2)
hold on
plot(xhat0(1,:),'m','LineWidth',2)
hold on
plot(xreal(:,1),'b','LineWidth',2)
xlabel('Hari Ke-')
ylabel('Jumlah Individu')
title('Hasil Estimasi Variabel Susceptibel (S)')
legend('Nilai Estimasi Susceptible')
grid on;

figure(2)
plot(xhat0(2,:),'r','LineWidth',2)
hold on
plot(xreal(:,2),'b','LineWidth',2)
xlabel('Hari Ke-')
ylabel('Jumlah Individu')
title('Hasil Estimasi Variabel Infected (I)')
legend('Nilai Estimasi Infected')
grid on;

figure(3)
plot(xhat0(3,:),'b','LineWidth',2)
hold on
plot(xreal(:,3),'r','LineWidth',2)
xlabel('Hari Ke-')
ylabel('Jumlah Individu')
title('Hasil Estimasi Variabel Recovered (R)')
legend('Nilai Estimasi Recovered')
grid on;

% 
% for j=1:k
% E1(j)=x_r1(j)-x_est1(j);
% SE1(j)=E1(j).^2;
% RMSE1(j)=sqrt(sum(SE1(j))/k)/100000;
% E2(j)=x_r2(j)-x_est2(j);
% SE2(j)=E2(j).^2;
% RMSE2(j)=sqrt(sum(SE2(j))/k)/100;
% end
% 
% rataRMSE1=mean(RMSE1);
% rataRMSE2=mean(RMSE2);
% figure(1)
% plot((2:k),x_est2(2:k),'-*b'),title('Perbandingan Jumlah Penderita COVID-19');
% xlabel('Tahun');
% ylabel('Jumlah Penderita COVID-19');
% legend('Nilai Extended Kalman Filter');
% figure(2)
% plot((2:k),x_est1(2:k),'-*b'),title('Perbandingan Jumlah Penderita COVID-19 yang Sembuh');
% xlabel('Tahun');
% ylabel('Jumlah Penderita COVID-19 yang Sembuh');
% legend('Nilai Extended Kalman Filter');
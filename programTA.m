clc,clear

%Data x1=Susceptible (S) x2=Infected (I), x3=Recovered (R)

%Pendefinisian Parameter
A    = 0.0261/365;
m    = 0.12;
b    = 0.75;
ro   = 0.075;
myu  = 0.0031;
gamma= 0.0023;
dt   = 1;
T    = 244;
R    = 0.001;
Q    = 0.001;
G    = 1;
p    = 100;

%pendefinisian matriks parameter
H  = [1,0,0;0,1,0;0,0,1];
Rk = eye(3)*R;
Qk = eye(3)*Q;
Gk = eye(3)*G;
vk = [normrnd(0,sqrt(R),1,1)*dt;normrnd(0,sqrt(R),1,1)*dt;normrnd(0,sqrt(R),1,1)*dt];
wk = [normrnd(0,sqrt(Q),1,1)*dt;normrnd(0,sqrt(Q),1,1)*dt;normrnd(0,sqrt(Q),1,1)*dt];

%pengalokasian variabel
f = zeros(3,T);
x = zeros(3,T);
z = zeros(3,T);
h = zeros(3,T);
hhat = zeros(3,T);
xhat = zeros(3,T);
xhat0= zeros(3,T);

%Pendefinisian nilai awal
p0 = [p,0,0;0,p,0;0,0,p];
x(:,1)= [667220,14,0];
xhat(:,1)= [667220,14,0];
xhat0(:,1)= [667220,14,0];
xreal1 = xlsread("xreal1");
xreal = xlsread("DataMatlab");

%Model Sistem & Model Pengukuran
for i=1:T-1
    N = sum(x(:,i));
%%%Model Sistem
%x(k+1)=f(xk,uk)+wk
    f(:,i+1) = [(A-m*b/N*x(1,i)*x(2,i)-myu*x(1,i))*dt+x(1,i);...
       (m*b/N*x(1,i)*x(2,i)-ro*x(2,i)-gamma*x(2,i))*dt+x(2,i);...
        (ro*x(2,i)-myu*x(3,i))*dt+x(3,i)];
    x(:,i+1) = f(:,i+1)+wk;
    
%%%Model Pengukuran
%z(k)=h(xk+1)+vk
       h(:,i) = [(A-m*b/N*x(1,i)*x(2,i)-myu*x(1,i))*dt+x(1,i);...
       (m*b/N*x(1,i)*x(2,i)-ro*x(2,i)-gamma*x(2,i))*dt+x(2,i);...
        (ro*x(2,i)-myu*x(3,i))*dt+x(3,i)];
       z(:,i) = h(:,i)+vk;

%%%Tahap Prediksi (Time Update)
%Estimasi
    A2 = [(-m*b/N*xhat0(2,i)-myu)*dt+1,-m*b/N*xhat0(1,i)*dt,0;...
        m*b/N*xhat0(2,i)*dt,(m*b/N*xhat0(1,i)-ro-gamma)*dt+1,0;...
        0,ro*dt,-myu*dt+1];
%x_hat(k+1)=f(xhatk,uk)
    xhat(:,i) = f(:,i);
     
%Kovariansi Error
%pMin(k+1)=A2kPkA2k'+PkA2k+GkQkGk'
    p0 = A2*p0+p0*A2'+Gk*Qk*Gk';

%%%Tahap Koreksi (Measurement Update)
%Kalman Gain`1z
%Kmin(k+1)=pMin(k+1)H'(Hp0H'+Rk+1)^-1
    K = p0*H'*(inv(H*p0*H'+Rk));

%Estimasi
%xhat(k+1)=x_hat(k+1)+Kmin(k+1)(z(k+1)-h(xhatk+1))
    hhat(:,i+1) = [(A-m*b/N*xhat(1,i+1)*xhat(2,i+1)-myu*xhat(1,i+1))*dt+xhat(1,i+1);...
        (m*b/N*xhat(1,i+1)*xhat(2,i+1)-ro*xhat(2,i+1)-gamma*xhat(2,i+1))*dt+xhat(2,i+1);...
        (ro*xhat(2,i+1)-myu*xhat(3,i+1))*dt+xhat(3,i+1)];
    xhat0(:,i+1)=xhat(:,i+1)+K*(xreal1(:,i)-hhat(:,i+1));

%KovarKiansi Error
%[I-K(:,i)*H]pMin
    KE = (eye(3)-K*H)*p0;
end
xhat0;

%RMSE
RMSE_x1 = sqrt(sum((xreal1(1,:)-xhat0(1,:)).^2))/T
RMSE_x2 = sqrt(sum((xreal1(2,:)-xhat0(2,:)).^2))/T
RMSE_x3 = sqrt(sum((xreal1(3,:)-xhat0(3,:)).^2))/T

E1 = abs((xreal1(1,:)-xhat0(1,:)));
E2 = abs((xreal1(2,:)-xhat0(2,:)));
E3 = abs((xreal1(3,:)-xhat0(3,:)));

figure(1)
plot(xhat0(1,:),'r')
hold on
plot(xreal(:,1),'b')
xlabel('Hari Ke-')
ylabel('Jumlah Individu')
legend('Nilai Estimasi Susceptible','Nilai Real Susceptible')
grid on;

figure(2)
plot(xhat0(2,:),'r')
hold on
plot(xreal(:,2),'b')
xlabel('Hari Ke-')
ylabel('Jumlah Individu')
legend('Nilai Estimasi Infected','Nilai Real Infected')
grid on;

figure(3)
plot(xhat0(3,:),'r')
hold on
plot(xreal(:,3),'b')
xlabel('Hari Ke-')
ylabel('Jumlah Individu')
legend('Nilai Estimasi Recovered','Nilai Real Recovered')
grid on;

figure(4)
plot(E1,'r')
xlabel('Hari Ke-')
ylabel('Nilai Error X1')

figure(5)
plot(E2,'r')
xlabel('Hari Ke-')
ylabel('Nilai Error X2')

figure(6)
plot(E3,'r')
xlabel('Hari Ke-')
ylabel('Nilai Error X3')

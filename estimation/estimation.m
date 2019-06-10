clear all;
close all;
 
M = 1.3608;
m = 0.096;
R = 0.0592;
I = 2.175*10^-4;
k = 186.3;
 
%pembagi di persamaan
h = (M+m)*(I+0.5*m*R^2);
 
%nilai di matriks A
a11 = 0;
a12 = 1;
a13 = 0;
a14 = 0;
a21 = (-k*(I+m*R^2))/h;
a22 = 0;
a23 = -(((m^2)*(R^3)*(0.555*((M+m)*I+(M*m*R^2))+m*(0.616*I+0.555*m*R^2)))/(h^2));
a24 = (1.57*I*m*R+1.11*(m^2)*(R^3))/h;
 
a31 = 0;
a32 = 0;
a33 = 0;
a34 = 1;
 
a41 = (0.707*m*R*k)/h;
a42 = 0;
a43 = 0;
a44 = -((0.785*(m^2)*(R^2))/h);
 
%nilai di matriks B
b11 = 0;
b12 = 0;
b21 = (I+m*R^2)/h;
b22 = (-0.707*m*R/h);
b31 = 0;
b32 = 0;
b41 = (-0.707*m*R/h);
b42 = (M+m)/h;
 


%matriks A,B,C, dan D
A = [a11 a12 a13 a14; a21 a22 a23 a24; a31 a32 a33 a34; a41 a42 a43 a44];
B = [0 0; b21 b22; 0 0; b41 b42];
C = eye(4);
D = zeros(4,2);
 
sys_ss = ss(A,B,C,D); % representasi state space
 
%Diskretisasi Sistem
fs = 250; %Hz
Ts = 1/fs; % periode sampling (sekon)
 
sys_ss_d = c2d(sys_ss, Ts)
[Ad, Bd, Cd, Dd] = ssdata(sys_ss_d)
Ob=obsv(Ad,Cd); 
[n,m]=size(Ob); 
unob = m-rank(Ob); %observable
if(unob==0)
   disp('Given System is Observable.');
else
   disp('Given System is Unobservable');
end 
 
Co = ctrb(Ad,Bd);
[n,m]=size(Co);
unco=n-rank(Co); %controllable
if(unco==0)
   disp('Given System is Controllable.');
else
   disp('Given System is Uncontrollable');
end
 
Tf = 10;
t = 0:Ts:Tf;
u = awgn(heaviside(t), 0.1); %Input Additive White Gaussian Noise (Derau Putih)
x0 = [0; 0; 0; 0]; %Nilai awal state
[y,t,x]=lsim(sys_ss_d,[u;u],t,x0);
 
%Non-Recursive Least Square
sy=size(y); %ukuran runtun output
sx=size(x); %ukuran runtun state
uT = transpose([u; u]); %ditranspose agar menjadi vektor kolom
su=size(uT); %ukuran runtun input, 
 
%Bacaan Output y[k+1]
for i=1:(sy(1)-1)
    for j=1:sy(2)
     Y(i,j)=y(i+1,j);
    end
end
 



%Bacaan State x[k]
for i=1:(sx(1)-1)
    for j=1:sx(2)
        psiX(i,j) = x(i,j);
    end
end
 
%Bacaan Input u[k]
for i=1:(su(1)-1)
    for j=1:su(2)
        psiU(i,j) = uT(i,j);
    end
end

%Identification Data Matrix Psi dengan menggabungkan State x[k] dan u[k]
psi = horzcat(psiX, psiU);
%Gabungan transpose matriks A dan B hasil estimasi ThetaHat
thetaHat = inv(transpose(psi)*psi)*transpose(psi)*Y
 
%Mendapatkan matriks A transpose
sth = size(thetaHat);
for(i=1:(sth(1)-su(2)))
    for(j=1:sth(2))
       At(i,j)=thetaHat(i,j);
    end
end
%Mendapatkan matriks B transpose
for(i=(sth(1)-su(2)+1):sth(1))
    for(j=1:sth(2))
       Bt(i-(sth(1)-su(2)),j)=thetaHat(i,j);
    end
end
 
%Hasil Estimasi State Space
Aest = transpose(At)
Best = transpose(Bt)
Cest = eye(sy(2),sx(2))
Dest = zeros(sy(2),su(2))
sysEst = ss(Aest, Best, Cest, Dest)

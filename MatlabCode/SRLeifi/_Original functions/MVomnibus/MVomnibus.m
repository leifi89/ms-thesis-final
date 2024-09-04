function [MVomnibus] = MVomnibus(Y)

%MV omnibus test proposed by Doornik and Hansen (1994).

% Syntax: function MVomnibus(Y) 
%      
% Inputs:
%      Y - data matrix (Size of matrix must be n-by-p; data=rows,
%          indepent variable=columns) 
%   significance level (default = 0.05)
%
% Output:
%        - MV Omnibus normality test  (significance level by default = 0.05)


% Example: From the Table 11.5 (Iris data) of Johnson and Wichern (1992,
%          p. 562), we take onlt the Isis setosa data to test if it has a 
%          multivariate normality distribution using the Doornik-Hansen
%          omnibus test. Data has 50 observations on 4 independent 
%          variables (Var1 (x1) = sepal length; Var2 (x2) = sepal width; 
%          Var3 (x3) = petal length; Var4 (x4) = petal width. 

% The data is:

Y=[5.1	3.5	1.4	0.2
4.9	3.0	1.4	0.2
4.7	3.2	1.3	0.2
4.6	3.1	1.5	0.2
5.0	3.6	1.4	0.2
5.4	3.9	1.7	0.4
4.6	3.4	1.4	0.3
5.0	3.4	1.5	0.2
4.4	2.9	1.4	0.2
4.9	3.1	1.5	0.1
5.4	3.7	1.5	0.2
4.8	3.4	1.6	0.2
4.8	3.0	1.4	0.1
4.3	3.0	1.1	0.1
5.8	4.0	1.2	0.2
5.7	4.4	1.5	0.4
5.4	3.9	1.3	0.4
5.1	3.5	1.4	0.3
5.7	3.8	1.7	0.3
5.1	3.8	1.5	0.3
5.4	3.4	1.7	0.2
5.1	3.7	1.5	0.4
4.6	3.6	1.0	0.2
5.1	3.3	1.7	0.5
4.8	3.4	1.9	0.2
5.0	3.0	1.6	0.2
5.0	3.4	1.6	0.4
5.2	3.5	1.5	0.2
5.2	3.4	1.4	0.2
4.7	3.2	1.6	0.2
4.8	3.1	1.6	0.2
5.4	3.4	1.5	0.4
5.2	4.1	1.5	0.1
5.5	4.2	1.4	0.2
4.9	3.1	1.5	0.2
5.0	3.2	1.2	0.2
5.5	3.5	1.3	0.2
4.9	3.6	1.4	0.1
4.4	3.0	1.3	0.2
5.1	3.4	1.5	0.2
5.0	3.5	1.3	0.3
4.5	2.3	1.3	0.3
4.4	3.2	1.3	0.2
5.0	3.5	1.6	0.6
5.1	3.8	1.9	0.4
4.8	3.0	1.4	0.3
5.1	3.8	1.6	0.2
4.6	3.2	1.4	0.2
5.3	3.7	1.5	0.2
5.0	3.3	1.4	0.2
];
% 

% Answer is (compare with Doornik and Hansen (1994)):
% 
% MV omnibus normality test
% -------------------------------------------------------------------
% Number of variables: 4
% Sample size: 50
% -------------------------------------------------------------------
% MV omnibus test statistic: 24.414494
% Equivalent degrees of freedom: 8.000000
% P-value associated to the Royston's statistic: 0.001952
% Lower critical value associated to the Royston's statistic: 2.179731
% Upper critical value associated to the Royston's statistic: 17.534546
% With a given significance = 0.050
% Data analyzed do not have a normal distribution.
% -------------------------------------------------------------------

% Created by Norli Anida Abdullah 
%             Institute of Mathematical 
%             Faculty of Science
%             University of Malaya
%            Kuala Lumpur
%             Malaysia
%             norlienih@yahoo.com / norlie@um.edu.my
%
% Copyright. July 2, 2013.
%
% To cite this file, this would be an appropriate format:
%   NA Abdullah (2013). MVOmnibus:Multivariate Omnibus Normality Test.   
%   A MATLAB file. [WWW document]. URL http:
%
% References:
% Johnson, R.A. and Wichern, D. W. (1992). Applied Multivariate Statistical
%      Analysis. 3rd. ed. New-Jersey:Prentice Hall.
% Doornik, JA and Hansen H (1994). An omnibus test for univariate and
%      multivariate normality. Working Paper, Nu±eld College, Oxford.
% Patrick J. Farrell, Matias Salibian-Barrera, Katarzyna NaczkOn (2007). Tests for 
%      multivariate normality and associated simulation studies.
%     Journal of Statistical Computation and Simulation - J STAT COMPUT SIM 01/2007; 77(12):1065-1080.     
%

%% calculating B1 and B2

X=Y';
[p,n]=size(X);
S=cov(X',1);
V=zeros(p,p);
V(eye(size(V))~=0)=[diag(S).^(-1/2)]' ;%diagonal V=diagonal S^(-1/2)
%check
I1=V*inv(V);

C=V*S*V;

[H,deltaEig]=eig(C); %H is eigenvectors, %delta is eigenvalues

%check
I=H'*H;
delStar=H'*C*H; %compare with delta;

Xhat=zeros(p,n);
for j=1:p
    Xhat(j,:)=X(j,:)-mean(X(j,:));
end
deltaStar=(deltaEig).^(-1/2);
delta=zeros(p,p);
delta(eye(size(delta))~=0)=[diag(deltaStar)]' ;%replace diagonal of deltaStar in delta


Rstar=H*(delta)*H'*V*Xhat;
R=Rstar';
 
m2star=zeros(n,p);
m2=zeros(1,p);
m3star=zeros(n,p);
m3=zeros(1,p);
m4star=zeros(n,p);
m4=zeros(1,p);
b1=zeros(1,p);
b2=zeros(1,p);
for j=1:p
    m2star(:,j)=(R(:,j)-mean(R(:,j))).^2;
    m2(j)=sum(m2star(:,j))/n;
    
    m3star(:,j)=(R(:,j)-mean(R(:,j))).^3;
    m3(j)=sum(m3star(:,j))/n;
    
    m4star(:,j)=(R(:,j)-mean(R(:,j))).^4;
    m4(j)=sum(m4star(:,j))/n;
    
    b1(j)=m3(j)/(m2(j)^(3/2));
    b2(j)=m4(j)/(m2(j)^2);
end

B1=b1;
B2=b2;

%% calculating Z1 and Z2

beta=(3*(n^2+27*n-70)*(n+1)*(n+3))/((n-2)*(n+5)*(n+7)*(n+9));
omega2=-1+ (2*(beta-1))^(1/2);
DELTA=1/((log(sqrt(omega2)))^(1/2));
y=B1.* (((omega2-1)/2)*(((n+1)*(n+3))/(6*(n-2))))^(1/2);
for j=1:p
    Z1(j)=DELTA*log((y(j)+(y(j)^2+1)^(1/2)));
end


DELTAz2=(n-3)*(n+1)*(n^2+15*n-4);
a=((n-2)*(n+5)*(n+7)*(n^2+27*n-70))/(6*DELTAz2);
c=((n-7)*(n+5)*(n+7)*(n^2+2*n-5))/(6*DELTAz2);
k=((n+5)*(n+7)*(n^3+37*n^2+11*n-313))/(12*DELTAz2);
alpha=a+((B1.^2)*c);
chi=(B2-1-(B1.^2))*2*k;
for j=1:p
    Z2(j)=( ((chi(j)./(2.*alpha(j)))^(1/3)) -1 + (1./(9.*alpha(j))) )*((9.*alpha(j))^(1/2));
end

%% test statistics for MV omnibus test, E

Z1=Z1';
Z2=Z2';
E=Z1'*Z1+Z2'*Z2; %omnibus statistic

df = 2*p; %equivalent degrees of freedom
P = 1 - chi2cdf(E,df); %P-value

cvL=chi2inv(0.025,df);
cvU=chi2inv(0.975,df);

disp(' ')
disp('MV omnibus normality test')
disp('-------------------------------------------------------------------')
fprintf('Number of variables: %i\n', p);
fprintf('Sample size: %i\n', n);
disp('-------------------------------------------------------------------')
fprintf('MV omnibus test statistic: %3.6f\n', E);
fprintf('Equivalent degrees of freedom: %3.6f\n', df);
fprintf('P-value associated to the Royston''s statistic: %3.6f\n', P);
fprintf('Lower critical value associated to the Royston''s statistic: %3.6f\n', cvL);
fprintf('Upper critical value associated to the Royston''s statistic: %3.6f\n', cvU);
fprintf('With a given significance = %3.3f\n', 0.05);
if P >= 0.05;
    disp('Data analyzed have a normal distribution.');
else
    disp('Data analyzed do not have a normal distribution.');
end
disp('-------------------------------------------------------------------')

end
% this is for China data

%least squares model linear regression
clear;
clc;
% x here represents the days among the 1-month interval
x= [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31];
y= [323 470 655 889 1128 1701 2036 2502 3089 3858 4636 5883 7375 9172 10149 12462 15113 17660 21157 24747 27980 31506 35713 41035 47021 53578 59138 63927 69176 74386 80589]; %y here as the total cases
%y=[11 12 17 21 29 41 52 79 107 148 197 233 366 463 631 827 1016 1266 1441 1809 2158 2503 2978 3405 4032 4825 5476 6077 6820 7503 8215]; %y here as the total deaths
n= 1;
%n here represents the polynomial order
X= ones(length(x),n+1);
%first I created a matrix X of length of x values rows and n+1 columns to
%match with the equation for deriving {a} in the proof explained in
%attached report.
for i=1:n
    X(:,i+1)=(x.^1)';
end
Y=y';
A=inv(X'*X)*(X'*Y);
%the A matrix is for getting the curve fitting function coefficients and
%the equation for deriving it is proved in the attached report.
for i=2:n+1
    yy(i,:)=A(i,1).*(x.^(i-1));
end
yy(1,:)=A(1,1).*ones(1,length(x));
yy=sum(yy);
plot(x,y,'r*',x,yy,'g');
%I plotted the real x and y values with the predicted y values on the
%fiting curve
%%
%least squares model quadratic polynomial 
clear;
clc;
% x here represents the days among the 1-month interval
x= [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31];
%y= [323 470 655 889 1128 1701 2036 2502 3089 3858 4636 5883 7375 9172 10149 12462 15113 17660 21157 24747 27980 31506 35713 41035 47021 53578 59138 63927 69176 74386 80589]; %y here as the total cases
y=[11 12 17 21 29 41 52 79 107 148 197 233 366 463 631 827 1016 1266 1441 1809 2158 2503 2978 3405 4032 4825 5476 6077 6820 7503 8215]; %y here as the total deaths
n= 2;
%n here represents the polynomial order
X= ones(length(x),n+1);
%first I created a matrix X of length of x values rows and n+1 columns to
%match with the equation for deriving {a} in the proof explained in
%attached report.
for i=1:n
    X(:,i+1)=(x.^1)';
end
Y=y';
A=inv(X'*X)*(X'*Y);
%the A matrix is for getting the curve fitting function coefficients and
%the equation for deriving it is proved in the attached report.
for i=2:n+1
    yy(i,:)=A(i,1).*(x.^(i-1));
end
yy(1,:)=A(1,1).*ones(1,length(x));
yy=sum(yy);
plot(x,y,'r*',x,yy,'g');
%I plotted the real x and y values with the predicted y values on the
%fiting curve
%%
%least squares model cubic polynomial 
clear;
clc;
% x here represents the days among the 1-month interval
x= [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31];
%y= [323 470 655 889 1128 1701 2036 2502 3089 3858 4636 5883 7375 9172 10149 12462 15113 17660 21157 24747 27980 31506 35713 41035 47021 53578 59138 63927 69176 74386 80589]; %y here as the total cases
y=[11 12 17 21 29 41 52 79 107 148 197 233 366 463 631 827 1016 1266 1441 1809 2158 2503 2978 3405 4032 4825 5476 6077 6820 7503 8215]; %y here as the total deaths
n= 3;
%n here represents the polynomial order
X= ones(length(x),n+1);
%first I created a matrix X of length of x values rows and n+1 columns to
%match with the equation for deriving {a} in the proof explained in
%attached report.
for i=1:n
    X(:,i+1)=(x.^1)';
end
Y=y';
A=inv(X'*X)*(X'*Y);
%the A matrix is for getting the curve fitting function coefficients and
%the equation for deriving it is proved in the attached report.
for i=2:n+1
    yy(i,:)=A(i,1).*(x.^(i-1));
end
yy(1,:)=A(1,1).*ones(1,length(x));
yy=sum(yy);
plot(x,y,'r*',x,yy,'g');
%I plotted the real x and y values with the predicted y values on the
%fiting curve
%%
%Interpolation model Newton's divided difference 
clear;
clc;
% x here represents the days among the 1-month interval
x1= [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31];
%y1= [323 470 655 889 1128 1701 2036 2502 3089 3858 4636 5883 7375 9172 10149 12462 15113 17660 21157 24747 27980 31506 35713 41035 47021 53578 59138 63927 69176 74386 80589]; %y here as the total cases
y1=[11 12 17 21 29 41 52 79 107 148 197 233 366 463 631 827 1016 1266 1441 1809 2158 2503 2978 3405 4032 4825 5476 6077 6820 7503 8215]; %y here as the total deaths
%TDD = divdiff( x1, y1); %divided difference table 
n=length(x1); %length of the x1
a=zeros(n,n); %create a square zero matrix that has the same length as x 
a(:,1)=x1'; %add the value of x1 in the first column
for i=2:n %loope over the rest of the zro matrix
    a(:,i)=circshift(x1',(i-1)); %shift the data of x1 by i-1 positions 
end
a=a';
syms x p l_k %create variables
for k=1:n %loop over the length for x1 (y1)
    xk=x1(1,k);
    yk=y1(1,k);
    xm=a(2:n,k); %the next element after the xk
    l_k(1,1)=(x-xm(1,1))/(xk-xm(1,1)); %calculate the first divided difference notation
    for m=2:n-1 %loop to calculate n-1th divided difference notation
        M=(x-xm(m,1))/(xk-xm(m,1));
        l_k(1,m)=M;
    end
    p(k,1)=yk.*(prod(l_k)); %calulate the newton polynomial interpolation formula for these data
end
pp=simplify(sum(p))%simplify the formula
[bx,lx]=sort(x1); %sort the data of x1 incase they are not ordered
for i=1:n %loop to sort the values of y to sutisfy their respective x value after sorting
    by(1,i)=y1(1,lx(1,i));
end
xx=[1:35]; % the range where we will test the formula against
syms x %creating x variable 
c = sym2poly(pp) %extract the coefficient of the newton polynomial equation with the highest degree stated first
plot(xx,polyval(c,xx),'r--',bx,by,'b*')
%plot the predicted values from the newton interpolationpolynomial equation
%in red and the original set of point in blue dots
val=[1:31];
residuals=norm(polyval(c,val)-y1,2) %calculate the l2 norm for the raw residuals
%%
%Lagrange Interpolation method
clear;
clc;
% x here represents the days among the 1-month interval
x= [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31];
y= [323 470 655 889 1128 1701 2036 2502 3089 3858 4636 5883 7375 9172 10149 12462 15113 17660 21157 24747 27980 31506 35713 41035 47021 53578 59138 63927 69176 74386 80589]; %y here as the total cases
%y=[11 12 17 21 29 41 52 79 107 148 197 233 366 463 631 827 1016 1266 1441 1809 2158 2503 2978 3405 4032 4825 5476 6077 6820 7503 8215]; %y here as the total deaths
[P,R,S] = lagrangepoly(x,y);% Lagrange interpolation polynomial fitting a set of points funtion where X and Y are row vectors defining a set of N points uses Lagrange's method to find 
%   the N-1th order polynomial in X that passes through these points. 
%P returns the N coefficients defining the polynomial(highest order first).
% R returns the x-coordinates of the N-1 extrema of the resulting polynomial (roots of its derivative),
% and S returns the y-values  at those extrema. 
xx = 1:35;
plot(xx,polyval(P,xx),x,y,'or');%plot the lagrange polynomial equation in blue and the original point in red circles
grid;
residuals=norm(S-y,2) %calculate the l2 norm for the raw residuals


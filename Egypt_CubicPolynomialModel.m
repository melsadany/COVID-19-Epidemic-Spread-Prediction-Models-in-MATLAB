%Egypt 

%least squares model cubic polynomial 
clear;
clc;
% x here represents the days among the 1-month interval
x= [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 61];
y= [1070 1173 1322 1450 1560 1699 1794 1939 2065 2190 2350 2505 2673 2844 3032 3144 3333 3490 3659 3891 4092 4319 4534 4782 5042 5268 5537 5895 6193 6465 6813 28615]; %y here as the total cases
%y=[71 78 85 94 103 118 135 146 159 164 178 183 196 205 224 239 250 264 276 287 294 307 317 337 359 380 392 406 415 429 436]; %y here as the total deaths
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
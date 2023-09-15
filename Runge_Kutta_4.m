function [X,Y1,Y2]=Runge_Kutta_4(h,limit,zhi)
%f是方程组[y' y"]'，h是步长，limit是区间[a b]，zhi是初值。
%k是y(y1)的斜率y2，l是y'（y2）的斜率f。
f=@(x,y1,y2) 2*(1-y1^2)*y2-y1;
n=(limit(2)-limit(1))/h;
x=limit(1);
y1=zhi(1);
y2=zhi(2);
X=[x];
Y1=[y1];
Y2=[y2];
for i=1:n
    k1=y2;
    l1=f(x,y1,y2);
    k2=y2+l1*h/2;
    l2=f(x+h/2,y1+k1*h/2,y2+l1*h/2);
    k3=y2+l2*h/2;
    l3=f(x+h/2,y1+k2*h/2,y2+l2*h/2);
    k4=y2+l3*h;
    l4=f(x+h,y1+k3*h,y2+l3*h);
    y1=y1+(k1+2*k2+2*k3+k4)*h/6;
    y2=y2+(l1+2*l2+2*l3+l4)*h/6;
    x=x+h;
    X=[X,x];
    Y1=[Y1,y1];
    Y2=[Y2,y2];
end
end




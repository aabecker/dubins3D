

syms theta1 theta2 theta4 theta5 
syms lx ly lz mx my mz nx ny nz rhox rhoy rhoz 




A1v = [cos(theta1), -sin(theta1), 0, 0; sin(theta1), cos(theta1), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
A2v = [cos(theta2), -sin(theta2), 0, 0; sin(theta2), cos(theta2), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
A3v = [cos(theta3), -sin(theta3), 0, 0; sin(theta3), cos(theta3), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
A4v = [cos(theta4), -sin(theta4), 0, 0; sin(theta4), cos(theta4), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
A5v = [cos(theta5), -sin(theta5), 0, 0; sin(theta5), cos(theta5), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];

% A1v=[cos(theta1),-sin(theta1),0,0; sin(theta1),cos(theta1)]

theta3=0;


lambdas=[0,0,1,0,0,1];

mus=[1,1,0,1,1,0];

a=[1,1,0,1,1,0];

d=[0,0,0,0,0];

A1s = [1,0,0,a(1);0,lambdas(1),-mus(1),0;0,mus(1),lambdas(1),d(1);0,0,0,1];

A2s = [1,0,0,a(2);0,lambdas(2),-mus(2),0;0,mus(2),lambdas(2),d(2);0,0,0,1];

A3s = [1,0,0,a(3);0,lambdas(3),-mus(3),0;0,mus(3),lambdas(3),d(3);0,0,0,1];

A4s = [1,0,0,a(4);0,lambdas(4),-mus(4),0;0,mus(4),lambdas(4),d(4);0,0,0,1];

A5s = [1,0,0,a(5);0,lambdas(5),-mus(5),0;0,mus(5),lambdas(5),d(5);0,0,0,1];




Ahand = [lx,mx,nx,rhox; ly, my, ny, rhoy; lz, mz, nz, rhoz; 0,0,0,1]

A1=A1v*A1s
A2=A2v*A2s
A3=A3v*A3s
A4=A4v*A4s
A5=A5v*A5s

A2vInv=inv(A2v)
A2vInv=simplify(inv(A2v)) % simplify works 

A1Inv=simplify(inv(A1)) % simplify works 

ModLHS=A2s*A3*A4*A5










syms theta1 theta2 d3 theta4 theta5 
syms lx ly lz mx my mz nx ny nz rhox rhoy rhoz 

assume([theta1,theta2,d3,theta4,theta5],'real')

assume([lx ly lz mx my mz nx ny nz rhox rhoy rhoz],'real')


theta3=0;

A1v = [cos(theta1), -sin(theta1), 0, 0; sin(theta1), cos(theta1), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
A2v = [cos(theta2), -sin(theta2), 0, 0; sin(theta2), cos(theta2), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
A3v = [cos(theta3), -sin(theta3), 0, 0; sin(theta3), cos(theta3), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
A4v = [cos(theta4), -sin(theta4), 0, 0; sin(theta4), cos(theta4), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
A5v = [cos(theta5), -sin(theta5), 0, 0; sin(theta5), cos(theta5), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];


lambdas=[0,0,1,0,0,1];

mus=[1,1,0,1,1,0];

a=[1,1,0,1,1,0];

d=[0,0,d3,0,0];

A1s = [1,0,0,a(1);0,lambdas(1),-mus(1),0;0,mus(1),lambdas(1),d(1);0,0,0,1];

A2s = [1,0,0,a(2);0,lambdas(2),-mus(2),0;0,mus(2),lambdas(2),d(2);0,0,0,1];

A3s = [1,0,0,a(3);0,lambdas(3),-mus(3),0;0,mus(3),lambdas(3),d(3);0,0,0,1];

A4s = [1,0,0,a(4);0,lambdas(4),-mus(4),0;0,mus(4),lambdas(4),d(4);0,0,0,1];

A5s = [1,0,0,a(5);0,lambdas(5),-mus(5),0;0,mus(5),lambdas(5),d(5);0,0,0,1];




Ahand = [lx,mx,nx,rhox; ly, my, ny, rhoy; lz, mz, nz, rhoz; 0,0,0,1];

A1=A1v*A1s;
A2=A2v*A2s;
A3=A3v*A3s;
A4=A4v*A4s;
A5=A5v*A5s;

A2vInv=simplify(inv(A2v)); % simplify works 

A1Inv=simplify(inv(A1)); % simplify works 


% Eq4RHS and Eq4LHS in the Mathematica 

ModLHS=A2s*A3*A4*A5;

ModRHS = A2vInv*A1Inv*Ahand;

ITLHS = ModLHS(1:3,3);

ITRHS = ModRHS(1:3,3);

PTLHS=ModLHS(1:3,4);

PTRHS=ModRHS(1:3,4);

% ---- generating Ideals ----

% so matlab assumes imaginary symbolic by default?

% PTilde Dot PTilde
PTLHSdotPTLHS=simplify(expand(dot(PTLHS,PTLHS)));

PTRHSdotPTRHS=simplify(expand(dot(PTRHS,PTRHS)));

% PTilde Dot ITilde

PTLHSdotITLHS = simplify(expand(dot(PTLHS,ITLHS)));

PTRHSdotITRHS = simplify(expand(dot(PTRHS,ITRHS)));


% PTilde Cross ITilde 

PTLHScrossITLHS = simplify(expand(cross(ITLHS, PTLHS)));


PTRHScrossITRHS = simplify(expand(cross(ITRHS, PTRHS)));

% (PTilde dot Ptilde)*ITilde - (2*PTilde dot PTilde)*PTilde

vec6LHS = (PTLHSdotPTLHS)*ITLHS - (2*PTLHSdotITLHS)*PTLHS;

vec6RHS = (PTRHSdotPTRHS)*ITRHS - (2*PTRHSdotITRHS)*PTRHS;


% testfunc = sin(theta1)+cos(theta1)
% 
% 
% testfunc2=convertTrigToVar(testfunc)


% this works
% vec6LHSVar=convertTrigToVar(vec6LHS)

% testArray=[1 1 1; 2 2 2; 3 3 3]
% fprintf('length')
% length(testArray)
% 
% testArrayp=scratchFunc(testArray)


% ---- 14 equations ---- 
% PTLHS;
eq1LHS=PTLHS(1);
eq2LHS=PTLHS(2);
eq3LHS=PTLHS(3);
% PTRHS;
eq1RHS=PTRHS(1);
eq2RHS=PTRHS(2);
eq3RHS=PTRHS(3);
% ITLHS;
eq4LHS=ITLHS(1);
eq5LHS=ITLHS(2);
eq6LHS=ITLHS(3);
% ITRHS
eq4RHS=ITRHS(1);
eq5RHS=ITRHS(2);
eq6RHS=ITRHS(3);
%
eq7LHS = PTLHSdotPTLHS;
eq7RHS = PTRHSdotPTRHS;
%
eq8LHS = PTLHSdotITLHS;
eq8RHS = PTRHSdotITRHS;
%
eq9LHS = PTLHScrossITLHS(1);
eq10LHS= PTLHScrossITLHS(2);
eq11LHS= PTLHScrossITLHS(3);
%
eq9RHS = PTRHScrossITRHS(1);
eq10RHS = PTRHScrossITRHS(2);
eq11RHS = PTRHScrossITRHS(3);
%
eq12LHS=vec6LHS(1);
eq13LHS=vec6LHS(2);
eq14LHS=vec6LHS(3);
%
eq12RHS=vec6RHS(1);
eq13RHS=vec6RHS(2);
eq14RHS=vec6RHS(3);
%

eqsLHSAll=simplify(expand([eq1LHS;eq2LHS;eq3LHS;eq4LHS;eq5LHS;eq6LHS;
    eq7LHS;eq8LHS;eq9LHS;eq10LHS;eq11LHS;eq12LHS;eq13LHS;eq14LHS]));

eqsRHSAll=simplify(expand([eq1RHS;eq2RHS;eq3RHS;eq4RHS;eq5RHS;eq6RHS;
    eq7RHS;eq8RHS;eq9RHS;eq10RHS;eq11RHS;eq12RHS;eq13RHS;eq14RHS]));

% --- converting to uni ---- % 


% % conversion functions seem to be working
% fprintf('vec6LHS')
% vec6LHS
% vec6LHSVar=convertTrigToVar(vec6LHS)
% 
% vec6LHSVarToTrig=convertVarToTrig(vec6LHSVar)

eqsLHSAllVar=convertTrigToVar(expand(eqsLHSAll));
eqsLHSAllUni=convertLHSToUni(eqsLHSAllVar);


eqsLHSAll;
eqsLHSAllVar;
eqsLHSAllUni;

% next needs to convert the RHS to var then Uni the var conversion can use
% the same function but the uni needs a new one 


eqsRHSAllVar=convertTrigToVar(expand(eqsRHSAll));
eqsRHSAllUni=convertRHSToUni(eqsRHSAllVar);

eqsRHSAllVar;
eqsRHSAllUni;

% fprintf('eqsLHSAllUni(7)')
% eqsLHSAllUni(7)


% looks like this is working 
% needs the QMat side, and then it still needs the constants moved over to
% the PMat side 
% PMatNotAdjusted=get9PMatTerms(eqsLHSAllUni)

eqsRHSAllUni
display('length eqsRHSAllUni')
length(eqsRHSAllUni)

QMatNotAdjusted=get9QMatTerms(eqsRHSAllUni)






















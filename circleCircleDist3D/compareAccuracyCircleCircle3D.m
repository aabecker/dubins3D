function [dists,timeReq] = compareAccuracyCircleCircle3D
% tests the functions  distanceCircleCircle3D(N0, r0, C0, N1, r1, C1 )
%
% based on https://www.geometrictools.com/Documentation/DistanceToCircle3.pdf
%
% Aaron Becker, atbecker@uh.edu
%  draws the circles, their closest points, and (optionally) the minimum
%  radius tori that touch.

format compact
showTorus = true;

nSamples = 100; %number of circles to randomly generate
%generate nSamples of circles
N0 = 10*(rand(nSamples, 3)-1/2);
N1 = 10*(rand(nSamples, 3)-1/2);
C0 = 10*(rand(nSamples, 3)-1/2);
C1 = 10*(rand(nSamples, 3)-1/2);
r0 = 3*rand(nSamples, 1);
r1 = 3*rand(nSamples, 1);
for c = 1:nSamples
    N0(c,:) = N0(c,:)/norm(N0(c,:));
    N1(c,:) = N1(c,:)/norm(N1(c,:));
end


dists = zeros(2,nSamples);
timeReq = [0,0];

%solve using each method

tic
for c = 1:nSamples
    res = distanceCircleCircle3D(N0(c,:), r0(c), C0(c,:), N1(c,:), r1(c), C1(c,:) );
    dists(1,c) = res.distance;
end
timeReq(1) = toc;

tic
for c = 1:nSamples
    res = fminConDistanceCircleCircle3D(N0(c,:), r0(c), C0(c,:), N1(c,:), r1(c), C1(c,:) );
    dists(2,c) = res.distance;
end
timeReq(2) = toc;
plot( [0,max(dists(:))], [0,max(dists(:))], 'r-', ...
    dists(1,:),dists(2,:),'b.')
title(['8th order solver took ',num2str(timeReq(1)),'s, fmincon ',num2str(timeReq(2)),'s, for ',num2str(nSamples),' trials'] )
%fmincon 76x slower, doesn't get the best answer, and only gets one answer
axis tight
axis equal
xlabel('8th order solver');
ylabel('fmincon')




function computeAnddrawSolutions( figNum, N0, r0, C0, N1, r1, C1  )


if numel(figNum) == 1
    figure(figNum); clf;
    view(45,45)
else
    subplot(figNum(1),figNum(2),figNum(3));cla;
end
plotCircle(N0, C0,r0, 'r', 1)
plotCircle(N1, C1,r1, 'b', 1)
hold on

tic
result = distanceCircleCircle3D(N0, r0, C0, N1, r1, C1 )
toc
tic
result = fminConDistanceCircleCircle3D(N0, r0, C0, N1, r1, C1 )
toc
%result
if showTorus
    drawTorus(N0, r0, C0, result.distance/2,'r')
    drawTorus(N1, r1, C1, result.distance/2,'b')
end
axis equal

title(strcat('Min dist = ',num2str(result.distance, 3),", ",num2str(result.numClosestPairs), ' pairs, equidistant =  ',num2str(result.equidistant)))
for i = 1:result.numClosestPairs
    plot3([result.circle0Closest(i,1),result.circle1Closest(i,1)],...
        [result.circle0Closest(i,2),result.circle1Closest(i,2)],...
        [result.circle0Closest(i,3),result.circle1Closest(i,3)],'.-','color','k','LineWidth',1,'MarkerSize',12)

end

hold off
end

function plotCircle(n_vec, PtCenter,r, colorLine, lineWidth)
% plots an arc of length theta about the circle with normal vector
% n_vec, radius r, and center PtCenter, starting at position PtLeave, in
% direction dir (-1 for CW, 1 for CCW).  The line is drawn in color
% colorLine with width lineWidth
n_vec = n_vec/norm(n_vec); %ensure it is a unit vector
theta = linspace(0,2*pi,40)';  % Draw the arc path
[U,V,~] = ComputeOrthogonalComplement(n_vec);
arc2Pts = PtCenter + U*r.*cos(theta) + V*r.*sin(theta);
hold on
harc = plot3( arc2Pts(:,1),arc2Pts(:,2),arc2Pts(:,3)); set(harc,'color' ,colorLine,'linewidth',lineWidth)
plot3(PtCenter(1)+[0,n_vec(1)],PtCenter(2)+[0,n_vec(2)],PtCenter(3)+[0,n_vec(3)],'.-','color',colorLine,'MarkerSize',10,'linewidth',lineWidth);
xlabel('x')
ylabel('y')
zlabel('z')
hold off
end

function drawTorus(N, r, C, R,color)
% draws a torus!
N = N/norm(N); %normalize the matrix
[theta,phi] = meshgrid(linspace(0,2*pi,20));

x = (r + R*cos(theta)).*cos(phi);
y = (r + R*cos(theta)).*sin(phi);
z = R*sin(theta);

if N(2)==0 && N(1)==0 % the longitude
    angRz = 0;
else
    angRz = atan2(N(2), N(1));
end
angRz = angRz-pi/2;

if N(3)==1  % the latitude
    angRx = 0;
elseif N(3)==-1
    angRx = pi;
else
    angRx = atan2(N(3),sqrt(N(1)^2+N(2)^2) )-pi/2;
end
% compute the rotation matrices.
Rx = [1 0 0; 0 cos(angRx) -sin(angRx); 0 sin(angRx) cos(angRx)];
Rz = [cos(angRz) -sin(angRz) 0; sin(angRz) cos(angRz) 0; 0 0 1];

xF = x; yF = y; zF = z;
for i = 1:numel(x)
    nCoord = Rz*Rx*[x(i);y(i);z(i)];
    xF(i) = nCoord(1); yF(i)= nCoord(2); zF(i)= nCoord(3);
end
%xyzR = Rx*Rz*[x;y;z];
%coordFinal = subs(xyzR, t, pi/4);

sPlot = surf(C(1)+ xF,C(2)+yF,C(3)+zF);
%shading interp
set(sPlot, 'facecolor',color,'EdgeColor', 'none')
alpha(sPlot,0.1);
end
function [U,V,W] = ComputeOrthogonalComplement(W)
% Robustly compute a right-handed orthonormal set { U, V, W }.
% The vector W is guaranteed to be unit-length, in which case
% there is no need to worry about a division by zero when
% computing invLength.
if (abs(W(1)) > abs(W(2)))
    % The component of maximum absolute value is either W[1]
    % or W[3].
    invLength = 1 / sqrt(W(1) * W(1) + W(3) * W(3));
    U = [ -W(3) * invLength, 0, +W(1) * invLength ];
else
    % The component of maximum absolute value is either W[1]
    % or W[2].
    invLength = 1 / sqrt(W(2) * W(2) + W(3) * W(3));
    U = [0, +W(3) * invLength, -W(2) * invLength ];
end
V = cross(W, U);
end


function testTorus( N0, r0, C0  ) %#ok<DEFNU>
% used to test plotting of the circle and torus in 3D.
figure(1)
clf;

plotCircle(N0, C0,r0, 'r', 1)
hold on
drawTorus(N0, r0, C0, r0/2,'r')

axis equal
hold off
end
end
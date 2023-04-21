function result = fminConDistanceCircleCircle3D(N0, r0, C0, N1, r1, C1 )
%  Numerically computes the shortest distances in 3D between two circles using fmincon.
% Requires the Optimization Toolbox
%
%  Inputs:  two circles circle 0, and circle 1, each described by 3
%  parameters:
%       N0: circle 0 normal vector (3x1)
%       r0: circle 0 radius
%       C0: circle 0 center position (3x1)
%       N1: circle 1 normal vector (3x1)
%       r1: circle 1 radius
%       C1: circle 1 center position (3x1)
%
%  Outputs:
%      result, a struct that includes the elements:
%       distance: shortest distance (scalar)
%       sqrDistance: squared shortest distance (scalar)
%       numClosestPairs: how many closest pairs on the circles exist ,
%       circle0Closest: coordinates of closest point on circle 0, 3xnumClosestPairs
%       circle1Closest: coordinates of closest point on circle 1 , 3xnumClosestPairs
%       equidistant: if infinite solutions exist
%
% Written by Aaron Becker, atbecker@uh.edu
N0 = N0/norm(N0);  % ENSURE N0 and N1 are unit vectors (Aaron added)
N1 = N1/norm(N1);

% X = fmincon(FUN,X0,A,B) starts at X0 and finds a minimum X to the
% function FUN, subject to the linear inequalities A*X <= B. FUN accepts
% input X and returns a scalar function value F evaluated at X. X0 may be
% a scalar, vector, or matrix.
X0 = [0;0];  %initial guess for the two angles around circle 0 and circle 1
A = [1,0;
    0,1;
    -1,0;
    0,-1];
B = [pi;pi;pi;pi];
[U0,V0,~] = ComputeOrthogonalComplement(N0);
[U1,V1,~] = ComputeOrthogonalComplement(N1);


X = fmincon(@(x) distSq(x,U0,V0, r0, C0, U1,V1, r1, C1  ),X0,A,B);

result = distResult(X,U0,V0, r0, C0, U1,V1, r1, C1  );

%%% DEBUGGING plots and outputs
% disp(X)
%
% figure(40); clf
% plotCircle(N0, C0,r0, 'r', 2)
% hold on
% plot3(result.circle0Closest(1), result.circle0Closest(2), result.circle0Closest(3), 'm.','markersize',14)
% axis equal
%
% plotCircle(N1, C1,r1, 'b', 2)
% hold on
% plot3(result.circle1Closest(1), result.circle1Closest(2), result.circle1Closest(3), 'c.','markersize',14)


    function dSq = distSq(x,U0,V0, r0, C0, U1,V1, r1, C1  )
        %x = [angCirc1, angCirc2]
        x0 = x(1);  x1 = x(2);
        P0 = circPt(x0,U0,V0,r0,C0);
        P1 = circPt(x1,U1,V1,r1,C1);
        dSq = sum((P1-P0).^2);
    end

    function result = distResult(x,U0,V0, r0, C0, U1,V1, r1, C1  )
        %x = [angCirc1, angCirc2]
        x0 = x(1);  x1 = x(2);
        P0 = circPt(x0,U0,V0,r0,C0);
        P1 = circPt(x1,U1,V1,r1,C1);
        result.sqrDistance = sum((P1-P0).^2);
        result.circle0Closest = P0;
        result.circle1Closest = P1;
        result.numClosestPairs = 1;
        result.distance = sqrt(result.sqrDistance);
        result.equidistant = -1;
    end

    function P = circPt(x,U,V,r,C)
        P = C + r*(U*cos(x) + V*sin(x));
    end

    function arc2Pts = circlePts(n_vec, PtCenter,r, nPts)
        % plots an arc of length theta about the circle with normal vector
        % n_vec, radius r, and center PtCenter, starting at position PtLeave, in
        % direction dir (-1 for CW, 1 for CCW).  The line is drawn in color
        % colorLine with width lineWidth
        n_vec = n_vec/norm(n_vec); %ensure it is a unit vector
        theta = linspace(0,2*pi,nPts)';  % Draw the arc path
        [U,V,~] = ComputeOrthogonalComplement(n_vec);
        arc2Pts = PtCenter + U*r.*cos(theta) + V*r.*sin(theta);
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

end
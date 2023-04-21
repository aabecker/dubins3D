function result = distanceCircleCircle3D(N0, r0, C0, N1, r1, C1 )
%  Computes the shortest distances in 3D between two circles.
%  Inputs:  two circles circle 0, and circle 1, each described by 3
%  parameters.
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
%  April 10, 2023:
%  Transcribed from C++ to Matlab by Aaron Becker atbecker@uh.edu, 
%  and Victor Montano, based on
%  the 3D circle-circle distance algorithm is described in
%  https://www.geometrictools.com/Documentation/DistanceToCircle3.pdf
%  and available at https://www.geometrictools.com/GTE/Mathematics/DistCircle3Circle3.h
%  The notation used in the code matches that of the document.
%  David Eberly, Geometric Tools, Redmond WA 98052
%  Copyright (c) 1998-2023
%  Distributed under the Boost Software License, Version 1.0.
%  https://www.boost.org/LICENSE_1_0.txt
%  https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
%  Version: 6.0.2022.01.06

PRECISION = 1e-6;  %values less than PRECISION are ignored (numerical error).
% set up the struct result
result.distance = NaN;
result.sqrDistance = NaN;
result.numClosestPairs = 0;
result.circle0Closest = [];
result.circle1Closest=[];
result.equidistant = false;

zero = 0;   % our value for 0.  Kept to match source.
vzero = [0,0,0]; % 3x1 zero vector.
D = C1 - C0;  % distance vector from 0 to 1.
N0 = N0/norm(N0);  % ENSURE N0 and N1 are unit vectors 
N1 = N1/norm(N1);

N0xN1 = cross(N0, N1);

if ~isequal(N0xN1,vzero)
    one = 1;
    two = 2;
    % Get parameters for constructing the degree-8 polynomial phi.
    r0sqr = r0 * r0;
    r1sqr = r1 * r1;

    %Compute U1 and V1 for the plane of circle1.
    [U1,V1] = ComputeOrthogonalComplement(N1);
    % Construct the polynomial phi(cos(theta)).
    N0xD  = cross(N0, D);
    N0xU1 = cross(N0, U1);
    N0xV1 = cross(N0, V1);
    a0 = r1 * dot(D, U1);
    a1 = r1 * dot(D, V1);
    a2 = dot(N0xD, N0xD);
    a3 = r1 * dot(N0xD, N0xU1);
    a4 = r1 * dot(N0xD, N0xV1);
    a5 = r1sqr * dot(N0xU1, N0xU1);
    a6 = r1sqr * dot(N0xU1, N0xV1);
    a7 = r1sqr * dot(N0xV1, N0xV1);
    % polynomials are defined in https://www.geometrictools.com/GTE/Mathematics/Polynomial1.h
    % This ordering is backwards of Matlab's, so fliplr() is used to make
    % them the same.
    % Matlab puts the highest power of x first, so p = [3,2,1] means 3*x^2 + 2*x+1
    p0 = fliplr([a2 + a7, two * a3, a5 - a7]) ; % coefficients of a polygon.  In Matlab, poly([r]) takes as input the roots of the poly and returns the coefficients
    p1 = fliplr([two * a4, two * a6 ] );
    p2 = fliplr([ zero, a1]) ;
    p3 =  -a0 ;
    p4 = fliplr([ -a6, a4, two * a6]) ;
    p5 = fliplr([ -a3, a7 - a5 ]);
    tmp0 = fliplr( [ one, zero, -one ]);
    % w = conv(u,v) returns the convolution of vectors u and v. If u and v are vectors of polynomial coefficients, convolving them is equivalent to multiplying the two polynomials.
    tmp1 = conv(p2, p2) + conv(conv(tmp0, p3),  p3);
    tmp2 = two * conv(p2, p3);
    tmp3 = conv(p4, p4) + conv(conv(tmp0, p5), p5);
    tmp4 = two * conv( p4, p5);
    p6 = conv(p0, tmp1) + conv(conv(tmp0 , p1), tmp2) -  r0sqr * tmp3;
    p7 = conv(p0, tmp2) + conv(p1, tmp1) - r0sqr * tmp4;

    %  Parameters for polynomial root finding. The roots[] array
    %  stores the roots. We need only the unique ones, which is
    %  the responsibility of the set uniqueRoots. The pairs[]
    %  array stores the (cosine,sine) information mentioned in the
    %  PDF. 
    % The following variables are unused since we use Matlab's built
    % in roots() to find the roots of the polynomial, and it doesn't take
    % any arguments.
    % maxIterations = 128;  
    % degree = 0;
    % numRoots = 0;
    % roots =[];
    % uniqueRoots = [];
    % temp = zero;
    % sn = zero;
    numPairs = 0;
    pairs = [];

    % if  p7.GetDegree() > 0 || p7(1) ~= zero    %https://www.geometrictools.com/GTE/Mathematics/Polynomial1.h
    % GetDegree is mCoefficient.size() - 1
    if  sum(abs(p7))>0 || p7(1) ~= zero % Matlab indexes from 1, not 0
        % H(cs,sn) = p6(cs) + sn * p7(cs)
        phi = conv(p6, p6) - conv(conv(tmp0, p7), p7);
        %degree = phi.GetDegree();  %return static_cast<uint32_t>(mCoefficient.size() - 1);
        %LogAssert(degree > 0, "Unexpected degree for phi.");
        uniqueRoots = unique([removeImaginary(roots(phi));0]); %remove duplicates, remove small imaginary components, add trivial x=0 root
        % unique is a bit slow.  Could be faster using diff https://stackoverflow.com/questions/8174578/faster-way-to-achieve-unique-in-matlab-if-assumed-1d-pre-sorted-vector
        % in 1000 calls, unique took 0.074 seconds, removeImaginary 0.028,
        % and roots() 0.193 seconds.
        %numRoots = numel(uniqueRoots);
        % for  i = 0: numRoots
        %      uniqueRoots.insert(roots[i]);
        % end
        % if numel(uniqueRoots) == 0
        %     uniqueRoots = 0;  %add the trivial root.
        % end

        for  csi = 1:numel(uniqueRoots)
            cs = uniqueRoots(csi);
            if cs >1 && cs < 1+PRECISION% numerical error
                cs = 1;
            elseif cs <-1 && cs > -1-PRECISION% numerical error
                cs = -1;
            end
            if abs(cs) <= one
                temp = polyval( p7,cs);
                if (temp ~= zero)
                    sn = -polyval(p6,cs) / temp;
                    numPairs = numPairs+1;
                    pairs(numPairs,:) = [cs, sn]; %#ok<*AGROW>
                else
                    temp = max(one - cs * cs, zero);
                    sn = sqrt(temp);
                    numPairs = numPairs+1;
                    pairs(numPairs,:) = [cs, sn];
                    if sn ~= zero
                        numPairs = numPairs+1;
                        pairs(numPairs,:) = [cs, -sn];
                    end
                end
            end
        end
    else % Circles in Nonparallel Planes but Centers on Normal Line  (p7 = 0)
        %%% the next 7 lines from the C are replaced by unique([removeImaginary(roots(p6));0]);
        % H(cs,sn) = p6(cs)
        % degree = static_cast<int32_t>(p6.GetDegree());
        % LogAssert(degree > 0, "Unexpected degree for p6.");
        % numRoots = RootsPolynomial<T>::Find(degree, &p6[0], maxIterations, roots.data());
        % for (size_t i = 0; i < numRoots; ++i)
        %     uniqueRoots.insert(roots[i]);
        % end
        uniqueRoots = unique([removeImaginary(roots(p6));0]);  % 5.3 Circles in Nonparallel Planes but Centers on Normal Line

        for csi = 1:numel(uniqueRoots)
            cs = uniqueRoots(csi);
            if abs(cs) <= one
                temp = max(one - cs * cs, zero);
                sn = sqrt(temp);
                numPairs = numPairs+1;
                pairs(numPairs,:) = [cs, sn];
                if sn ~= zero
                    numPairs = numPairs+1;
                    pairs(numPairs,:) = [cs, -sn];
                end
            end
        end
    end

    candidates = cell(numPairs,1);  % list of closest info
    for i = 1:numPairs
        %info = candidates(i);
        delta = D + r1 * (pairs(i,1) * U1 + pairs(i,2) * V1);
        info.circle1Closest = C0 + delta;
        N0dDelta = dot(N0, delta);
        lenN0xDelta = norm(cross(N0, delta));  % the length of this?
        if lenN0xDelta > 0
            diff = lenN0xDelta - r0;
            info.sqrDistance = N0dDelta * N0dDelta + diff * diff;
            delta = delta - N0dDelta * N0;
            delta = delta/norm(delta);  %Note: NORMALIZE in Matlab is not a unit norm vector. 
            info.circle0Closest = C0 + r0 * delta;
            info.equidistant = false;
        else
            r0U0 = r0 * GetOrthogonal(N0);
            diff = delta - r0U0;
            info.sqrDistance = dot(diff, diff);
            info.circle0Closest = C0 + r0U0;
            info.equidistant = true;
        end
        candidates{i} = info;
    end
    sqrDists = zeros(size(candidates));
    for i = 1:numel(candidates)
        sqrDists(i) = candidates{i}.sqrDistance;
    end
    [~,idx]=sort(sqrDists);
    candidates=candidates(idx);

    result.numClosestPairs = 1;
    result.sqrDistance = candidates{1}.sqrDistance;
    result.circle0Closest = candidates{1}.circle0Closest;
    result.circle1Closest = candidates{1}.circle1Closest;
    result.equidistant = candidates{1}.equidistant;
    if (numel(uniqueRoots) > 1 && abs(candidates{2}.sqrDistance - candidates{1}.sqrDistance) < PRECISION )
        result.numClosestPairs = 2;
        result.circle0Closest(2,:) = candidates{2}.circle0Closest;
        result.circle1Closest(2,:) = candidates{2}.circle1Closest;
    end
else
    % The planes of the circles are parallel.  Whether the planes
    % are the same or different, the problem reduces to
    % determining how two circles in the same plane are
    % separated, tangent with one circle outside the other,
    % overlapping, or one circle contained inside the other
    % circle.
    result = DoQueryParallelPlanes(N0, r0, C0, r1, C1, D);
end
result.distance = sqrt(result.sqrDistance);


% The two circles are in parallel planes where D = C1 - C0, the
% difference of circle centers.
    function result = DoQueryParallelPlanes(N0, r0, C0, r1, C1, D)  % don't need N1 since circles are parallel

        N0dD = dot(N0, D);
        normProj = N0dD * N0;
        compProj = D - normProj;
        U = compProj;
        d = norm(U);  
        U = U/d; %normalize the U  (NORMALIZE in Matlab for vectors)
        %The configuration is determined by the relative location of the
        %intervals of projection of the circles on to the D-line.
        %Circle0 projects to [-r0,r0] and circle1 projects to
        %[d-r1,d+r1].
        dmr1 = d - r1;
        if (dmr1 >= r0)  % d >= r0 + r1
            % The circles are separated (d > r0 + r1) or tangent with one
            % outside the other (d = r0 + r1).
            distance = dmr1 - r0;
            result.numClosestPairs = 1;
            result.circle0Closest = C0 + r0 * U;
            result.circle1Closest = C1 - r1 * U;
            result.equidistant = false;
        else % d < r0 + r1
            % The cases implicitly use the knowledge that d >= 0.
            dpr1 = d + r1;
            if (dpr1 <= r0)
                % Circle1 is inside circle0.
                distance = r0 - dpr1;
                result.numClosestPairs = 1;
                if (d > 0)
                    result.circle0Closest = C0 + r0 * U;
                    result.circle1Closest = C1 + r1 * U;
                    result.equidistant = false;
                else
                    % The circles are concentric, so U = (0,0,0).
                    % Construct a vector perpendicular to N0 to use for
                    % closest points.
                    U = GetOrthogonal(N0);
                    result.circle0Closest = C0 + r0 * U;
                    result.circle1Closest = C1 + r1 * U;
                    result.equidistant = true;
                end
            elseif (dmr1 <= -r0)
                % Circle0 is inside circle1.
                distance = -r0 - dmr1;
                result.numClosestPairs = 1;
                if (d > 0)
                    result.circle0Closest = C0 - r0 * U;
                    result.circle1Closest = C1 - r1 * U;
                    result.equidistant = false;
                else
                    % The circles are concentric, so U = (0,0,0).
                    % Construct a vector perpendicular to N0 to use for
                    % closest points.
                    U = GetOrthogonal(N0);
                    result.circle0Closest = C0 + r0 * U;
                    result.circle1Closest = C1 + r1 * U;
                    result.equidistant = true;
                end
            else
                % The circles are overlapping.  The two points of
                % intersection are C0 + s*(C1-C0) +/- h*Cross(N,U), where
                % s = (1 + (r0^2 - r1^2)/d^2)/2 and
                % h = sqrt(r0^2 - s^2 * d^2).
                dsqr = d * d;
                s = (1 + (r0 * r0 - r1 * r1) / dsqr) / 2;
                arg = max(r0 * r0 - dsqr * s * s, 0);
                h = sqrt(arg);
                midpoint = C0+ s * compProj;
                NxU = cross(N0, U);
                hNxU = h * NxU/norm(NxU);  % Aaron added the norm here.
                distance = 0;
                result.numClosestPairs = 2;
                result.circle0Closest(1,:) = midpoint + hNxU;
                result.circle0Closest(2,:) = midpoint - hNxU;
                result.circle1Closest(1,:) = result.circle0Closest(1,:) + normProj;
                result.circle1Closest(2,:) = result.circle0Closest(2,:) + normProj;
                result.equidistant = false;
            end
        end
        result.sqrDistance = distance * distance + N0dD * N0dD;
    end

    function U = GetOrthogonal(W)
        % generates a vector orthogonal to W in a principled manner
        [U,~,~] = ComputeOrthogonalComplement(W);
    end

    function rootsOut = removeImaginary( rootsIn)
        % due to numerical instability, sometimes we get small imaginary
        % components. This code removes the imaginary component if it is
        % less than PRECISION, and only saves real roots.
        %rootsOut = rootsIn;
        %indices = imag(rootsIn) < PRECISION;
        %rootsOut(indices) = real(rootsIn(indices));

        rootsOut = real(rootsIn(imag(rootsIn)<PRECISION));
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

end
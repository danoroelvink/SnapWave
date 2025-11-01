function[dist_to_line,xcr,ycr,scr]=get_distance_to_polyline(Xc,Yc,Xr,Yr);
% DIST_TO_LINE compute shortest distance and projected points on reference
% polyline Xr,Yr from polyline Xc,Yc

sr=[0 cumsum(hypot(diff(Xr),diff(Yr)))];
for i=1:length(Xc)
    % compute distances to all line sections of Xr, Yr
    P=[Xc(i),Yc(i)];
    B=[Xr(1),Yr(1)];
    for j=1:length(Xr)-1
        A=B;
        B=[Xr(j+1),Yr(j+1)];
        [~,~,dist(j),t] = pointToSegmentDist(P,A,B);
    end
    % select line segment with shortest absolute distance
    [distance_to_line(i),jmin] = min (abs(dist));
    % store distance (including sign) and projected points on reference
    % line
    A=[Xr(jmin),Yr(jmin)];
    B=[Xr(jmin+1),Yr(jmin+1)];
    [xcr(i),ycr(i),dist_to_line(i),t] = pointToSegmentDist (P,A,B);
    scr(i)=sr(jmin)+t*(sr(jmin+1)-sr(jmin));
end
end

function [xcr,ycr,dist,t] = pointToSegmentDist(P, A, B)
    % P: The point [Px, Py]
    % A: One endpoint of the line segment [Ax, Ay]
    % B: The other endpoint of the line segment [Bx, By]

    % Convert points to vectors
    AP = P - A;    % Vector from A to P
    AB = B - A;    % Vector from A to B

    % Squared length of segment AB
    AB_squared = dot(AB, AB);

    if AB_squared == 0
        % A and B are the same point, distance is from P to A
        dist = norm(P - A);
        xcr  = A(1);
        ycr  = A(2);
        t    = 0;
        return;
    end

    % Projection of AP onto AB, normalized by the squared length of AB
    t = dot(AP, AB) / AB_squared;

    % Clamp t to the range [0, 1] to restrict the closest point to the segment
    t = max(0, min(1, t));

    % Closest point on the segment to P
    closestPoint = A + t * AB;
    xcr=closestPoint(1);
    ycr=closestPoint(2);

    % Compute the distance from P to the closest point
    dist = norm(P - closestPoint);

    % Determine the side of the segment using the cross product
    % Cross product in 2D: z = (ABx * APy) - (ABy * APx)
    crossProductZ = (AB(1) * AP(2)) - (AB(2) * AP(1));

    % Make the distance negative if the point is on one side (e.g., left side)
    if crossProductZ < 0
        dist = -dist;
    end

    
end

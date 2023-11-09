function [a1, a2] = computeTimeParallel(a1_old, q0, q)
% a1_old is a matrix of size nex3: It contains the
% time-parallel reference director (a1) at each of
% the ne edges at the old time step
% q0 is the DOF vector at the old time step
% q is the DOF vector at the new time step
% Outputs:
% a1 is the time-parallel reference director (a1) at
% each of the ne edges at the new time step. Size of
% a1 is ne x 3.
% a2 is the second reference director. It's size is
% ne x 3.
nv = (length(q)+1)/4;
ne = nv - 1;
tangent0 = computeTangent(q0); % Tangents at old step
tangent = computeTangent(q); % Tangents at new step
a1 = zeros(ne, 3);
a2 = zeros(ne, 3);
for c=1:ne % loop over edges
t0 = tangent0(c,:); % tangent on c-th edge at old step
t = tangent(c,:); % tangent on c-th edge at new step
a1_local_old = a1_old(c,:); % a1 director on c-th edge
% at old step
a1_local = parallel_transport( a1_local_old, t0, t);
% a1_local is the first reference director on c-th
% edge at new step
% Just to be careful: enforce a1 and t are perp.
a1_local = a1_local - dot(a1_local, t) * t;
a1_local = a1_local / norm(a1_local);
a1(c,:) = a1_local; % store
a2(c,:) = cross(t, a1_local);
end
end

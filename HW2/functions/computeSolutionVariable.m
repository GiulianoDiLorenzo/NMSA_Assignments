function [sol] = computeSolutionVariable(NX, NT, dx, dt, S, u_0, u_1, lambda, f)

% This function solves the Webster's problem
%     S utt = gamma^2 (S ux)x + f
% Where S(x) = (1+2x)^2
% 
% Args:   NX              number of elements in space
%         NT              number of elements in time
%         dx              space step
%         dt              time step
%         u_0             u(x,0)
%         u_1             ut(x,0)
%         lambda          gamma dt/dx
%         f               force
% 
% Output: sol             numerical solution

disp("Computing solution with variable profile, NX = " + NX + ", NT = " + NT);

tic

% initializing the solution
sol = zeros(NX+1, NT+1);

% applying the boundary conditions
sol(2:end, 2) = u_0;                        % 3rd condition
sol(2:end, 1) = sol(2:end, 2) - dt * u_1;   % 4th condition
sol(end,:) = u_0(end);                      % 2nd condition
sol(1,:) = sol(2,:);                        % 1st condition

% extending S(x) to be coherent with the virtual lines
S_ext = [(1-2*dx)^2; S];

for n = 2:NT      % evaluating n+1, up to NT to consider virtual line

    for k = 2:NX
        sol(k,n+1) = 2*sol(k,n) - sol(k,n-1) + lambda^2/S_ext(k) * (((S_ext(k+1)-S_ext(k-1))/4 + S_ext(k)) * sol(k+1,n) - 2*S_ext(k) * sol(k,n) + ((S_ext(k-1)-S_ext(k+1))/4 + S_ext(k)) * sol(k-1,n)) + dt^2/S_ext(k) * f(k-1,n);
    end

    sol(1,n+1) = sol(2,n+1);                % 1st condition
end

% discarding virtual lines
sol = sol(2:end, 2:end);

toc

end

function dt = dt_CFL(u,v,mu,T,rho,gamma,R,Pr,dx,dy)
    % Credit: Ryan Dunn
    sos = sqrt(gamma*R*T);
    v_prime = max((4*gamma*mu.^2)./(3*Pr*rho),[],'all');
    temp = abs(u)/dx + abs(v)/dy + sqrt(1/(dx^2)+1/(dy^2))*sos + 2*v_prime*(1/(dx^2)+1/(dy^2));
    dt = min(1./temp,[],'all');
end
function [rho,u,v,T,p,e,Et] = cons2prim(U,R,Cv)
    % Credit: Ryan Dunn
    rho = squeeze(U(1,:,:));
    u = squeeze(U(2,:,:))./rho;
    v = squeeze(U(3,:,:))./rho;
    Et = squeeze(U(4,:,:));
    T = (Et./rho - 0.5*(u.^2 + v.^2))./Cv;
    p = R*rho.*T;
    e = Cv*T;
end
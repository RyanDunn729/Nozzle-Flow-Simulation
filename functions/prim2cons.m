function U = prim2cons(rho,u,v,T,Cv)

[nx,ny] = size(rho);
U = zeros(4,nx,ny);
U(1,:,:) = rho;
U(2,:,:) = rho.*u;
U(3,:,:) = rho.*v;
U(4,:,:) = rho.*(Cv.*T + 0.5*(u.^2 + v.^2));

end
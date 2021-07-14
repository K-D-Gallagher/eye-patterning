function [Nx, Ny, Nz] = GaussNormal(X,Y,Z,sigma)
%GaussNormal computes smoothened normal using Gaussian derivative
%
% GaussNormal(X,Y,Z,sigma)

    [Xz, Xphi] = GaussGradient(X, sigma);
    [Yz, Yphi] = GaussGradient(Y, sigma);
    [Zz, Zphi] = GaussGradient(Z, sigma);

    Nx = Yz.*Zphi - Zz.*Yphi;
    Ny = Zz.*Xphi - Xz.*Zphi;
    Nz = Xz.*Yphi - Yz.*Xphi;
    Nnorm = sqrt(Nx.^2 + Ny.^2 + Nz.^2);
    Nx = Nx./Nnorm;
    Ny = Ny./Nnorm;
    Nz = Nz./Nnorm;

end


% This function generates the dielectric and permeability tensors for a system that is half
% filled with a Haldane and Raghu photonic Chern insulator, and half filled with air.

function [epsMuSysTM, epsMuSysTE, epsMuSysTM_InvSqrt, epsMuSysTE_InvSqrt, Nx, Ny, dx, dy] = eps_mu_gen(ax,ay,nUCx,nUCy,nPts,rodR,epsAir,epsRod,epsAirMO,muAir,muRod,epsTrivInsulBkg,radTrivInsulBkg,phcBoundary,rodAbsorb,noise_strength,rng_seed,CHECK_EPS_PLOTS)

    %CHECK_EPS_PLOTS = 1;

    assert(muAir == 1);
    assert(muRod == 1);

    rng(rng_seed);

    Nx = round(nUCx*ax*nPts);
    Ny = round(nUCy*ay*nPts);

    % Ensure Nx is even
    Nx = 2*floor(Nx/2);
    Ny = 2*floor(Ny/2);
    dx = (nUCx*ax)/Nx;
    dy = (nUCy*ay)/Ny;

    y_uc = (-(Ny-1)/2:(Ny-1)/2)*dy;
    assert(length(y_uc)==Ny);
    x_uc = (-(Nx-1)/2:(Nx-1)/2)*dx;
    assert(length(x_uc)==Nx);
    %[X,Y] = meshgrid(x_uc,y_uc);
    [Ycoords,Xcoords] = meshgrid(y_uc,x_uc);

    % rod centers:
    xCens = zeros(2*nUCx*nUCy+nUCx,1);
    yCens = zeros(2*nUCx*nUCy+nUCx,1);
    for yy=1:nUCy
        for xx=1:nUCx
            idx1 = 2*nUCx*(yy-1) + 2*(xx-1) + 1;
            idx2 = 2*nUCx*(yy-1) + 2*(xx-1) + 2;

            xCens(idx1) = ax*(xx-1) + ax/4 - nUCx*ax/2 + ay*noise_strength*(rand(1)-0.5); % disorder scales with ay for every dim
            yCens(idx1) = ay*(yy-1) + ay/2 - nUCy*ay/2 + ay*noise_strength*(rand(1)-0.5);
            xCens(idx2) = ax*(xx-1) + 3*ax/4 - nUCx*ax/2 + ay*noise_strength*(rand(1)-0.5);
            yCens(idx2) = ay*(yy-1) - nUCy*ay/2 + ay*noise_strength*(rand(1)-0.5);
            
            if (yy==nUCy)
                idx3 = 2*nUCx*nUCy + xx;
                xCens(idx3) = ax*(xx-1) + 3*ax/4 - nUCx*ax/2 + ay*noise_strength*(rand(1)-0.5);
                yCens(idx3) = ay*(yy-1) + ay - nUCy*ay/2 + ay*noise_strength*(rand(1)-0.5);
            end
        end
    end

    % matrices for M = diag(mu, eps)
    epsDiag = epsAir*ones(Nx,Ny);
    epsOffDiag = 1i*epsAirMO*ones(Nx,Ny);
    muDiag = muAir*ones(Nx,Ny);
    muOffDiag = 1i*zeros(Nx,Ny);

    % matrices for M^-1, M^(1/2), and M^(-1/2)
    epsRodMO = 0;

    function [retVal] = approxIsEqual(mat1,mat2)
        errorTol = 1e-9;
        checkXX = mat1-mat2;
        retVal = (sum(sum(abs(checkXX)))<errorTol);
    end

    epsTensorAirPixel = [epsAir, 1i*epsAirMO; -1i*epsAirMO, epsAir];
    [VV,EE] = eig(epsTensorAirPixel);
    epsTensorAirPixelInverse = VV*diag(1./diag(EE))*inv(VV);
    epsTensorAirPixelSqrt = VV*diag(sqrt(diag(EE)))*inv(VV);
    epsTensorAirPixelInvSqrt = VV*diag(1./sqrt(diag(EE)))*inv(VV);

    assert(approxIsEqual(epsTensorAirPixel,VV*EE*inv(VV)));
    assert(approxIsEqual(epsTensorAirPixel*epsTensorAirPixelInverse, eye(2)));
    assert(approxIsEqual(epsTensorAirPixelSqrt*epsTensorAirPixelInvSqrt, eye(2)));
    assert(approxIsEqual(epsTensorAirPixelSqrt*epsTensorAirPixelSqrt,epsTensorAirPixel));
    assert(approxIsEqual(epsTensorAirPixelInvSqrt*epsTensorAirPixel*epsTensorAirPixelInvSqrt, eye(2)));

    assert(ishermitian(epsTensorAirPixelInverse));
    assert(ishermitian(epsTensorAirPixelSqrt));
    assert(ishermitian(epsTensorAirPixelInvSqrt));

    %epsTensorMatPixel = [epsRod, 1i*epsRodMO; -1i*epsRodMO, epsRod];
    %[VV,EE] = eig(epsTensorMatPixel);
    %epsTensorMatPixelInverse = VV*diag(1./diag(EE))*inv(VV);
    %epsTensorMatPixelSqrt = VV*diag(sqrt(diag(EE)))*inv(VV);
    %epsTensorMatPixelInvSqrt = VV*diag(1./sqrt(diag(EE)))*inv(VV);
    assert(epsRodMO==0);
    epsTensorMatPixel = [epsRod + 1i*rodAbsorb, 0; 0, epsRod + 1i*rodAbsorb];
    epsTensorMatPixelInverse = [1/(epsRod + 1i*rodAbsorb), 0; 0, 1/(epsRod + 1i*rodAbsorb)];
    epsTensorMatPixelSqrt = [sqrt(epsRod + 1i*rodAbsorb), 0; 0, sqrt(epsRod + 1i*rodAbsorb)];
    epsTensorMatPixelInvSqrt = [1/sqrt(epsRod + 1i*rodAbsorb), 0; 0, 1/sqrt(epsRod + 1i*rodAbsorb)];

    %assert(approxIsEqual(epsTensorMatPixel,VV*EE*inv(VV)));
    assert(approxIsEqual(epsTensorMatPixel*epsTensorMatPixelInverse, eye(2)));
    assert(approxIsEqual(epsTensorMatPixelSqrt*epsTensorMatPixelInvSqrt, eye(2)));
    assert(approxIsEqual(epsTensorMatPixelSqrt*epsTensorMatPixelSqrt,epsTensorMatPixel));
    assert(approxIsEqual(epsTensorMatPixelInvSqrt*epsTensorMatPixel*epsTensorMatPixelInvSqrt, eye(2)));

    %assert(ishermitian(epsTensorMatPixelInverse));
    %assert(ishermitian(epsTensorMatPixelSqrt));
    %assert(ishermitian(epsTensorMatPixelInvSqrt));

    epsInvDiag = epsTensorAirPixelInverse(1,1)*ones(Nx,Ny);
    epsInvOffDiag = epsTensorAirPixelInverse(1,2)*ones(Nx,Ny);
    muInvDiag = (1/muAir)*ones(Nx,Ny);
    muInvOffDiag = zeros(Nx,Ny);

    epsSqrtDiag = epsTensorAirPixelSqrt(1,1)*ones(Nx,Ny);
    epsSqrtOffDiag = epsTensorAirPixelSqrt(1,2)*ones(Nx,Ny);
    muSqrtDiag = sqrt(muAir)*ones(Nx,Ny);
    muSqrtOffDiag = zeros(Nx,Ny);

    epsInvSqrtDiag = epsTensorAirPixelInvSqrt(1,1)*ones(Nx,Ny);
    epsInvSqrtOffDiag = epsTensorAirPixelInvSqrt(1,2)*ones(Nx,Ny);
    muInvSqrtDiag = (1/sqrt(muAir))*ones(Nx,Ny);
    muInvSqrtOffDiag = zeros(Nx,Ny);

    % resetting the air:
    phcBoundary_y = phcBoundary + 2*ay*noise_strength;
    secondRegionCriteria = (Xcoords < -ax*phcBoundary/4) | (Xcoords > ax*phcBoundary/4) | (Ycoords < -ay*phcBoundary_y/2) | (Ycoords > ay*phcBoundary_y/2);

    epsDiag(secondRegionCriteria) = epsTrivInsulBkg;
    epsOffDiag(secondRegionCriteria) = 0;

    epsInvDiag(secondRegionCriteria) = 1/epsTrivInsulBkg;
    epsInvOffDiag(secondRegionCriteria) = 0;

    epsSqrtDiag(secondRegionCriteria) = sqrt(epsTrivInsulBkg);
    epsSqrtOffDiag(secondRegionCriteria) = 0;

    epsInvSqrtDiag(secondRegionCriteria) = 1/sqrt(epsTrivInsulBkg);
    epsInvSqrtOffDiag(secondRegionCriteria) = 0;

    % filling in:
    rod_noise_scaling_factor = 1/2;
    nsr = rod_noise_scaling_factor*noise_strength;

    for rr=1:length(xCens)
        if (xCens(rr) >= -ax*phcBoundary/4) && (xCens(rr) <= ax*phcBoundary/4) && (yCens(rr) >= -ay*phcBoundary_y/2) && (yCens(rr) <= ay*phcBoundary_y/2)
            rodMat = ((Xcoords-xCens(rr))/(rodR+ (nsr*(rand(1)-0.5)))).^2 + ((Ycoords-yCens(rr))/(rodR+ (nsr*(rand(1)-0.5)))).^2;

            % filling in M:
            epsDiag(rodMat<=1) = epsTensorMatPixel(1,1);
            muDiag(rodMat<=1) = muRod;
            epsOffDiag(rodMat<=1) = epsTensorMatPixel(1,2);

            % filling in M^-1
            epsInvDiag(rodMat<=1) = epsTensorMatPixelInverse(1,1);
            muInvDiag(rodMat<=1) = 1/muRod;
            epsInvOffDiag(rodMat<=1) = epsTensorMatPixelInverse(1,2);

            % filling in M^(1/2)
            epsSqrtDiag(rodMat<=1) = epsTensorMatPixelSqrt(1,1);
            muSqrtDiag(rodMat<=1) = sqrt(muRod);
            epsSqrtOffDiag(rodMat<=1) = epsTensorMatPixelSqrt(1,2);

            % filling in M^(-1/2)
            epsInvSqrtDiag(rodMat<=1) = epsTensorMatPixelInvSqrt(1,1);
            muInvSqrtDiag(rodMat<=1) = 1/sqrt(muRod);
            epsInvSqrtOffDiag(rodMat<=1) = epsTensorMatPixelInvSqrt(1,2);
        else
            rodMat = ((Xcoords-xCens(rr))/radTrivInsulBkg).^2 + ((Ycoords-yCens(rr))/radTrivInsulBkg).^2;

            % filling in M:
            epsDiag(rodMat<=1) = 1;
            muDiag(rodMat<=1) = 1;
            epsOffDiag(rodMat<=1) = 0;

            % filling in M^-1
            epsInvDiag(rodMat<=1) = 1;
            muInvDiag(rodMat<=1) = 1;
            epsInvOffDiag(rodMat<=1) = 0;

            % filling in M^(1/2)
            epsSqrtDiag(rodMat<=1) = 1;
            muSqrtDiag(rodMat<=1) = 1;
            epsSqrtOffDiag(rodMat<=1) = 0;

            % filling in M^(-1/2)
            epsInvSqrtDiag(rodMat<=1) = 1;
            muInvSqrtDiag(rodMat<=1) = 1;
            epsInvSqrtOffDiag(rodMat<=1) = 0;
        end
    end

    %% checking:
    if CHECK_EPS_PLOTS
        figure(11); imagesc(real(epsDiag.')); % can plot here.
        axis equal;
        figure(12); imagesc(imag(epsOffDiag.'));
        axis equal;
        figure(13); imagesc(imag(epsDiag.')); % can plot here.
        axis equal;
        assert(false);
    end

    %% output formatting:
    nVals = Nx*Ny;
    z0 = sparse(nVals,nVals);

    % constructing and checking M
    epsDiag = diag(sparse(reshape(epsDiag,Nx*Ny,1)));
    epsOffDiag = diag(sparse(reshape(epsOffDiag,Nx*Ny,1)));
    muDiag = diag(sparse(reshape(muDiag,Nx*Ny,1)));
    muOffDiag = diag(sparse(reshape(muOffDiag,Nx*Ny,1)));
    
    epsMuSysTM = [muDiag, muOffDiag, z0;...
                  -1*muOffDiag, muDiag, z0;...
                  z0, z0, epsDiag];
    epsMuSysTE = [muDiag, z0, z0;...
                  z0, epsDiag, epsOffDiag;...
                  z0, -1*epsOffDiag, epsDiag];

    %assert(ishermitian(epsMuSysTM));
    %assert(ishermitian(epsMuSysTE));
    assert(issparse(epsMuSysTM));
    assert(issparse(epsMuSysTE));

    % constructing and checking M^-1
    epsInvDiag = diag(sparse(reshape(epsInvDiag,Nx*Ny,1)));
    epsInvOffDiag = diag(sparse(reshape(epsInvOffDiag,Nx*Ny,1)));
    muInvDiag = diag(sparse(reshape(muInvDiag,Nx*Ny,1)));
    muInvOffDiag = diag(sparse(reshape(muInvOffDiag,Nx*Ny,1)));
    
    epsMuSysTM_Inv = [muInvDiag, muInvOffDiag, z0;...
                  -1*muInvOffDiag, muInvDiag, z0;...
                  z0, z0, diag(sparse(1./diag(epsDiag)))];
    epsMuSysTE_Inv = [muInvDiag, z0, z0;...
                  z0, epsInvDiag, epsInvOffDiag;...
                  z0, -1*epsInvOffDiag, epsInvDiag];

    assert(approxIsEqual(epsMuSysTM*epsMuSysTM_Inv, speye(3*Nx*Ny)));
    assert(approxIsEqual(epsMuSysTE*epsMuSysTE_Inv, speye(3*Nx*Ny)));
    %assert(ishermitian(epsMuSysTM_Inv));
    %assert(ishermitian(epsMuSysTE_Inv));
    assert(issparse(epsMuSysTM_Inv));
    assert(issparse(epsMuSysTE_Inv));

    % constructing and checking M^(1/2)
    epsSqrtDiag = diag(sparse(reshape(epsSqrtDiag,Nx*Ny,1)));
    epsSqrtOffDiag = diag(sparse(reshape(epsSqrtOffDiag,Nx*Ny,1)));
    muSqrtDiag = diag(sparse(reshape(muSqrtDiag,Nx*Ny,1)));
    muSqrtOffDiag = diag(sparse(reshape(muSqrtOffDiag,Nx*Ny,1)));
    
    epsMuSysTM_Sqrt = [muSqrtDiag, muSqrtOffDiag, z0;...
                  -1*muSqrtOffDiag, muSqrtDiag, z0;...
                  z0, z0, diag(sqrt(diag(epsDiag)))];
    epsMuSysTE_Sqrt = [muSqrtDiag, z0, z0;...
                  z0, epsSqrtDiag, epsSqrtOffDiag;...
                  z0, -1*epsSqrtOffDiag, epsSqrtDiag];

    assert(approxIsEqual(epsMuSysTM_Sqrt*epsMuSysTM_Sqrt, epsMuSysTM));
    %assert(approxIsEqual(epsMuSysTE_Sqrt*epsMuSysTE_Sqrt, epsMuSysTE));
    %assert(ishermitian(epsMuSysTM_Sqrt));
    %assert(ishermitian(epsMuSysTE_Sqrt));
    assert(issparse(epsMuSysTM_Sqrt));
    assert(issparse(epsMuSysTE_Sqrt));

    % constructing and checking M^(-1/2)
    epsInvSqrtDiag = diag(sparse(reshape(epsInvSqrtDiag,Nx*Ny,1)));
    epsInvSqrtOffDiag = diag(sparse(reshape(epsInvSqrtOffDiag,Nx*Ny,1)));
    muInvSqrtDiag = diag(sparse(reshape(muInvSqrtDiag,Nx*Ny,1)));
    muInvSqrtOffDiag = diag(sparse(reshape(muInvSqrtOffDiag,Nx*Ny,1)));
    
    epsMuSysTM_InvSqrt = [muInvSqrtDiag, muInvSqrtOffDiag, z0;...
                  -1*muInvSqrtOffDiag, muInvSqrtDiag, z0;...
                  z0, z0, diag(sparse(1./sqrt(diag(epsDiag))))];
    epsMuSysTE_InvSqrt = [muInvSqrtDiag, z0, z0;...
                  z0, epsInvSqrtDiag, epsInvSqrtOffDiag;...
                  z0, -1*epsInvSqrtOffDiag, epsInvSqrtDiag];

    assert(approxIsEqual(epsMuSysTM_Sqrt*epsMuSysTM_InvSqrt, speye(3*Nx*Ny)));
    assert(approxIsEqual(epsMuSysTE_Sqrt*epsMuSysTE_InvSqrt, speye(3*Nx*Ny)));
    assert(approxIsEqual(epsMuSysTM_InvSqrt*epsMuSysTM_InvSqrt, epsMuSysTM_Inv));
    assert(approxIsEqual(epsMuSysTE_InvSqrt*epsMuSysTE_InvSqrt, epsMuSysTE_Inv));
    assert(approxIsEqual(epsMuSysTM_InvSqrt*epsMuSysTM*epsMuSysTM_InvSqrt, speye(3*Nx*Ny)));
    assert(approxIsEqual(epsMuSysTE_InvSqrt*epsMuSysTE*epsMuSysTE_InvSqrt, speye(3*Nx*Ny)));
    %assert(ishermitian(epsMuSysTM_InvSqrt));
    %assert(ishermitian(epsMuSysTE_InvSqrt));
    assert(issparse(epsMuSysTM_InvSqrt));
    assert(issparse(epsMuSysTE_InvSqrt));

end
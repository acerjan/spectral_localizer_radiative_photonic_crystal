% gen_2d_yee_derivs_coords.m
% Alex Cerjan
% 2.21.2022

% Input Arguments
% =================
% Nx, Ny 1X grid size
% dx, dy 1X grid resolution
% BC    [xhi yhi xlo ylo] boundary conditions
%       -2: pseudo-periodic (requires kinc)
%       -1: periodic
%        0: Dirichlet
% kx, ky incident wave vector
%       This argument is only needed for pseudo-periodic boundaries.
%
% Note: For normalized grids, use dx=k0*dx and kinc=kinc/k0 ????
%
% Output Arguments
% =================
%
%          E(i+1,j) - E(i,j)                     E(i,j+1) - E(i,j)
% dEdx*E = -----------------           dEdy*E = -------------------
%                  dx                                   dy
%
%          H(i,j) - H(i-1,j)                   H(i,j) - H(i,j-1)
% dHdx*H = -----------------         dHdy*H = -------------------
%                  dx                                  dy
%
% Thus, if we have site (i,j), we need to implement the boundary conditions
% on the "upper" sides for the E-fields, and on the "lower" sides for the H-fields

% XX, YY are the x and y coordinates of the field vectors.
% for TM waves, these are arranged in the order (Hx, Hy, Ez)
% for TE waves, these are arranged in the order (Hz, Ex, Ey)

% In 2D the grids look like:
% TM:
%   |
%   |
% Hx(x_i,y_j+1/2)
%   |
%   |
% Ez(x_i,y_j) -------- Hy(x_i+1/2,y_j)

% TE:
%   |
%   |
% Ey(x_i,y_j+1/2)      Hz(x_i+1/2,y_j+1/2)
%   |
%   |
%   ------------------ Ex(x_i+1/2,y_j)

% To realize an absorbing boundary condition, we will utilize a stretched coordinate PML
% this means that, inside the PML
% d/dx --> (1/(1+ i sigma_x(x)/omega)) * (d/dx)

function [WW,XX,YY] = gen_2d_yee_derivs_coords(nPtsX,nPtsY,dx,dy,kx,ky,NxAbsorb,NyAbsorb,sigma0,omegaIn,BC,polz,CHECK_ABS_REGION)

    assert(length(BC)==2);

    nVals = nPtsX*nPtsY;

    sigmaPow = 3;

    function [sXout] = sX(xCur)
        
        sXout = 1;
        if (xCur <= NxAbsorb) %|| (xCur > (nPtsX-NxAbsorb))
            normL = (NxAbsorb - (xCur-1))/NxAbsorb;
            sXout = (1 + 1i*(sigma0/omegaIn)*(normL)^sigmaPow);
        end
        if (xCur > (nPtsX-NxAbsorb))
            normL = (NxAbsorb - (nPtsX-xCur))/NxAbsorb;
            sXout = (1 + 1i*(sigma0/omegaIn)*(normL)^sigmaPow);
        end

    end

    function [sYout] = sY(yCur)
        
        sYout = 1;
        if (yCur <= NyAbsorb) %|| (yCur > (nPtsY-NyAbsorb))
            normL = (NyAbsorb - (yCur-1))/NyAbsorb;
            sYout = (1 + 1i*(sigma0/omegaIn)*(normL)^sigmaPow);
        end
        if (yCur > (nPtsY-NyAbsorb))
            normL = (NyAbsorb - (nPtsY-yCur))/NyAbsorb;
            sYout = (1 + 1i*(sigma0/omegaIn)*(normL)^sigmaPow);
        end

    end

    %% generate wave operator:
    x1vec = zeros(2*nVals,1);
    y1vec = zeros(2*nVals,1);
    dedxVec = zeros(2*nVals,1);
    ii1=1;

    x2vec = zeros(2*nVals,1);
    y2vec = zeros(2*nVals,1);
    dedyVec = zeros(2*nVals,1);
    ii2=1;

    x3vec = zeros(2*nVals,1);
    y3vec = zeros(2*nVals,1);
    dhdxVec = zeros(2*nVals,1);
    ii3=1;

    x4vec = zeros(2*nVals,1);
    y4vec = zeros(2*nVals,1);
    dhdyVec = zeros(2*nVals,1);
    ii4=1;

    for yy=1:nPtsY
        for xx=1:nPtsX
            idxCur = nPtsX*(yy-1) + xx;
            idxNextX = nPtsX*(yy-1) + xx + 1;
            idxPrevX = nPtsX*(yy-1) + xx - 1;
            idxNextY = nPtsX*(yy-0) + xx;
            idxPrevY = nPtsX*(yy-2) + xx;

            % cur positions for PML:
            switch polz
            case 'TM' % f1 = Hx -- f2 = Hy -- f3 = Ez
                f1x = (xx-1); f1y = (yy-1) + 1/2;
                f2x = (xx-1) + 1/2; f2y = (yy-1);
                f3x = (xx-1); f3y = (yy-1);
            case 'TE' % f1 = Hz -- f2 = Ex -- f3 = Ey
                f1x = (xx-1) + 1/2; f1y = (yy-1) + 1/2;
                f2x = (xx-1) + 1/2; f2y = (yy-1);
                f3x = (xx-1); f3y = (yy-1) + 1/2;
            end

            % on-site:
            switch polz
            case 'TM'
                x1vec(ii1)=idxCur; y1vec(ii1)=idxCur; dedxVec(ii1)=-1/(sX(f2x)); ii1=ii1+1;
                x2vec(ii2)=idxCur; y2vec(ii2)=idxCur; dedyVec(ii2)=-1/(sY(f1y)); ii2=ii2+1;
                x3vec(ii3)=idxCur; y3vec(ii3)=idxCur; dhdxVec(ii3)=+1/(sX(f3x)); ii3=ii3+1;
                x4vec(ii4)=idxCur; y4vec(ii4)=idxCur; dhdyVec(ii4)=+1/(sY(f3y)); ii4=ii4+1;
            case 'TE'
                x1vec(ii1)=idxCur; y1vec(ii1)=idxCur; dedxVec(ii1)=-1/(sX(f1x)); ii1=ii1+1;
                x2vec(ii2)=idxCur; y2vec(ii2)=idxCur; dedyVec(ii2)=-1/(sY(f1y)); ii2=ii2+1;
                x3vec(ii3)=idxCur; y3vec(ii3)=idxCur; dhdxVec(ii3)=+1/(sX(f3x)); ii3=ii3+1;
                x4vec(ii4)=idxCur; y4vec(ii4)=idxCur; dhdyVec(ii4)=+1/(sY(f2y)); ii4=ii4+1;
            end

            % low-x:
            if xx==1
                switch BC(1)
                case -2
                    x3vec(ii3)=idxCur; y3vec(ii3)=nPtsX*(yy-1)+nPtsX; dhdxVec(ii3)=-exp(-1i*kx*nPtsX*dx); ii3=ii3+1;
                case -1
                    x3vec(ii3)=idxCur; y3vec(ii3)=nPtsX*(yy-1)+nPtsX; dhdxVec(ii3)=-1; ii3=ii3+1;
                case 0
                    % do nothing. = 0.
                otherwise
                    error('Unrecognized x-low boundary condition.');
                end
            else
                switch polz
                case 'TM'
                    x3vec(ii3)=idxCur; y3vec(ii3)=idxPrevX; dhdxVec(ii3)=-1/(sX(f3x)); ii3=ii3+1;
                case 'TE'
                    x3vec(ii3)=idxCur; y3vec(ii3)=idxPrevX; dhdxVec(ii3)=-1/(sX(f3x)); ii3=ii3+1;
                end
            end

            % low-y:
            if yy==1
                switch BC(2)
                case -2
                    x4vec(ii4)=idxCur; y4vec(ii4)=nPtsX*(nPtsY-1)+xx; dhdyVec(ii4)=-exp(-1i*ky*nPtsY*dy); ii4=ii4+1;
                case -1
                    x4vec(ii4)=idxCur; y4vec(ii4)=nPtsX*(nPtsY-1)+xx; dhdyVec(ii4)=-1; ii4=ii4+1;
                case 0
                    % do nothing. = 0.
                otherwise
                    error('Unrecognized y-low boundary condition.');
                end
            else
                switch polz
                case 'TM'
                    x4vec(ii4)=idxCur; y4vec(ii4)=idxPrevY; dhdyVec(ii4)=-1/(sY(f3y)); ii4=ii4+1;
                case 'TE'
                    x4vec(ii4)=idxCur; y4vec(ii4)=idxPrevY; dhdyVec(ii4)=-1/(sY(f2y)); ii4=ii4+1;
                end
            end

            % high-x:
            if xx==nPtsX
                switch BC(1)
                case -2
                    x1vec(ii1)=idxCur; y1vec(ii1)=nPtsX*(yy-1)+1; dedxVec(ii1)=+exp(+1i*kx*nPtsX*dx); ii1=ii1+1;
                case -1
                    x1vec(ii1)=idxCur; y1vec(ii1)=nPtsX*(yy-1)+1; dedxVec(ii1)=+1; ii1=ii1+1;
                case 0
                    % do nothing. = 0.
                otherwise
                    error('Unrecognized x-high boundary condition.');
                end
            else
                switch polz
                case 'TM'
                    x1vec(ii1)=idxCur; y1vec(ii1)=idxNextX; dedxVec(ii1)=+1/(sX(f2x)); ii1=ii1+1;
                case 'TE'
                    x1vec(ii1)=idxCur; y1vec(ii1)=idxNextX; dedxVec(ii1)=+1/(sX(f1x)); ii1=ii1+1;
                end
            end

            % high-y:
            if yy==nPtsY
                switch BC(2)
                case -2
                    x2vec(ii2)=idxCur; y2vec(ii2)=xx; dedyVec(ii2)=+exp(+1i*ky*nPtsY*dy); ii2=ii2+1;
                case -1
                    x2vec(ii2)=idxCur; y2vec(ii2)=xx; dedyVec(ii2)=+1; ii2=ii2+1;
                case 0
                    % do nothing. = 0.
                otherwise
                    error('Unrecognized x-high boundary condition.');
                end
            else
                switch polz
                case 'TM'
                    x2vec(ii2)=idxCur; y2vec(ii2)=idxNextY; dedyVec(ii2)=+1/(sY(f1y)); ii2=ii2+1;
                case 'TE'
                    x2vec(ii2)=idxCur; y2vec(ii2)=idxNextY; dedyVec(ii2)=+1/(sY(f1y)); ii2=ii2+1;
                end
            end
        end
    end

    dedxVec(x1vec==0) = [];
    y1vec(x1vec==0) = [];
    x1vec(x1vec==0) = [];

    dedyVec(x2vec==0) = [];
    y2vec(x2vec==0) = [];
    x2vec(x2vec==0) = [];

    dhdxVec(x3vec==0) = [];
    y3vec(x3vec==0) = [];
    x3vec(x3vec==0) = [];

    dhdyVec(x4vec==0) = [];
    y4vec(x4vec==0) = [];
    x4vec(x4vec==0) = [];

    dEdx = sparse(x1vec,y1vec,dedxVec/dx,nVals,nVals);
    dEdy = sparse(x2vec,y2vec,dedyVec/dy,nVals,nVals);
    dHdx = sparse(x3vec,y3vec,dhdxVec/dx,nVals,nVals);
    dHdy = sparse(x4vec,y4vec,dhdyVec/dy,nVals,nVals);

    if CHECK_ABS_REGION
        fieldToCheck = dHdx;
        fieldToCheck = reshape(diag(fieldToCheck),nPtsX,nPtsY);
        figure(21); imagesc(real(fieldToCheck.')); % can plot here.
        axis equal;
        figure(22); imagesc((imag(fieldToCheck.') < 0));
        axis equal;
        figure(23); imagesc(abs(fieldToCheck.')); % can plot here.
        axis equal;

        fieldToCheck = dEdy;
        fieldToCheck = reshape(diag(fieldToCheck),nPtsX,nPtsY);
        figure(24); imagesc(real(fieldToCheck.')); % can plot here.
        axis equal;
        figure(25); imagesc((imag(fieldToCheck.') > 0));
        axis equal;
        assert(false);
    end

    z0 = sparse(nVals,nVals);

    switch polz
    case 'TM'
        line1 = [z0, z0, +dEdy];
        line2 = [z0, z0, -dEdx];
        line3 = [-dHdy, +dHdx, z0];
        WW = [line1; line2; line3];

    case 'TE'
        line1 = [z0, -dEdy, +dEdx];
        line2 = [+dHdy, z0, z0];
        line3 = [-dHdx, z0, z0];
        WW = [line1; line2; line3];

    otherwise
        error('Unrecognized Polarization Type');
    end

    %figure(1);
    %imagesc(real(WW));
    %figure(2);
    %imagesc(imag(WW));

    if sigma0 == 0
        assert(ishermitian(WW));
    else
        %assert(issymmetric(WW));
    end

    %% generate coords:
    xvec = zeros(nVals,1);
    yvec = zeros(nVals,1);
    f1xvec = zeros(nVals,1);
    f2xvec = zeros(nVals,1);
    f3xvec = zeros(nVals,1);
    f1yvec = zeros(nVals,1);
    f2yvec = zeros(nVals,1);
    f3yvec = zeros(nVals,1);
    ii=1;

    for yy=1:nPtsY
        for xx=1:nPtsX
            idxCur = nPtsX*(yy-1) + xx;

            xvec(ii)=idxCur; yvec(ii)=idxCur;
            switch polz
            case 'TM' % f1 = Hx -- f2 = Hy -- f3 = Ez
                f1xvec(ii) = (xx-1)*dx; f1yvec(ii) = (yy-1)*dy + dy/2;
                f2xvec(ii) = (xx-1)*dx + dx/2; f2yvec(ii) = (yy-1)*dy;
                f3xvec(ii) = (xx-1)*dx; f3yvec(ii) = (yy-1)*dy;

            case 'TE' % f1 = Hz -- f2 = Ex -- f3 = Ey
                f1xvec(ii) = (xx-1)*dx + dx/2; f1yvec(ii) = (yy-1)*dy + dy/2;
                f2xvec(ii) = (xx-1)*dx + dx/2; f2yvec(ii) = (yy-1)*dy;
                f3xvec(ii) = (xx-1)*dx; f3yvec(ii) = (yy-1)*dy + dy/2;

            end
            ii = ii+1;
        end
    end

    F1X = sparse(xvec,yvec,f1xvec,nVals,nVals);
    F1Y = sparse(xvec,yvec,f1yvec,nVals,nVals);
    F2X = sparse(xvec,yvec,f2xvec,nVals,nVals);
    F2Y = sparse(xvec,yvec,f2yvec,nVals,nVals);
    F3X = sparse(xvec,yvec,f3xvec,nVals,nVals);
    F3Y = sparse(xvec,yvec,f3yvec,nVals,nVals);

    XX = [F1X, z0, z0; z0, F2X, z0; z0, z0, F3X];
    YY = [F1Y, z0, z0; z0, F2Y, z0; z0, z0, F3Y];

    switch polz
    case 'TM'
        entry_offset = 2*nPtsX*nPtsY;
    case 'TE'
        entry_offset = 0;
    end

    if (mod(nPtsX,2)==0)
        idxMidX1 = round(nPtsX/2) + entry_offset;
        idxMidX2 = round(nPtsX/2)+1 + entry_offset;
        xMidPoint = (XX(idxMidX1,idxMidX1) + XX(idxMidX2,idxMidX2))/2;
    else
        idxMidX = ceil(nPtsX/2) + entry_offset;
        xMidPoint = XX(idxMidX,idxMidX);
    end
    if (mod(nPtsY,2)==0)
        idxMidY1= (round(nPtsY/2)-1)*nPtsX+1 + entry_offset;
        idxMidY2= (round(nPtsY/2))*nPtsX+1 + entry_offset;
        yMidPoint = (YY(idxMidY1,idxMidY1) + YY(idxMidY2,idxMidY2))/2;
    else
        idxMidY = nPtsX*(ceil(nPtsY/2)-1)+1 + entry_offset;
        yMidPoint = YY(idxMidY,idxMidY);
    end
    if strcmp(polz,'TM')
        xMidPoint = xMidPoint + dx/2;
        yMidPoint = yMidPoint + dy/2;
    end
    XX = XX - xMidPoint*speye(size(XX));
    YY = YY - yMidPoint*speye(size(YY));

    %figure(3);
    %imagesc(YY);
    %assert(false);

end

function [] = haldane_2008_wAir_localizer(nuIn,omegaIn,noise_strength)

    %% visualization flags for debugging:
    CHECK_EPS_PLOTS = 0;
    CHECK_ABS_REGION = 0;
    STOP_AFTER_NOISE_NORM = 0;

    %% INPUTS:
    ax = sqrt(3);
    ay = 1;
    omegaLook = omegaIn*(2*pi/ay); %omegaLook = 0.37*(2*pi/ay);

    kappa = 0.04 * (2*pi/(ay^2));

    nUCxAbsorb = 4;
    nUCyAbsorb = 4;
    absorbMax = 4;

    numPtsPerHalfWavelength = 30;
    nUCx = 12;
    nUCy = 22;

    phcBoundary = 8;

    nkeep = 8;

    % phc parameters:
    rodRadius = 0.37; %(3.9mm / 40mm)

    epsAir = 1;
    epsRod = 14;
    epsAirMO = nuIn;
    rodAbsorb = 0;

    muAir = 1;
    muRod = 1;
    %muRodMO = nuIn;

    epsBackground = 1;

    rng_seed = 1;

    polz = 'TE';

    %% Run it:
    GEN_DATA = 1;
    SAVE_BOOL = 0;

    if GEN_DATA
        %parpool(3);

        BC = [0, 0];
        NxAbsorb = round(nUCxAbsorb*ay*numPtsPerHalfWavelength);
        NyAbsorb = round(nUCyAbsorb*ay*numPtsPerHalfWavelength);

        [~,~, epsMuSysTM_InvSqrt, epsMuSysTE_InvSqrt, nPtsX, nPtsY, dx, dy] = eps_mu_gen(ax,ay,nUCx,nUCy,numPtsPerHalfWavelength,rodRadius,epsAir,epsRod,epsAirMO,muAir,muRod,epsBackground,0,phcBoundary,rodAbsorb,noise_strength,rng_seed,CHECK_EPS_PLOTS);
        [WW,XX,YY] = gen_2d_yee_derivs_coords(nPtsX,nPtsY,dx,dy,0,0,NxAbsorb,NyAbsorb,absorbMax,omegaLook,BC,polz,CHECK_ABS_REGION);
        [~,mirrorY] = gen_mirrors(nPtsX,nPtsY,BC,polz);
        switch polz
        case 'TM'
            [WW,epsMuSys_InvSqrt,XX,YY] = trim_to_PEC(WW,epsMuSysTM_InvSqrt,XX,YY,BC,polz);
            [XX,YY,epsMuSys_InvSqrt] = check_and_sym(XX,YY,epsMuSys_InvSqrt,WW,mirrorY,nPtsX,nPtsY,BC,polz);
        case 'TE'
            [WW,epsMuSys_InvSqrt,XX,YY] = trim_to_PEC(WW,epsMuSysTE_InvSqrt,XX,YY,BC,polz);
            [XX,YY,epsMuSys_InvSqrt] = check_and_sym(XX,YY,epsMuSys_InvSqrt,WW,mirrorY,nPtsX,nPtsY,BC,polz);
        end

        II = speye(size(WW));
        HH = epsMuSys_InvSqrt*WW*epsMuSys_InvSqrt;

        % check strength of the disorder:
        if abs(noise_strength) > 0

            [~,~, epsMuSysTM_InvSqrt, epsMuSysTE_InvSqrt, nPtsX, nPtsY, ~, ~] = eps_mu_gen(ax,ay,nUCx,nUCy,numPtsPerHalfWavelength,rodRadius,epsAir,epsRod,epsAirMO,muAir,muRod,epsBackground,0,phcBoundary,rodAbsorb,0,rng_seed,CHECK_EPS_PLOTS);
            [WW,XX,YY] = gen_2d_yee_derivs_coords(nPtsX,nPtsY,dx,dy,0,0,NxAbsorb,NyAbsorb,absorbMax,omegaLook,BC,polz,CHECK_ABS_REGION);
            switch polz
            case 'TM'
                [WW,epsMuSys_InvSqrt,XX,YY] = trim_to_PEC(WW,epsMuSysTM_InvSqrt,XX,YY,BC,polz);
                [XX,YY,epsMuSys_InvSqrt] = check_and_sym(XX,YY,epsMuSys_InvSqrt,WW,mirrorY,nPtsX,nPtsY,BC,polz);
            case 'TE'
                [WW,epsMuSys_InvSqrt,XX,YY] = trim_to_PEC(WW,epsMuSysTE_InvSqrt,XX,YY,BC,polz);
                [XX,YY,epsMuSys_InvSqrt] = check_and_sym(XX,YY,epsMuSys_InvSqrt,WW,mirrorY,nPtsX,nPtsY,BC,polz);
            end
            HHclean = epsMuSys_InvSqrt*WW*epsMuSys_InvSqrt;

            normVal = normest(HH - HHclean,1e-3);

            disp('norm of perturbation on full Hilbert space in units of (2 pi c/a):')
            disp(normVal*(ay/(2*pi)));

            nKeepCheck = 50;
            opts.tol = 1e-3;
            opts.maxit = 100;
            [eVecs,~] = eigs(HHclean,nKeepCheck,omegaLook,opts);
            reducedPertHH = eVecs.' * (HH - HHclean) * eVecs;
            normValR = normest(reducedPertHH,1e-3);
            
            disp('norm of perturbation on reduced Hilbert space in units of (2 pi c/a):')
            disp(normValR*(ay/(2*pi)));

            % projected localizer gap:
            rHH = eVecs.' * HH * eVecs;
            rXX = eVecs.' * XX * eVecs;
            rYY = eVecs.' * YY * eVecs;

            rII = speye(size(rHH));

            rLL = [(rHH - omegaLook*rII), kappa*(rXX - 0*rII) - 1i*kappa*(rYY - 0*rII);...
                    kappa*(rXX - 0*rII) + 1i*kappa*(rYY - 0*rII), -(rHH - omegaLook*rII)'];

            rGap = eigs(rLL,1,'sm',opts);

            disp('Gap of the reduced localizer in units of (2 pi c/a):')
            disp(abs(rGap)*(ay/(2*pi)));
            
            disp('ratio of disorder to gap:');
            disp(normValR/abs(rGap));

            if STOP_AFTER_NOISE_NORM
                assert(false);
            end
        end

        % set up where to calculate the localizer gap
        yFixed = 0;

        xVec = (-(nUCx/2+1/2)*ax):(3*dx):((nUCx/2+1/2)*ax);
        %xVec = -(nUCx/2)*ax*(2/3);
        %xVec = [0, 1, 2, 3, 4, 5, 6, 7]*ax;

        gapVec = zeros(length(xVec),1);
        flowMat = zeros(nkeep,length(xVec));

        % run it:
        %opts.disp = 0;
        opts.tol = 1e-3;
        opts.maxit = 100;
        for xx=1:length(xVec)
            disp(xx/length(xVec));
            
            LL = [(HH - omegaLook*II), kappa*(XX - xVec(xx)*II) - 1i*kappa*(YY - yFixed*II);...
                kappa*(XX - xVec(xx)*II) + 1i*kappa*(YY - yFixed*II), -(HH - omegaLook*II)'];

            flowMat(:,xx) = eigs(LL,nkeep,'sm',opts);
            gapVec(xx) = min(abs(real(flowMat(:,xx))));
        end

        %% save it:
        if SAVE_BOOL
            save(['data_haldane/haldane_2008_twoPh_local_kappa_',sprintf('%02ip%02i',floor(kappa),round(100*(kappa-floor(kappa)))),'_nu_',sprintf('%02ip%02i',floor(abs(nuIn)),round(100*(abs(nuIn)-floor(abs(nuIn))))),'.mat'],'xVec','gapVec','flowMat','absorbMax','yFixed','numPtsPerHalfWavelength','nUCx','nUCy','ax','ay','rodRadius','epsAir','epsRod','epsAirMO','muAir','muRod','kappa');
        end

        if length(xVec)==1
            disp(gapVec*(ay/(2*pi)));
            disp(flowMat*(ay/(2*pi)));
            assert(false);
        end
    else
        load(['data_haldane/haldane_2008_twoPh_local_kappa_',sprintf('%02ip%02i',floor(kappa),round(100*(kappa-floor(kappa)))),'_nu_',sprintf('%02ip%02i',floor(abs(nuIn)),round(100*(abs(nuIn)-floor(abs(nuIn))))),'.mat'],'xVec','gapVec','flowMat','absorbMax','yFixed','numPtsPerHalfWavelength','nUCx','nUCy','ax','ay','rodRadius','epsAir','epsRod','epsAirMO','muAir','muRod');
    end

    %% Plot it:
    figure(1);
    plot(xVec/ax,gapVec*(ay/(2*pi)),'Marker','x','Color',[0.5, 0.5, 0.5],'LineWidth',2);
    AX = gca; % get current axes for figure.
    AX.FontSize = 16;
    AX.TickLabelInterpreter = 'latex';
    AX.LineWidth = 1.5;
    xlabel('Position, $x$ $(a_x)$','Interpreter','latex','FontSize',18);
    ylabel('Localizer gap $(2\pi c/a_y)$','Interpreter','latex','FontSize',18);
    %axis([min(xVec)/ax, max(xVec)/ax, 0, 0.04]);
    grid on;

    figure(2);
    for nn=1:length(flowMat(:,1))
        plot(xVec/ax,real(flowMat(nn,:))*(ay/(2*pi)),'Marker','x','Color',[0.5, 0.5, 0.5],'LineStyle','none','LineWidth',2);
        hold on;
    end
    AX = gca; % get current axes for figure.
    AX.FontSize = 16;
    AX.TickLabelInterpreter = 'latex';
    AX.LineWidth = 1.5;
    xlabel('Position, $x$ $(a_x)$','Interpreter','latex','FontSize',18);
    ylabel('Localizer spectral flow','Interpreter','latex','FontSize',18);
    %axis([min(xVec)/ax, max(xVec)/ax, -1.1, 1.1]);
    grid on;
    hold off;

end

function [symX,symY,symEpsMu] = check_and_sym(XX,YY,epsMu,WW,mirrorY,nPtsX,nPtsY,BC,polz)

    errorTol = 1e-8;
    checkXX = mirrorY*XX-XX*mirrorY;
    %checkXX(abs(checkXX) < errorTol) = 0;
    checkYY = mirrorY*YY+YY*mirrorY;
    %checkYY(abs(checkYY) < errorTol) = 0;
    assert(sum(sum(abs(checkXX)))<errorTol);
    assert(sum(sum(abs(checkYY)))<errorTol);

    symX = (XX + mirrorY*XX*mirrorY)/2;
    symY = (YY - mirrorY*YY*mirrorY)/2;
    assert(issparse(symX));
    assert(issparse(symY));

    assert(isequal(mirrorY*symX,symX*mirrorY));
    assert(isequal(mirrorY*symY,-symY*mirrorY));

    symEpsMu = real(epsMu + mirrorY*epsMu*mirrorY)/2 + 1i*imag(epsMu);
    %assert(isequal(mirrorY*symEpsMu,symEpsMu*mirrorY));

    %checkEE = mirrorX*mirrorY*symEpsMu-symEpsMu*mirrorY*mirrorX;
    %checkEE(abs(checkEE) < errorTol) = 0;
    %assert(sum(sum(abs(checkEE)))==0);

    %figure(1);
    %imagesc(WW);
    %figure(2);
    %imagesc(mirrorX*mirrorY);
    %figure(3);
    %imagesc(WW - mirrorX*mirrorY*WW*mirrorY*mirrorX);

    %assert(isequal(mirrorX*mirrorY*WW,WW*mirrorY*mirrorX));

    CHECK_INV_PLOTS=0;
    if CHECK_INV_PLOTS
        switch polz
        case 'TM'
            f1Len = (nPtsX-1)*nPtsY;
            f2Len = nPtsX*(nPtsY-1);
            f3Len = (nPtsX-1)*(nPtsY-1);
        case 'TE'
            f1Len = nPtsX*nPtsY;
            f2Len = nPtsX*(nPtsY-1);
            f3Len = (nPtsX-1)*nPtsY;
        end
        XXm = mirrorY*XX*mirrorY;
        YYm = mirrorY*YY*mirrorY;
        figure(10);
        for ii=1:f1Len
            plot(XX(ii,ii),YY(ii,ii),'Marker','x','Color','r');
            hold on;
            plot(XXm(ii,ii),YYm(ii,ii),'Marker','o','Color','g');
        end
        for ii=(f1Len+1):(f1Len+f2Len)
            plot(XX(ii,ii),YY(ii,ii),'Marker','>','Color','r');
            plot(XXm(ii,ii),YYm(ii,ii),'Marker','<','Color','g');
        end
        for ii=(f1Len+f2Len+1):(f1Len+f2Len+f3Len)
        plot(XX(ii,ii),YY(ii,ii),'Marker','^','Color','r');
        plot(XXm(ii,ii),YYm(ii,ii),'Marker','+','Color','g');
        end
        figure(10); hold off;
        assert(false);
    end

    CHECK_EPS_AVE_PLOTS=0;
    if CHECK_EPS_AVE_PLOTS
        switch polz
        case 'TM'
            epsMuSysTM = symEpsMu;
            epsZZ = diag(epsMuSysTM((nPtsX-1)*nPtsY+nPtsX*(nPtsY-1)+1:(nPtsX-1)*nPtsY+nPtsX*(nPtsY-1)+(nPtsX-1)*(nPtsY-1),(nPtsX-1)*nPtsY+nPtsX*(nPtsY-1)+1:(nPtsX-1)*nPtsY+nPtsX*(nPtsY-1)+(nPtsX-1)*(nPtsY-1)));
            epsZZ = reshape(epsZZ,nPtsX-1,nPtsY-1);
            figure(1);
            imagesc(epsZZ);
            axis equal;
        case 'TE'
            epsMuSysTE = symEpsMu;
            epsXX = diag(epsMuSysTE((nPtsX*nPtsY+1):(nPtsX*nPtsY+nPtsX*(nPtsY-1)),(nPtsX*nPtsY+1):(nPtsX*nPtsY+nPtsX*(nPtsY-1))));
            epsYY = diag(epsMuSysTE((nPtsX*nPtsY+nPtsX*(nPtsY-1)+1):end,(nPtsX*nPtsY+nPtsX*(nPtsY-1)+1):end));
            epsXX = reshape(epsXX,nPtsX,nPtsY-1);
            epsYY = reshape(epsYY,nPtsX-1,nPtsY);
            figure(1);
            imagesc(epsXX);
            axis equal;
            figure(2);
            imagesc(epsYY);
            axis equal;

            epsXY = epsMuSysTE((nPtsX*nPtsY+1):(nPtsX*nPtsY+nPtsX*(nPtsY-1)),(nPtsX*nPtsY+nPtsX*(nPtsY-1)+1):end);
            figure(3);
            imagesc(imag(epsXY));
        end

        assert(false);
    end
end

function [mirrorX,mirrorY] = gen_mirrors(nPtsX,nPtsY,BC,polz)

    assert(BC(1)==0);
    assert(BC(2)==0);

    switch polz
    case 'TM'
        nVals = (nPtsX-1)*nPtsY + nPtsX*(nPtsY-1) + (nPtsX-1)*(nPtsY-1);
        xXvec = zeros(nVals,1);
        yXvec = zeros(nVals,1);
        mXvec = zeros(nVals,1);
        xYvec = zeros(nVals,1);
        yYvec = zeros(nVals,1);
        mYvec = zeros(nVals,1);

        iiX = 1;
        iiY = 1;

        for yy=1:nPtsY
            for xx=1:(nPtsX-1)
                idxCur = (nPtsX-1)*(yy-1) + xx;
                idxMirrorX = (nPtsX-1)*(yy-1) + ((nPtsX-1) - (xx-1));
                idxMirrorY = (nPtsX-1)*(nPtsY-1 - (yy-1)) + xx;

                xXvec(iiX) = idxMirrorX; yXvec(iiX) = idxCur; mXvec(iiX) = 1; iiX=iiX+1;
                xYvec(iiY) = idxMirrorY; yYvec(iiY) = idxCur; mYvec(iiY) = 1; iiY=iiY+1;
            end
        end
        for yy=1:(nPtsY-1)
            for xx=1:nPtsX
                idxCur = (nPtsX-1)*nPtsY + nPtsX*(yy-1) + xx;
                idxMirrorX = (nPtsX-1)*nPtsY + nPtsX*(yy-1) + (nPtsX - (xx-1));
                idxMirrorY = (nPtsX-1)*nPtsY + nPtsX*(nPtsY-2 - (yy-1)) + xx;

                xXvec(iiX) = idxMirrorX; yXvec(iiX) = idxCur; mXvec(iiX) = 1; iiX=iiX+1;
                xYvec(iiY) = idxMirrorY; yYvec(iiY) = idxCur; mYvec(iiY) = 1; iiY=iiY+1;
            end
        end
        for yy=1:(nPtsY-1)
            for xx=1:(nPtsX-1)
                idxCur = (nPtsX-1)*nPtsY + nPtsX*(nPtsY-1) + (nPtsX-1)*(yy-1) + xx;
                idxMirrorX = (nPtsX-1)*nPtsY + nPtsX*(nPtsY-1) + (nPtsX-1)*(yy-1) + (nPtsX-1 - (xx-1));
                idxMirrorY = (nPtsX-1)*nPtsY + nPtsX*(nPtsY-1) + (nPtsX-1)*(nPtsY-2 - (yy-1)) + xx;

                xXvec(iiX) = idxMirrorX; yXvec(iiX) = idxCur; mXvec(iiX) = 1; iiX=iiX+1;
                xYvec(iiY) = idxMirrorY; yYvec(iiY) = idxCur; mYvec(iiY) = 1; iiY=iiY+1;
            end
        end

    case 'TE'
        nVals = nPtsX*nPtsY + (nPtsX-1)*nPtsY + nPtsX*(nPtsY-1);
        xXvec = zeros(nVals,1);
        yXvec = zeros(nVals,1);
        mXvec = zeros(nVals,1);
        xYvec = zeros(nVals,1);
        yYvec = zeros(nVals,1);
        mYvec = zeros(nVals,1);

        iiX = 1;
        iiY = 1;

        for yy=1:nPtsY
            for xx=1:nPtsX
                idxCur = nPtsX*(yy-1) + xx;
                idxMirrorX = nPtsX*(yy-1) + (nPtsX - (xx-1));
                idxMirrorY = nPtsX*(nPtsY-1 - (yy-1)) + xx;

                xXvec(iiX) = idxMirrorX; yXvec(iiX) = idxCur; mXvec(iiX) = 1; iiX=iiX+1;
                xYvec(iiY) = idxMirrorY; yYvec(iiY) = idxCur; mYvec(iiY) = 1; iiY=iiY+1;
            end
        end
        for yy=1:(nPtsY-1)
            for xx=1:nPtsX
                idxCur = nPtsX*nPtsY + nPtsX*(yy-1) + xx;
                idxMirrorX = nPtsX*nPtsY + nPtsX*(yy-1) + (nPtsX - (xx-1));
                idxMirrorY = nPtsX*nPtsY + nPtsX*(nPtsY-2 - (yy-1)) + xx;

                xXvec(iiX) = idxMirrorX; yXvec(iiX) = idxCur; mXvec(iiX) = 1; iiX=iiX+1;
                xYvec(iiY) = idxMirrorY; yYvec(iiY) = idxCur; mYvec(iiY) = 1; iiY=iiY+1;
            end
        end
        for yy=1:nPtsY
            for xx=1:(nPtsX-1)
                idxCur = nPtsX*nPtsY + nPtsX*(nPtsY-1) + (nPtsX-1)*(yy-1) + xx;
                idxMirrorX = nPtsX*nPtsY + nPtsX*(nPtsY-1) + (nPtsX-1)*(yy-1) + (nPtsX-1 - (xx-1));
                idxMirrorY = nPtsX*nPtsY + nPtsX*(nPtsY-1) + (nPtsX-1)*(nPtsY-1 - (yy-1)) + xx;

                xXvec(iiX) = idxMirrorX; yXvec(iiX) = idxCur; mXvec(iiX) = 1; iiX=iiX+1;
                xYvec(iiY) = idxMirrorY; yYvec(iiY) = idxCur; mYvec(iiY) = 1; iiY=iiY+1;
            end
        end
    end

    mirrorX = sparse(xXvec,yXvec,mXvec,nVals,nVals);
    mirrorY = sparse(xYvec,yYvec,mYvec,nVals,nVals);

    assert(isequal(mirrorX*mirrorX,speye(nVals,nVals)));
    assert(isequal(mirrorY*mirrorY,speye(nVals,nVals)));
    assert(isequal(mirrorX*mirrorY,mirrorY*mirrorX));

end
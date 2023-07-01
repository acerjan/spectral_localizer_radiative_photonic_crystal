function [] = haldane_2008_wAir_LDOS(nuIn,omegaIn,noise_strength)

    %% visualization flags for debugging:
    CHECK_EPS_PLOTS = 0;
    CHECK_ABS_REGION = 0;

    %% INPUTS:
    ax = sqrt(3);
    ay = 1;
    omegaLook = omegaIn*(2*pi/ay);

    nUCxAbsorb = 4;
    nUCyAbsorb = 4;
    absorbMax = 4;

    nKeep = 50;

    numPtsPerHalfWavelength = 30;
    nUCx = 12;
    nUCy = 22;

    phcBoundary = 8;

    % use source?
    srcLocX = 1.5;
    srcLocY = 1;
    sourceStr = 0;

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
        BC = [0, 0];
        NxAbsorb = round(nUCxAbsorb*ay*numPtsPerHalfWavelength);
        NyAbsorb = round(nUCyAbsorb*ay*numPtsPerHalfWavelength);

        [epsMuSysTM, epsMuSysTE, ~,~, nPtsX, nPtsY, dx, dy] = eps_mu_gen(ax,ay,nUCx,nUCy,numPtsPerHalfWavelength,rodRadius,epsAir,epsRod,epsAirMO,muAir,muRod,epsBackground,0,phcBoundary,rodAbsorb,noise_strength,rng_seed,CHECK_EPS_PLOTS);
        [WW,XX,YY] = gen_2d_yee_derivs_coords(nPtsX,nPtsY,dx,dy,0,0,NxAbsorb,NyAbsorb,absorbMax,omegaLook,BC,polz,CHECK_ABS_REGION);

        if sourceStr == 0 % calculate LDOS
            [~,mirrorY] = gen_mirrors(nPtsX,nPtsY,BC,polz);
            switch polz
            case 'TM'
                [WW,epsMuSysTM,XX,YY] = trim_to_PEC(WW,epsMuSysTM,XX,YY,BC,polz);
                [XX,YY,epsMuSys] = check_and_sym(XX,YY,epsMuSysTM,WW,mirrorY,nPtsX,nPtsY,BC,polz);
                assert(false);
            case 'TE'
                [WW,epsMuSysTE,XX,YY] = trim_to_PEC(WW,epsMuSysTE,XX,YY,BC,polz);
                [XX,YY,epsMuSys] = check_and_sym(XX,YY,epsMuSysTE,WW,mirrorY,nPtsX,nPtsY,BC,polz);
            end

            opts.tol = 1e-3;
            opts.maxit = 100;
            [eVecs,omegas] = eigs(WW,epsMuSys,nKeep,omegaLook,opts);
            FzLDOS = zeros(nPtsX,nPtsY);
            omegas = diag(omegas);

            for ii=1:nKeep
                omegaDist = (1/(abs(imag(omegas(ii)))*sqrt(2*pi)))*exp(-(real(omegas(ii))-omegaLook)^2/(2*(imag(omegas(ii)))^2));
                FzCur = eVecs(1:(nPtsX*nPtsY),ii); %eVecs(2*(nPtsX*nPtsY)+1:3*(nPtsX*nPtsY),ii);%eVecs(1:(nPtsX*nPtsY),ii);
                FzCur = omegaDist * FzCur/(sqrt(FzCur' * FzCur));
                FzCur = reshape(FzCur,nPtsX,nPtsY);
                FzLDOS = FzLDOS + abs(FzCur).^2;
            end

            Fz = reshape(eVecs(1:(nPtsX*nPtsY),1),nPtsX,nPtsY);

            %HzLDOS(:,(nPtsY)/2+2:end) = fliplr(HzLDOS(:,1:(nPtsY)/2));

            %% save it:
            if SAVE_BOOL
                save(['data_haldane/haldane_2008_wAir_LDOS_nu_',sprintf('%02ip%02i',floor(abs(nuScale)),round(100*(abs(nuScale)-floor(abs(nuScale))))),'_wIn_',sprintf('%02ip%04i',floor(abs(omegaIn)),round(10000*(abs(omegaIn)-floor(abs(omegaIn))))),'.mat'],'numPtsPerHalfWavelength','nUCx','nUCy','ax','ay','rodR','epsAir','epsMat','epsAirMO','epsMatMO','muAir','muRod');
            end

            omegas*(ay/(2*pi))

        else % or, calculate the response to the source:

            SidxX = round((nUCx/2+srcLocX)*ax*numPtsPerHalfWavelength);
            SidxY = round((nUCy/2+srcLocY)*ay*numPtsPerHalfWavelength);

            switch polz
            case 'TM'
                epsMuSys = epsMuSysTM;
                sourceJ = sparse(2*nPtsX*nPtsY+ nPtsX*(SidxY-1)+SidxX,1,1i*sourceStr,3*nPtsX*nPtsY,1);
            case 'TE'
                epsMuSys = epsMuSysTE;
                sourceJ = sparse(0*nPtsX*nPtsY+ nPtsX*(SidxY-1)+SidxX,1,1i*sourceStr,3*nPtsX*nPtsY,1);
            end

            MEmat = (WW - omegaLook*epsMuSys);
            outVec = MEmat\sourceJ;

            switch polz
            case 'TM'
                Fz = reshape(outVec((2*nPtsX*nPtsY+1):(3*nPtsX*nPtsY),1),nPtsX,nPtsY);
            case 'TE'
                Fz = reshape(outVec(1:(nPtsX*nPtsY),1),nPtsX,nPtsY);
            end

            FzLDOS = abs(Fz).^2;

        end

    else
        load(['data_haldane/haldane_2008_wAir_LDOS_nu_',sprintf('%02ip%02i',floor(abs(nuScale)),round(100*(abs(nuScale)-floor(abs(nuScale))))),'_wIn_',sprintf('%02ip%04i',floor(abs(omegaIn)),round(10000*(abs(omegaIn)-floor(abs(omegaIn))))),'.mat'],'numPtsPerHalfWavelength','nUCx','nUCy','ax','ay','rodR','epsAir','epsMat','epsAirMO','epsMatMO','muAir','muRod');
    end

    %% Plot it:
    xMin = -(ax*nUCx/2 + ax/2);
    xMax = (ax*nUCx/2 + ax/2);
    yMin = -(ay*nUCy/2 + ay/2);
    yMax = (ay*nUCy/2 + ay/2);
    buf = 0;

    % plot the structure:
    %[epsMuSysTM, epsMuSysTE, ~, ~, nPtsX, nPtsY, dx, dy] = eps_mu_gen_HR_wAir(ax,ay,nUCx,nUCy,numPtsPerHalfWavelength,rodRadius,epsAir,epsRod,epsAirMO,muAir,muRod,epsBackground,0,phcBoundary);
    epsDiag = reshape(diag(epsMuSys(nPtsX*nPtsY+1:2*nPtsX*nPtsY,nPtsX*nPtsY+1:2*nPtsX*nPtsY)),nPtsX,nPtsY);

    XX = full(diag(XX));
    YY = full(diag(YY));
    XX = reshape(XX(1:nPtsX*nPtsY),nPtsX,nPtsY);
    YY = reshape(YY(1:nPtsX*nPtsY),nPtsX,nPtsY);

    figure(1);
    pp = pcolor(XX,YY,real(epsDiag));
    pp.EdgeColor = 'none';
    AX = gca; % get current axes for figure.
    AX.FontSize = 16;
    AX.TickLabelInterpreter = 'latex';
    AX.LineWidth = 1.5;
    xlabel('$x_0$','Interpreter','latex','FontSize',18);
    ylabel('$y_0$','Interpreter','latex','FontSize',18);
    colorbar('LineWidth',1.5,'TickLabelInterpreter','latex','FontSize',16);
    %caxis([0, 0.02]);
    %colormap(bone);
    hold on;
    ms = 1;
    cc = [0,0,0];
    plot([xMin-buf],[yMin-buf],'Marker','.','Color',cc,'MarkerSize',ms);
    plot([xMax+buf],[yMin-buf],'Marker','.','Color',cc,'MarkerSize',ms);
    plot([xMin-buf],[yMax+buf],'Marker','.','Color',cc,'MarkerSize',ms);
    plot([xMax+buf],[yMax+buf],'Marker','.','Color',cc,'MarkerSize',ms);
    hold off;
    axis([xMin-buf, xMax+buf, yMin-buf, yMax+buf]);
    axis equal;

    % plot the LDOS:
    figure(2);
    pp = pcolor(XX,YY,FzLDOS);
    pp.EdgeColor = 'none';
    %xticks(-2:2:10);
    %yticks(-2:2:10);
    AX = gca; % get current axes for figure.
    AX.FontSize = 16;
    AX.TickLabelInterpreter = 'latex';
    AX.LineWidth = 1.5;
    xlabel('$x_0$','Interpreter','latex','FontSize',18);
    ylabel('$y_0$','Interpreter','latex','FontSize',18);
    colorbar('LineWidth',1.5,'TickLabelInterpreter','latex','FontSize',16);
    %caxis([0, 0.02]);
    colormap(magma);
    hold on;
    ms = 1;
    cc = [0,0,0];
    plot([xMin-buf],[yMin-buf],'Marker','.','Color',cc,'MarkerSize',ms);
    plot([xMax+buf],[yMin-buf],'Marker','.','Color',cc,'MarkerSize',ms);
    plot([xMin-buf],[yMax+buf],'Marker','.','Color',cc,'MarkerSize',ms);
    plot([xMax+buf],[yMax+buf],'Marker','.','Color',cc,'MarkerSize',ms);
    hold off;
    axis([xMin-buf, xMax+buf, yMin-buf, yMax+buf]);
    axis equal;

    figure(3);
    pp = pcolor(XX,YY,real(Fz));
    pp.EdgeColor = 'none';
    %xticks(-2:2:10);
    %yticks(-2:2:10);
    AX = gca; % get current axes for figure.
    AX.FontSize = 16;
    AX.TickLabelInterpreter = 'latex';
    AX.LineWidth = 1.5;
    xlabel('$x_0$','Interpreter','latex','FontSize',18);
    ylabel('$y_0$','Interpreter','latex','FontSize',18);
    colorbar('LineWidth',1.5,'TickLabelInterpreter','latex','FontSize',16);
    %caxis([0, 0.02]);
    colormap(redblue);
    hold on;
    ms = 1;
    cc = [0,0,0];
    plot([xMin-buf],[yMin-buf],'Marker','.','Color',cc,'MarkerSize',ms);
    plot([xMax+buf],[yMin-buf],'Marker','.','Color',cc,'MarkerSize',ms);
    plot([xMin-buf],[yMax+buf],'Marker','.','Color',cc,'MarkerSize',ms);
    plot([xMax+buf],[yMax+buf],'Marker','.','Color',cc,'MarkerSize',ms);
    hold off;
    axis([xMin-buf, xMax+buf, yMin-buf, yMax+buf]);
    axis equal;

    figure(4);
    pp = pcolor(XX,YY,imag(Fz));
    pp.EdgeColor = 'none';
    %xticks(-2:2:10);
    %yticks(-2:2:10);
    AX = gca; % get current axes for figure.
    AX.FontSize = 16;
    AX.TickLabelInterpreter = 'latex';
    AX.LineWidth = 1.5;
    xlabel('$x_0$','Interpreter','latex','FontSize',18);
    ylabel('$y_0$','Interpreter','latex','FontSize',18);
    colorbar('LineWidth',1.5,'TickLabelInterpreter','latex','FontSize',16);
    %caxis([0, 0.02]);
    colormap(redblue);
    hold on;
    ms = 1;
    cc = [0,0,0];
    plot([xMin-buf],[yMin-buf],'Marker','.','Color',cc,'MarkerSize',ms);
    plot([xMax+buf],[yMin-buf],'Marker','.','Color',cc,'MarkerSize',ms);
    plot([xMin-buf],[yMax+buf],'Marker','.','Color',cc,'MarkerSize',ms);
    plot([xMax+buf],[yMax+buf],'Marker','.','Color',cc,'MarkerSize',ms);
    hold off;
    axis([xMin-buf, xMax+buf, yMin-buf, yMax+buf]);
    axis equal;

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

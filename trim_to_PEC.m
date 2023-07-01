function [WW,epsMu,XX,YY] = trim_to_PEC(WW,epsMu,XX,YY,BC,polz)

    XX = full(diag(XX));
    YY = full(diag(YY));

    minX = min(XX);
    idxMinX = find(XX==minX);

    minY = min(YY);
    idxMinY = find(YY==minY);

    nPtsY = length(idxMinX);
    nPtsX = length(idxMinY);

    if strcmp(polz,'TM')
        nPtsY = nPtsY/2;
        nPtsX = nPtsX/2;
    end
    assert(3*nPtsX*nPtsY == length(XX));

    CHECK_TRIM = 0;
    if CHECK_TRIM
        figure(10);
        for ii=1:nPtsX*nPtsY
            plot(XX(ii),YY(ii),'Marker','x','Color','r');
            hold on;
        end
        for ii=(nPtsX*nPtsY+1):2*nPtsX*nPtsY
            plot(XX(ii),YY(ii),'Marker','>','Color','r');
        end
        for ii=(2*nPtsX*nPtsY+1):3*nPtsX*nPtsY %(nPtsX-1)*(nPtsY)
            plot(XX(ii),YY(ii),'Marker','^','Color','r');
        end
    end

    idxRemove = [];
    if BC(1)==0
        idxRemove = [idxRemove; idxMinX];
    end
    if BC(2)==0
        idxRemove = [idxRemove; idxMinY];
    end

    XX(idxRemove) = [];
    YY(idxRemove) = [];
    WW(idxRemove,:) = [];
    WW(:,idxRemove) = [];
    epsMu(idxRemove,:) = [];
    epsMu(:,idxRemove) = [];

    if (BC(1)==0) && (BC(2)==0)
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
    elseif (BC(1)==0)
        switch polz
        case 'TM'
            f1Len = (nPtsX-1)*nPtsY;
            f2Len = nPtsX*(nPtsY);
            f3Len = (nPtsX-1)*(nPtsY);
        case 'TE'
            f1Len = nPtsX*nPtsY;
            f2Len = (nPtsX)*nPtsY;
            f3Len = (nPtsX-1)*(nPtsY);
        end
    elseif (BC(2)==0)
        switch polz
        case 'TM'
            f1Len = (nPtsX)*nPtsY;
            f2Len = nPtsX*(nPtsY-1);
            f3Len = (nPtsX)*(nPtsY-1);
        case 'TE'
            f1Len = nPtsX*nPtsY;
            f2Len = (nPtsX)*(nPtsY-1);
            f3Len = nPtsX*nPtsY;
        end
    else
        disp('fully periodic. Why trim?');
        assert(false);
    end
    assert(length(XX) == f1Len + f2Len + f3Len);

    XX = diag(sparse(XX));
    YY = diag(sparse(YY));

    %assert(issymmetric(WW));
    %assert(ishermitian(epsMu));

    if CHECK_TRIM
        for ii=1:f1Len
            %plot(XX(ii,ii),YY(ii,ii),'Marker','x','Color','r');
            hold on;
            plot(XX(ii,ii),YY(ii,ii),'Marker','o','Color','g');
        end
        for ii=(f1Len+1):(f1Len+f2Len) %(nPtsX)*(nPtsY-1)
            %plot(XX(nPtsX*nPtsY+ii,nPtsX*nPtsY+ii),YY(nPtsX*nPtsY+ii,nPtsX*nPtsY+ii),'Marker','>','Color','r');
            plot(XX(ii,ii),YY(ii,ii),'Marker','<','Color','g');
        end
        for ii=(f1Len+f2Len+1):(f1Len+f2Len+f3Len) %(nPtsX-1)*(nPtsY)
            %plot(XX(2*nPtsX*nPtsY+ii,2*nPtsX*nPtsY+ii),YY(2*nPtsX*nPtsY+ii,2*nPtsX*nPtsY+ii),'Marker','^','Color','r');
            plot(XX(ii,ii),YY(ii,ii),'Marker','+','Color','g');
        end
        figure(10);
        hold off;
        assert(false);
    end
end
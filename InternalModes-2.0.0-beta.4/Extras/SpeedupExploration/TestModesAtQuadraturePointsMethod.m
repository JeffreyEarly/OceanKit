N0 = 3*2*pi/3600;
L_gm = 1300;
N2 = @(z) N0*N0*exp(2*z/L_gm);
Lz = 4000;
zDomain = [-Lz 0];

[N2,zDomain] = stratificationProfileWithMixedLayerModified();

im = InternalModesWKBSpectral(N2=N2,zIn=zDomain,latitude=31,rho0=1025);


%% We need to establish just how many points we need *for a given k*. So we
% will run a simple convergence test.
cachePath = fullfile(fileparts(mfilename("fullpath")), "k_convergence_ml.mat");
if exist(cachePath,"file")
    load(cachePath);
else
    lambda = [1e5;1e4;1e3;3.1e2;1e2;3.1e1;1e1;3.1]; % ;1e-1 % 10 cm spirals out of control
    nModes = sort([20; 40; 80; 160],'descend');

    k = 2*pi./lambda;
    [K,MODES] = ndgrid(k,nModes);
    minEVP = nan(size(K));

    dEVP = 5;

    for i=1:numel(minEVP)
        isConverged = false;

        nModes = MODES(i);
        im.nEVP = MODES(i) + dEVP;

        [A,B] = im.eigenmatricesForWavenumber(K(i));
        [~,D] = eig( A, B );
        h_lb = sort(real(im.hFromLambda(diag(D))),'descend');

        while (~isConverged)
            im.nEVP = im.nEVP + dEVP;
            [A,B] = im.eigenmatricesForWavenumber(K(i));
            [~,D] = eig( A, B );
            h_ub = sort(real(im.hFromLambda(diag(D))),'descend');

            rel_error = abs((h_ub(1:nModes)-h_lb(1:nModes))./h_ub(1:nModes));
            isConverged = max(rel_error) < 1e-8;
            h_lb = h_ub;
        end
        im.nEVP = im.nEVP - dEVP;
        minEVP(i) = im.nEVP;
    end

    % for iK = 1:size(K,1)
    %     iMode = 1;
    %     i = sub2ind(size(K),iK,iMode);
    %
    %     % first we get the high resolved modes to converge
    %     isConverged = false;
    %
    %     nModes = MODES(i);
    %     im.nEVP = MODES(i) + dEVP;
    %
    %     [A,B] = im.eigenmatricesForWavenumber(K(i));
    %     [~,D] = eig( A, B );
    %     h_lb = sort(real(im.hFromLambda(diag(D))),'descend');
    %
    %     while (~isConverged)
    %         im.nEVP = im.nEVP + dEVP;
    %         [A,B] = im.eigenmatricesForWavenumber(K(i));
    %         [~,D] = eig( A, B );
    %         h_ub = sort(real(im.hFromLambda(diag(D))),'descend');
    %
    %         rel_error = abs((h_ub(1:nModes)-h_lb(1:nModes))./h_ub(1:nModes));
    %         isConverged = max(rel_error) < 1e-8;
    %         h_lb = h_ub;
    %     end
    %     h_high_res = h_lb;
    %     im.nEVP = im.nEVP - dEVP;
    %     minEVP(i) = im.nEVP;
    %
    %
    %     for iMode = 2:size(MODES,2)
    %         i = sub2ind(size(K),iK,iMode);
    %         disp(i)
    %         isConverged = false;
    %
    %         nModes = MODES(i);
    %         im.nEVP = MODES(i) + dEVP;
    %
    %         [A,B] = im.eigenmatricesForWavenumber(K(i));
    %         [~,D] = eig( A, B );
    %         h_lb = sort(real(im.hFromLambda(diag(D))),'descend');
    %
    %         while (~isConverged)
    %             im.nEVP = im.nEVP + dEVP;
    %             [A,B] = im.eigenmatricesForWavenumber(K(i));
    %             [~,D] = eig( A, B );
    %             h_ub = sort(real(im.hFromLambda(diag(D))),'descend');
    %
    %             rel_error = abs((h_ub(1:nModes)-h_high_res(1:nModes))./h_high_res(1:nModes));
    %             isConverged = max(rel_error) < 1e-8;
    %         end
    %         im.nEVP = im.nEVP;
    %         minEVP(i) = im.nEVP;
    %     end
    % end

    save(cachePath,"k","nModes","K","MODES","minEVP");
end

%%
fit = fit_fkn(K(:), MODES(:), minEVP(:));
fit.theta

%%
shouldShowModelFit = false;

nModes = reshape(MODES(1,:),[],1);
k = K(:,1);
k_cpm = k/(2*pi);

model = @(N,k) 3.07 * N.^0.88 .*(1 + (k./(0.00062*N)) ).^0.3058;
% model = @(N,k) 3 * N .*(1 + (k./(0.00062*N)) ).^0.33;
model = @(N,k) fit.model(k,N);

figure
tiledlayout(2,1)

nexttile
labels = cell(size(MODES,2),1);
for iMode = 1:size(MODES,2)
    plot(K(:,iMode)/2/pi,minEVP(:,iMode)), hold on
    labels{iMode} = string(MODES(1,iMode));
end
xlog, ylog
xlabel("k (cycles/m)")
ylabel("minimum EVP points")
legend(labels,Location="northwest")
title("minimum number of EVP points vs k, for a required number of resolved modes")

if shouldShowModelFit
    for iMode = 1:size(MODES,2)
        plot(k_cpm,model(nModes(iMode),k),Color=0*[1 1 1], LineWidth=2)
    end
end

% for iMode = 1:size(MODES,2)
%     fit = fitBrokenPowerLaw_a0(K(:,iMode),minEVP(:,iMode));
%     fit.theta
%     plot(k_cpm,fit.model(k),Color=0*[1 1 1], LineWidth=2)
% end

% iK = size(K,1);
% iL = iK-3;
% iMode = 1;
% slope_k = ( log(minEVP(iL,iMode))-log(minEVP(iK,iMode)) )/( log(K(iL,iMode))-log(K(iK,iMode)) );
% k_star = 2*pi/3.1e2;
% delta = 1e-3;
% m = @(k) (slope_k/2)*(1+tanh((k-k_star)/delta));
%
%
% plot(k_cpm,100*k.^m(k),Color=0*[1 1 1], LineWidth=2)

nexttile
labels = cell(size(K,1),1);
for iK = 1:size(K,1)
    plot(MODES(iK,:),minEVP(iK,:)), hold on
    labels{iK} = string(round(2*pi/K(iK,1)));
end
xlog, ylog
xlabel("resolved modes")
ylabel("minimum EVP points")
legend(labels,Location="northwest")
title("minimum number of EVP points vs required number of resolved modes, for a given wavelength")

if shouldShowModelFit
    for iK = 1:size(K,1)
        plot(nModes,model(nModes,k(iK)),Color=0*[1 1 1], LineWidth=2)
    end
end


% iK = size(K,1);
% slope_mode = ( log(minEVP(iK,1))-log(minEVP(iK,end)) )/( log(MODES(iK,1))-log(MODES(iK,end)) );
% C = minEVP(iK,1)/(MODES(iK,1)^slope_mode);
% plot(nModes,1.1*C*nModes.^slope_mode,Color=0*[1 1 1], LineWidth=2)
%
% k0 = [0.0138; 0.0237; 0.0481; 0.0790];

% profile on; [F,G,h] = im.modesAtQuadraturePoints(nPoints=64,k=2*pi/1e3); profile viewer;

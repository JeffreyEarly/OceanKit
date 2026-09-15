%% Cosine half-window eigenmodes
% This standalone analytic benchmark compares three ways to form modes from
% a partial-depth window.
%
% The full-depth orthonormal cosine basis e_n is used for numerical clarity.
% To mimic hydrostatic F modes with full-depth norms h_j, we also define
% named modes phi_j = sqrt(h_j) e_j with h_j = 1/(j+1)^2.
%
% The three constructions are:
%
% 1. Visibility eigenmodes: diagonalize the full-depth-normalized window
%    matrix M. Eigenvalues are window-energy fractions.
% 2. Visible subspace, h-ordered: keep the visible eigenvectors of M, then
%    diagonalize H = diag(h_j) inside that visible subspace.
% 3. Raw Gamma_win modes: diagonalize Gamma_win for the named modes without
%    normalizing by H. Eigenvalues are absolute partial-window h-like norms.

D = 1;
nModes = 32;
nPlot = 2001;
nShow = 6;
visibilityThreshold = 0.5;

%% Build full-depth and window matrices
M = cosineHalfWindowMatrix(nModes);
M = (M + M.')/2;

modeNumber = (0:(nModes-1)).';
hProxy = 1./(modeNumber + 1).^2;
HProxy = diag(hProxy);
sqrtHProxy = diag(sqrt(hProxy));
GammaWinRaw = sqrtHProxy*M*sqrtHProxy;

%% Construction 1: visibility eigenmodes
[QVisibility,LambdaMatrix] = eig(M);
lambda = real(diag(LambdaMatrix));
[lambda,sortIndex] = sort(lambda,'descend');
QVisibility = orientColumns(real(QVisibility(:,sortIndex)));

%% Construction 2: visible subspace ordered by proxy equivalent depth
visibleIndex = lambda >= visibilityThreshold;
UVisible = QVisibility(:,visibleIndex);

HVisible = UVisible.'*HProxy*UVisible;
[SVisible,XiMatrix] = eig(HVisible);
xi = real(diag(XiMatrix));
[xi,xiSortIndex] = sort(xi,'descend');
SVisible = real(SVisible(:,xiSortIndex));

QHOrdered = orientColumns(UVisible*SVisible);
effectiveVisibilityH = diag(QHOrdered.'*M*QHOrdered);
windowGramHOrdered = QHOrdered.'*M*QHOrdered;

%% Construction 3: raw Gamma_win eigenmodes
% VRaw contains coefficients in the named basis phi_j = sqrt(h_j) e_j.
% The plotted modes are scaled to unit full-depth norm after diagonalizing
% Gamma_win, because the raw eigenvectors are not full-depth orthonormal.
[VRaw,MuMatrix] = eig(GammaWinRaw);
muRaw = real(diag(MuMatrix));
[muRaw,muSortIndex] = sort(muRaw,'descend');
VRaw = orientColumns(real(VRaw(:,muSortIndex)));

rawFullDepthNorm = diag(VRaw.'*HProxy*VRaw);
QRawNormalized = sqrtHProxy*VRaw*diag(1./sqrt(rawFullDepthNorm));
rawVisibility = diag(QRawNormalized.'*M*QRawNormalized);
rawFullDepthGram = QRawNormalized.'*QRawNormalized;
rawWindowGram = VRaw.'*GammaWinRaw*VRaw;

%% Validate the finite-band constructions
x = linspace(0,D,nPlot).';
E = cosineBasis(x,D,nModes);
nShowH = min(nShow,size(QHOrdered,2));
nShowRaw = min(nShow,size(QRawNormalized,2));

PsiVisibility = E*QVisibility(:,1:nShow);
PsiHOrdered = E*QHOrdered(:,1:nShowH);
PsiRaw = E*QRawNormalized(:,1:nShowRaw);

windowIndex = x >= D/2;
xWindow = x(windowIndex);
PsiVisibilityWindow = PsiVisibility(windowIndex,:);
PsiHOrderedWindow = PsiHOrdered(windowIndex,:);
PsiRawWindow = PsiRaw(windowIndex,:);

numericalVisibilityGram = zeros(nShow,nShow);
numericalHVisibility = zeros(nShowH,1);
numericalRawVisibility = zeros(nShowRaw,1);
for i = 1:nShow
    for j = 1:nShow
        numericalVisibilityGram(i,j) = trapz(xWindow,PsiVisibilityWindow(:,i).*PsiVisibilityWindow(:,j));
    end
end
for i = 1:nShowH
    numericalHVisibility(i) = trapz(xWindow,PsiHOrderedWindow(:,i).*PsiHOrderedWindow(:,i));
end
for i = 1:nShowRaw
    numericalRawVisibility(i) = trapz(xWindow,PsiRawWindow(:,i).*PsiRawWindow(:,i));
end

eigenvalueBoundError = max([max(lambda - 1), max(-lambda), 0]);
visibilityOrthogonalityError = norm(QVisibility.'*QVisibility - eye(nModes),2);
hOrderedOrthogonalityError = norm(QHOrdered.'*QHOrdered - eye(size(QHOrdered,2)),2);
visibilityQuadratureError = norm(numericalVisibilityGram - diag(lambda(1:nShow)),2);
hVisibilityQuadratureError = norm(numericalHVisibility - effectiveVisibilityH(1:nShowH),2);
rawVisibilityQuadratureError = norm(numericalRawVisibility - rawVisibility(1:nShowRaw),2);
hDiagonalizationError = norm(QHOrdered.'*HProxy*QHOrdered - diag(xi),2);
hOrderedWindowCross = offDiagonalNorm(windowGramHOrdered);
rawWindowDiagonalizationError = norm(rawWindowGram - diag(muRaw),2);
rawFullDepthCross = offDiagonalNorm(rawFullDepthGram);

if mod(nModes,2) == 0
    pairingError = max(abs(lambda + flipud(lambda) - 1));
else
    pairingError = NaN;
end

[SFull,XiFull] = eig(HProxy);
[~,fullSortIndex] = sort(diag(XiFull),'descend');
SFull = SFull(:,fullSortIndex);
fullDepthRecoveryError = norm(abs(SFull) - eye(nModes),2);

fprintf('Cosine half-window comparison, nModes = %d\n', nModes);
fprintf('Visible subspace threshold: %.3f, visible modes: %d\n', visibilityThreshold, size(QHOrdered,2));
fprintf('Visibility eigenvalue bounds: %.3e\n', eigenvalueBoundError);
fprintf('Visibility basis orthogonality: %.3e\n', visibilityOrthogonalityError);
fprintf('H-ordered basis orthogonality:  %.3e\n', hOrderedOrthogonalityError);
fprintf('Visibility quadrature check:    %.3e\n', visibilityQuadratureError);
fprintf('H-ordered visibility check:     %.3e\n', hVisibilityQuadratureError);
fprintf('Raw Gamma visibility check:     %.3e\n', rawVisibilityQuadratureError);
fprintf('H diagonalization error:        %.3e\n', hDiagonalizationError);
fprintf('H-ordered window cross term:    %.3e\n', hOrderedWindowCross);
fprintf('Raw Gamma_win diagonalization:  %.3e\n', rawWindowDiagonalizationError);
fprintf('Raw full-depth cross term:      %.3e\n', rawFullDepthCross);
fprintf('Full-depth H recovery error:    %.3e\n', fullDepthRecoveryError);
fprintf('Half-window pairing error:      %.3e\n', pairingError);
fprintf('Leading visibility lambda:      %s\n', mat2str(lambda(1:nShow).',6));
fprintf('Leading h-ordered xi:           %s\n', mat2str(xi(1:nShowH).',6));
fprintf('Leading h-ordered visibility:   %s\n', mat2str(effectiveVisibilityH(1:nShowH).',6));
fprintf('Leading raw Gamma_win mu:       %s\n', mat2str(muRaw(1:nShowRaw).',6));
fprintf('Leading raw Gamma visibility:   %s\n', mat2str(rawVisibility(1:nShowRaw).',6));

if eigenvalueBoundError > 1e-12
    error("CosineHalfWindowEigenmodes:EigenvalueBounds", "Window visibility eigenvalues should lie in [0,1] up to roundoff.");
end
if visibilityOrthogonalityError > 1e-12
    error("CosineHalfWindowEigenmodes:VisibilityOrthogonality", "Visibility eigenvectors are not orthonormal.");
end
if visibilityQuadratureError > 1e-4
    error("CosineHalfWindowEigenmodes:VisibilityQuadrature", "Numerical quadrature does not match the analytic visibility eigenvalues.");
end
if hOrderedOrthogonalityError > 1e-12
    error("CosineHalfWindowEigenmodes:HOrderedOrthogonality", "Equivalent-depth ordered visible modes are not full-depth orthonormal.");
end
if hDiagonalizationError > 1e-12
    error("CosineHalfWindowEigenmodes:HDiagonalization", "Visible modes do not diagonalize the proxy equivalent-depth operator.");
end
if rawWindowDiagonalizationError > 1e-12
    error("CosineHalfWindowEigenmodes:RawWindowDiagonalization", "Raw Gamma_win modes do not diagonalize Gamma_win.");
end

%% Plot the three constructions
figure('Name','Cosine half-window mode comparison')
tiledlayout(3,2,'TileSpacing','compact','Padding','compact')

nexttile
stem(1:nModes,lambda,'filled','DisplayName','visibility fraction'), hold on
stem(1:nModes,muRaw/max(hProxy),'DisplayName','raw \Gamma_{win} \mu / max(h)')
yline(visibilityThreshold,'r--','visibility threshold','DisplayName','threshold')
xlabel('mode index')
ylabel('normalized value')
title('Visibility fractions and raw partial h')
legend('Location','southwest')
ylim([-0.05 1.05])
grid on

nexttile
semilogy(0:(nModes-1),hProxy,'k-o','DisplayName','full-depth h'), hold on
semilogy(1:nModes,max(muRaw,realmin),'r-o','DisplayName','raw window \mu')
xlabel('mode index')
ylabel('h-like value')
title('Raw \Gamma_{win} recovers h ordering at full depth')
legend('Location','southwest')
grid on

nexttile
plotModesWithWindow(x,PsiVisibility,D,"Visibility-ordered window modes",compose('\\alpha=%d, \\lambda=%.3f',(1:nShow).',lambda(1:nShow)))

nexttile
plotModesWithWindow(x,PsiHOrdered,D,"Visible modes ordered by h",compose('\\alpha=%d, \\xi=%.3f',(1:nShowH).',xi(1:nShowH)))

nexttile
plotModesWithWindow(x,PsiRaw,D,"Raw \Gamma_{win} modes",compose('\\alpha=%d, \\mu=%.3g',(1:nShowRaw).',muRaw(1:nShowRaw)))

nexttile
imagesc(1:nShowRaw,0:(nModes-1),QRawNormalized(:,1:nShowRaw))
set(gca,'YDir','normal')
colorbar
xlabel('raw \Gamma_{win} mode index')
ylabel('full-depth cosine mode n')
title('Raw mode coefficients after full-depth normalization')

drawnow

%% Local helpers
function M = cosineHalfWindowMatrix(nModes)
M = zeros(nModes,nModes);
M(1,1) = 1/2;

for n = 1:(nModes-1)
    value = -sqrt(2)*sin(n*pi/2)/(n*pi);
    M(1,n+1) = value;
    M(n+1,1) = value;
end

for m = 1:(nModes-1)
    for n = 1:(nModes-1)
        if m == n
            M(m+1,n+1) = 1/2;
        else
            differenceTerm = sin((m-n)*pi/2)/(m-n);
            sumTerm = sin((m+n)*pi/2)/(m+n);
            M(m+1,n+1) = -(differenceTerm + sumTerm)/pi;
        end
    end
end
end

function E = cosineBasis(x,D,nModes)
E = zeros(length(x),nModes);
E(:,1) = 1/sqrt(D);

for n = 1:(nModes-1)
    E(:,n+1) = sqrt(2/D)*cos(n*pi*x/D);
end
end

function Q = orientColumns(Q)
for iMode = 1:size(Q,2)
    [~,iMax] = max(abs(Q(:,iMode)));
    Q(:,iMode) = sign(Q(iMax,iMode))*Q(:,iMode);
end
end

function value = offDiagonalNorm(A)
value = norm(A - diag(diag(A)),'fro');
end

function plotModesWithWindow(x,Psi,D,plotTitle,legendLabels)
windowColor = [0.92 0.92 0.92];
yLimits = [min(Psi(:)) max(Psi(:))];
patch([D/2 D D D/2],[yLimits(1) yLimits(1) yLimits(2) yLimits(2)],windowColor,'EdgeColor','none')
hold on
plot(x,Psi,'LineWidth',1.2)
xlabel('x')
ylabel('\psi(x)')
title(plotTitle)
legend(legendLabels,'Location','eastoutside')
grid on
end

function transform = discreteTransform(self,grid,options)
    % Build a fixed-grid F/G transform from the retained rigid hydrostatic basis.
    arguments
        self IMAnalyticalInternalModesBasis
        grid (1,1) IMVerticalGrid
    options.recoveryF = []
    options.recoveryG = []
    options.forwardOperator (1,1) string = ""
    end
    if ~isequal(grid.zDomain,self.zDomain)
        error("IMAnalyticalInternalModesBasis:IncompatibleDomain","The basis and grid must have the same zDomain.");
    end
    rigid = true;
    for boundary = [self.evp.surfaceBoundary,self.evp.bottomBoundary]
        if self.evp.formulation == "F"
            rigid = rigid && boundary.a == 0 && boundary.b ~= 0 && boundary.c == 0 && boundary.d == 0;
        else
            rigid = rigid && boundary.a ~= 0 && boundary.b == 0 && boundary.c == 0 && boundary.d == 0;
        end
    end
    if self.evp.modeFamily ~= "hydrostatic" || ~rigid || any(self.h < 0)
        error("IMAnalyticalInternalModesBasis:UnavailableDiscreteTransform","The fixed-grid analytical bridge requires a rigid hydrostatic F/G basis including its barotropic sector.");
    end
    F = self.F(grid.z); G = self.G(grid.z);
    HF = self.gramMatrix(variable="F");
    HG = self.gramMatrix(variable="G");
    majorantF = self.majorantGramMatrix(variable="F");
    majorantG = self.majorantGramMatrix(variable="G");
    N2Values = self.N2(grid.z);
    channel = @(interior,target,majorant,active) struct(available=true,reason="",activeModeMask=active,interiorWeight=interior,endpointMetricMatrix=zeros(grid.nPoints),metricMatrix=diag(interior.*grid.weights),targetGramMatrix=target,targetMajorantGramMatrix=majorant);
    channels = struct(F=channel(ones(size(grid.z)),HF,majorantF,true(1,self.nColumns)),G=channel(N2Values/self.evp.g,HG,majorantG,isfinite(self.h)));
    if ~isempty(options.recoveryF), channels.F.recovery = options.recoveryF; end
if ~isempty(options.recoveryG), channels.G.recovery = options.recoveryG; end
if options.forwardOperator ~= ""
    channels.F.forwardOperator = options.forwardOperator;
    channels.G.forwardOperator = options.forwardOperator;
end
transformMetadata = struct(construction="IMAnalyticalInternalModesBasis.discreteTransform",columnLabels=self.columnLabels,relation="F=h Gz with F0=1 and G0=0");
    transform = IMInternalModesDiscreteTransform(z=grid.z,weights=grid.weights,modeNumber=self.modeNumber,h=self.h,normalization=string(self.normalization),inverseF=F,inverseG=G,endpointF=self.F(self.zDomain([2 1]).'),endpointG=self.G(self.zDomain([2 1]).'),channelData=channels,primaryVariable="F",zDomain=self.zDomain,g=self.evp.g,modeFamily="hydrostatic",N2Values=N2Values,problemMetadata=transformMetadata,basis=self);
end

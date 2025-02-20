
function  SDF1 = OptSDF(splineint,moments,interiorknots )
% INPUTS:   splineint=matrix of integrated spline bases
%           moments=number of moment conditions used
%           interiorknots = set of interior knots used to create spline
% OUTPUTS:  f=spline coefficients with identity weights 
%           g=spline coefficients with 2nd stage weight matrix 
%
knots=augknt(interiorknots,4); %augmented knot sequence which contains the first and last knot 
                           % with exact multiplicity k (k=4 for cubic-b-splines)
numbases=size(interiorknots,2)+2; %number of basis functions (note interiorknots actually contains 2 non-interior points)  
Mspline=spmak(knots,eye(numbases));


randn('seed',1234);
rand('seed',0);
%--------------------------------------------------------------------------
%  first stage gmm estimate
%--------------------------------------------------------------------------
startvec=ones(numbases,1);
funkie=@(x) gmmmoments(splineint,x,eye(moments),moments); %number of moments


A = [];
b = [];
Aeq = [];
beq = [];
lb = -100*ones(numbases,1);
ub = 100*ones(numbases,1);

xmin1 = fmincon(funkie,startvec,A,b,Aeq,beq,lb,ub);


range=min(interiorknots):.001:max(interiorknots);
splinefullrange=(spval(Mspline,range))';
invSDF1=splinefullrange*xmin1;
SDF1=invSDF1.^(-1);


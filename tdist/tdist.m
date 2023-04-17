function [yfun,xfun,result] = tdist(X,df,lambda,funtype,options)
%TDIST  Computes the distribution (PDF, CDF or QF - quantiles) 
%       of a linear combination of independent SYMMETRIC 
%       (zero-mean) random variables (RVs) with specific
%       distributions:
%       - STUDENT's t distribution with 0 < df < Inf,
%       - NORMAL distribution N(0,1), 
%       - symmetric RECTANGULAR (uniform) distribution R(-1,1), 
%       - symmetric TRIANGULAR distribution T(-1,1), 
%       - symmetric ARCSINE distribution (U-distribution) U(-1,1),
%
%       TDIST is based on numerical inversion of the
%       characteristic function (Gil-Pelaez method). 
%       The required integration is performed by the 14-points
%       Gaussian quadrature over N subintervals of [0, 10*pi].
%           
% SYNTAX: 
% [yfun,xfun,results] = TDIST(X,df,lambda,funtype,options)
%
% INPUT: 
% X       - vector of appropriate (funtype) x values
%           if X = [] then xfun is generated automatically
% df      - vector of degrees of freedom of independent
%           STUDENT's t RVs, t_df, for 0 < df < Inf.
%           Further, 
%           set df = Inf for the NORMAL RVs, N(0,1),
%           set df = -1  for the RECTANGULAR RVs, R(-1,1),
%           set df = -2  for the symmetric TRIANGULAR RVs, T(-1,1),
%           set df = -3  for the symmetric ARCSINE RVs, U(-1,1),
%           set df < -10  for the symmetric mixture of CHI2 and -CHI2 RVs
%           with nu = abs(df+10) degrees of freedom
% lambda  - vector of coefficients of the linear combination 
% funtype - default value is 1 (calculates the CDF)
%           The following funtypes are legible:
%           0: TDIST calculates CDF and PDF, yfun = [CDF,PDF].
%           1: TDIST calculates the cumulative distribution
%              function, CDF at X, yfun = CDF.
%           2: TDIST calculates the probability density function,
%              PDF, at X, yfun = PDF.
%           3: TDIST calculates the quantile function, 
%              QF at X, yfun = QF.
%           4: TDIST calculates the characteristic function,
%              CHF at X, yfun = CHF.
% options - structure with the following parameters:
% Tmax    - Upper limit for integration interval [0, Tmax]. The default 
%           value is Tmax = 10 * pi;
% SixSigmaCoef - coefficient for setting the (approximate) support 
%           [-SixSigmaCoef * STD, -SixSigmaCoef * STD], where
%           STD is the calculate standard deviation of the distribution 
%           (if it exists). The default value is SixSigmaCoef = 6.
% N       - Number of subintervals of the interval [0, Tmax] for 
%           integration by the 14-points Gaussian quadrature.
%           The default value is N = 2^7;
% n       - Number of automatically generated X values (when the input 
%           is X = []). The default value is n = 2^7;
% isPlot  - logical indicator for plotting the PDF, default
%           value is isPlot = true.
% isVerbose - logical indicator for presenting structure with further 
%           detailed results, default value is isVerbose = true.
% isChebfun - logical indicator for generating X values (when the input 
%           is X = []) in the Chebyshev points. This allows to 
%           present the resulted function alternatively as a chebfun 
%           (for more details see the Chebfun Guide). The default value 
%           is isChebfun = false.
%
% OUTPUT:
% yfun    - column vector with calculated function values, 
%           the result depends on funtype.
%           If funtype = 0, yfun has two columns (CDF and PDF).
% xfun    - column vector of function input values
%           Typically xfun = X.
%           If X = [], xfun is generated automatically.
%
% EXAMPLE 1: (CDF of a linear combination of RVs defined by df)
% % (Normal, Student's t_1, Rectangular, Triangular, and U-distribution)
% df      = [Inf 1 -1 -2 -3];
% lambda  = [1 1 5 1 10];
% funtype = 1;
% [cdf,x] = tdist([],df,lambda,funtype);
%
% EXAMPLE 2: (PDF of a linear combination of RVs defined by df)
% df      = [Inf 1 -1 -2 -3];
% lambda  = [1 1 5 1 10];
% funtype = 2;
% [pdf,x] = tdist([],df,lambda,funtype);
%
% EXAMPLE 3: (QF of a linear combination of RVs defined by df)
% options.isPlot = false;
% df      = [Inf 1 -1 -2 -3];
% lambda  = [1 1 5 1 10];
% funtype = 3;
% prob    = [0.9 0.95 0.99]';
% qf      = tdist(prob,df,lambda,funtype,options);
% disp([prob qf]);
%
% EXAMPLE 4: (CHF of a linear combination of RVs defined by df)
% options.isPlot = true;
% df      = [Inf 1 -1 -2 -3];
% lambda  = [1 1 5 1 10];
% funtype = 4;
% t = linspace(0,pi);
% chf = tdist(t,df,lambda,funtype,options);
%
% EXAMPLE 5 (Create PDF as a CHEBFUN function and use it to compute CDF)
% df_true =  [1 2 3 10];
% df      = -10 - df_true;    % symmetric chi2-mixture distributions
% lambda  = [1 1 1 1];
% funtype = 2;
% N = 2^10;
% xmax = 130;
% x = -xmax * cos((0:N)*pi/N);
% f = tdist(x,df,lambda,funtype);
% pdf = chebfun(f,[-xmax,xmax]);
% integrate = sum(pdf);
% cdf = cumsum(pdf);
% xnew = linspace(-50,50);
% figure
% plot(xnew,cdf(xnew))
%
% EXAMPLE 6: (Generate a CHEBFUN from the QF)
% options.isChebfun = true;
% options.n = 2^8;
% options.N = 2^9;
% df      = [Inf 1 -1 -2 -3];
% lambda  = [1 1 5 1 10];
% funtype = 3;
% [qf,prob,results] = tdist([],df,lambda,funtype,options);
% QF = results.chebfun;
% disp(QF([0.9 0.95 0.99]'))
%
% The algorithm requires evaluation of the BesselK and BesselJ functions.
%  
% REFERENCES:
%
% 1. GIL-PELAEZ, J. Note on the inversion theorem. Biometrika 38
%    (1951), 481482.
% 2. WITKOVSKY , V. On the exact computation of the density
%    and of the quantiles of linear combinations of t and F
%    random variables. Journal of Statistical Planning and
%    Inference 94 (2001), 113.
% 3. WITKOVSKY , V. Matlab algorithm TDIST: The distribution of a
%    linear combination of Students t random variables. In
%    COMPSTAT 2004 Symposium (2004), J. Antoch, Ed.,
%    Physica-Verlag/Springer 2004, Heidelberg, Germany,
%    19952002.
% 4. DRISCOLL, T. A., HALE, N., TREFETHEN, L. N. Chebfun Guide. Pafnuty 
%    Publications, Oxford, 2014.

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 24-Jan-2017 12:54:29

%% CHECK THE INPUT PARAMETERS
narginchk(2, 5);

if nargin < 5, options = []; end
if nargin < 4, funtype = []; end
if nargin < 3, lambda  = []; end

if ~isfield(options, 'Tmax')
    options.Tmax = 10 * pi;
end

if ~isfield(options, 'SixSigmaCoef')
    options.SixSigmaCoef = 6;
end

if ~isfield(options, 'n')
    options.n = 2^7;
end

if ~isfield(options, 'N')
    options.N = 2^7;
end

if ~isfield(options, 'minprob')
    options.minprob = 0.001;
end

if ~isfield(options, 'isPlot')
    options.isPlot = true;
end

if ~isfield(options, 'isVerbose')
    options.isVerbose = true;
end

if ~isfield(options, 'isChebfun')
    options.isChebfun = false;
end

if isempty(funtype), funtype = 1; end
if isempty(lambda), lambda = 1; end
if isempty(df), warning('DF should be specified'); end

switch lower(funtype)
    case 0
        funtype = 0;
    case {'cdf', 'c', 1}
        funtype = 1;
    case {'pdf', 'p', 2}
        funtype = 2;
    case {'qf', 'quantile', 'q', 'inv', 3}
        funtype = 3;
    case {'chf', 'cf', 'characteristic', 4}
        funtype = 4;
    otherwise
        warning('Unexpected funtype.');
end

df = df(:); lambda = lambda(:); X = X(:);

if length(lambda) == 1
    lambda = lambda * ones(size(df));
elseif length(lambda) ~= length(df)
    error('Dimension mismatch.');
end

% Exclude RVs with zero coefficients or wrong df
lambdaOK   = (lambda ~= 0);
lambda = lambda(lambdaOK);
df     = df(lambdaOK);

dfOK   = (df == Inf | df > 0);
dfOK   = dfOK | (df == -1 | df == -2 | df == -3);
dfOK   = dfOK | (df < 10);
df     = df(dfOK);
lambda = lambda(dfOK);

% If df > 100 set to be standard normal random variables
dfLarge   = (df > 100 & df ~= Inf);
if any(dfLarge)
    df(dfLarge) = Inf;
    warning('DF > 100 was changed to DF = Inf.');
end

% Estimate the (approximate) range
condt1 = (df > 0 & df <= 1);
condt2 = (df > 1 & df <= 2);
condt3 = (df > 2 & df <= 3);
condt4 = (df > 3 & df < Inf);
condN  = (df == Inf);
condR  = (df == -1);
condT  = (df == -2);
condU  = (df == -3);
condChi= (df < -10);

var = 0;
var = var + sum(15*lambda(condt1).^2);    % surrogate for var = 15
var = var + sum(7.5*lambda(condt2).^2);   % surrogate for var = 7.5
var = var + sum(3.5*lambda(condt3).^2);   % surrogate for var = 3.5 
var = var + sum(lambda(condt4).^2 .* ...  % var = df/(df-2)
    (df(condt4)./(df(condt4)-2)));
var = var + sum(lambda(condN).^2);        % var = 1
var = var + sum(lambda(condR).^2/3);      % var = 1/3
var = var + sum(lambda(condT).^2/6);      % var = 1/6
var = var + sum(lambda(condU).^2/2);      % var = 1/2
var = var + sum(lambda(condChi).^2 .* ... % var = 2nd raw-moment of Chi2
    (4 * gamma(abs(df(condChi)+10) + 2) ./ gamma(abs(df(condChi)+10))));

xmax = options.SixSigmaCoef * sqrt(var);

symmetry = 0;
if isempty(X)
    switch funtype
        case 3
            if options.isChebfun
                minprob = options.minprob;
                X = 0.5 - (0.5-minprob) * ...
                    cos((0:options.n)'*pi/(options.n));
                low = 0 + minprob;
                upp = 1 - minprob;
            else
                X = linspace(0,1,options.n)';
                low = 0;
                upp = 1;
            end
        case 4
            if options.isChebfun
                X = -(3*pi) * cos((0:2*options.n)'*pi/(2*options.n));
                low = -3*pi;
                upp = 3*pi;
            else
                X = linspace(0,3*pi,options.n)';
                low = 0;
                upp = 3*pi;
            end
        otherwise
            if options.isChebfun
                symmetry = 1;
                X = -xmax * cos((0:2*options.n)*pi/(2*options.n));
                X = X((options.n+1):(2*options.n+1))';
                low = -xmax;
                upp = xmax;
            else
                symmetry = 1;
                X = linspace(0,xmax,options.n)';
                low = -xmax;
                upp = xmax;
            end
    end
end
xfun = X;

norml = sqrt(lambda'*lambda);
yfun  = [];

result.options = options;
if options.isVerbose && nargout == 3
    result.fun     = [];
    result.x       = [];
    result.chebfun = [];
    result.xmax    = xmax;
    result.df      = df;
    result.lambda  = lambda;
    result.N       = options.N;
    result.n       = options.n;
    result.minprob = options.minprob;
    result.var     = var;
    result.norml   = norml;
    result.xlow    = low;
    result.xupp    = upp;
end

%% ALGORITHM 
% Set the Gaussian-Quadrature nodes and weights
N       = options.N;
dt      = (10*pi) / N;
Tmax    = options.Tmax;
dt2     = pi / xmax;
Tmax2   = 2^9 * dt2;

% Tmax    = Tmax2;
% dt      = dt2;
% For better precision of the results, (for large x values),
% the first subinterval [0, dt] is further splitted to 9 subintervals
limits  = [0; dt/100; dt/10; dt; Tmax];
subInts = [3; 3; 3; N-1];

% Alternatively, use the simple equidistant splitting to N subintervals
%limits  = [0; 10*pi];
%subInts = N;

% Evaluate the characteristic function (CHF) at X     
if funtype == 4
    yfun  = tchfvw(X,df,lambda);
    weightedChf = [];
    t           = [];
    weights     = [];
else
    [t,weights] = gweights(limits,subInts);
    chf = tchfvw(t,df,lambda/norml);
    weightedChf  = weights .* chf;
end

% Evaluate the required function (CDF, PDF or QF) at X
if (funtype == 0 || funtype == 1 || funtype == 2)
    yfun = tcdfpdf(X,t,weightedChf,funtype,norml);
    xfun = X;
    if symmetry == 1
        xfun  = [-flip(X);X(2:end)];              
        if funtype == 0
            cdf   = yfun(:,1);
            pdf   = yfun(:,2);
            yfun  = zeros(length(xfun),2);
            yfun(:,1) = [1 - flip(cdf) ; cdf(2:end)];
            yfun(:,2) = [flip(pdf) ; pdf(2:end)];
        elseif funtype == 1
            cdf   = yfun;
            yfun  = [1 - flip(cdf) ; cdf(2:end)];
        elseif funtype == 2
            pdf   = yfun;
            yfun  = [flip(pdf) ; pdf(2:end)];
        end
    end
elseif funtype == 3
    [yfun,xfun] = tinvvw(X,t,weightedChf,norml);
end

if options.isVerbose  && nargout == 3
    if options.isChebfun
        FUN = chebfun(yfun,[low,upp]);
    else
        FUN = [];
    end
    result.Tmax      = Tmax;
    result.Tmax2     = Tmax2;
    result.dt        = dt;
    result.dt2       = dt2;
    result.limits    = limits;
    result.subInts   = subInts;
    result.t         = t;
    result.weights   = weights;
    result.weightChf = weightedChf;
    result.fun       = yfun;
    result.x         = xfun;
    result.chebfun   = FUN;
end

%% PLOT
if options.isPlot
    plot(xfun,yfun)
    grid on
    xlabel('x')
    switch funtype
        case 0
            ylabel('CDF and PDF')
        case 1
            ylabel('CDF')
        case 2
            ylabel('PDF')
        case 3
            xlabel('prob')
            ylabel('Quantile function (QF)')
        case 4
            xlabel('t')
            ylabel('Characteristic function (CHF)')
    end
end

%% FUNCTION TCHFVW
function  [chfout,tout] = tchfvw(t,nu,l)
%TCHFVW   Evaluates the characteristic function (CHF) for linear 
%         combination of independent random variables, RVs, with
%         specific symmetrical zero mean distributions:
%         STUDENT'S T DISTRIBUTION with 0 < df < Inf
%         NORMAL DISTRIBUTION, N(0,1)
%         RECTANGULAR DISTRIBUTION, R(-1,1)
%         symmetric TRIANGULAR DISTRIBUTION, T(-1,1)
%         symmetric ARCSINE DISTRIBUTION (U-distribution), U(-1,1)
           
% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 29-Nov-2014 18:20:52

chf = ones(size(t));
tout = t;
chfout = chf;

idOK = abs(t) > 50*eps;
t = t(idOK);
chf = chf(idOK);

% STUDENT'S T DISTRIBUTION
idxt = find(nu > 0 & nu < Inf); 
lt = l(idxt);
df = nu(idxt);
if ~isempty(idxt)
    for k = 1:length(lt)
        chf = chf .* besselk(df(k)/2, abs(lt(k)*t).*sqrt(df(k)),1) ...
          .* exp(-abs(lt(k)*t).*sqrt(df(k))) ...
          .* (sqrt(df(k)).*abs(lt(k)*t)).^(df(k)/2) ...
          / 2^(df(k)/2-1)/gamma(df(k)/2);
    end
end

% NORMAL DISTRIBUTION N(0,1)
idxN = find(nu == Inf); 
lN = l(idxN);
if ~isempty(idxN)
    for k = 1:length(lN)
        chf = chf .* exp(-(lN(k) * t).^2 / 2);
    end
end

% RECTANGULAR DISTRIBUTION over [-1,1]
idxU = find(nu == -1); 
lU = l(idxU);
if ~isempty(idxU)
    for k = 1:length(lU)
        chf = chf .* sin(lU(k) * t) ./ (lU(k) * t);
    end
end

% SYMMETRIC TRIANGULAR DISTRIBUTION over [-1,1]
idxT = find(nu == -2); 
lT = l(idxT);
if ~isempty(idxT)
    for k = 1:length(lT)
        chf = chf .* (2-2*cos((lT(k)*t))) ./ (lT(k)*t).^2;
    end
end

% SYMMETRIC ARCSINE DISTRIBUTION over [-1,1]
idxT = find(nu == -3); 
lT = l(idxT);
if ~isempty(idxT)
    for k = 1:length(lT)
        chf = chf .* besselj(0,lT(k)*t);
    end
end

% SYMMETRIC MIXTURE OF CHI2 AND -CHI2 DISTRIBUTION WITH NU DFs
idxChi = find(nu < -10); 
dfChi  = abs(nu(idxChi)+10);
lChi   = l(idxChi);
if ~isempty(idxChi)
    for k = 1:length(lChi)
        chf = chf .* real((1-2*1i*lChi(k)*t).^(-dfChi(k)/2));
    end
end

chfout(idOK) = chf;

%% FUNCTION TCDFPDF 
function  yfun=tcdfpdf(x,t,weightedChf,funtype,normalize)
%TCDFPDF  Computes the CDF and/or PDF by numerical inversion of
%         the characteristic function

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 29-Nov-2014 18:20:52

% Integration by the GAUSS QUADRATURE
x = x / normalize;
xsize = length(x);

pdf = zeros(xsize,1);
cdf = zeros(xsize,1);
wft = weightedChf ./ t;
for i = 1:xsize
    pdf(i) = weightedChf'  * cos(x(i)*t);
    cdf(i) = wft' * sin(x(i)*t);
end
pdf  = max(0,pdf(:))/pi / normalize;
cdf  = min(1,max(0,1/2+cdf(:)/pi));

switch funtype
    case 1
        yfun = cdf;
    case 2
        yfun = pdf;
    otherwise
        yfun = [cdf pdf];
end
 
%% FUNCTION TINVVW
function [yfun,xfun] = tinvvw(pr,t,weightedChf,normalize)
%TINVVW Computes the quantile function by the Newton's method

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 29-Nov-2014 18:20:52

q = zeros(size(pr));
q(pr == 0) = -Inf;
q(pr == 1) = Inf;
q(pr < 0 | pr > 1) = NaN;

ok = find(pr > 0 & pr < 1);
if isempty(ok)
    yfun = q;
    return
end
prob = pr(ok);

% Newton's iterative method
maxiter   = 100;
count     = 0;
crit      = 1e-12;
criterion = true;
quantile  = zeros(size(prob));
while criterion
    count  = count + 1;
    CdfPdf = tcdfpdf(quantile,t,weightedChf,0,normalize);
    correction  = (CdfPdf(:,1) - prob) ./ CdfPdf(:,2);
    quantile = quantile - correction;
    criterion = any(abs(correction) > crit * abs(quantile)) ...
        && max(abs(correction)) > crit && count < maxiter; 
end

if count == maxiter
    warning('TINVVW did not converge');
end

q(ok) = quantile;
yfun = q;
xfun = pr;

%% FUNCTION  GWEIGHTS
function [nodes,weights] = gweights(limits,subInts)
%GWEIGHTS Gaussian quadrature nodes and weights for all
%         subintervals defined by their limits and by the 
%         number of their subintervals (subInts)

GQ14rule = [ 0.006858095651593830579201, 0.017559730165875931516;
             0.035782558168213241331804, 0.040079043579880104903;
             0.086399342465117503405103, 0.060759285343951592345;
             0.156353547594157264925990, 0.078601583579096767285;
             0.242375681820922954017355, 0.092769198738968906871;
             0.340443815536055119782164, 0.102599231860647801983;
             0.445972525646328168966877, 0.107631926731578895098;
             0.554027474353671831033122, 0.107631926731578895098;
             0.659556184463944880217836, 0.102599231860647801983;
             0.757624318179077045982645, 0.092769198738968906871;
             0.843646452405842735074010, 0.078601583579096767285;
             0.913600657534882496594897, 0.060759285343951592340;
             0.964217441831786758668196, 0.040079043579880104900;
             0.993141904348406169420799, 0.017559730165875931516];

nquad   = 14;
on      = ones(nquad,1);
Nints   = sum(subInts);
nodes       = zeros(14*Nints,1);
weights       = zeros(14*Nints,1);
Nlimits = length(limits);
startid = 0; 
for i  = 1:Nlimits-1
    d  = (limits(i+1) - limits(i)) / subInts(i);
    om = ones(1,subInts(i));
    shift = 0:(subInts(i)-1);
    ti  = d * (GQ14rule(:,1) * om + on * shift) + limits(i);
    wi  = d * GQ14rule(:,2) * om;
    id  = startid + (1:14*subInts(i));
    nodes(id) = ti(:);
    weights(id) = wi(:);
    startid = startid + 14 * subInts(i);
end
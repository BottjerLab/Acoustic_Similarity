function hs = plotDerivGram(spec, params, varargin)

if nargin < 2 || isempty(params)
    params = defaultParams;
end
params = processArgs(params, varargin{:});
% arrange into log space
deps = params.dgram.minContrast;
lDeriv = spec.deriv;

lDeriv(abs(spec.deriv) < deps) = 0;

lDeriv(spec.deriv > deps) = -log(lDeriv(spec.deriv > deps)/deps);
lDeriv(spec.deriv < -deps) = log(-lDeriv(spec.deriv < -deps)/deps);

% convert to direct color mapping, for efficiency's sake
nColors = 256;
lmin = min(lDeriv(:)); lmax = max(lDeriv(:));
%lDerivMap = flipud(uint8(fix((lDeriv - lmin)/(lmax-lmin) * nColors) + 1));

hs = imagesc([min(spec.times) max(spec.times)],...
             [min(spec.freqs) max(spec.freqs)],... 
             lDeriv);
colormap(gray(nColors));

% flip the y-axis
set(gca,'YDir','normal');
% a much slower way to do it, but imagesc may not be supported in all matlab
% versions
%hs = surf(spec.times, spec.freqs, lDerivMap,'EdgeColor','none',...
%    'CDataMapping','direct','facecolor','texturemap');%,...
%     'xlimmode','manual','ylimmode','manual','zlimmode','manual',...
%     'climmode','manual','alimmode','manual');
%view(0,90)

end

function ans = roundN(val,N)
    round(val / N) * N;
end
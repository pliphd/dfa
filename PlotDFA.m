function [haxes, varargout] = PlotDFA(n, Fn, varargin)
%PLOTDFA plot fluctuation function Fn as a function of scale n
% 
% Inputs:  n           - scale, a column vector, point
%          Fn          - fluctuation function, a column vector of the same
%                        length of n, or a matrix that has the same row
%                        number to n. Each column will be treated as one
%                        fluctuation function.
%          haxes       - optional
%          SpecStruc   - plot optioins, containing the following fields
%                Drift   * y shift for each Fn curve
%                Region  * fit region, a vector or matrix with two columns,
%                          the first column indicates the start scale and
%                          the end scale of fitting, point
%                Alpha   * scaling exponent, a matrix, size(Alpha, 1) ==
%                          size(Region, 1), size(Alpha, 2) == size(Fn, 2)
%                Interse * intersection of each fitting line, has the same
%                          size as Alpha
%                XLim    * XLim of axes
%                YLim    * YLim of axes
%                Epoch   * epoch length
%          Type        - F(n) or F(n)/n
% 
% Outputs: haxes       - the handle of axes
%          alpha       - output alpha in case plot F(n)/n
% 
% $Author:  Peng Li, Ph.D.
%           MBP, Div Sleep Med, BWH & HMS
% $Date:    Apr 20, 2016
% $Modif.:  Apr 21, 2016
%               Add options to plot F(n)/n
%           May 27, 2016
%               Change scale and range's unit to point (corresponding to
%               DFA fitting results)
%           Jun 10, 2016
%               Add Epoch info to SpecStruc structure
%           Aug 16, 2017
%               stop extend x axis when plotting fit line
%           Jul 08, 2021
%               no longer force xtick to be hour since the grid becomes
%               difficult to read in log scale
% 

% parse inputs
narginchk(2, 7);

% set properties
ColorSpec = load('ColorSpec.txt');
if size(ColorSpec, 1) < size(Fn, 2)
    disp('Number of curves were beyond the number of predefined colors. Using default color (black) for all remains.');
    ColorSpec = [ColoSpec; zeros(size(Fn, 2) - size(ColorSpec, 1), 3)];
end
% following properties to be changed to another way to define
MarkerSpec = {'o', 's', 'd', '^'};
FillSpec   = [0 0 0 0 0 0 0 0 0];
% plot options
if nargin >= 3
    haxes = varargin{1};
    SpecStruc = varargin{2};
    alpha     = zeros(size(SpecStruc.Alpha));
    interse   = zeros(size(SpecStruc.Interse));
    if nargin >= 4
        Type  = varargin{3};
    else
        Type  = 'Fn';
    end
    if nargin >= 5
        saveit = varargin{4};
    else
        saveit = 0;
    end
    if saveit
        outdir = varargin{5};
    else
        outdir = '';
    end
else
    SpecStruc.Drift = zeros(size(Fn, 2), 1);
    Type = 'Fn';
end

% plot
if isempty(haxes)
    figWidth = 5; figHeight = 4;
    screenUn = get(0, 'Units');
    set(0, 'Units', 'inches');
    screensz = get(0, 'ScreenSize');
    figLeft  = (screensz(3) - figWidth)/2;
    figBot   = (screensz(4) - figHeight)/2;
    set(0, 'Units', screenUn);
    
    figure('Unit', 'inches', 'Position', [figLeft figBot figWidth figHeight]);
    haxes = axes('Parent', gcf, 'Units', 'normalized', 'Position', [.12 .14 .85 .82]);
end
haxes.NextPlot = 'add';

for iC = 1:size(Fn, 2)
    switch Type
        case 'Fn'
            hc = plot(haxes, n(1:2:end), Fn(1:2:end, iC) / SpecStruc.Drift(iC), 'Color', ColorSpec(iC, :) ./ 255, ...
                'Marker', MarkerSpec{iC}, 'MarkerEdgeColor', ColorSpec(iC, :) ./ 255, 'MarkerSize', 6, ...
                'LineStyle', 'none');
        case 'FnRatio'
            hc = plot(haxes, n(1:2:end), Fn(1:2:end, iC) ./ n(1:2:end) / SpecStruc.Drift(iC), 'Color', ColorSpec(iC, :) ./ 255, ...
                'Marker', MarkerSpec{iC}, 'MarkerEdgeColor', ColorSpec(iC, :) ./ 255, 'MarkerSize', 6, ...
                'LineStyle', 'none');
    end
    if FillSpec(iC)
        set(hc, 'MarkerFaceColor', ColorSpec(iC, :) ./ 255);
    else
        set(hc, 'MarkerFaceColor', 'w');
    end
    
    if nargin >= 3
        for iF   = 1:size(SpecStruc.Region, 1)
            fitx = [log(SpecStruc.Region(iF, 1)) - .0, log(SpecStruc.Region(iF, 2)) + .0];
            
            % specify start point and end point for plotting alpha line
            stx  = log(exp(fitx(1))) + diff(log(exp(fitx)))/5;
            edx  = log(exp(fitx(2))) - diff(log(exp(fitx)))/5;
            
            switch Type
                case 'Fn'
                    plot(haxes, exp(fitx), exp(SpecStruc.Alpha(iF, iC).*fitx + SpecStruc.Interse(iF, iC)) / SpecStruc.Drift(iC), ...
                        'Color', ColorSpec(iC, :) ./ 255, 'LineWidth', 2);
                    
                    sty = exp(SpecStruc.Alpha(iF, iC).*stx + SpecStruc.Interse(iF, iC)) / SpecStruc.Drift(iC);
                    edy = exp(SpecStruc.Alpha(iF, iC).*edx + SpecStruc.Interse(iF, iC)) / SpecStruc.Drift(iC);
                    
                    if mod(iC, 2) == 1
                        plot(haxes, [exp(stx) exp(edx)], [sty sty], 'Color', ColorSpec(iC, :) ./ 255);
                        plot(haxes, [exp(edx) exp(edx)], [sty edy], 'Color', ColorSpec(iC, :) ./ 255);
                        text(haxes, exp(edx), exp(log(sty) - log(1.8)), sprintf('%.2f', SpecStruc.Alpha(iF, iC)), 'FontSize', 12, 'Color', ColorSpec(iC, :) ./ 255);
                    else
                        plot(haxes, [exp(stx) exp(stx)], [sty edy], 'Color', ColorSpec(iC, :) ./ 255);
                        plot(haxes, [exp(stx) exp(edx)], [edy edy], 'Color', ColorSpec(iC, :) ./ 255);
                        text(haxes, exp(log(exp(stx)) - .5), exp(log(edy) + log(1.8)), sprintf('%.2f', SpecStruc.Alpha(iF, iC)), 'FontSize', 12, 'Color', ColorSpec(iC, :) ./ 255);
                    end
                    
                case 'FnRatio'
                    plot(haxes, exp(fitx), exp(SpecStruc.Alpha(iF, iC).*fitx + SpecStruc.Interse(iF, iC)) ./ exp(fitx) / SpecStruc.Drift(iC), ...
                        'Color', ColorSpec(iC, :) ./ 255, 'LineWidth', 2);
                    
                    alpha(iF, iC)   = diff(log(exp(SpecStruc.Alpha(iF, iC).*fitx + SpecStruc.Interse(iF, iC)) ./ exp(fitx))) ./ diff(fitx);
                    interse(iF, iC) = log(exp(SpecStruc.Alpha(iF, iC).*fitx(1) + SpecStruc.Interse(iF, iC)) ./ exp(fitx(1))) - alpha(iF, iC)*fitx(1);
                    
                    sty = exp(alpha(iF, iC).*stx + interse(iF, iC)) / SpecStruc.Drift(iC);
                    edy = exp(alpha(iF, iC).*edx + interse(iF, iC)) / SpecStruc.Drift(iC);
                    
                    if mod(iC, 2) == 1
                        if alpha(iF, iC) > 0
                            plot(haxes, [exp(stx) exp(edx)], [sty sty], 'Color', ColorSpec(iC, :) ./ 255);
                            plot(haxes, [exp(edx) exp(edx)], [sty edy], 'Color', ColorSpec(iC, :) ./ 255);
                            text(haxes, exp(edx), exp(log(sty) - log(1.2)), sprintf('%.2f', alpha(iF, iC)), 'FontSize', 12, 'Color', ColorSpec(iC, :) ./ 255);
                        else
                            plot(haxes, [exp(stx) exp(edx)], [edy edy], 'Color', ColorSpec(iC, :) ./ 255);
                            plot(haxes, [exp(stx) exp(stx)], [sty edy], 'Color', ColorSpec(iC, :) ./ 255);
                            text(haxes, exp(log(exp(stx)) - .5), exp(log(edy) - log(1.2)), sprintf('%.2f', alpha(iF, iC)), 'FontSize', 12, 'Color', ColorSpec(iC, :) ./ 255);
                        end
                    else
                        if alpha(iF, iC) > 0
                            plot(haxes, [exp(stx) exp(stx)], [sty edy], 'Color', ColorSpec(iC, :) ./ 255);
                            plot(haxes, [exp(stx) exp(edx)], [edy edy], 'Color', ColorSpec(iC, :) ./ 255);
                            text(haxes, exp(log(exp(stx)) - .5), exp(log(edy) + log(1.2)), sprintf('%.2f', alpha(iF, iC)), 'FontSize', 12, 'Color', ColorSpec(iC, :) ./ 255);
                        else
                            plot(haxes, [exp(edx) exp(edx)], [sty edy], 'Color', ColorSpec(iC, :) ./ 255);
                            plot(haxes, [exp(stx) exp(edx)], [sty sty], 'Color', ColorSpec(iC, :) ./ 255);
                            text(haxes, exp(edx), exp(log(edy) + log(1.2)), sprintf('%.2f', alpha(iF, iC)), 'FontSize', 12, 'Color', ColorSpec(iC, :) ./ 255);
                        end
                    end
            end
        end
    end
end

set(haxes, 'XScale', 'log', 'YScale', 'log', 'Box', 'off');
%    'XTick', [.1 1 10 100] * 60 * 60 / SpecStruc.Epoch, 'XTickLabel', {'10^{-1}', '10^0', '10^1', '10^2'});

if nargin >= 3
    if isfield(SpecStruc, 'XLim')
        haxes.XLim = SpecStruc.XLim;
    end
    if isfield(SpecStruc, 'YLim')
        haxes.YXLim = SpecStruc.YLim;
    end
    
    switch Type
        case 'Fn'
            varargout{1} = SpecStruc.Alpha;
        case 'FnRatio'
            varargout{1} = alpha;
    end
end

xlabel(haxes, 'n');
switch Type
    case 'Fn'
        ylabel(haxes, 'F(n)');
    case 'FnRatio'
        ylabel(haxes, 'F(n)/n');
end

if saveit
    print(haxes.Parent, outdir, '-djpeg', '-r0');
end
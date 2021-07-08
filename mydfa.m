%MYDFA class for detrended fluctuation analysis
%
% $Author:  Peng Li
% $Date:    Arp 1, 2020
%

classdef mydfa < handle
    properties
        pts           % a physiological time series object, in unit: sec
        
        order = 2
        minNumMotif = 4
        
        windowLength  % nan-empty for sliding window analysis, in unit: h
        fitRegion     % M by 2 matrix, each row corresponding to one region, in unit: min
    end
    
    properties (SetAccess = protected)
        fluctuationFunction
    end
    
    properties (SetAccess = protected, Hidden = true)
        timescale     % in unit: points
        numberMotifs
    end
    
    properties (SetAccess = protected)
        fitResult
    end
    
    properties (SetAccess = protected, Hidden = true)
        tsSplitIntegrate
        tsReconstruct
    end
    
    % construct
    methods
        function this = mydfa
        end
    end
    
    % public, dfa, fit, and plot
    methods
        function dfa(this)
            if isempty(this.windowLength)
                this.doDfa;
            else
                [windows, awindows] = ClipWindow(this.pts.TrueTime(1), this.pts.TrueTime(end), this.pts.UserData.Epoch/60, this.windowLength);
                PTS = this.pts; % back up as it will be overwritten
                RES = [];
                for iW = 1:size(windows, 1)
                    if windows(iW, 2) - windows(iW, 1) < .9 * this.windowLength * 3600 / this.pts.UserData.Epoch
                        continue;
                    end
                    this.pts = getinterestedsample(PTS, awindows(iW, 1), this.windowLength * 3600 / this.pts.UserData.Epoch);
                    this.doDfa;
                    
                    res = [table(string(PTS.Name) + "_" + datestr(awindows(iW, 1), 'yyyymmddHHMMss') + ".dfa" + num2str(this.order), ...
                        'VariableNames', {'fileName'}) this.fluctuationFunction(:, 2:end)];
                    RES = [RES; res];
                end
                this.pts = PTS;
                this.fluctuationFunction = RES;
            end
        end
        
        function fit(this)
            if isempty(this.fitRegion)
               return;
            end
            
            res = splitapply(@(r) {fitThisRegion(this.fluctuationFunction, r, this.pts.UserData.Epoch)}, ...
                this.fitRegion, (1:size(this.fitRegion, 1))');
            
            this.fitResult = [];
            for iR = 1:length(res)
                curRes  = res{iR};
                if isempty(curRes)
                    continue;
                end
                
                fileTbl = table(this.fluctuationFunction.fileName(curRes(:, 1)), 'VariableNames', {'filename'});
                resTbl  = array2table(curRes(:, 2:end), 'VariableNames', {'minscale', 'maxscale', 'alpha', 'goodness', 'intercept'});
                
                this.fitResult = [this.fitResult; [fileTbl resTbl]];
            end
        end
        
        function save(this, varargin)
            % parse inputs
            ni = nargin - 1;
            
            if ni > 0
                % parse PV pairs
                pvStart = 1;
                
                if mod(ni - pvStart + 1, 2) ~= 0
                    error('Missing property names or values.');
                end
                for iP = pvStart:2:ni
                    value = varargin{iP+1};
                    switch lower(char(varargin{iP}))
                        case 'option'
                            saveOption = value;
                        case 'outdir'
                            outdir = value;
                            if ~(exist(outdir, 'dir') == 7)
                                mkdir(outdir);
                            end
                    end
                end
            end
            
            % save option
            switch lower(char(saveOption))
                case {'fluctuationfunction', 'fn', 'fluctuation'}
                    saveFluc(this.fluctuationFunction, outdir);
                case {'result', 'fit', 'results'}
                    saveFit(this.fitResult, outdir);
                case 'all'
                    saveFluc(this.fluctuationFunction, outdir);
                    saveFit(this.fitResult, outdir);
            end
        end
        
        function ha = plot(this, varargin)
            % to be swiched, using old function PlotDFA temporarily
            SpecStruc.Epoch = this.pts.UserData.Epoch;
            plotType = 'Fn';
            outdir   = '';
            saveit   = 0;
            
            % parse inputs
            ni = nargin - 1;
            
            if ni > 0
                % parse PV pairs
                if string(class(varargin{1})) == "dfaPlot"
                    hf      = varargin{1};
                    pvStart = 2;
                else
                    hf      = dfaPlot;
                    pvStart = 1;
                end
                
                if mod(ni - pvStart + 1, 2) ~= 0
                    error('Missing property names or values.');
                end
                for iP = pvStart:2:ni
                    value = varargin{iP+1};
                    switch lower(char(varargin{iP}))
                        case 'drift'
                            SpecStruc.Drift = value;
                        case 'xlim'
                            SpecStruc.XLim  = value;
                        case 'ylim'
                            SpecStruc.YLim  = value;
                        case 'type'
                            plotType  = value;
                        case 'outdir'
                            outdir = value;
                            if ~(exist(outdir, 'dir') == 7)
                                mkdir(outdir);
                            end
                            saveit = 1;
                    end
                end
            end
            
            plotName = this.fluctuationFunction.fileName;
            for iP = 1:length(plotName)
                if saveit
                    savedir = fullfile(outdir, plotName(iP)+".jpg");
                else
                    savedir = outdir;
                end
                
                cla(hf.actiAxis);
                cla(hf.shallowAxis);
                cla(hf.dfaAxis);
                
                % show data, copied from ezActi.dataShow, possibly to make it static
                y = this.fluctuationFunction.timeSeries{iP};
                x = 1:length(y);
                
                nameparts = strsplit(this.fluctuationFunction.fileName(iP), {'_', '.'});
                startTime = datetime(nameparts(end-1), 'InputFormat', 'yyyyMMddHHmmss');
                
                yrange = [min(y) max(y)];
                yrange = [yrange(1) - .1*diff(yrange), yrange(2) + .1*diff(yrange)];
                
                % signal end bar
                plot(hf.actiAxis, [x(end), x(end)], [-1e6 1e6], 'LineWidth', 3, 'Color', 'c');
                set(hf.actiAxis, 'NextPlot', 'add');
                
                plot(hf.actiAxis, x, y, 'Color', 'k');
                
                % plot day seperator
                PointPerDay = 24*60*60 / SpecStruc.Epoch;
                TotalDays   = ceil(length(x) / PointPerDay);
                StartBorder = (1:TotalDays-1) .* PointPerDay + 1;
                Seperator   = repmat([-1e6; 1e6], 1, length(StartBorder));
                plot(hf.actiAxis, [StartBorder; StartBorder], Seperator, 'Color', [.15 .15 .65], 'LineStyle', '--', 'LineWidth', 2);
                xlabel(hf.actiAxis, 'Points');
                
                % range set
                try
                    hf.actiAxis.YLim    = yrange;
                end
                hf.actiAxis.XLim    = [0 length(x)];
                hf.shallowAxis.XLim = [0 length(x)];
                
                % display time
                hf.actiAxis.Box        = 'off';
                hf.shallowAxis.Visible = 'on';
                
                peerXTick = hf.actiAxis.XTick;
                xTick     = startTime + seconds(SpecStruc.Epoch).*peerXTick - x(1).*seconds(SpecStruc.Epoch);
                hf.shallowAxis.XTick = peerXTick;
                hf.shallowAxis.XTickLabel = datestr(xTick, 'dd-mmm-yyyy HH:MM');
                hf.shallowAxis.XTickLabelRotation = 10;
                
                % show toolbar
                hf.shallowAxis.Toolbar.Visible = 'off';
                hf.actiAxis.Toolbar.Visible    = 'on';
                
                % add quality info
                gap = detConstantOne(~this.fluctuationFunction.isNotGap{iP});
                for iG = 1:size(gap, 1)
                    item = gap(iG, :);
                    patch(hf.actiAxis, 'XData', ...
                        [item(1) item(1) item(2) item(2)], ...
                        'YData', [yrange(1) yrange(2) yrange(2) yrange(1)], ...
                        'FaceColor', 'm', 'EdgeColor', 'm', 'FaceAlpha', 0.2);
                    
                    % switch layers
                    child = hf.actiAxis.Children;
                    child = [child(end); child(1:end-1)];
                    hf.actiAxis.Children = child;
                end
                
                hf.actiAxis.Layer = 'top';
                hf.shallowAxis.Layer = 'top';
                
                % show dfa
                if isempty(this.fitResult)
                    drawnow;
                    continue;
                end
                
                currentResult = this.fitResult(this.fitResult.filename == plotName(iP), :);
                if isempty(currentResult)
                    drawnow;
                    continue;
                end
                
                SpecStruc.Region  = [currentResult.minscale currentResult.maxscale] * 60 / SpecStruc.Epoch;
                SpecStruc.Alpha   = currentResult.alpha;
                SpecStruc.Interse = currentResult.intercept;
                
                % show DFA
                PlotDFA(this.fluctuationFunction.timeScale{iP}, this.fluctuationFunction.fluctuationFunction{iP}, hf.dfaAxis, SpecStruc, plotType, saveit, savedir);
                drawnow;
                
                % link dfa and dfaShallow
                linkaxes([hf.dfaShallowAxis, hf.dfaAxis], 'x');
                
                % copy obj to shallow since log scale sometimes make the
                % XLim different visually (0 is way left but if with obj, 0
                % is not shown
                copyobj(allchild(hf.dfaAxis), hf.dfaShallowAxis);
                
                % color none
                allc = allchild(hf.dfaShallowAxis);
                set(allc, 'Color', 'none');
                set(allc(strcmpi(get(allc, 'Type'), 'line')), 'Marker', 'none');
                
                hf.dfaShallowAxis.XScale = hf.dfaAxis.XScale;
                hf.dfaShallowAxis.XTick  = hf.dfaAxis.XTick;
                hf.dfaShallowAxis.XMinorTick = hf.dfaAxis.XMinorTick;
                hf.dfaShallowAxis.XTickLabel = string(hf.dfaAxis.XTick * SpecStruc.Epoch / 60);
                hf.dfaShallowAxis.XLim   = hf.dfaAxis.XLim;
                
                % labels
                xlabel(hf.dfaAxis, 'time scale (points)');
                xlabel(hf.dfaShallowAxis, 'time scale (min)');
            end
        end
    end
    
    % support functions
    methods (Access = protected)
        function doDfa(this)
            if all(this.pts.Quality == 0 | isnan(this.pts.Data))
                this.fluctuationFunction = table(string(this.pts.Name) + "_" + datestr(this.pts.TrueTime(1), 'yyyymmddHHMMss') + ".dfa" + num2str(this.order), ...
                    {this.pts.Data}, {this.pts.Quality}, {[]}, {[]}, {[]}, ...
                    'VariableNames', {'fileName', 'timeSeries', 'isNotGap', 'timeScale', 'fluctuationFunction', 'numberMotifs'});
                return;
            end
            
            this.splitIntegrate;
            this.calculateTimescale;
            this.reconstructSplit;
            this.fitReconstruct;
        end
    end
    
    methods (Access = protected)
        function splitIntegrate(this)
            indGap             = this.pts.Quality == 0 | isnan(this.pts.Data);
            splitLineIndicator = 1 + cumsum(indGap);
            origLineIndicator  = 1:this.pts.Length;
            
            tempSplit = accumarray( ...
                splitLineIndicator(~indGap), ...
                origLineIndicator(~indGap), ...
                [], ...
                @(x) {cumsum(this.pts.Data(x, :))});
            this.tsSplitIntegrate = tempSplit(~cellfun(@isempty, tempSplit));
        end
        
        function calculateTimescale(this)
            splitLength    = cellfun(@(x) size(x, 1), this.tsSplitIntegrate);
            logscale       = log(this.order+3):log(2)/10:log(max(splitLength));
            this.timescale = unique(floor(exp(logscale)))';
        end
        
        function fitReconstruct(this)
            fn       = cellfun(@(m) doFitReconst(m, this.order), this.tsReconstruct);
            numMotif = cellfun(@(m) size(m, 2), this.tsReconstruct);
            
            Fn     = fn(numMotif >= this.minNumMotif);
            nMotif = numMotif(numMotif >= this.minNumMotif);
            tscale           = this.timescale(numMotif >= this.minNumMotif);
            
            this.fluctuationFunction = table(string(this.pts.Name) + "_" + datestr(this.pts.TrueTime(1), 'yyyymmddHHMMss') + ".dfa" + num2str(this.order), ...
                {this.pts.Data}, {this.pts.Quality}, {tscale}, {Fn}, {nMotif}, ...
                'VariableNames', {'fileName', 'timeSeries', 'isNotGap', 'timeScale', 'fluctuationFunction', 'numberMotifs'});
        end
        
        function reconstructSplit(this)
            this.tsReconstruct = arrayfun(@(n) reconstructSplitAtScale(this, n), this.timescale, 'UniformOutput', 0);
        end
        
        function tsReconstructAtScale = reconstructSplitAtScale(this, currentScale)
            tsReconstCell = cellfun(@(x) doReconstruct(x, currentScale), ...
                this.tsSplitIntegrate, 'UniformOutput', 0); % suppose to be a cell with each element's row = currentScale
            tsReconstructAtScale = cell2mat(tsReconstCell');
        end
    end
end

function tsReconstruct = doReconstruct(ts, n)
%DORECONSTRUCT reconstruct continuous time-series TS at time scale N
%

if length(ts) < n
    tsReconstruct = [];
    return;
end

intLength     = floor(size(ts, 1) / n) * n;
tsReconstruct = reshape(ts(1:intLength), n, []);
end

function fn = doFitReconst(tsReconstruct, order)
x = zscore(1:size(tsReconstruct, 1))';

% The built-in polyfit automatically treat a matrix input as a column vec,
% thus need splitapply to apply it to each column, might be slower
% p = splitapply(@(sts) {polyfit(x, sts, order)}, tsReconstruct, 1:size(tsReconstruct, 2));

% use a customized polyfit instead which is quicker
p    = polyfitM(x, tsReconstruct, order);

yhat = splitapply(@(e) polyval(e, x), p, 1:size(p, 2));

fn   = sqrt(nanmean((tsReconstruct(:) - yhat(:))...
    .*(tsReconstruct(:) - yhat(:))));
end

function P = polyfitM(x, Y, n)
%POLYFITM do polyfit along each column of Y, using the same x vector
%
% Construct Vandermonde matrix:
x = x(:);
V = ones(length(x), n + 1);
for j = n:-1:1
    V(:, j) = V(:, j + 1) .* x;
end
% Solve least squares problem:
[Q, R] = qr(V, 0);
P      = R \ (transpose(Q) * Y);  % equivalent to (V \ Y)
end

function res = fitThisRegion(fn, r, epoch)
%FITTHISREGION fit fluctuation function
%
res = splitapply(@(n, f) {doFit(n{1}, f{1}, r, epoch)}, fn.timeScale, fn.fluctuationFunction, (1:height(fn))');
index    = 1:height(fn);
validInd = index(cellfun(@(x) ~isempty(x), res));
res = [validInd(:), cell2mat(res)];
end


function res = doFit(n, f, r, epoch)
%DOFIT fit f agains n using log-log scale in region r
%
if isempty(n)
    res = [];
    return;
end

minscale  = max(n(1),   r(1) * 60 / epoch);
maxscale  = min(n(end), r(2) * 60 / epoch);

if minscale >= maxscale
    res = [];
    return;
end

rangeind  = n >= minscale & n <= maxscale;
LogFn     = log(f(rangeind));
LogTimesc = log(n(rangeind));

p         = polyfit(LogTimesc, LogFn, 1);
alpha     = p(1);
intercept = p(2);
fnhat     = polyval(p, LogTimesc);
goodness  = 1 - ((fnhat - LogFn)'*(fnhat - LogFn)) / ((LogFn - mean(LogFn))'*(LogFn - mean(LogFn)));

res = [minscale*epoch/60, maxscale*epoch/60, alpha, goodness, intercept];
end

function saveFluc(fn, outdir)
fn.fileName = fullfile(outdir, fn.fileName);
rowfun(@doSaveFluc, fn, 'InputVariables', {'fileName', 'timeScale', 'fluctuationFunction', 'numberMotifs'}, 'NumOutputs', 0);
end

function doSaveFluc(fileName, timeScale, fluctuationFunction, numberMotifs)
if ~isempty(timeScale{1})
    fid = fopen(fileName, 'w');
    fprintf(fid, '%d,%f,%d\r\n', [timeScale{1}, fluctuationFunction{1}, numberMotifs{1}]');
    fclose(fid);
end
end

function saveFit(fitResult, outdir)
if isempty(fitResult)
    return;
end
filename = fullfile(outdir, "fitResults.csv");
if exist(filename, 'file') == 2
    fid = fopen(filename, 'a');
    for iL = 1:height(fitResult)
        fprintf(fid, '%s,%.2f,%.2f,%f,%f,%f\r\n', ...
            fitResult.filename(iL), ...
            fitResult.minscale(iL), ...
            fitResult.maxscale(iL), ...
            fitResult.alpha(iL), ...
            fitResult.goodness(iL), ...
            fitResult.intercept(iL));
    end
else
    fid = fopen(filename, 'w');
    fprintf(fid, '%s,%s,%s,%s,%s,%s\r\n', 'filename', 'minscale', 'maxscale', 'alpha', 'goodness', 'intercept');
    for iL = 1:height(fitResult)
        fprintf(fid, '%s,%.2f,%.2f,%f,%f,%f\r\n', ...
            fitResult.filename(iL), ...
            fitResult.minscale(iL), ...
            fitResult.maxscale(iL), ...
            fitResult.alpha(iL), ...
            fitResult.goodness(iL), ...
            fitResult.intercept(iL));
    end
end
fclose(fid);
end
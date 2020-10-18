function handleOut = stackplot(Waveform, Range, div, gridopt, iTitle, iXlabel)
%STACKPLOT plot signal with stacking method.
%Syntax
%  HANDLEOUT = STACKPLOT(WAVEFORM, RANGE, DIV, GRIDOPT, TITLE, XLABEL)
%    Plot signal in stack
%    WAVEFORM 是嵌套的2层cell: {{time, data, Unit}, ...}
%    其中Unit是可选的
%    RANGE is the time axis range, default is using MATLAB's default value.
%    DIV: div value for stacking axes, default =0.02. if div=0, the
%         adjacent axes will be connected.
%    handleOut  :   return the handle object of each figure. 
%    title: the title for the figure
%    xlabel: the xlabel in the last axe
%
% Example:
%
%  stackplot({{t1, data1,'Volt (V)'},{t2, data2, 'Is (A)'}}, [], [], [], 'I-V plot', 'Time (ms)');
%
% See also:
% PLOT, SUBPLOT, AXES, XLIM, YLIM, SUBLABEL              

if nargin < 1
    display(['>> help ' mfilename])
    eval(['help ' mfilename]);
    error('Too less input parameters!');
end

Range_default = [];  % default plot range
div_default = 0.01; % divide default value
gridopt_default = 'gridoff'; % grid option default
iTitle_default = [];  % default title

iXlabel_default = 'Times(s)';  % xlabel

if nargin < 2
    Range = Range_default;
elseif isempty(Range)
    Range = Range_default;
end

if nargin < 3
    div = div_default;
elseif isempty(div)
    div = div_default;
end

if nargin < 4
    gridopt = gridopt_default;
elseif isempty(gridopt)
    gridopt = gridopt_default;
end

if nargin < 5
    iTitle = iTitle_default;
elseif isempty(iTitle)
    iTitle = iTitle_default;
end % default

if nargin < 6
    iXlabel = iXlabel_default;
end

if ~iscell(Waveform)
    error('Input waveform must be cell.')
end

wavenum = length(Waveform);
stackplottag = 'stackplot';

if wavenum ~= 1
    clf
end

XLim = Range;
for i = 1:wavenum
    position = [0 0 0 0];
    position(1) = 0.2; % position is defined as [left bottom width height].
    position(2) = 0.95 - i*0.82/wavenum;
    position(3) = 0.7;
    position(4) = 0.82/wavenum - div;
    if wavenum ~= 1
        subplot('Position', position);
    end
    temp = Waveform{i};
    
    if ischar(temp{end})
        hg{i} = plot(temp{1:end-1}, 'Tag', stackplottag,'linewidth',2);
    else
        hg{i} = plot(temp{1:end}, 'Tag', stackplottag,'linewidth',2);
    end
    
    allaxes(i) = gca;
    
    set(gca, 'FontWeight', 'normal', 'FontSize', 16, 'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on','ticklength',[0.025 0.01],'Xgrid','off')
    if strcmpi(gridopt,'grid')
        grid on;
    end
    
    %     box off;
    
    if ~isempty(XLim)
        xlim(XLim)
    end
    if length(temp) > 2
        if ~ischar(temp{end})
            ylab = ['Sig' num2str(i)];
        else
            if strncmpi(temp{end}, 'dx', 2) || strncmpi(temp{end}, 'dy', 2)
                ylim([-1 1]);
            end
            ylab = temp{end};
        end
    else
        ylab = ['Sig' num2str(i)];
    end
    set(gca,'FontName','Times New Roman')
    ax = gca;
    h= get(ax, 'ylabel');
    set(h, 'FontAngle',  'normal', ...
           'FontName',   'Times New Roman', ...
           'FontUnits',  'points',...
            'FontSize',  16, ...
           'FontWeight', 'normal', ...
           'string', ylab,'Interpreter','LaTex');
    %ylabel(ylab);
    
    if i ~= wavenum
        set(gca, 'XTickLabel',[])
    end
end


if ~isempty(iTitle)
    temp = hg{1};
    title(get(temp(1), 'parent'), iTitle)
end


ax = gca;
h = get(ax, 'xlabel');
      set(h, 'FontAngle',  'normal', ...
           'FontName',   'Times New Roman', ...
           'FontUnits',  'points',...
            'FontSize',  16, ...
           'FontWeight', 'normal', ...
           'string', iXlabel,'Interpreter','LaTex');
           


% link all x axes
nn = 1;
for ii = 1:length(hg)
    for jj = length(hg{ii})
        all_axes(nn) = get(hg{ii}(jj), 'parent');
        nn = nn + 1;
    end
end
if length(all_axes)>1
    linkaxes(all_axes, 'x');    % link all x-axes
end

% link all x-axes' scale
hlinks = linkprop(all_axes, 'xscale');
KEY = 'graphics_linkaxes_scale';
for ii = 1:length(all_axes)
    setappdata(all_axes(ii), KEY, hlinks);
end

for i = 1:wavenum
    % place subfigure label, such as (a), (b) ...
    if wavenum ~= 1
        labeltext = ['('  char('a' + i - 1)  ') '];
        ylimtemp = ylim(allaxes(i));
        xlimtemp = xlim(allaxes(i));
        th = text(xlimtemp(2)-0.08*(xlimtemp(2)-xlimtemp(1)), ylimtemp(2)-0.08*(ylimtemp(2)-ylimtemp(1)), labeltext, 'Color', 'k', 'BackgroundColor', 'w', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
            'FontWeight', 'normal', 'FontSize', 15, ...
            'handlevisibility','off', ...
            'parent', allaxes(i), ...
            'Tag', 'sublabel');
%         uistack(th, 'top');
    end
end   


% output handles
if nargout >= 1
    handleOut = hg;
end
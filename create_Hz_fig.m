function create_Hz_fig(cdata1)
%CREATEFIGURE(CDATA1)
%  CDATA1:  image cdata

%  Auto-generated by MATLAB on 03-May-2017 15:46:35

% Create figure
figure;

% Create axes
axes1 = axes;
hold(axes1,'on');

% Create image
image(cdata1,'CDataMapping','scaled');

% Create xlabel
xlabel('12 Frequencies (Hz) each with 6 dB Levels');

% Create ylabel
ylabel('12 Frequencies (Hz) each with 6 dB Levels');

% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0.5 72.5]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0.5 72.5]);
box(axes1,'on');
axis(axes1,'ij');
% Set the remaining axes properties
set(axes1,'FontSize',24,'Layer','top','XTick',...
    [0.5 3 6.5 9 12.5 15 18.5 21.5 24.5 27.5 30.5 33 36.5 39.5 42.5 45.5 48.5 51.5 54.5 57.5 60.5 64 66.5 70],...
    'XTickLabel',...
    {'','250','','500','','630','','800','','1000','','1250','','1600','','2000','','2500','','3150','','5040','','8000'},...
    'YTick',...
    [0.5 3 6.5 9.5 12.5 15 18.5 21 24.5 27.5 30.5 34 36.5 39.5 42.5 45.5 48.5 51.5 54.5 57.5 60.5 63.5 66.5 70],...
    'YTickLabel',...
    {'','250','','500','','630','','800','','1000','','1250','','1600','','2000','','2500','','3150','','5040','','8000'});
% Create colorbar
colorbar('peer',axes1,'FontSize',21.6);


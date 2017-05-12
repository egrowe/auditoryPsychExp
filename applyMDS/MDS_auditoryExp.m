%% Applying MDS to the auditory dissimilarity data
clear all

filename = ['INSERT FILENAME HERE'] %insert filename for results
load(filename, 'resultsMatrix')

%Checks whether frequency variable exists
if exist('nFreq', 'var')
    noFreq = nFreq
end

%Checks whether results are for dB or wave types
if exist('ndB', 'var')
    noOtherVariable = ndB
elseif exist('nWaves', 'var')
    noOtherVariable = nWaves
end

%% Setup the dissimilarity matrix from our experimental results
%Shift the values so '4' in our results similarity matrix now equals '0'
%(req. for MDS that most similar values equal zero)
addMatrix = ones(length(resultsMatrix),length(resultsMatrix))*4;
disimMatrix = (abs(resultsMatrix-addMatrix))/8

% MDS will not work without diagonal equally zero
for ii = 1:length(resultsMatrix)
    disimMatrix(ii,ii) = 0;
end

% Run MDS on the dissimilarity matrix data
[Y,eigvals] = cmdscale(disimMatrix)


%% PLOT RESULTS (different colours and sizes for data points)
% Setup the colours for each of the 12 Frequencies
C_250 = [0.9 0.7 0.9] %light pink
C_500 = [1 0 1] %pink
C_630 = [0.7 0 0.7] %purple
C_800 = [1 0 0] %red
C_1000 = [1 0.5 0] %orange
C_1250 = [1 1 0] %yellow
C_1600 = [0.7 0.7 0] %green/yellow
C_2000 = [0 1 0] %green
C_2500 = [0 0.7 0.7] %teal
C_3150 = [0 0 1] %blue
C_5040 = [0 0 0.5] %blue/black
C_8000 = [0 0 0] %black

colours = vertcat(C_250,C_500,C_630,C_800,C_1000,C_1250,C_1600, ...
    C_2000, C_2500, C_3150, C_5040, C_8000)

colours(noFreq+1:end,:) = [] %remove colours if less than 12 freqs used

colourMatrix = []
for ii = 1:length(colours)
    Hz_use = colours(ii,:)
    colourMatrix = vertcat(colourMatrix, repmat(Hz_use,noOtherVariable,1)) %change ndB to nWaves
end

% Set the size of the markers from small to large
sizes = [1000;3000;5000;7000;9000;11000]; %use code below is less than 6 other factors
sizes(noOtherVariable+1:end) = [] %this removes some sizes to suit nWaves

useSizes = repmat(sizes,noFreq,1)
%% Plot MDS results
h = scatter(Y(:,1),Y(:,2), useSizes, '.')
set(gca,'FontSize', 20)
c = h.CData; c = colourMatrix; h.CData = c;
axis(max(max(abs(Y))) * [-1.1,1.1,-1.1,1.1]); axis('square'); 
line([-1,1],[0 0],'XLimInclude','off','Color',[.7 .7 .7])
line([0 0],[-1,1],'YLimInclude','off','Color',[.7 .7 .7])
xlabel('Dimension 1','FontSize',22);
ylabel('Dimension 2','FontSize',20);

%text(Y(:,1),Y(:,2),labels,'HorizontalAlignment','left');
% DID HAVE LABELS BUT TOO MANY VARIABLE = MESSY - use if you like

%labels = {' 250',' 500',' 630',' 800',' 1000',' 1250',' 1600',' 2000','
%2500','3150', '5040', '8000'}; % uncomment this line if only Hz was tested
%(no dB changes) -- will also need to comment out the next lines (31 to 42)
% 
% labels = {'250 11', '250 21', '250 31', '250 41', '250 51', '250 61', ...
%     '500 11', '500 21', '500 31', '500 41', '500 51', '500 61', ...
%     '630 11', '630 21', '630 31', '630 41', '630 51', '630 61', ...
%     '800 11', '800 21', '800 31', '800 41', '800 51', '800 61', ...
%     '1000 11', '1000 21', '1000 31', '1000 41', '1000 51', '1000 61', ...
%     '1250 11', '1250 21', '1250 31', '1250 41', '1250 51', '1250 61', ...
%     '1600 11', '1600 21', '1600 31', '1600 41', '1600 51', '1600 61', ...
%     '2000 11', '2000 21', '2000 31', '2000 41', '2000 51', '2000 61', ...
%     '2500 11', '2500 21', '2500 31', '2500 41', '2500 51', '2500 61', ...
%     '3150 11', '3150 21', '3150 31', '3150 41', '3150 51', '3150 61', ...
%     '5040 11', '5040 21', '5040 31', '5040 41', '5040 51', '5040 61', ...
%     '8000 11', '8000 21', '8000 31', '8000 41', '8000 51', '8000 61'}

%% Plot the eigenvalues
plot(1:length(eigvals),eigvals,'bo-');
line([1,length(eigvals)],[0 0],'LineStyle',':','XLimInclude','off',...
     'Color',[.7 .7 .7])
axis([1,length(eigvals),min(eigvals),max(eigvals)*1.1]);
xlabel('Eigenvalue number');
ylabel('Eigenvalue');


%% Plot legend of colour scale
points = [1:nFreq]'
pointsB = ones(1,nFreq)
h = scatter(pointsB,points, 5000, '.')
c = h.CData; c = colours; h.CData = c;

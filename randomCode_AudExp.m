
%Combine results matrices from 2 separate sessions of trials
totalSize = length(resultsMatrix)*length(resultsMatrix)
combMatrix = resultsMatrix;

for RR = 1:length(resultsMatrix)
    
    for CC = 1:length(resultsMatrix)
        
        if resultsMatrix(RR,CC) == 0
            
            combMatrix(RR,CC) = resultsMatrix2000(RR,CC);
            
        else
            continue
        end
        
    end
end


%Create stimulus matrix (trials to use) after some have already been
%partially completed earlier

for RR = 1:length(resultsMatrix)
    
    for CC = 1:length(resultsMatrix)
        
        if resultsMatrix(RR,CC) < 0 || resultsMatrix(RR,CC) > 0
            
            stimMatrix(RR,CC) = 0;
        end
        
    end
end

%% CHANGE SCALE OF DATA FROM POSITIVE TO ZERO (incorrect rating on 8AFC)
scaledMatrix = ones(length(resultsMatrix),length(resultsMatrix))

for RR = 1:length(resultsMatrix)
    
    for CC = 1:length(resultsMatrix)
        
        if resultsMatrix(RR,CC) > 0
            scaledMatrix(RR,CC) = 0;
        else
            scaledMatrix(RR,CC) =  resultsMatrix(RR,CC)
        end
        
    end
end


%% Change png to jpg

letterArray = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', ...
    'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'Z'};

targetArray = {'X', 'Y'};


for ii = 1:length(targetArray)
    useLetter = targetArray{ii}
    S = imread([num2str(useLetter) '.png']);
    imwrite(S,[num2str(useLetter) '.jpg'])
end

%% NEW WAY TO INDEX - should WORK!!!

takeColfrom = [1, 13, 25, 37, 49, 61]
idx2 = [1, 13, 25, 37, 49, 61]

for ii = 1:11
    idx2 = horzcat(idx2,takeColfrom+ii)
end

matrix2 = repmat(idx2,72,1)

(ind2,ind2)

for nn = 1:72
    
    newMat = {idx2(nn),idx2(nn)}
    
end

%idx2 = [1,7,13,19,25,31,37,43,49,55,61,67]

for ii = 1:5
    idx2 = horzcat(idx2,takeColfrom+ii)
end

% Matrix idxs to extract values


%% MDS plot using different dependencies

% USING A DIFFERENT DEPENDENCY
[eigvals eigvals./max(abs(eigvals))]

plot(Y(:,1),Y(:,2),'bx');
axis(max(max(abs(Y))) * [-1.1,1.1,-1.1,1.1]); axis('square');
text(Y(:,1),Y(:,2),labels,'HorizontalAlignment','left');
line([0 0],[-1,1],'XLimInclude','off','Color',[.7 .7 .7])
line([-1,1],[0 0],'YLimInclude','off','Color',[.7 .7 .7])



%% NEXT ATTEMPT ---- WORKS PERFECTLY!!!!
takeColfrom = [1,7,13,19,25,31,37,43,49,55,61,67]
takeRowfrom = [1,7,13,19,25,31,37,43,49,55,61,67]
putColIn = [1, 13, 25, 37, 49, 61, 72]
putRowIn = putColIn;

%Has to extract 12 numbers then start a new row?
% DO THIS FOR THE FIRST dB ONLY
reorderedMatrix = zeros(72,72);


for nextSubset = 1:6
    
    takeRowfrom = [1,7,13,19,25,31,37,43,49,55,61,67]+(nextSubset-1)
    useRowAddition = putRowIn(nextSubset)-1
    
    for withindBPair = 1:6
        
        takeColfrom = [1,7,13,19,25,31,37,43,49,55,61,67]+(withindBPair-1)
        useColAddition = putColIn(withindBPair)-1
        
        for newMatRow = 1:12
            
            for newMatCol = 1:12
                
                resultsRow = takeRowfrom(newMatRow)
                
                resultsCol = takeColfrom(newMatCol)
                
                reorderedMatrix((newMatRow+useRowAddition),(newMatCol+useColAddition)) = resultsMatrix(resultsRow,resultsCol);
                
            end
        end
        
    end
    
    
end





%% NEXT ATTEMPT ---- works for ONE row
takeColfrom = [1,7,13,19,25,31,37,43,49,55,61,67]
takeRowfrom = [1,7,13,19,25,31,37,43,49,55,61,67]
putColIn = [1, 13, 25, 37, 49, 61, 72]

%Has to extract 12 numbers then start a new row?
% DO THIS FOR THE FIRST dB ONLY
reorderedMatrix = zeros(72,72)

for withindBPair = 1:6
    
    takeRowfrom = [1,7,13,19,25,31,37,43,49,55,61,67]+(withindBPair-1)
    useColAddition = putColIn(withindBPair)-1
    
    for newMatRow = 1:12
        
        for newMatCol = 1:12
            
            resultsRow = takeRowfrom(newMatRow)
            
            resultsCol = takeColfrom(newMatCol)
            
            reorderedMatrix(newMatRow,(newMatCol+useColAddition)) = resultsMatrix(resultsRow,resultsCol);
            
        end
    end
    
end


%% Apply FFT to differently shaped soundwaves
% % 
% % [myWavfile, Fs] = audioread('250Hz_1.wav')
% % 
% % figure(8)
% % plot(abs(fft((myWavfile))))
% % 
% % %that plots you a the freq-domain-coefficients without
% % %relation to samplefreq and so on. Now you have to take care
% % %to generate a frequency-axes according to your specific wavfile:
% % 
% % %say N is the choosen length for fft (i.E. 1024). If you want
% % %the spec of the entire wavfile,
% % N=length(audioread('250Hz_1.wav'))
% % 
% % f = Fs/N .* (0:N-1);
% % 
% % % transformation -> returns complex coefficients
% % Y = fft(myWavfile, N);
% % % compute abs for plot and normalize to N
% % Y = abs( Y(1:N) ) ./ (N/2);
% % Y(Y==0) = eps; % avoiding log(0) --> >> help eps
% % Y = 20 * log10(Y); % turn spec logarithmic
% % 
% % figure(2)
% % plot(f,Y);
% % grid on;
% % 
% % %you can limit the spec to nyquist-frequency by calculating
% % %up to N/2 points
% % 
% % 
% % 
% % % %% TRY AGAIN
% 
% [y,Fs] = audioread('500Hz_1.wav');
% sound(y, Fs);
% F = fftshift(abs(fft(y)));
% f = linspace(-Fs/2, Fs/2, numel(y)+1);
% f(end) = [];
% figure(3); plot(f, F);
% 
% %% AND AGAIN
% 
% [data, Fs] = audioread('500Hz_1.wav');
% freqdata = fft(data);
% plot(freqdata)

%% FFT working

freqs = [250.00,281.25,312.50,343.75,375.00,406.25,437.50,468.75,500.00]
waveType = [1,2,3,4];
waveNames = {'Sine', 'Square', 'Sawtooth', 'Triangle'}
figCount = 1;

%Fourier Transform of Sound File
for chooseFreq = 1:length(freqs)
    
    for chooseType = 1:length(waveType)
        %for ii = 1:length(
        %Load File
        
        currFreq = freqs(chooseFreq)
        currType = waveType(chooseType)
        typeName = waveNames{chooseType}
        
        titlename = [num2str(currFreq) ' Hz ' num2str(typeName) ' Wave FFT']
        filename = [num2str(currFreq) 'Hz_Type_' num2str(currType) '_FFT.jpg']

        
        %[y,Fs] = audioread('281.25Hz_4.wav');
        
        [y,Fs] = audioread(sprintf('%gHz_%d.wav',currFreq,currType))
        
        Nsamps = length(y);
        t = (1/Fs)*(1:Nsamps)          %Prepare time data for plot
        
        %Do Fourier Transform
        y_fft = abs(fft(y));            %Retain Magnitude
        y_fft = y_fft(1:Nsamps/2);      %Discard Half of Points
        f = Fs*(0:Nsamps/2-1)/Nsamps;   %Prepare freq data for plot
        
        %Plot Sound File in Time Domain
        % figure
        % plot(t, y)
        % xlabel('Time (s)')
        % ylabel('Amplitude')
        % title('Tuning Fork A4 in Time Domain')
        
        %Plot Sound File in Frequency Domain
        figure(figCount)
        plot(f, y_fft)
        xlim([0 1000])
        xlabel('Frequency (Hz)')
        ylabel('Amplitude')
        title(['Frequency Response of Tuning Fork: ' num2str(titlename)])
        
        saveas(gcf, filename)
        
        figCount = figCount+1;
        clf; close all;
        
    end
end
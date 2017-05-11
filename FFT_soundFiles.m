%% FFT on sounds at different frequencies - and different wave shapes
%Determining where harmonics lie

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
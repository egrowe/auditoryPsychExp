%function runExp(pNo,block)
%% -----------------------------------------------------------------------%
%------------------ Auditory Perception Task (Qualia) ------------------%
%-------------- Testing Likeness of Auditory Experience ----------------%
%-----------------------------------------------------------------------%
clear all;
%Written by Elise Rowe May 2017
%PhD Project in Tsuchiya t-lab Monash Neuroscience of Consciousness
clc
%----------------------------------------------------------------------%
%-------------- PARTICIPANT DETAILS & BLOCK SELECTION ------------------%
%-----------------------------------------------------------------------%
Priority(2); %set priority in Matlab to ensure timing precision
%TR.subjNo = 1%pNo;       %number of the participant
ListenChar(2); % prevent output from keypresses displayed in command window

%% ---------------------------------------------------------------------%
%----------------- INITIALISE SCREEN & ADD PARAMETERS ------------------%
%-----------------------------------------------------------------------%
Screen('Preference', 'SkipSyncTests', 1);
[windowPtr, rect] = Screen('OpenWindow', 0, 1);%, [0 0 500 500]); %uncomment parentheses for editing
Cfg.windowPtr = windowPtr; %window pointer number
Cfg.winRect = rect; Cfg.width = rect(3); Cfg.height = rect(4); %determine screen size
Cfg.xCentre = rect(3)/2; Cfg.yCentre = rect(4)/2; Cfg.black = [0 0 0]; %setup colours
Cfg.WinColor = [0 0 0]; Cfg.white = [255 255 255]; Cfg.gray = [130 130 130]; %setup colours
HideCursor; %hide cursor from screen

t = clock;
rng(t(3) * t(4) * t(5),'twister'); %setup the random number generator accoriding to clock
Cfg.aux_buffer = 0;
Screen('TextSize',Cfg.windowPtr,20); %%Set the size of the text displayed to participants

%-----------------------------------------------------------------------%
%% --------------------- SETUP THE STIMULUS TABLE ----------------------%
%-----------------------------------------------------------------------%

%freqs = [250,500,630,800,1000,1250,1600,2000,2500,3150,5040,8000];
freqs = [250.00,281.25,312.50,343.75,375.00,406.25,437.50,468.75,500.00]; %better setup/coding of this!
%,531.25,562.50,593.75,625.00,656.25,687.50,718.75,750.00]
nFreq = length(freqs);

waveType = [1,2,3,4];
nWaves = length(waveType);

ntotalStim = nFreq*nWaves;

stimMatrix = triu(ones(ntotalStim,ntotalStim),0); %select freqs out of here!
%load('newStimMatrix.mat')

resultsMatrix = zeros(ntotalStim,ntotalStim); %setup to record responses to freqs
%load('Elise_5068_TrialsComplete.mat','resultsMatrix')
%resultsMatrix = combMatrix;

%The above sets some freqs to be the same as the order is not important
[row,col] = find(stimMatrix);
nPairs = length(row);

for ii = 1:nPairs
    stimPairings{ii} = [row(ii),col(ii)];
end

orderStim = randperm(nPairs);
noTrials = nPairs;


%% ---------------------------------------------------------------------%
%----------------------- PRESENT INSTRUCTIONS --------------------------%
%-----------------------------------------------------------------------%

%Blank screen (clear previous text)
Screen('FillRect', Cfg.windowPtr, 0);
for m = 1 : 45 % In Frames! (eg 60hz = 15 = 250ms)
    Screen('Flip', Cfg.windowPtr, [], Cfg.aux_buffer);
end

%Give instructions to participants about which ear to attend depending on which trial type is selected
DrawFormattedText(Cfg.windowPtr, 'Pay attention to the TWO tones and rate their similarity', 'center', Cfg.yCentre, [255 255 255]);
%Screen('DrawText', Cfg.windowPtr, 'Pay attention to the TWO tones and rate their similarity', ...
%    Cfg.xCentre-370, Cfg.yCentre, [255 255 255]);

Screen('Flip', Cfg.windowPtr, [], Cfg.aux_buffer);
KbWait();

%Blank screen (clear previous text)
Screen('FillRect', Cfg.windowPtr, 0);
for m = 1 : 60 % In Frames! (eg 60hz = 15 = 250ms)
    Screen('Flip', Cfg.windowPtr, [], Cfg.aux_buffer);
end

%Ready screen
DrawFormattedText(Cfg.windowPtr, 'Ready... (Press Spacebar to Continue)', 'center', Cfg.yCentre, [255 255 255]);
Screen('Flip', Cfg.windowPtr, [], Cfg.aux_buffer);
KbWait();

%Blank screen
Screen('FillRect', Cfg.windowPtr, 0);
for m = 1 : 30 % In Frames! (eg 60hz = 15 = 250ms)
    Screen('Flip', Cfg.windowPtr, [], Cfg.aux_buffer);
end


%% ---------------------------------------------------------------------%
%----------------------- START OF EXPERIMENT ---------------------------%
%--------------------- Start of Auditory Trials ------------------------%
%-----------------------------------------------------------------------%
cd('/Users/EliseRowe/Desktop/AuditoryPsych_Exp/audStim_waveDiff')
InitializePsychSound;

for trialNo = 1:noTrials
    
    % Different frequency rates for different sound files! IMPORTANT for playback timing
    sampRateTones = 44100; %sampling rate of the pure tones (kHz)
    pahandle = PsychPortAudio('Open', [], [], 2, sampRateTones, 2);
    
    % Stim table that randomly cycles through tone pairs
    toneCoords = stimPairings{orderStim(trialNo)};
    
    % SETUP FOR TONE A
    freqToneA = freqs(ceil(toneCoords(1)/nWaves)); %Frequency
    wave_ToneA = waveType((toneCoords(1)-((ceil(toneCoords(1)/nWaves))*nWaves))+nWaves);    
    
    % SETUP FOR TONE B
    freqToneB = freqs(ceil(toneCoords(2)/nWaves)); %Frequency
    wave_ToneB = waveType((toneCoords(2)-((ceil(toneCoords(2)/nWaves))*nWaves))+nWaves);    

    %Load the appropriate wav file for the block depending on trialBlock
    [toneA, FsA] = audioread(sprintf('%gHz_%d.wav',freqToneA,wave_ToneA));
    toneA = [toneA, toneA];
    [toneB, FsB] = audioread(sprintf('%gHz_%d.wav',freqToneB,wave_ToneB));
    toneB = [toneB, toneB];
    
    %% ---------------------------------------------------------------------%
    %--- Setup the Timing Structure to Record Responses and Send Triggers --%
    %-----------------------------------------------------------------------%
    
    %Fill audio buffer with the trial sound to be played
    PsychPortAudio('FillBuffer', pahandle, toneA'); % TONE A PLAY
    PsychPortAudio('Start', pahandle, [], [], 1);
    WaitSecs(0.5);
    
    PsychPortAudio('DeleteBuffer');
    PsychPortAudio('FillBuffer', pahandle, toneB'); % TONE B PLAY
    PsychPortAudio('Start', pahandle, [], [], 1);
    WaitSecs(0.5);
    %     %% ---------------------------------------------------------------------%
    %     %----------- RECORD RESPONSE TO 8AFC AFTER EACH TRIAL ------------------%
    %     %-----------------------------------------------------------------------%
    %
    %     %present the 8AFC screen
    question_string = sprintf('How similar or different were the two sounds?');
    
    ShowCursor;
    
    Cfg = response_screen(Cfg,question_string,'SAME','DIFFERENT');
    
    WaitSecs(.3);
    
    %% COLLECT A RESPONSE
    clicks = 0;
    respStartTime = GetSecs;
    
    % Wait until subject has given a response
    while clicks == 0
        
        [clicks, x,y] = GetClicks();
        
        % Check whether the click went inside a box area
        for m = 1 : size(Cfg.polyL, 1)
            idxs_left(m) = inpolygon(x,y,squeeze(Cfg.polyL(m,1,:)),squeeze(Cfg.polyL(m,2,:))); %#ok<*AGROW>
            
            idxs_right(m) = inpolygon(x,y,squeeze(Cfg.polyR(m,1,:)),squeeze(Cfg.polyR(m,2,:)));
        end
        
        idx_pos_left = find(idxs_left == 1);
        idx_pos_right = find(idxs_right == 1);
        
        if sum(idx_pos_left) > 0
            [kk, valuePos] = find(idxs_left == 1);
            possResp = [1 2 3 4];
            partResponse = possResp(valuePos);
        elseif sum(idx_pos_right) > 0
            [kk,valuePos] = find(idxs_right == 1);
            possResp = [-1 -2 -3 -4];
            partResponse = possResp(valuePos);
        end
        
        % Left boxes click . %% NEED TO ASSIGN 1-4 rating!!
        if length(idx_pos_left) == 1 %~isempty(idx_pos_left)
            keyid = -1;
            keyid2 = idx_pos_left;
            
            clicks = 1;
            
            % Paint selected box blue
            Screen('FillPoly', Cfg.windowPtr, [0 0 255], squeeze(Cfg.polyL(idx_pos_left,:,:))',1);
            for wait = 1:10
                Screen('Flip', Cfg.windowPtr,  [], Cfg.aux_buffer);
            end
            
        end
        
        if length(idx_pos_right) == 1 %~isempty(idx_pos_right)
            keyid = 1;
            keyid2 = idx_pos_right;
            
            clicks= 1;
            
            % Paint selected box blue
            Screen('FillPoly', Cfg.windowPtr, [0 0 255], squeeze(Cfg.polyR(idx_pos_right,:,:))',1);
            for wait = 1:10
                Screen('Flip', Cfg.windowPtr,  [], Cfg.aux_buffer);
            end
            
        end
    end
    
    respEndTime = GetSecs;
    reactionTime = respEndTime - respStartTime; % In seconds
    
    
    %setup structures to save the responses
    % Define response to TR
    if keyid == -1
        % Response on left: 'present'
        TR(trialNo).mouseResponsesPer = 'SAME';
        
    elseif keyid == 1
        % Response on right: 'absent'
        TR(trialNo).mouseResponsesPer = 'DIFFERENT';
    end
    
    %% EXAMINE RESPONSE
    TR(trialNo).response = partResponse;
    TR(trialNo).peripheral_confid = partResponse;
    TR(trialNo).peripheral_mouse_pos = [x y];
    %here, we want to take the coordinates of the selected Tone A and Tone B
    %and then place the value assigned in the 8AFC step to these coords in the
    %matrix (one-by-one) til we have filled all of the matrix (well one half on
    %the diagonal)
    
    % UPDATE MATRIX
    resultsMatrix(toneCoords(1),toneCoords(2)) = partResponse; %exact tone seq.
    resultsMatrix(toneCoords(2),toneCoords(1)) = partResponse; %opposite tone seq.
    
    clear toneA
    clear toneB
    PsychPortAudio('DeleteBuffer');
    
    WaitSecs(1);
    
end

%end

%Close the audio driver
PsychPortAudio('Close', pahandle);

ListenChar(0); % Turn responses back on at command window

ShowCursor; %Show the cursor again

%Blank screen at end of trial
Screen('FillRect', Cfg.windowPtr, 0);
for m = 1 : 60 % In Frames! (eg 60hz = 15 = 250ms)
    Screen('Flip', Cfg.windowPtr, [], Cfg.aux_buffer);
end


%% ---------------------------------------------------------------------%
%---------------- END EXPERIMENT: Clean up and go home -----------------%
%-----------------------------------------------------------------------%
%Thank you for participating screen
DrawFormattedText(Cfg.windowPtr, 'End of Experiment / Thanks for Participating', 'center', Cfg.yCentre, [255 255 255]);
Screen('Flip', Cfg.windowPtr, [], Cfg.aux_buffer);
WaitSecs(2);
%% ---------------------------------------------------------------------%%
% ------------------------- SAVE ALL VARIABLES --------------------------%
% -------------------- Using SubjNo and TrialBlock ----------------------%
% ----------------------------------------------------------------------%%
%UPDATE the Stimulus Matrix with the trials we just completed and save
for RR = 1:length(ntotalStim)
    
    for CC = 1:length(ntotalStim)
        
        if resultsMatrix(RR,CC) ~= 0 %|| resultsMatrix(RR,CC) > 0
            
            stimMatrix(RR,CC) = 0;
        end
        
    end
end

trialsComplete = length(find(stimMatrix == 0))

filename = [num2str(participantNo) 'noFreqs_' num2str(nFreq) 'noWaves_' num2str(nWaves) '_' (trialsComplete) '_TrialsComplete.mat']
%save(datafilename, 'TR', 'stimPairings', 'stimMatrix', 'resultsMatrix', 'freqs', 'orderStim'); %save all files for the subject
save(filename)%, 'resultsMatrix')
% Clear all - viola! All done! (test for accuracy using getAccuracy_x.m)
sca

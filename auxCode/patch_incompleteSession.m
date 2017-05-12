%% 2 Scripts if testing session is disrupted/ends before complete

% Use these scripts to either (1) combine 2 separate results sessions OR
% (2) reupdate the stimulus matrix to not replay these tones (will also
% need to load this stimMatrix and resultsMatrix file before starting next
% testing session

filename = [' LOAD YOUR FILENAME FOR DATA HERE ']
load(filename, 'resultsMatrix', 'stimMatrix')
resultsMatrixB = resultsMatrix

filenameSecond = [' LOAD YOUR FILENAME FOR SECOND SESSION DATA HERE ']
load(filename, 'resultsMatrix', 'stimMatrix')

%Combine results matrices from 2 separate sessions of trials
%totalSize = length(resultsMatrix)*length(resultsMatrix)
combMatrix = resultsMatrix;

for RR = 1:length(resultsMatrix)
    
    for CC = 1:length(resultsMatrix)
        
        if resultsMatrix(RR,CC) == 0 %if no response recorded at this location..
            combMatrix(RR,CC) = resultsMatrixB(RR,CC); %replace with value in 2nd matrix (even if it is empty too)
        else
            continue
        end
        
    end
end


%% Alter the stimulus matrix for the next session (no duplicate pairings)

filename = [' LOAD YOUR FILENAME FOR DATA HERE ']
load(filename, 'resultsMatrix', 'stimMatrix')

%Create stimulus matrix (trials to use) after some have already been
%partially completed earlier

for RR = 1:length(resultsMatrix)
    
    for CC = 1:length(resultsMatrix)
        
        if resultsMatrix(RR,CC) ~= 0 %|| resultsMatrix(RR,CC) > 0 
            % if results is recorded (i.e not zero), delete from stim selection matrix
            stimMatrix(RR,CC) = 0;
        end
        
    end
end

%This STIMULUS MATRIX and the previous RESULTS matrix will need to be added
%to the runExp script before testing a new session :)

%% -----------------------------------------------------------------------%
%------------------ Rearrange Results Matrix (Hz/Wave) -------------------%
%--------- Change format of Matrix when in Hz and Wave format ------------%
%-------------------------------------------------------------------------%

%% Plot results using imagesc according to Hz
clear all
load('FINAL_Elise_RESULTS_allTrials.mat', 'resultsMatrix')

create_Hz_fig(resultsMatrix)

%% Plot results using imagesc according to wavetype
clear all
load('FINAL_Elise_RESULTS_allTrials.mat', 'resultsMatrix', 'nFreq', 'nWaves')

%Setup indexes for the 2 categories (dB and Hz)
useLength = length(resultsMatrix)

%Indexes for the dB segmentation
start = 1;
for ii = 1:nFreq
    wave_Idx(ii) = start + ((ii-1)*nWaves)
end

%Indexes for the Hz segmentation
start = 1;
for rr = 1:nWaves
    Hz_Idx(rr) = start + ((rr-1)*nFreq)
end
Hz_Idx(rr+1) = useLength; %add on the final point

%Setup new empty results matrix
resultsMatrix_waves = zeros(length(resultsMatrix),length(resultsMatrix));

for nextSubset = 1:length(Hz_Idx)-1
    
    takeRowfrom = (wave_Idx)+(nextSubset-1);
    useRowAddition = Hz_Idx(nextSubset)-1;

    for withindBPair = 1:length(Hz_Idx)-1
        
        takeColfrom = (wave_Idx)+(withindBPair-1);
        useColAddition = Hz_Idx(withindBPair)-1;
        
        for newMatRow = 1:length(wave_Idx)
            
            for newMatCol = 1:length(wave_Idx)
                
                resultsRow = takeRowfrom(newMatRow);
                resultsCol = takeColfrom(newMatCol);
                
                resultsMatrix_waves((newMatRow+useRowAddition),(newMatCol+useColAddition)) = resultsMatrix(resultsRow,resultsCol);
            end
        end
    end
end

imagesc(resultsMatrix_waves);

create_waves_fig(resultsMatrix_waves);
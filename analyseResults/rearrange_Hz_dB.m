%% -----------------------------------------------------------------------%
%-------------------- Rearrange Results Matrix (Hz/dB) -------------------%
%----------- Change format of Matrix when in Hz and dB format ------------%
%-------------------------------------------------------------------------%

%% Plot results using imagesc according to Hz
clear all
load('FINAL_Elise_RESULTS_allTrials.mat', 'resultsMatrix')

create_Hz_fig(resultsMatrix)


%% Plot results using imagesc according to dB
clear all
load('FINAL_Elise_RESULTS_allTrials.mat', 'resultsMatrix', 'nFreq', 'ndB')

%Setup indexes for the 2 categories (dB and Hz)
useLength = length(resultsMatrix)

%Indexes for the dB segmentation
start = 1;
for ii = 1:nFreq
    dB_Idx(ii) = start + ((ii-1)*ndB)
end

%Indexes for the Hz segmentation
start = 1;
for rr = 1:ndB
    Hz_Idx(rr) = start + ((rr-1)*nFreq)
end
Hz_Idx(rr+1) = useLength; %add on the final point

%Setup new empty results matrix
resultsMatrix_dB = zeros(length(resultsMatrix),length(resultsMatrix));

for nextSubset = 1:length(Hz_Idx)-1
    
    takeRowfrom = (dB_Idx)+(nextSubset-1);
    useRowAddition = Hz_Idx(nextSubset)-1;

    for withindBPair = 1:length(Hz_Idx)-1
        
        takeColfrom = (dB_Idx)+(withindBPair-1);
        useColAddition = Hz_Idx(withindBPair)-1;
        
        for newMatRow = 1:length(dB_Idx)
            
            for newMatCol = 1:length(dB_Idx)
                
                resultsRow = takeRowfrom(newMatRow);
                resultsCol = takeColfrom(newMatCol);
                
                resultsMatrix_dB((newMatRow+useRowAddition),(newMatCol+useColAddition)) = resultsMatrix(resultsRow,resultsCol);
            end
        end
    end
end

imagesc(resultsMatrix_dB);

create_dB_fig(resultsMatrix_dB);

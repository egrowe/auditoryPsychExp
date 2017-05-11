%% Plot results using imagesc according to Hz
clear all
load('FINAL_Elise_RESULTS_allTrials.mat', 'resultsMatrix')

create_Hz_fig(resultsMatrix)


%% Plot results using imagesc according to dB
clear all
load('FINAL_Elise_RESULTS_allTrials.mat', 'resultsMatrix')

%Setup indexes for the 2 categories (dB and Hz)
dB_Idx = [1,7,13,19,25,31,37,43,49,55,61,67]
Hz_Idx = [1, 13, 25, 37, 49, 61, 72]

%Setup new empty results matrix
resultsMatrix_dB = zeros(72,72);

for nextSubset = 1:6
    
    takeRowfrom = (dB_Idx)+(nextSubset-1)
    useRowAddition = Hz_Idx(nextSubset)-1

    for withindBPair = 1:6
        
        takeColfrom = (dB_Idx)+(withindBPair-1)
        useColAddition = Hz_Idx(withindBPair)-1
        
        for newMatRow = 1:12
            
            for newMatCol = 1:12
                
                resultsRow = takeRowfrom(newMatRow)
                resultsCol = takeColfrom(newMatCol)
                
                resultsMatrix_dB((newMatRow+useRowAddition),(newMatCol+useColAddition)) = resultsMatrix(resultsRow,resultsCol);
            end
        end
    end
end

%imagesc(resultsMatrix_dB)

create_dB_fig(resultsMatrix_dB);

%imagesc(resultsMatrix)
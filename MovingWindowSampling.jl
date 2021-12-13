using Statistics
using Combinatorics

function findNanRows(data)
    #filter out rows with any nan occurances
    nanRows = falses(size(data,1))
    for i = 1:size(data,2)
        nanRows[isnan.(data[:,i])] .= 1
    end
    return nanRows
end

function multipleRegression(depVar,predVars)
    #calculate regression coefficients using linear inversion
    coeffs = ((predVars' * predVars)^-1 * predVars') * depVar
    r2 = cor(depVar,predVars*coeffs)
    return coeffs,r2
end

function multipleCorrelation(depVar,predVars)
    #filter out rows with any nan occurances
    predVarsSize = size(predVars)
    #calculate correlation coefficients for each dependent/predictor variable pair
    coeffs = Array{Float64}(undef,predVarsSize[2])
    for i = 1:predVarsSize[2]
        coeffs[i] = cor(depVar,predVars[:,i])
    end
    return coeffs
end

function dominanceRecursion(depVar,predVars, parentCombination,parentR2,dominanceResultsMatrix::Matrix{Vector{Float64}})
    #part of the function dominanceAnalysis
    #calculate the difference in r2 values between regressions with a larger number of
    #predictor variables (parent) and regressions with combination of predictor
    #variables one less than the first regrssion (child). Do this recursively to test
    #the regressions of every combination of predictor variables.


    #base case with a single predictor variable, use the r2 of the single regression
    if length(parentCombination) <= 1
        ~,childR2 =multipleRegression(depVar, predVars[:,parentCombination])
        missingIdx = parentCombination[1]
        rowIdx = 1
        differenceList = copy(dominanceResultsMatrix[rowIdx,missingIdx])
        append!(differenceList,childR2)
        dominanceResultsMatrix[rowIdx,missingIdx] = differenceList

    else
        #find combinations of the parent combination with one less variable
        childCombinations = collect(combinations(parentCombination,length(parentCombination)-1))

        #for each combination find the r2 value, calculate the difference between it
        #and the parent r2 value, and record it in the right location of the results matrix
        for i = 1:length(childCombinations)

            ~,childR2 = multipleRegression(depVar, predVars[:,childCombinations[i]])
            missingIdx = sum(parentCombination)-sum(childCombinations[i])
            rowIdx = length(parentCombination)

            differenceList = copy(dominanceResultsMatrix[rowIdx,missingIdx])
            append!(differenceList,parentR2-childR2)
            dominanceResultsMatrix[rowIdx,missingIdx] = differenceList

            dominanceRecursion(depVar, predVars, childCombinations[i],childR2,dominanceResultsMatrix)
        end
    end
    return dominanceResultsMatrix
end

function dominanceAnalysis(depVar,predVars)
    #calculate dominance coefficients by determining the indiviual contribution of
    #each predictor variable for regression of every combination of predictor variables

    predVarsSize = size(predVars)

    #enter the full set of combinations into function dominanceRecursion
    allPredVarsCombinations = collect(combinations(1:predVarsSize[2]))[end]
    ~,allPredVarsR2 = multipleRegression(depVar,predVars)
    dominanceResultsMatrix = fill(Float64[], predVarsSize[2], predVarsSize[2])
    dominanceResultsMatrix = dominanceRecursion(depVar,predVars, allPredVarsCombinations,allPredVarsR2,dominanceResultsMatrix)

    #average the relative contribution of each predictor variable for all regressions
    #with the same number of predictor variables
    dominanceAvgMatrix = Array{Float64,2}(undef,predVarsSize[2], predVarsSize[2])
    for i = 1:predVarsSize[2]
        for j = 1:predVarsSize[2]
            dominanceAvgMatrix[i,j] = mean(dominanceResultsMatrix[i,j])
        end
    end

    #average the average contributions of each predictor variable for all regressions
    #regardless of the number of predictor variables to determine raw coefficients
    dominanceCoeffsRaw = vec(mean(dominanceAvgMatrix, dims = 1))

    #normalize raw coefficients to determine relative contribution
    dominanceCoeffsNorm = copy(dominanceCoeffsRaw)
    dominanceCoeffsNorm[dominanceCoeffsNorm .< 0] .= 0
    dominanceCoeffsNorm = dominanceCoeffsNorm./(sum(dominanceCoeffsNorm))

    return dominanceCoeffsNorm,dominanceCoeffsRaw,dominanceResultsMatrix,dominanceAvgMatrix
end

function movingWindowSampling(depVar, predVars, heights, windowSize::Int64,stepSize::Int64,onlyReturnWindowNum,onlyReturnCoeffs;analysis)

    #remove rows with nan occurances
    nanIdx = findNanRows([depVar predVars heights]);
    depVar = depVar[nanIdx .== 0]
    predVars = predVars[nanIdx .== 0,:]
    heights = heights[nanIdx .== 0]

    #produce StepRange of starting values for window iterations
    startingWindowVals = floor(Int,minimum(heights)):stepSize:ceil(Int64,(maximum(heights)-windowSize));
    #if the stepRange does not cover all height values, increase by the stepSize until it does
    if startingWindowVals[end] + windowSize < maximum(heights)
        startingWindowVals = StepRange(startingWindowVals[1], stepSize, startingWindowVals[end] + stepSize)
    end
    #calculate ending values for the window iterations
    endingWindowVals = startingWindowVals .+ windowSize;

    #find the number of window iterations and hosts
    windowNum = length(startingWindowVals);
    predVarNum = size(predVars,2);

    #if the fucntion is being used only to find the number of window iterations
    if onlyReturnWindowNum
    coeffItr = NaN ;
    allCoeffDataItr = NaN;
    return coeffItr, allCoeffDataItr, windowNum #need way to end function here
    end

    #if only coefficients are needed
    if onlyReturnCoeffs
        allCoeffDataItr = NaN;

    #if all regression data is to be collected
    elseif analysis =="MRA"
        allCoeffDataItr = Array{Any,2}(undef,windowNum,3);
    elseif analysis =="cor"
        allCoeffDataItr = Array{Any,2}(undef,windowNum,2);
    elseif analysis =="dom"
        allCoeffDataItr = Array{Any,2}(undef,windowNum,5);
    end
    coeffItr = Array{Float64,2}(undef,windowNum,predVarNum);


    #begin moving window up dataset, collecting data at each iteration
    for windowItr = 1:windowNum
        #find indexes of heights between window edges
        paneIdx = heights .>= startingWindowVals[windowItr] .&& heights .<= endingWindowVals[windowItr];

        #run chosen analysis
        if analysis == "MRA"
            coeffs,r2 = multipleRegression(depVar[paneIdx],predVars[paneIdx,:])
            if !onlyReturnCoeffs
                #record coeffs, relevant heights,and r2 values
                allCoeffDataItr[windowItr,1] = coeffs
                allCoeffDataItr[windowItr,2] = heights[paneIdx]
                allCoeffDataItr[windowItr,3] = r2
            end
            coeffItr[windowItr,:] .= coeffs
        elseif analysis == "cor"
            coeffs = cor(depVar[paneIdx],predVars[paneIdx,:], dims = 1)
            if !onlyReturnCoeffs
                #record coeffs and relevant heights i
                allCoeffDataItr[windowItr,1] = coeffs
                allCoeffDataItr[windowItr,2] = heights[paneIdx]
            end
            coeffItr[windowItr,:] = coeffs
        elseif analysis == "dom"
            coeffs,coeffsRaw,resultsMatrix,avgMatrix = dominanceAnalysis(depVar[paneIdx],predVars[paneIdx,:])
            if !onlyReturnCoeffs
                #record coeffs, relevant heights, raw dominance coeffs, the full
                #results matrix, and the averages matrix
                allCoeffDataItr[windowItr,1] = coeffs
                allCoeffDataItr[windowItr,2] = heights[paneIdx]
                allCoeffDataItr[windowItr,3] = coeffsRaw
                allCoeffDataItr[windowItr,4] = resultsMatrix
                allCoeffDataItr[windowItr,5] = avgMatrix
            end
            coeffItr[windowItr,:] .= coeffs
        end

    end
    return coeffItr,allCoeffDataItr,windowNum
end

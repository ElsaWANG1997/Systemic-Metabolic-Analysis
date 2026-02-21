%% Individual metabolic network construction via third-order polynomial TAC fitting
% Reads harmonized SUL data, normalizes by liver TAC, and saves per-patient Z-score connectivity maps.

clear; clc;

% Input and output configuration
excelPath = ' ';
outputDir = ' ';

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Load the Excel table preserving original column labels
rawTbl = readtable(excelPath, 'PreserveVariableNames', true);
patientIDs = string(rawTbl.Subject);
numPatients = height(rawTbl);

varNames = rawTbl.Properties.VariableNames;

smokingVarIdx = find(strcmpi(varNames, 'Smoking'), 1);
diseaseVarIdx = find(strcmpi(varNames, 'Disease'), 1);
if isempty(smokingVarIdx) || isempty(diseaseVarIdx)
    error('Expected Smoking and Disease columns in %s.', excelPath);
end

smokingVals = rawTbl{:, smokingVarIdx};
diseaseVals = rawTbl{:, diseaseVarIdx};

if isa(smokingVals, 'categorical')
    smokingVals = str2double(string(smokingVals));
elseif ~isnumeric(smokingVals)
    smokingVals = str2double(smokingVals);
end

if isa(diseaseVals, 'categorical')
    diseaseVals = str2double(string(diseaseVals));
elseif ~isnumeric(diseaseVals)
    diseaseVals = str2double(diseaseVals);
end

smokingVals = double(smokingVals(:));
diseaseVals = double(diseaseVals(:));

if any(isnan(smokingVals)) || any(isnan(diseaseVals))
    error('Smoking or Disease columns contain missing values in %s.', excelPath);
end

smokingUnique = unique(smokingVals);
diseaseUnique = unique(diseaseVals);
if any(~ismember(smokingUnique, [0, 1]))
    warning('Smoking column in %s contains unexpected codes.', excelPath);
end
if any(~ismember(diseaseUnique, 1:4))
    warning('Disease column in %s contains unexpected codes.', excelPath);
end

% Determine cohort sizes per patient (Smoking 0/1 x Disease 1-4)
patientGroupMatrix = [smokingVals, diseaseVals];
[uniqueGroups, ~, patientGroupIdx] = unique(patientGroupMatrix, 'rows', 'stable');
patientGroupCounts = accumarray(patientGroupIdx, 1);
N_pat_values = patientGroupCounts(patientGroupIdx);

if numel(patientGroupCounts) > 8
    warning('Detected %d Smoking/Disease cohorts in %s (expected at most 8).', numel(patientGroupCounts), excelPath);
end

roiStartCol = 8; % columns 8:end contain ROI TAC information
roiVarNames = varNames(roiStartCol:end);

roiMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
roiOrder = {};
numTimePts = 30;
timeAxis = (1:numTimePts)';

for i = 1:numel(roiVarNames)
    varName = roiVarNames{i};
    token = regexp(varName, '^(?<roi>.+?)_T(?<time>\d+)(?:_(?<rep>\d+))?$', 'names');
    if isempty(token)
        warning('Skipping unrecognized column name: %s', varName);
        continue;
    end

    roiName = token.roi;
    timeIdx = str2double(token.time);
    if timeIdx < 1 || timeIdx > numTimePts
        error('Time index out of range for column %s', varName);
    end

    if ~isKey(roiMap, roiName)
        roiMap(roiName) = nan(1, numTimePts);
        roiOrder{end + 1} = roiName; 
    end

    colIndices = roiMap(roiName);
    colIndices(timeIdx) = roiStartCol - 1 + i;
    roiMap(roiName) = colIndices;
end

roiNames = roiOrder;

% Remove Aorta from analysis ROI set while keeping original mapping intact
analysisROIs = roiNames;
aortaMask = strcmpi(analysisROIs, 'Aorta');
analysisROIs(aortaMask) = [];

numROIs = numel(analysisROIs);

% Confirm that every ROI being analysed has a full 30-frame TAC
for r = 1:numROIs
    assignedCols = roiMap(analysisROIs{r});
    if any(isnan(assignedCols))
        error('ROI %s is missing one or more time points.', analysisROIs{r});
    end
end

if ~isKey(roiMap, 'Liver')
    error('Liver ROI was not found in the dataset.');
end

liverCols = roiMap('Liver');
liverTAC = rawTbl{:, liverCols};

% Guard against zero liver activity values
zeroMask = abs(liverTAC) < eps;
liverTAC(zeroMask) = eps;

relSUL = zeros(numPatients, numTimePts, numROIs);
for r = 1:numROIs
    currentCols = roiMap(analysisROIs{r});
    roiTAC = rawTbl{:, currentCols};
    relSUL(:, :, r) = roiTAC ./ liverTAC;
end

% Ensure invalid ratios stay marked and clip any negative values to zero
relSUL(~isfinite(relSUL)) = NaN;
relSUL(relSUL < 0) = 0;


for p = 1:numPatients
    N_pat = N_pat_values(p);
    if ~isfinite(N_pat)
        error('Unable to determine N_pat for patient %d (Smoking=%g, Disease=%g).', p, smokingVals(p), diseaseVals(p));
    elseif N_pat < 2
        warning('Group size for patient %d (Smoking=%g, Disease=%g) is less than 2. Results may be unstable.', p, smokingVals(p), diseaseVals(p));
    end

    residuals = nan(numTimePts, numROIs);
    for r = 1:numROIs
        y = reshape(relSUL(p, :, r), [], 1);
        validIdx = isfinite(y);
        if sum(validIdx) >= 4
            coeffs = polyfit(timeAxis(validIdx), y(validIdx), 3);
            fitted = polyval(coeffs, timeAxis);
            residualVec = y - fitted;
            residualVec(~validIdx) = NaN;
            residuals(:, r) = residualVec;
        else
            residuals(:, r) = NaN;
        end
    end

    corrMat = corrcoef(residuals, 'Rows', 'pairwise');
    corrMat(~isfinite(corrMat)) = 0;
    corrMat = max(min(corrMat, 0.999999), -0.999999);
    corrMat(1:numROIs + 1:end) = 0;

    fisherZ = 0.5 * log((1 + corrMat) ./ (1 - corrMat));
    fisherZ(~isfinite(fisherZ)) = 0;

    denom = 1 - fisherZ .^ 2;
    nearZero = abs(denom) < eps;
    denom(nearZero) = eps .* sign(denom(nearZero));
    ZCC = fisherZ .* (N_pat - 1) ./ denom;
    ZCC(~isfinite(ZCC)) = 0;
    ZCC(1:numROIs + 1:end) = 0;

    patientLabel = strtrim(patientIDs(p));
    if strlength(patientLabel) == 0
        patientLabel = "Patient" + p;
    end
    cleanID = regexprep(char(patientLabel), '[^0-9A-Za-z_-]', '_');

    outPath = fullfile(outputDir, [cleanID '.mat']);
    roiNamesOut = analysisROIs(:); 
    timeAxisOut = timeAxis; 
    save(outPath, 'ZCC', 'fisherZ', 'corrMat', 'roiNamesOut', 'timeAxisOut', 'patientLabel');
end

fprintf('Processed %d patients. Outputs saved to %s\n', numPatients, outputDir);

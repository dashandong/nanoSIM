% Image Reconstruction for Structured Illumination Microscopy(SIM)
% Copyright (C) 2025 by Dr. Dashan Dong (dongdashan@pku.edu.cn)
%
% This program implements the image reconstruction algorithm for SIM
% using a simalar approach as the one described in:
% Perez, V., Chang, BJ. & Stelzer, E. "Optimal 2D-SIM reconstruction
% by two filtering steps with Richardson-Lucy deconvolution"
% Sci Rep 6, 37149 (2016).
%
% The Energy conservation is carefully considered in this
% implementation of the algorithm. To enable a fair comparison
% between the SIM and wide-field images, the energy of the images
% is normalized to ensure consistency across different imaging modalities.
% This normalization process helps to mitigate the effects of varying
% illumination and detector sensitivity, providing a more accurate
% representation of the sample being imaged.
%
% This program is free software: you can redistribute it and / or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.

% last edited @ 2025/5/1 by Dashan Dong

function nanoSIM
    clc;

    fprintf('\tdsSIM  Copyright (C) 2025  by Dashan Dong\n');
    fprintf('\tThis program comes with ABSOLUTELY NO WARRANTY.\n');
    fprintf('\tThis is free software, and you are welcome to redistribute it\n');
    fprintf('\tunder certain conditions.\n\n');

    %% Basic Parameters
    % Number of SIM phases
    numPhase = 3;
    % Number of SIM directions
    numDirection = 3;
    % Iteration time for Lucy-Richardson Deconvolution
    paraLRIter = 5;
    % Magnification
    paraMag = 100;
    % Numerical aperture
    paraNA = 1.45;
    % Pixel size
    paraPixelSz = 6.5;
    % Gain for high frequency
    paraAttAmp = 0;
    % High gain range for high frequency
    paraAttSigma = 2 * pi * 1.5;
    % Wavelength for emission
    paraWavelength = 0.520;

    % ODT Dampening factor
    paraOTFDampening = 0.95;

    %% Select SIM raw data (standard Tiff format <4G)
    [filename_raw, pathname_raw] = uigetfile( ...
        '*.tif', 'Select SIM Raw Tiff Stack');

    if isequal(filename_raw, 0)
        error('SIM raw tiff selection cancelled');
    else
        disp('Open SIM raw tiff file:');
        disp(fullfile(pathname_raw, filename_raw));
    end

    %% Reading size of raw data
    warning off;
    tiff_raw = Tiff([pathname_raw, filename_raw], 'r');
    szWidth = tiff_raw.getTag('ImageWidth');
    szHeight = tiff_raw.getTag('ImageLength');
    szStack = 1;

    while (~tiff_raw.lastDirectory())
        tiff_raw.nextDirectory();
        szStack = szStack + 1;
    end

    warning on;
    disp(['Tiff Size: ', num2str(szWidth), ' x ', num2str(szHeight), ' x ', num2str(szStack)]);

    if (szStack < numPhase * numDirection)
        error('Raw Tiff stack is too short to be SIM rawdata!');
    end

    %% Stack and average raw data
    disp('Start to solve the parameters...');
    numStackGroup = input([ ...
                               'How many time points to stack together?\n', ...
                               '(Range: 1-', ...
                               num2str(floor(szStack / (numPhase * numDirection))), ...
                           ', default 10):']);

    if isempty(numStackGroup)
        numStackGroup = 10;
    end

    numStackGroup = round(numStackGroup);

    if (numStackGroup < 1) || (numStackGroup > floor(szStack / (numPhase * numDirection)))
        error('Wrong number input!');
    end

    disp(['Stack ', num2str(numStackGroup), ' time points to stack.']);

    numStackStart = 1;

    if (szStack > numStackGroup * numPhase * numDirection)
        numStackStart = input([ ...
                                   'Start stacking from frame #?\n', ...
                                   '(Range: 1-', ...
                                   num2str(szStack - numStackGroup * numPhase * numDirection + 1), ...
                               ', default 1):']);

        if isempty(numStackStart)
            numStackStart = 1;
        end

        numStackStart = round(numStackStart);

        if (numStackStart < 1) || (numStackStart > szStack - numStackGroup * numPhase * numDirection + 1)
            error('Wrong number input!');
        end

    end

    imgRaw = zeros(szWidth, szHeight, numPhase * numDirection, 'single');
    warning off;

    for i = 1:numStackGroup %#ok<FXUP>

        for k = 1:(numPhase * numDirection)
            tiff_raw.setDirectory(numStackStart + (i - 1) * (numPhase * numDirection) + (k - 1));
            imgRaw(:, :, k) = imgRaw(:, :, k) + single(tiff_raw.read());
        end

    end

    %% Average the stack
    imgRaw = imgRaw ./ numStackGroup;
    warning on;
    tiff_raw.close();

    imgWF = sum(imgRaw, 3);

    %% Read OTF
    [filename_otf, pathname_otf] = uigetfile( ...
        '*.tif', 'Select OTF');

    if isequal(filename_otf, 0)
        error('OTF Selection Cancelled!');
    end

    OTF = im2double(imread([pathname_otf, filename_otf]));
    OTF = imresize(OTF, [szWidth, szHeight], 'bilinear');
    OTF = OTF / max(OTF(:));
    OTFx2 = zeros(size(OTF) * 2);
    OTFx2((floor((szWidth - 1) / 2) + 2):(floor((szWidth - 1) / 2) + 1 + szWidth), ...
        (floor((szHeight - 1) / 2) + 2):(floor((szHeight - 1) / 2) + 1 + szHeight)) = OTF;
    OTF = ifftshift(OTF);
    OTFx2 = ifftshift(OTFx2);

    % PSF = abs(fftshift(ifft2(OTF)));

    %% Calculate intermediate variables
    paraKm = 2 * pi * paraNA / paraWavelength;
    paraKWSz = 2 * pi / paraPixelSz / szWidth;
    paraKHSz = 2 * pi / paraPixelSz / szHeight;
    paraPhaseVector = ((1:numPhase) - 1) * 2 * pi / numPhase;
    paraPhaseVectorFull = ones(numPhase, numDirection);
    paraPhaseMatrix = ones(numPhase, 3);
    paraPhaseMatrixFull = ones(numPhase, 3, numDirection);
    paraIMatrix = ones(3, numDirection); % suppose a 100 % modulation depth

    % phase matrix
    for p = 1:numPhase %#ok<FXUP>
        paraPhaseVectorFull(p, :) = paraPhaseVector(p);
        paraPhaseMatrix(p, :) = paraPhaseMatrix(p, :) .* ...
            exp(1i * [-paraPhaseVector(p), 0, paraPhaseVector(p)]);
    end

    % intensity matrix
    for d = 1:numDirection %#ok<FXUP>
        % paraIMatrix(:, d) = [1, 0.5, 1];
        paraIMatrix(:, d) = [2, 1, 2];
    end

    % generate meshgrid
    axisW = ((1:szWidth) - 1 - floor(szWidth / 2)) / szWidth;
    axisH = ((1:szHeight) - 1 - floor(szHeight / 2)) / szHeight;
    % [meshW, meshH] = meshgrid(axisW, axisH);

    axisKW = ((1:szWidth) - 1 - floor(szWidth / 2)) * paraKWSz;
    axisKH = ((1:szHeight) - 1 - floor(szHeight / 2)) * paraKHSz;

    axisWx2 = ((1:(szWidth * 2)) - 1 - szWidth) / (szWidth * 2);
    axisHx2 = ((1:(szHeight * 2)) - 1 - szHeight) / (szHeight * 2);
    [meshWx2, meshHx2] = meshgrid(axisWx2, axisHx2);

    axisKWx2 = ((1:(szWidth * 2)) - 1 - szWidth) * paraKWSz;
    axisKHx2 = ((1:(szHeight * 2)) - 1 - szHeight) * paraKHSz;

    [meshKW, meshKH] = meshgrid (axisKW, axisKH);
    [~, meshKR] = cart2pol(meshKW, meshKH);
    meshKR_shift = ifftshift(meshKR);

    [meshKWx2, meshKHx2] = meshgrid(axisKWx2, axisKHx2);
    [~, meshKRx2] = cart2pol(meshKWx2, meshKHx2);
    meshKRx2_shift = ifftshift(meshKRx2);

    % OTF dampening mask
    OTFdampening = paraOTFDampening .^ meshKR_shift;
    OTFdampeningx2 = paraOTFDampening .^ meshKRx2_shift;

    % Gain function for high frequency
    attFun2x = 1 - paraAttAmp * exp(-meshKRx2_shift .^ 2 / (2 * paraAttSigma ^ 2));

    % Filter mask
    OTFMask = (meshKR_shift < 2 * paraKm);
    % OTFRingMask = (meshKR_shift < 2 * paraKm) & (meshKR_shift > 1.0 * paraKm);
    OTFMaskx2 = (meshKRx2_shift < 2 * paraKm);
    OTFRingMaskx2 = (meshKRx2_shift < 2 * paraKm) & (meshKRx2_shift > 0.5 * paraKm);

    paraIllVector = zeros(numDirection, 2);

    %% Lucy-Richardson Deconvolution over Raw Images
    % Default iteration number is 10
    for i = 1:(numPhase * numDirection) %#ok<FXUP>
        % imgRaw(:, :, i) = deconvlucy(imgRaw(:, :, i), PSF, paraLRIter);
        imgRaw(:, :, i) = mydeconvlucy(imgRaw(:, :, i), OTF .* OTFdampening, paraLRIter);
    end

    % % Low-pass filter on RAW images
    % for i = 1:(numPhase * numDirection) %#ok<FXUP>
    %     imgRaw(:, :, i) = fftshift(ifft2(fft2(ifftshift(imgRaw(:, :, i))) .* OTFMask));
    % end

    imgWFdeconv = sum(imgRaw, 3);

    %% Background substraction with optical sectioned SIM
    osImg = zeros(szWidth, szHeight, numDirection, 'single');
    osBG = zeros(szWidth, szHeight, numDirection, 'single');

    for d = 1:numDirection %#ok<FXUP>

        for p = 1:numPhase %#ok<FXUP>

            for p2 = 1:numPhase %#ok<FXUP>

                if p2 ~= p
                    osImg(:, :, d) = osImg(:, :, d) + (imgRaw(:, :, p + (d - 1) * numPhase) - imgRaw(:, :, p2 + (d - 1) * numPhase)) .^ 2;
                end

            end

            osBG(:, :, d) = osBG(:, :, d) + imgRaw(:, :, p + (d - 1) * numPhase) / numPhase;
        end

        osImg(:, :, d) = sqrt(osImg(:, :, d) / (numPhase * (numPhase - 1)));
        osBG(:, :, d) = osBG(:, :, d) - osImg(:, :, d);
    end

    imgBG = sum(osBG, 3) .* numPhase;
    imgOS = imgWFdeconv - imgBG;

    %% Notch filter
    notchFilter = 1 - exp(- (meshKR_shift) .^ 4 ./ (2 * (0.5 * paraKm) .^ 4));
    notchFilterx2 = 1 - exp(- (meshKRx2_shift) .^ 4 ./ (2 * (0.5 * paraKm) .^ 4));
    % for i = 1:(numPhase * numDirection) %#ok<FXUP>
    %     imgRaw(:, :, i) = fftshift(ifft2(fft2(ifftshift(imgRaw(:, :, i))) .* OTFMask .* notchFilter));
    % end

    %% FFT of raw data
    spectrumRaw = complex(zeros(szHeight, szWidth, ...
        numPhase * numDirection, 'single'));

    % FFT of raw data
    for i = 1:(numPhase * numDirection)
        spectrumRaw(:, :, i) = fftshift(fft2(ifftshift(imgRaw(:, :, i))));
    end

    %% Separation of three frequency components in each direction
    spectrumSep = complex(zeros(3, szHeight, szWidth, numDirection));
    spectrumSep2x = complex(zeros(3, szHeight * 2, szWidth * 2, ...
        numDirection, 'single'));

    for d = 1:numDirection %#ok<FXUP>
        spectrumSepEachP = complex(zeros(numPhase, szHeight * szWidth));

        for p = 1:numPhase %#ok<FXUP>
            spectrumTemp = spectrumRaw(:, :, p + (d - 1) * numPhase);
            spectrumSepEachP(p, :) = spectrumTemp(:);
        end

        spectrumSepEach = paraPhaseMatrix \ spectrumSepEachP;
        spectrumSep((((d - 1) * 3 * szHeight * szWidth) + 1) ...
            :(d * 3 * szHeight * szWidth)) = spectrumSepEach(:);

        for i = 1:3
            spectrumSep(i, :, :, d) = paraIMatrix(i, d) .* ...
                spectrumSep(i, :, :, d) .* numPhase;
            % x numPhase because of the speration average over numPhase images
        end

    end

    % Padding the spectrum 2x
    for d = 1:numDirection

        for i = 1:3
            spectrumSep2x(i, (floor((szHeight + 1) / 2) + 1):(floor((szHeight + 1) / 2) + szHeight), ...
                (floor((szWidth + 1) / 2) + 1):(floor((szWidth + 1) / 2) + szWidth), d) = ...
                spectrumSep(i, :, :, d) .* 4; % x4 because of the padding
        end

    end

    % Shift zero frequency back to corner
    for d = 1:numDirection

        for i = 1:3
            spectrumSep2x(i, :, :, d) = ifftshift(spectrumSep2x(i, :, :, d));
        end

    end

    % Plot the spectrum and roughly locate the vectors
    f = figure( ...
        'MenuBar', 'none', ...
        'ToolBar', 'none', ...
        'Name', ['SIM parameters for ', filename_raw], ...
        'Units', 'normalized', ...
        'Position', [0.1 0.1 0.8 0.8]);
    tiledlayout(f, 3, 4);
    spectrumCrossCorr = zeros(szHeight * 2, szWidth * 2, numDirection);
    disp('Locating vectors...');

    for d = 1:numDirection %#ok<FXUP>
        spectrum1 = squeeze(spectrumSep2x(1, :, :, d)) .* OTFMaskx2 .* notchFilterx2;
        spectrum0 = squeeze(spectrumSep2x(2, :, :, d)) .* OTFMaskx2 .* notchFilterx2;
        spectrum2 = squeeze(spectrumSep2x(3, :, :, d)) .* OTFMaskx2 .* notchFilterx2;

        % Cross correlation of three frequency components
        spectrumCrossCorr1 = abs(fft2(conj(ifft2(spectrum1)) .* ifft2(spectrum0)));
        spectrumCrossCorr2 = abs(fft2(conj(ifft2(spectrum2)) .* ifft2(spectrum0)));
        spectrumCrossCorr(:, :, d) = spectrumCrossCorr1 + spectrumCrossCorr2;

        nexttile((d - 1) * 4 + 1);
        pbg = imagesc( ...
            [-szWidth, szWidth - 1], ...
            [-szHeight, szHeight - 1], ...
            fftshift(spectrumCrossCorr(:, :, d) .* OTFRingMaskx2));
        colormap(hot);
        axis image off;
        pbg.Parent.CLim = [0, 0.3 * max(pbg.CData, [], "all")];
        xlim([-floor(szWidth / 2), floor(szWidth / 2) - 1]);
        ylim([-floor(szHeight / 2), floor(szHeight / 2) - 1]);
        title(['Direction #', num2str(d)]);
        hold on;

        spectrumLoc = fftshift(spectrumCrossCorr1 .* OTFRingMaskx2);
        [~, indMax] = max(spectrumLoc(:));
        [indH, indW] = ind2sub(size(spectrumLoc), indMax);
        indW = indW - szWidth - 1;
        indH = indH - szHeight - 1;
        plot(indW, indH, 'go');
        paraIllVector(d, :) = [indW, indH];

        spectrumLoc = fftshift(spectrumCrossCorr2 .* OTFRingMaskx2);
        [~, indMax] = max(spectrumLoc(:));
        [indH, indW] = ind2sub(size(spectrumLoc), indMax);
        indW = indW - szWidth - 1;
        indH = indH - szHeight - 1;
        plot(indW, indH, 'ro');

        paraIllVector(d, :) = (paraIllVector(d, :) - [indW, indH]) / 2;
    end

    % Plot the zoomed-in spectrum and locate the vectors
    for d = 1:numDirection %#ok<FXUP>
        nexttile((d - 1) * 4 + 2);
        imagesc( ...
            [-szWidth, szWidth - 1], ...
            [-szHeight, szHeight - 1], ...
            fftshift(spectrumCrossCorr(:, :, d) .* OTFRingMaskx2));
        axis image;
        xlim([-2 2] + paraIllVector(d, 1))
        ylim([-2 2] + paraIllVector(d, 2))
        title(['Direction #', num2str(d), ' Zoom In']);
        hold on;
        plot(paraIllVector(d, 1), paraIllVector(d, 2), 'go');
    end

    drawnow;

    %% Iteratively solve the vector precisely
    for d = 1:numDirection %#ok<FXUP>
        disp(['Vector iteration for direction #', num2str(d), ' ...']);
        iterGradStep = 1e-4;
        iterBeta = 0.8;
        vector = paraIllVector(d, :);
        corValw1 = getCorrVal(spectrumSep2x(:, :, :, d), ...
            vector - [0.5 * iterGradStep, 0]);
        corValh1 = getCorrVal(spectrumSep2x(:, :, :, d), ...
            vector - [0, 0.5 * iterGradStep]);
        corValw2 = getCorrVal(spectrumSep2x(:, :, :, d), ...
            vector + [0.5 * iterGradStep, 0]);
        corValh2 = getCorrVal(spectrumSep2x(:, :, :, d), ...
            vector + [0, 0.5 * iterGradStep]);
        corGrad = [corValw2 - corValw1, corValh2 - corValh1] ...
            / iterGradStep;
        iterStep = 0.1 / norm(corGrad);
        corStep = iterStep * corGrad;

        nexttile((d - 1) * 4 + 2);
        p = plot(vector(1), vector(2), 'rx');
        drawnow;

        while norm(corStep) > iterGradStep
            vector = vector + corStep;
            corValw1 = getCorrVal(spectrumSep2x(:, :, :, d), ...
                vector - [0.5 * iterGradStep, 0]);
            corValh1 = getCorrVal(spectrumSep2x(:, :, :, d), ...
                vector - [0, 0.5 * iterGradStep]);
            corValw2 = getCorrVal(spectrumSep2x(:, :, :, d), ...
                vector + [0.5 * iterGradStep, 0]);
            corValh2 = getCorrVal(spectrumSep2x(:, :, :, d), ...
                vector + [0, 0.5 * iterGradStep]);
            corGrad = [corValw2 - corValw1, corValh2 - corValh1] ...
                / iterGradStep;
            corStep = iterBeta * iterStep * corGrad + ...
                (1 - iterBeta) * corStep;

            p.XData = [p.XData, vector(1)];
            p.YData = [p.YData, vector(2)];
            xlabel(['k_x=', num2str(vector(1))]);
            ylabel(['k_y=', num2str(vector(2))]);
            drawnow;
            % disp(norm(corStep));
        end

        vector = vector + corStep;
        paraIllVector(d, :) = vector;
    end

    %% Solve the initial phase and modulation depth
    for d = 1:numDirection %#ok<FXUP>
        disp(['Solving initial phase and modulation depth for direction #', num2str(d), ' ...']);
        [gamma, phi] = getParameters(spectrumSep2x(:, :, :, d), ...
            paraIllVector(d, :), d);
        paraPhaseVectorFull(:, d) = wrapToPi( ...
            paraPhaseVectorFull(:, d) - phi);
        paraIMatrix(:, d) = [2 / gamma, 1, 2 / gamma];

        for p = 1:numPhase %#ok<FXUP>
            paraPhaseMatrixFull(p, :, d) = ...
                paraPhaseMatrixFull(p, :, d) .* exp(1i * ...
                [-paraPhaseVectorFull(p, d), ...
                 0, ...
                 paraPhaseVectorFull(p, d)]);
        end

    end

    %% Save a report of parameters
    exportgraphics(f, [pathname_raw, filesep, filename_raw, '.pdf'], 'ContentType', 'vector');

    %% Separation of three frequency components in each direction with precise vectors
    spectrumSep = complex(zeros(3, szHeight, szWidth, numDirection));
    spectrumSep2x = complex(zeros(3, szHeight * 2, szWidth * 2, ...
        numDirection, 'single'));

    for d = 1:numDirection %#ok<FXUP>
        spectrumSepEachP = complex(zeros(numPhase, szHeight * szWidth));

        for p = 1:numPhase %#ok<FXUP>
            spectrumTemp = spectrumRaw(:, :, p + (d - 1) * numPhase);
            spectrumSepEachP(p, :) = spectrumTemp(:);
        end

        spectrumSepEach = squeeze(paraPhaseMatrixFull(:, :, d)) ...
            \ spectrumSepEachP;
        spectrumSep((((d - 1) * 3 * szHeight * szWidth) + 1) ...
            :(d * 3 * szHeight * szWidth)) = spectrumSepEach(:);

        for i = 1:3
            spectrumSep(i, :, :, d) = paraIMatrix(i, d) .* ...
                spectrumSep(i, :, :, d) .* numPhase;
            % x numPhase because of the speration average over numPhase images
        end

    end

    % Padding the spectrum 2x
    for d = 1:numDirection

        for i = 1:3
            spectrumSep2x(i, (floor((szHeight + 1) / 2) + 1):(floor((szHeight + 1) / 2) + szHeight), ...
                (floor((szWidth + 1) / 2) + 1):(floor((szWidth + 1) / 2) + szWidth), d) = ...
                spectrumSep(i, :, :, d) .* 4; % x4 because of the 2x padding
        end

    end

    % Shift zero frequency back to corner
    for d = 1:numDirection

        for i = 1:3
            spectrumSep2x(i, :, :, d) = ifftshift(spectrumSep2x(i, :, :, d));
        end

    end

    %% Wiener inverse filtering
    spectrumSum = complex(zeros(szHeight * 2, szWidth * 2));
    otfSum = complex(zeros(szHeight * 2, szWidth * 2));
    maskSum = zeros(szHeight * 2, szWidth * 2);

    for d = 1:numDirection %#ok<FXUP>
        spectrumSum = spectrumSum + ...
            specCombine(spectrumSep2x(:, :, :, d), paraIllVector(d, :));
        [otftemp, masktemp] = otfCombine(paraIllVector(d, :), paraIMatrix(:, d));
        otfSum = otfSum + otftemp;
        maskSum = maskSum + masktemp;
    end

    % Normalize the combined spectrum
    % this step is essential to keep the energy conservation during the deconvolution
    % maskSum(maskSum == 0) = 1; % avoid division by zero
    % spectrumSum = spectrumSum ./ maskSum .* numPhase .* numDirection;
    % otfSum = otfSum ./ maskSum .* numPhase;

    combinedImg = real(fftshift(ifft2(spectrumSum)));
    combinedOTF = otfSum ./ numDirection;
    combinedPSF = real(fftshift(ifft2(combinedOTF)));

    %% Preview the result and fine-tune the Wiener parameter
    close(f);
    uiResult = uifigure('Name', 'SIM Result', ...
        'Units', 'normalized', ...
        'Position', [0.1, 0.1, 0.8, 0.8]);
    uiLayout = uigridlayout(uiResult, [2, 4]);
    uiLayout.RowHeight = {'1x', '1x', 50, 20};
    uiLayout.ColumnWidth = {'1x', '1x'};
    uiAxes1 = uiaxes(uiLayout);
    uiAxes1.Layout.Row = 1;
    uiAxes1.Layout.Column = 1;
    uiAxes2 = uiaxes(uiLayout);
    uiAxes2.Layout.Row = 2;
    uiAxes2.Layout.Column = 1;
    uiAxes3 = uiaxes(uiLayout);
    uiAxes3.Layout.Row = 1;
    uiAxes3.Layout.Column = 2;
    uiAxes4 = uiaxes(uiLayout);
    uiAxes4.Layout.Row = 2;
    uiAxes4.Layout.Column = 2;

    imagesc(uiAxes1, imgWF);
    title(uiAxes1, 'Wide Field');
    axis(uiAxes1, 'image');
    axis(uiAxes1, 'off');
    colormap(uiAxes1, gray);
    clim(uiAxes1, [0, 0.90 * max(imgWF(:))]);
    colorbar(uiAxes1);

    imagesc(uiAxes2, imgWFdeconv);
    title(uiAxes2, 'Wide Field Deconvolved');
    axis(uiAxes2, 'image');
    axis(uiAxes2, 'off');
    colormap(uiAxes2, gray);
    clim(uiAxes2, [0, 0.90 * max(imgWF(:))]);
    colorbar(uiAxes2);

    imagesc(uiAxes3, imgOS);
    title(uiAxes3, 'Optical Sectioned');
    axis(uiAxes3, 'image');
    axis(uiAxes3, 'off');
    colormap(uiAxes3, gray);
    clim(uiAxes3, [0, 0.90 * max(imgWF(:))]);
    colorbar(uiAxes3);

    simResult = mydeconvlucy(combinedImg, combinedOTF, paraLRIter);
    rimg = imagesc(uiAxes4, simResult);
    title(uiAxes4, 'SIM');
    axis(uiAxes4, 'image');
    axis(uiAxes4, 'off');
    colormap(uiAxes4, gray);
    clim(uiAxes4, [0, 0.90 * max(imgWF(:))]);
    colorbar(uiAxes4);

    uiIterSlider = uislider(uiLayout, 'Limits', [1, 20], 'Value', 10, 'ValueChangedFcn', @drawResult);
    uiIterSlider.Layout.Row = 3;
    uiIterSlider.Layout.Column = 2;
    uiProcessButton = uibutton(uiLayout, 'push', ...
        'Text', 'Start Batch Process', ...
        'ButtonPushedFcn', @batchProcess);
    uiProcessButton.Layout.Row = 4;
    uiProcessButton.Layout.Column = 2;

    %% Batch processing for all frames
    function batchProcess(~, ~)
        disp('Start batch processing ...');
        paraLRIter = round(uiIterSlider.Value);
        uiProcessButton.Enable = 'off';
        uiIterSlider.Enable = 'off';
        numStackStart = mod(numStackStart - 1, 9) + 1;
        numStackGroup = (szStack - numStackStart + 1) / (numPhase * numDirection);
        uiProgress = uiprogressdlg(uiResult, 'Title', 'Batch processing...', ...
            'Message', ['0 / ', num2str(numStackGroup)]);

        % 打开文件
        tiff_outDescription = {'ImageJ', ...
                                   'unit=\u00B5m', ...
                                   'loop=false', ...
                                   'min=0', ...
                                   ['max=', num2str(round(0.9 * max(imgOS(:))))]};
        warning off;
        tiff_raw = Tiff([pathname_raw, filename_raw], 'r');
        tiff_out = Tiff([pathname_raw, 'SIM_', filename_raw], 'w');
        tiff_outTag = struct('ImageWidth', szWidth * 2, ...
            'ImageLength', szHeight * 2, ...
            'Photometric', Tiff.Photometric.MinIsBlack, ...
            'BitsPerSample', 32, ...
            'SampleFormat', Tiff.SampleFormat.IEEEFP, ...
            'SamplesPerPixel', 1, ...
            'PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
        tiff_out.setTag('ImageDescription', sprintf('%s\n', tiff_outDescription{:}));
        tiff_out.setTag('XResolution', 2 / paraPixelSz);
        tiff_out.setTag('YResolution', 2 / paraPixelSz);
        tiff_out.setTag('ResolutionUnit', 1);

        tiff_WF = Tiff([pathname_raw, 'WF_', filename_raw], 'w');
        tiff_WFTag = struct('ImageWidth', szWidth, ...
            'ImageLength', szHeight, ...
            'Photometric', Tiff.Photometric.MinIsBlack, ...
            'BitsPerSample', 32, ...
            'SampleFormat', Tiff.SampleFormat.IEEEFP, ...
            'SamplesPerPixel', 1, ...
            'PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
        tiff_WF.setTag('ImageDescription', sprintf('%s\n', tiff_outDescription{:}));
        tiff_WF.setTag('XResolution', 1 / paraPixelSz);
        tiff_WF.setTag('YResolution', 1 / paraPixelSz);
        tiff_WF.setTag('ResolutionUnit', 1);

        tiff_WFdeconv = Tiff([pathname_raw, 'WFdeconv_', filename_raw], 'w');
        tiff_WFdeconvTag = struct('ImageWidth', szWidth, ...
            'ImageLength', szHeight, ...
            'Photometric', Tiff.Photometric.MinIsBlack, ...
            'BitsPerSample', 32, ...
            'SampleFormat', Tiff.SampleFormat.IEEEFP, ...
            'SamplesPerPixel', 1, ...
            'PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
        tiff_WFdeconv.setTag('ImageDescription', sprintf('%s\n', tiff_outDescription{:}));
        tiff_WFdeconv.setTag('XResolution', 1 / paraPixelSz);
        tiff_WFdeconv.setTag('YResolution', 1 / paraPixelSz);
        tiff_WFdeconv.setTag('ResolutionUnit', 1);

        tiff_OSBG = Tiff([pathname_raw, 'OSBG_', filename_raw], 'w');
        tiff_OSBGTag = struct('ImageWidth', szWidth, ...
            'ImageLength', szHeight, ...
            'Photometric', Tiff.Photometric.MinIsBlack, ...
            'BitsPerSample', 32, ...
            'SampleFormat', Tiff.SampleFormat.IEEEFP, ...
            'SamplesPerPixel', 1, ...
            'PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
        tiff_OSBG.setTag('ImageDescription', sprintf('%s\n', tiff_outDescription{:}));
        tiff_OSBG.setTag('XResolution', 1 / paraPixelSz);
        tiff_OSBG.setTag('YResolution', 1 / paraPixelSz);
        tiff_OSBG.setTag('ResolutionUnit', 1);

        tiff_OS = Tiff([pathname_raw, 'OS_', filename_raw], 'w');
        tiff_OSTag = struct('ImageWidth', szWidth, ...
            'ImageLength', szHeight, ...
            'Photometric', Tiff.Photometric.MinIsBlack, ...
            'BitsPerSample', 32, ...
            'SampleFormat', Tiff.SampleFormat.IEEEFP, ...
            'SamplesPerPixel', 1, ...
            'PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
        tiff_OS.setTag('ImageDescription', sprintf('%s\n', tiff_outDescription{:}));
        tiff_OS.setTag('XResolution', 1 / paraPixelSz);
        tiff_OS.setTag('YResolution', 1 / paraPixelSz);
        tiff_OS.setTag('ResolutionUnit', 1);
        warning on;

        for g = 1:numStackGroup

            for f = 1:(numPhase * numDirection) %#ok<FXUP>
                warning off;
                tiff_raw.setDirectory(numStackStart + (g - 1) * (numPhase * numDirection) + (f - 1));
                imgRaw(:, :, f) = single(tiff_raw.read());
                warning on;
            end

            imgWF = sum(imgRaw, 3);

            %% Lucy-Richardson Deconvolution over Raw Images
            % Default iteration number is 10
            for f = 1:(numPhase * numDirection) %#ok<FXUP>
                imgRaw(:, :, f) = mydeconvlucy(imgRaw(:, :, f), OTF, 10);
            end

            imgWFdeconv = sum(imgRaw, 3);

            %% Background substraction with optical sectioned SIM
            osImg(:) = 0;
            osBG(:) = 0;

            for d = 1:numDirection

                for p = 1:numPhase

                    for p2 = 1:numPhase

                        if p2 ~= p
                            osImg(:, :, d) = osImg(:, :, d) + (imgRaw(:, :, p + (d - 1) * numPhase) - imgRaw(:, :, p2 + (d - 1) * numPhase)) .^ 2;
                        end

                    end

                    osBG(:, :, d) = osBG(:, :, d) + imgRaw(:, :, p + (d - 1) * numPhase) / numPhase;
                end

                osImg(:, :, d) = sqrt(osImg(:, :, d) / (numPhase * (numPhase - 1)));
                osBG(:, :, d) = osBG(:, :, d) - osImg(:, :, d);

            end

            imgBG = sum(osBG, 3) .* numPhase;
            imgOS = imgWFdeconv - imgBG;

            for f = 1:(numPhase * numDirection) %#ok<FXUP>
                spectrumRaw(:, :, f) = fftshift(fft2(ifftshift(imgRaw(:, :, f))));
            end

            for d = 1:numDirection %#ok<FXUP>
                spectrumSepEachP = complex(zeros(numPhase, szHeight * szWidth));

                for p = 1:numPhase %#ok<FXUP>
                    spectrumTemp = spectrumRaw(:, :, p + (d - 1) * numPhase);
                    spectrumSepEachP(p, :) = spectrumTemp(:);
                end

                spectrumSepEach = squeeze(paraPhaseMatrixFull(:, :, d)) ...
                    \ spectrumSepEachP;
                spectrumSep((((d - 1) * 3 * szHeight * szWidth) + 1) ...
                    :(d * 3 * szHeight * szWidth)) = spectrumSepEach(:);

                for r = 1:3
                    spectrumSep(r, :, :, d) = paraIMatrix(r, d) .* ...
                        spectrumSep(r, :, :, d) .* numPhase;
                    % x numPhase because of the speration average over numPhase images
                end

            end

            spectrumSep2x(:) = 0;

            for d = 1:numDirection

                for i = 1:3
                    spectrumSep2x(i, (floor((szHeight + 1) / 2) + 1):(floor((szHeight + 1) / 2) + szHeight), ...
                        (floor((szWidth + 1) / 2) + 1):(floor((szWidth + 1) / 2) + szWidth), d) = ...
                        spectrumSep(i, :, :, d) .* 4;
                    % x4 because of the padding
                end

            end

            % Shift zero frequency back to corner
            for d = 1:numDirection

                for i = 1:3
                    spectrumSep2x(i, :, :, d) = ifftshift(spectrumSep2x(i, :, :, d));
                end

            end

            spectrumSum = complex(zeros(szHeight * 2, szWidth * 2));

            for d = 1:numDirection %#ok<FXUP>
                spectrumSum = spectrumSum + ...
                    specCombine(spectrumSep2x(:, :, :, d), paraIllVector(d, :));
            end

            simResult = mydeconvlucy(real(fftshift(ifft2(spectrumSum))), combinedOTF, paraLRIter);

            tiff_WF.setTag(tiff_WFTag);
            tiff_WF.write(single(imgWF));
            tiff_WF.writeDirectory();

            tiff_WFdeconv.setTag(tiff_WFdeconvTag);
            tiff_WFdeconv.write(single(imgWFdeconv));
            tiff_WFdeconv.writeDirectory();

            tiff_OSBG.setTag(tiff_OSBGTag);
            tiff_OSBG.write(single(imgBG));
            tiff_OSBG.writeDirectory();

            tiff_OS.setTag(tiff_OSTag);
            tiff_OS.write(single(imgOS));
            tiff_OS.writeDirectory();

            tiff_out.setTag(tiff_outTag);
            tiff_out.write(single(simResult));
            tiff_out.writeDirectory();

            uiProgress.Message = [num2str(g), ' / ', num2str(numStackGroup)];
        end

        tiff_raw.close();
        tiff_out.close();
        tiff_WF.close();
    end

    %% Callback function for tuning the Wiener parameter
    function drawResult(src, ~)
        paraLRIter = round(src.Value);
        simResult = mydeconvlucy(combinedImg, combinedOTF, paraLRIter);
        rimg.CData = simResult;
        clim(uiAxes5, [0, 0.90 * max(imgOS(:))]);
    end

    %% Get correlation value at specific wave vector
    function corrValve = getCorrVal(spectrumSep3, shiftVector)
        shiftMask = exp(1i * 2 * pi * ...
            (shiftVector(1) * meshWx2 + shiftVector(2) * meshHx2));
        mask0 = OTFMaskx2;
        mask1 = real(fft2(ifftshift(fftshift(ifft2(OTFMaskx2)) ...
            .* shiftMask)));
        mask2 = real(fft2(ifftshift(fftshift(ifft2(OTFMaskx2)) ...
            .* conj(shiftMask))));
        mask1(mask1 < 0.5) = 0;
        mask2(mask2 < 0.5) = 0;
        mask1(mask1 ~= 0) = 1;
        mask2(mask2 ~= 0) = 1;
        % OTF1 = real(fft2(ifftshift(fftshift(ifft2(OTFx2)) ...
        %     .* shiftMask)));
        % OTF2 = real(fft2(ifftshift(fftshift(ifft2(OTFx2)) ...
        %     .* conj(shiftMask))));
        spec0 = squeeze(spectrumSep3(2, :, :));
        spec1 = fft2(ifftshift(fftshift(ifft2( ...
            squeeze(spectrumSep3(1, :, :)))) .* shiftMask));
        spec2 = fft2(ifftshift(fftshift(ifft2( ...
            squeeze(spectrumSep3(3, :, :)))) .* conj(shiftMask)));
        maskCom0 = mask1 .* mask2;
        maskCom1 = mask0 .* mask1 .* (~maskCom0);
        maskCom2 = mask0 .* mask2 .* (~maskCom0);
        % specCom1 = maskCom1 .* OTFx2 .* OTF1 .* ...
        %     spec0 .* conj(spec1);
        % specCom2 = maskCom2 .* OTFx2 .* OTF2 .* ...
        %     spec0 .* conj(spec2);
        specCom1 = maskCom1 .* spec0 .* conj(spec1);
        specCom2 = maskCom2 .* spec0 .* conj(spec2);
        corrValve = (abs(sum(specCom1(:))) / sum(maskCom1(:)) + ...
            abs(sum(specCom2(:))) / sum(maskCom2(:))) / 2;
    end

    %% Joint the three frequency components at specific wave vector
    function spectrumNew = specCombine(spectrumSep3, shiftVector)
        shiftMask = exp(1i * 2 * pi * ...
            (shiftVector(1) * meshWx2 + shiftVector(2) * meshHx2));
        mask0 = OTFMaskx2;
        mask1val = real(fft2(ifftshift(fftshift(ifft2(OTFMaskx2)) ...
            .* shiftMask)));
        mask2val = real(fft2(ifftshift(fftshift(ifft2(OTFMaskx2)) ...
            .* conj(shiftMask))));
        mask1 = (mask1val > 0.5);
        mask2 = (mask2val > 0.5);
        spec0 = squeeze(spectrumSep3(2, :, :));
        spec1 = fft2(ifftshift(fftshift(ifft2( ...
            squeeze(spectrumSep3(1, :, :)))) .* shiftMask));
        spec2 = fft2(ifftshift(fftshift(ifft2( ...
            squeeze(spectrumSep3(3, :, :)))) .* conj(shiftMask)));
        attFun2x1 = fft2(ifftshift(fftshift(ifft2( ...
            attFun2x)) .* shiftMask));
        attFun2x2 = fft2(ifftshift(fftshift(ifft2( ...
            attFun2x)) .* conj(shiftMask)));
        sum0 = mask0 .* attFun2x .* spec0;
        sum1 = mask1 .* attFun2x1 .* spec1;
        sum2 = mask2 .* attFun2x2 .* spec2;
        spectrumNew = sum0 + sum1 +sum2;
    end

    %% Joint the OTF at specific wave vector
    function [otfNew, maskNew] = otfCombine(shiftVector, depthVector)
        shiftMask = exp(1i * 2 * pi * ...
            (shiftVector(1) * meshWx2 + shiftVector(2) * meshHx2));
        mask0 = OTFMaskx2;
        mask1val = real(fft2(ifftshift(fftshift(ifft2(OTFMaskx2)) ...
            .* shiftMask)));
        mask2val = real(fft2(ifftshift(fftshift(ifft2(OTFMaskx2)) ...
            .* conj(shiftMask))));
        mask1 = (mask1val > 0.5);
        mask2 = (mask2val > 0.5);
        OTF1 = real(fft2(ifftshift(fftshift(ifft2(OTFx2)) ...
            .* shiftMask)));
        OTF2 = real(fft2(ifftshift(fftshift(ifft2(OTFx2)) ...
            .* conj(shiftMask))));
        % attFun2x1 = fft2(ifftshift(fftshift(ifft2( ...
        %     attFun2x)) .* shiftMask));
        % attFun2x2 = fft2(ifftshift(fftshift(ifft2( ...
        %     attFun2x)) .* conj(shiftMask)));
        sum0 = mask0 .* OTFx2;
        sum1 = mask1 .* OTF1;
        sum2 = mask2 .* OTF2;
        otfNew = (1 / depthVector(1)) .* sum1 + (1 / depthVector(2)) .* sum0 + (1 / depthVector(3)) .* sum2;
        maskNew = mask0 + mask1 + mask2;
    end

    %% Solve the initial phase and modulation depth at specific wave vector
    function [gamma, phi] = getParameters(spectrumSep3, shiftVector, ind)
        shiftMask = exp(1i * 2 * pi * ...
            (shiftVector(1) * meshWx2 + shiftVector(2) * meshHx2));
        mask0 = OTFMaskx2;
        mask1 = real(fft2(ifftshift(fftshift(ifft2(OTFMaskx2)) ...
            .* shiftMask)));
        mask2 = real(fft2(ifftshift(fftshift(ifft2(OTFMaskx2)) ...
            .* conj(shiftMask))));
        mask1(mask1 < 0.5) = 0;
        mask2(mask2 < 0.5) = 0;
        mask1(mask1 ~= 0) = 1;
        mask2(mask2 ~= 0) = 1;
        OTF1 = real(fft2(ifftshift(fftshift(ifft2(OTFx2)) ...
            .* shiftMask)));
        OTF2 = real(fft2(ifftshift(fftshift(ifft2(OTFx2)) ...
            .* conj(shiftMask))));
        spec0 = squeeze(spectrumSep3(2, :, :));
        spec1 = fft2(ifftshift(fftshift(ifft2( ...
            squeeze(spectrumSep3(1, :, :)))) .* shiftMask));
        spec2 = fft2(ifftshift(fftshift(ifft2( ...
            squeeze(spectrumSep3(3, :, :)))) .* conj(shiftMask)));
        maskCom0 = mask1 .* mask2;
        maskcom1 = mask0 .* mask1 .* (~maskCom0) .* (OTF1 > 0.1) .* (OTFx2 > 0.1);
        maskcom2 = mask0 .* mask2 .* (~maskCom0) .* (OTF2 > 0.1) .* (OTFx2 > 0.1);

        % method from fairSIM
        test1 = maskcom1 .* conj(spec0) .* spec1;
        test10 = maskcom1 .* conj(spec0) .* spec0;
        test2 = maskcom2 .* conj(spec0) .* spec2;
        test20 = maskcom2 .* conj(spec0) .* spec0;
        gamma1 = abs(sum(test1(:)) ./ sum(test10(:)));
        gamma2 = abs(sum(test2(:)) ./ sum(test20(:)));
        gamma = (gamma1 + gamma2) ./ 2;

        phi1 = angle(sum(test1(:)) ./ sum(test10(:)));
        phi2 = angle(sum(test2(:)) ./ sum(test20(:)));
        phi = sign(phi1) * (abs(phi1) + abs(phi1)) / 2;

        disp(gamma);
        disp(rad2deg(phi));

    end

    %% Fitting an ellipse in polar coordinates
    function [alpha_rad, ellipse] = getElipseAlpha(Theta, Rho)
        [x, y] = pol2cart(Theta(:), Rho(:));
        mean_x = mean(x(:));
        mean_y = mean(y(:));
        x = x - mean_x;
        y = y - mean_y;
        X = [x .^ 2, x .* y, y .^ 2, x, y];
        cm = sum(X) / (X' * X);
        [A, B, C, D, E] = deal(cm(1), cm(2), cm(3), cm(4), cm(5));
        orientation_rad = 1/2 * atan(B / (C - A));
        cos_phi = cos(orientation_rad);
        sin_phi = sin(orientation_rad);
        [A, B, C, D, E] = deal( ...
            A * cos_phi ^ 2 - B * cos_phi * sin_phi + C * sin_phi ^ 2, ...
            0, ...
            A * sin_phi ^ 2 + B * cos_phi * sin_phi + C * cos_phi ^ 2, ...
            D * cos_phi - E * sin_phi, ...
            D * sin_phi + E * cos_phi); %#ok<ASGLU>
        [mean_x, mean_y] = deal( ...
            cos_phi * mean_x - sin_phi * mean_y, ...
            sin_phi * mean_x + cos_phi * mean_y);

        X0 = mean_x - D / 2 / A;
        Y0 = mean_y - E / 2 / C;
        F = 1 + (D ^ 2) / (4 * A) + (E ^ 2) / (4 * C);
        [A, B] = deal(sqrt(F / A), sqrt(F / C));

        R = [cos_phi sin_phi; -sin_phi cos_phi];

        ellipse_x_r = X0 + A * cos(Theta);
        ellipse_y_r = Y0 + B * sin(Theta);
        rotated_ellipse = R * [ellipse_x_r; ellipse_y_r];
        [ellipse(1, :), ellipse(2, :)] = cart2pol( ...
            rotated_ellipse(1, :), rotated_ellipse(2, :));
        alpha_rad = atan2(-sin_phi * X0 + cos_phi * Y0, ...
            cos_phi * X0 + sin_phi * Y0);
    end

end

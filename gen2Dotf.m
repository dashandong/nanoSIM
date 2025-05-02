function gen2Dotf()
    clc;

    %% -------------------- Parameters --------------------
    % Image size
    szWidth = 512;
    szHeight = 512;

    % Optical system parameters
    paraMag = 100; % Objective magnification
    paraImmersion = 1.518; % Refractive index of immersion medium
    paraNA = 1.45; % Numerical aperture
    paraWavelength = 0.520; % Emission wavelength (µm)
    paraPixelSz = 6.5 / paraMag; % Pixel size in object space (µm)

    % Z offsets in object space (µm)
    deltaZ = [0]; % For 2D OTF; extend this vector for 3D OTFs

    %% -------------------- Frequency Domain Setup --------------------
    % Wave number (k = 2πn / λ)
    paraK = 2 * pi * paraImmersion / paraWavelength;
    paraKm = 2 * pi * paraNA / paraWavelength;

    % Frequency sampling intervals
    paraKWSz = 2 * pi / (paraPixelSz * szWidth);
    paraKHSz = 2 * pi / (paraPixelSz * szHeight);

    % Spatial frequency axes (centered)
    axisKW = ((1:szWidth) - 1 - floor(szWidth / 2)) * paraKWSz;
    axisKH = ((1:szHeight) - 1 - floor(szHeight / 2)) * paraKHSz;
    [meshKW, meshKH] = meshgrid(axisKW, axisKH);
    [~, meshKR] = cart2pol(meshKW, meshKH); % Radial frequency magnitude

    %% -------------------- Initialize OTF Volume --------------------
    otf = complex(zeros(szWidth, szHeight, length(deltaZ)));

    % Calculate pupil amplitude normalization
    pupilArea = pi * (paraKm / (2 * pi)) ^ 2;
    powerTotal = 1;
    ampInput = sqrt(powerTotal / pupilArea); % Uniform input amplitude

    %% -------------------- Compute OTF for Each Z --------------------
    for currZ = 1:length(deltaZ)
        % Coherent pupil function with defocus phase
        pupil = ampInput .* exp(1i .* sqrt(paraK ^ 2 - meshKR .^ 2) .* deltaZ(currZ));
        pupil(meshKR > paraKm) = 0; % Apply aperture cutoff

        % Amplitude PSF from inverse FFT of pupil
        apsf = fftshift(ifft2(ifftshift(pupil))) * ...
            szWidth * szHeight * paraKWSz ^ 2 / (2 * pi) ^ 2;

        % Intensity PSF
        ipsf = abs(apsf) .^ 2;

        % Compute OTF via FFT of the IPSF
        otf(:, :, currZ) = fftshift(fft2(ifftshift(ipsf))) * paraPixelSz ^ 2;

        % Normalize OTF (so DC value = 1 if total input power is 1)
        otf(:, :, currZ) = otf(:, :, currZ) / powerTotal;
    end

    %% -------------------- Final OTF Post-Processing --------------------
    finalotf = real(otf); % Discard any imaginary residue

    % Hard cutoff: zero out values beyond 2× cutoff frequency
    for currZ = 1:length(deltaZ)
        tempotf = finalotf(:, :, currZ);
        tempotf(meshKR > 2 * paraKm) = 0;
        finalotf(:, :, currZ) = tempotf;
    end

    %% -------------------- Save as TIFF Stack --------------------
    % Prepare metadata for ImageJ
    tiff_outDescription = {
                           'ImageJ', ...
                               'unit=µm⁻¹', ...
                               'loop=false', ...
                               'min=0', ...
                           'max=1'};

    % Create TIFF file
    tiff_out = Tiff(['OTF_NA', num2str(paraNA, '%.2f'), '.tif'], 'w');
    warning off;

    % Set constant TIFF tag values
    tiff_outTag = struct('ImageWidth', szWidth, ...
        'ImageLength', szHeight, ...
        'Photometric', Tiff.Photometric.MinIsBlack, ...
        'BitsPerSample', 32, ...
        'SampleFormat', Tiff.SampleFormat.IEEEFP, ...
        'SamplesPerPixel', 1, ...
        'PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);

    tiff_out.setTag('ImageDescription', sprintf('%s\n', tiff_outDescription{:}));
    tiff_out.setTag('XResolution', 2 / paraKWSz); % For µm⁻¹ units
    tiff_out.setTag('YResolution', 2 / paraKHSz);
    tiff_out.setTag('ResolutionUnit', 1); % Unit = None (custom units)

    % Write each Z-slice to TIFF
    for currZ = 1:length(deltaZ)
        tiff_out.setTag(tiff_outTag);
        tiff_out.write(single(finalotf(:, :, currZ)));
        tiff_out.writeDirectory();
    end

    tiff_out.close();
    warning on;
end

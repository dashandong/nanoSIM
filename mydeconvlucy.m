function J = mydeconvlucy(I, OTF, NUMIT)
    %MYDECONVLUCY Custom Lucy-Richardson deconvolution using a frequency-domain OTF.
    %   modified from MATLAB's deconvlucy function.
    %   modify lambda to be in [-1, 1], disable non-negative constraint
    %
    %   J = mydeconvlucy(I, OTF, NUMIT) deconvolves image I using
    %   optical transfer function (OTF) and runs for NUMIT iterations.

    narginchk(3, 3);

    % Validate input image
    validateattributes(I, {'uint8', 'uint16', 'double', 'int16', 'single'}, ...
        {'real', 'nonempty', 'finite'}, mfilename, 'I', 1);

    % Default iteration count
    defaultNumIter = 10;

    if isempty(NUMIT)
        NUMIT = defaultNumIter;
    else
        validateattributes(NUMIT, {'double'}, ...
            {'scalar', 'positive', 'finite'}, mfilename, 'NUMIT', 3);
    end

    %% Initialization
    sizeI = size(I);
    numNSdim = find(sizeI ~= 1); % Non-singleton dimensions

    % History buffers for accelerated iteration (e.g., momentum)
    J = cell(4, 1);
    J{1} = I; % Original input
    J{2} = J{1}; % Current estimate
    J{3} = 0; % Previous estimate
    J{4}(prod(sizeI), 2) = 0; % Buffer for momentum/acceleration

    % Weighting image (all ones = no spatial weighting)
    WEIGHT = ones(sizeI);
    READOUT = 0; % Constant readout noise (set to 0)

    %% Prepare FFT scaling
    H = OTF; % Use provided OTF in frequency domain

    % Index slice for dimension alignment (in case of singleton dims)
    idx = repmat({':'}, [1, length(sizeI)]);

    for k = numNSdim
        idx{k} = 1:sizeI(k);
    end

    % Weighted input image
    wI = WEIGHT .* (READOUT + J{1});
    J{2} = J{2}(idx{:});

    % Compute scaling term to normalize updates
    scale = real(ifftn(conj(H) .* fftn(WEIGHT(idx{:})))) + sqrt(eps);
    clear WEIGHT;

    %% Iterative Lucy-Richardson Updates
    lambda = 2 * any(J{4}(:) ~= 0); % Momentum flag

    for k = lambda + (1:NUMIT)
        % === Step 1: Estimate with momentum ===
        if k > 2
            % Compute Nesterov-style acceleration (safe and bounded)
            num = J{4}(:, 1)' * J{4}(:, 2);
            denom = J{4}(:, 2)' * J{4}(:, 2) + eps;
            % modify lambda to be in [-1, 1], disable non-negative constraint
            lambda = max(min(num / denom, 1), -1);
        end

        Y = J{2} + lambda * (J{2} - J{3}); % Accelerated estimate

        % === Step 2: Blur the estimate and compute ratio ===
        Y_conv = real(ifftn(fftn(Y) .* H)); % Convolve estimate
        ratio = wI ./ (Y_conv + READOUT + eps); % Observed/predicted

        % === Step 3: Update estimate ===
        correction = real(ifftn(conj(H) .* fftn(ratio))); % Back-projection
        J{3} = J{2}; % Save previous estimate
        J{2} = Y .* correction ./ scale; % Update with positivity constraint

        % Update momentum buffer
        J{4} = [J{2}(:) - Y(:), J{4}(:, 1)];
    end

    %% Final output
    J = J{2};
end

function [X, Stats] = ccstack(method, A, varargin)

p = inputParser;
p.addOptional('ref', 0);
p.addOptional('ccthreshold', 0);
p.addOptional('eps', 1e-6);
p.addOptional('maxiter', 100);
p.addOptional('v', 0);

p.parse(varargin{:});

ref         = p.Results.ref;
ccthreshold = p.Results.ccthreshold;
eps         = p.Results.eps;
maxiter     = p.Results.maxiter;
v           = p.Results.v;

if p.Results.v == 1; fprintf("Stack method: %s\n", method); end

if method == "linear"

    X = mean(A, 2);
    Stats.Ntimelag = size(A, 1);
    Stats.Ntrace = size(A, 2);

elseif method == "selective"

    if length(ref) ~= size(A, 1)
        error("Length of reference is diferent from trace matrix.");
    end

    [X, Stats] = selectivestack(A, ref, ccthreshold, v);

elseif method == "robust"

    [X, Stats] = robuststack(A, eps, maxiter, v);

end
end

function [X, Stats] = selectivestack(A, ref, ccthreshold, v)
%% Process flow:
%   1. compute Pearson's correlation coefficient between each trace and
%   ref.
%   2. thresholding out them and stack

% Nan check
if any(any(isnan(A)))
    warning("Trace contain NaN value. Ignore the trace.");
    nancol = any(isnan(A));
    A = A(:, ~nancol);
end

cclist = zeros(size(A, 2), 1);
for i = 1:size(A, 2)
    temp = corrcoef(A(:, i), ref);
    cclist(i) = temp(1,2);
end
ind = cclist >= ccthreshold;
X = mean(A(:, ind), 2);

Stats.accceptance_ratio = sum(cclist >= ccthreshold) / size(A, 2) * 100;
temp_ref = corrcoef(X, ref);
Stats.cc_to_ref = temp_ref(1,2);

if v == 1
    fprintf("Acceptance ratio: %4.2f [%%]\n", Stats.accceptance_ratio);
    fprintf("CC between stacked and reference: %4.2f\n", Stats.cc_to_ref);
end
end

function [Bnew, Stats] = robuststack(A, eps, maxiter, v)
% Original code is from SeisNoise.jl
% (https://github.com/tclements/SeisNoise.jl/blob/master/src/stacking.jl)

    if any(any(isnan(A)))
        % Nan check
        warning("Traces contain NaN value. Ignore the trace.");
        nancol = any(isnan(A));
        A = A(:, ~nancol);
    elseif any(all(A==0))
        % All zero check
        warning("Traces contain all zero trace. Ignore the trace.");
        zerocol = all(A==0);
        A = A(:, ~zerocol);
    end

    N = size(A,2);
    Bold = median(A, 2, 'omitnan');
    Bold_norm = Bold ./ norm(Bold,2);
    w =zeros(N, 1);
    r = zeros(N, 1);
    d2 = zeros(N, 1);

    % do 2-norm for all columns in A
    for ii = 1:N
        d2(ii) = norm(A(:,ii),2);
    end

    BdotD = sum(A .* Bold_norm, 1);

    for ii = 1:N
        r(ii) = norm(A(:,ii) - (BdotD(ii) .* Bold_norm),2);
        w(ii) = abs(BdotD(ii)) ./ d2(ii) ./ r(ii);
    end

    Bnew = mean(w'.*A,2);
    Bnew_norm = Bnew ./ norm(Bnew,2);

    % check convergence
    epsN = norm(Bnew_norm - Bold_norm,1) / (norm(Bnew_norm,2) * N);
    Bold_norm = Bnew_norm;
    iter = 0;

    while (epsN > eps) && (iter <= maxiter)

        if v==1
            fprintf("Iteration %d\n", iter);
        end

        BdotD = sum(A .* Bold_norm, 1);

        for ii = 1:N
            r(ii) = norm(A(:,ii) - (BdotD(ii) .* Bold_norm),2);
            w(ii) = abs(BdotD(ii)) ./ d2(ii) ./ r(ii);
        end

        Bnew = mean(w'.*A, 2);
        Bnew_norm = Bnew ./ norm(Bnew, 2);
        % check convergence
        epsN = norm(Bnew_norm - Bold_norm,1) / (norm(Bnew_norm,2) * N);
        Bold_norm = Bnew_norm;
        iter = iter + 1;

        Stats.epsN(iter) = epsN;

        if v==1
            fprintf("epsN = %4.4e\n", epsN);
        end
    end

    if iter >= maxiter
        warning("Robutst stack is not converged.");
    end

    Stats.weight = w;
    Stats.iter = iter-1;

end



function r = slope_limiter(r, type)
    switch type
        case 'minmod'
            r = max(0, min(1, r));
        case 'superbee'
            r = max(0, max(min(1, 2*r), min(2, r)));
        case 'vanLeer'
            r = (r + abs(r)) ./ (1 + abs(r));
        case 'MC'
            % MC limiter (monotonized central)
            r = max(0, min(min(2, 2*r), (1+r)/2));
        case 'none'
            % No limiting (may lead to oscillations)
            r = ones(size(r));
        otherwise
            error('Unknown limiter type. Choose "minmod", "superbee", "vanLeer", "MC", or "none".');
    end
end

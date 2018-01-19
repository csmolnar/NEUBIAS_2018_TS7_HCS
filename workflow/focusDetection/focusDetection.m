function focusDetection(options)
% focus detection

fprintf('\n====================\nFocus detection started...\n');

% setting default focus detection method if it does not exist
if ~isfield(options, 'focusDetection') || isempty(options.focusDetection)
    options.focusDetection = 'adaptive';
end

switch options.focusDetection
    case 'adaptive'
        focusmethod = @adaptiveFocus;
    case 'bestplane'
        focusMethod = @bestPlaneFocus;
    otherwise
        error('HCS:focus','No focus method %s is implemented.',options.focusDetecion);
end

fprintf('\nFocus detection DONE.\n');


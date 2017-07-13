function calib_eight_lum(lums)
%%%%%%%%%%
%
% The program calib_eight_lums.m takes one input argument, a vector of eight
% luminances ranging between 0 and 1. Those luminances are normalized, meaning
% that 0 represents 'black' (the lowest luminance that the monitor can produces)
% and 1 represents 'white' (the highest producible lumaince).
%
% For example:
%   calib_eight_lums([0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7])
%
% The program renders a test pattern on the computer monitor -- eight square
% patches of luminance on a uniform (ie., 0.5) background. The program quits
% after 'j' or 'k' is pressed.
%
% v1.0 20150223 Luke Hallum
%
%%%%%%%%%%

% Check input is sensible...
if ((length(lums) ~= 8) | ...
    (max(lums) > 1) | ...
    (min(lums) < 0))
  error('Your input argument seems wacky. See ''help calib_eight_lums''.')
end

%%%
% Experiment parameters.
%%%%%%%%%%
RES_STIMULUS = [1000 1000];  %% Stimulus must be square!

%%%
% Display parameters and calculations you'll probably need.
%%%%%%%%%%
%% NB: the unit here is 'stimulus width'.
[x,y] = meshgrid(linspace(0,1,RES_STIMULUS(2)));
fieldMeanLum = 0.5*ones(size(x));
patternTest = 0*x;
% Inefficient, but easy to read...
for iilum = 1:8
  patternTest(x > ((iilum-1)/8)) = lums(iilum);
end

%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%

mglOpen();
%mglOpen(0,RES_STIMULUS(1),RES_STIMULUS(2),60,32);

%%%%%%%%%%
mglgrey = mglCreateTexture(255.0 * fieldMeanLum);
mglpattern = mglCreateTexture(255.0 * patternTest);
mglBltTexture(mglgrey,[0 0]); mglFlush; % Set both display buffers to grey.
mglBltTexture(mglgrey,[0 0]); mglFlush; %
%%%%%%%%%%

%%%%%%%%%%
%%%%%%%%%%
%%%
%%%
%%%
%%%%%%%%%%
%%%%%%%%%%

%%%%%%%%%%
% Pattern appears after a key press...
%%%%%%%%%%
mglBltTexture(mglgrey,[0 0]); mglFlush();
respKeys = 0 * mglGetKeys();
KEYS_RESPONSE = [39 41];
while (sum(respKeys([KEYS_RESPONSE])) < 1), respKeys = mglGetKeys(); end
%Debug..
disp('start')

%%%%%%%%%%
%%%%%%%%%%
% Render...
%%%%%%%%%%
%%%%%%%%%%
mglBltTexture(mglpattern,[0 0]); mglFlush();
pause(1);
%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%

respKeys = 0 * mglGetKeys();
while (sum(respKeys([KEYS_RESPONSE])) < 1), respKeys = mglGetKeys(); end
mglClose;

return;


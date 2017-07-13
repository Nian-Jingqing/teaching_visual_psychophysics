function [behaviour] = basic_disc(seedRandom,subjectID,oriStandard_deg,conStim,oriTestRange_deg,numRepeats,fullScreen)
%%%%%%%%%%
%
% The program basic_disc.m takes severn input arguments and
% outputs trial-by-trial data in a matrix. 
%
% For example:
% for the full experiment
%   behavior = basic_disc(1,'LH',45,0.05,20,10,1);
%
% for the initial characterization, horizontal grating, 
% three trials per condition
%   behavior = basic_disc(1,'LH',90,0.05,20,3,1);
% 
% Argument 1: A “random seed”, which enables you to run the
% exact same experiment over and over, if desired. For our
% purposes, use a different nonzero integer each time you
% run a new experiment. 
%
% Argument 2: A subject identifier which will be used to
% name automatically saved files. 
%
% Argument 3: The orientation of the standard in degrees.
% You can input any number here, but we’re interested in the
% cardinals (0 and 90) and the obliques (45 and 135).
% Vertical is 0 degrees.
% 
% Argument 4: The contrast of the test stimuli ranging
% between 0 and 1. For our purposes, use a low contrast,
% say, 0.05. 
%
% Argument 5: The range of orientations, in degrees, of the
% test stimuli used. 
%
% Argument 6: The number of trials of each stimulus type.
% For our purposes, we’ll use 3 for piloting, and 25 for the
% main event. 
%
% Argument 7: Flag indicating whether the program should 
% occupy the whole monitor (1), or just a window (0).
%
% The output of basic_disc.m is a N-by-3 matrix, where N is
% the total number of trials that comprised the experiment.
% So, each row of that matrix contains 3 numbers: The first
% number indicates the orientation of the test stimulus on
% that trial.  The second number is either 0 or 1, and
% indicates whether, on that trial, the subject responded
% “clockwise” (1) or “counter-clockwise” (0).  The third
% number is either 0 or 1, and indicates whether, on that
% trial, the subject’s response was correct (1) or incorrect
% (0).
%
% v1.0 20140130 Luke Hallum
% v1.1 20150209 Luke Hallum
%               Minor improvements: matrix 'behaviour' now contains the
%               absolute orientation tested, not an index from 1 to 8; and
%               'behaviour' is now also written to a text file.
% v1.2 20150220 Mike Hawken
%               Added flag for 'full screen'. Added the standard itself
%               to the set of test stimuli for direct estimation of
%               bias.
%%%%%%%%%%

% Feedback tones disabled.
%global MGL

%%%
% Experiment parameters.
%%%%%%%%%%
NUM_TRIALS = 100000;       %% Really big? Then number of trials is determined by the design vector (see below).
if fullScreen
    RES_STIMULUS = [1000 1000];%% Stimulus must be square!
else
    RES_STIMULUS = [500 500];%% Stimulus must be square!
end
RATE_RENDER_FPS = 40;      %% These numbers need to be chosen 'together'. E.g., a
DURATION_TOPUP_S = 0.1;    %% render rate of 10fps won't work with a probe
DURATION_BLANK_S = 0.1;    %% duration of 0.55 seconds.
DURATION_PROBE_S = 0.05;   %%
DURATION_RESPONSE_S = 1.5; %%
durationTrial_s = sum([DURATION_TOPUP_S DURATION_BLANK_S DURATION_PROBE_S DURATION_RESPONSE_S]);
numFramesInTrial = RATE_RENDER_FPS * (durationTrial_s - DURATION_RESPONSE_S);
ixFramesTopup = 1:(RATE_RENDER_FPS*DURATION_TOPUP_S);
ixFramesBlank = (ixFramesTopup(end)+1):(ixFramesTopup(end)+1+RATE_RENDER_FPS*DURATION_BLANK_S);
ixFramesProbe = (ixFramesBlank(end)+1):(ixFramesBlank(end)+1+RATE_RENDER_FPS*DURATION_PROBE_S);
%%%
% Display parameters and calculations you'll probably need.
%%%%%%%%%%
%% NB: the unit here is 'stimulus width'.
[x,y] = meshgrid(linspace(-0.5,0.5,RES_STIMULUS(2)));
%%%%%%%%%%

%%%
% Build the required rotated coordinate systems.
%%%%%%%%%%
NUM_TEST_ORI = 8; % This must be even. Like, 8.
XR = cell(NUM_TEST_ORI+1,1);
orisTest_deg = linspace(oriStandard_deg-oriTestRange_deg/2, oriStandard_deg+oriTestRange_deg/2, NUM_TEST_ORI);
orisTest_deg = [orisTest_deg(1:(NUM_TEST_ORI/2)) oriStandard_deg orisTest_deg(1+(NUM_TEST_ORI/2):NUM_TEST_ORI)]; 
for iiori = 1:(NUM_TEST_ORI+1)
  XR{iiori} = baRotateCoords(orisTest_deg(iiori),x,y);
end
%%%%%%%%%%

%%%
% Seed rand -- see 'help rand'.
%%%%%%%%%%
rand('twister',seedRandom);
%%%%%%%%%%

%%%
% Design issues.
%%%%%%%%%%
vectorDesignOri = repmat(1:(NUM_TEST_ORI+1),[1 numRepeats]);
vectorDesignOri = vectorDesignOri(randperm(length(vectorDesignOri)));
%%%%%%%%%%

%%%
% Some stimulus parameters...
% 
% Here, 1080 pixels occupies about 27cm on the screen. So,
% assuming a viewing distance of 57cm, the stimulus
% (probably 1000 pixels tall) occupies 25 deg.  That is, a 3
% cyc/deg grating cycles 75 times across the stimulus.
%%%%%%%%%%
F1_CPD = 4;
F1_CPSTIM = F1_CPD*100*(RES_STIMULUS(1)/1080 * 0.27);
%%%%%%%%%%

%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%

mglOpen();
%mglOpen(0,RES_STIMULUS(1),RES_STIMULUS(2),60,32);
%%%
% Gamma correction.
%%%%%%%%%%
%%% load 'calib_brownie_20100503';
%%% mglSetGammaTable(calib.linear_table);
%%%%%%%%%%
fieldMeanLum = 0.5 * ones(RES_STIMULUS);
mglgrey = mglCreateTexture(255.0 * fieldMeanLum);
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

%%%
% The format of matrix 'behaviour'.
%%%%%%%%%%
IX_BEH_LEVEL_MOD = 1;
IX_BEH_CW = 2;
IX_BEH_CORRECT = 3;
behaviour = -1*ones(min(NUM_TRIALS,length(vectorDesignOri)),3);
%%%%%%%%%%
% Trial 1 starts after a key press.
%%%%%%%%%%
mglBltTexture(mglgrey,[0 0]); mglFlush();
respKeys = 0 * mglGetKeys();
KEYS_RESPONSE = [39 41];
while (sum(respKeys([KEYS_RESPONSE])) < 1), respKeys = mglGetKeys(); end
for iiTrial = 1:min(NUM_TRIALS,length(vectorDesignOri))

  clear framesTrial mgltex

  mglBltTexture(mglgrey,[0 0]); mglFlush();
  % clk2 is used to keep the ITI about constant.
  clk2 = clock;

  % Beep signals trial's start.
  %fprintf(1,'\a')

  %%%
  % Here, build frames for this trial.
  %%%%%%%%%%
  framesTrial = zeros([RES_STIMULUS floor(numFramesInTrial)]);
  % Top-up frames.
  framesTrial(:,:,ixFramesTopup) = 0.5*zeros([size(x) length(ixFramesTopup)]);
  % Blank frames.
  framesTrial(:,:,ixFramesBlank) = 0.5*zeros([size(x) length(ixFramesBlank)]);
  % Probe frames.
  grating = baAnnulus(conStim*cos(2*pi*F1_CPSTIM*XR{vectorDesignOri(iiTrial)} - 1000*rand),x,y,0.0,0.1);
  framesTrial(:,:,ixFramesProbe) = repmat(grating,[1 1 length(ixFramesProbe)]);
  %%%%%%%%%%
  %%%%%%%%%%
%DEBUG_FILENAME = 'frames.mat';
%save(DEBUG_FILENAME,'framesTrial');
  for iiFrame = 1:size(framesTrial,3), mgltex(iiFrame) = mglCreateTexture(255.0 * (0.5 + 0.5*framesTrial(:,:,iiFrame))); end
  while (etime(clock,clk2) < 1.0), end % Keeps the ITI about constant.

  %%%%%%%%%%
  %%%%%%%%%%
  % Main rendering loop. Then delete all those MGL textures.
  %%%%%%%%%%
  %%%%%%%%%%
  clk1 = clock;
%%tic % Outer tic/toc for verifying total trial duration.
  for iiFrame = 1:size(framesTrial,3)

%%tic % Inner tic/toc for checking MGL render is very, very fast.
    mglBltTexture(mgltex(iiFrame),[0 0]); mglFlush();
%%toc

    while (etime(clock,clk1) < iiFrame * (1/RATE_RENDER_FPS)), end
  end
%%toc
  for iiFrame = 1:size(framesTrial,3), mglDeleteTexture(mgltex(iiFrame)); end
  %%%%%%%%%%
  %%%%%%%%%%
  %%%%%%%%%%

  %%%
  % Get behaviour. Feed back.
  %%%%%%%%%%
  mglBltTexture(mglgrey,[0 0]); mglFlush();
  respKeys = 0 * mglGetKeys();
  % These key codes correspond to 'j' and 'k'. We'll use 'j' for 'counter
  % clockwise' and 'k' for 'clockwise'.
  while (etime(clock,clk1) < durationTrial_s & sum(respKeys([KEYS_RESPONSE])) < 1), respKeys = mglGetKeys(); end
  if (sum(respKeys([KEYS_RESPONSE])) > 0)
    behaviour(iiTrial,IX_BEH_LEVEL_MOD) = orisTest_deg(vectorDesignOri(iiTrial));
    behaviour(iiTrial,IX_BEH_CW) = respKeys(KEYS_RESPONSE(2)) == 1;
    
    behaviour(iiTrial,IX_BEH_CORRECT) = vectorDesignOri(iiTrial) <= NUM_TEST_ORI/2 == respKeys(KEYS_RESPONSE(1));
    %behaviour(iiTrial,IX_BEH_CORRECT) = vectorDesignOri(iiTrial) <= NUM_TEST_ORI/2 == respKeys(KEYS_RESPONSE(1));
  end
%% Feedback disabled.
%%if(behaviour(iiTrial,IX_BEH_CORRECT) < 0.5), mglPlaySound(find(strcmp(MGL.soundNames,'Basso')));
%%else, mglPlaySound(find(strcmp(MGL.soundNames,'Blow')));
%%end
  while (etime(clock,clk1) < durationTrial_s), end
  %%%%%%%%%%
end

mglClose;

filename = sprintf('basic_disc_beh_%s_%s',subjectID,datestr(now,30));
save(filename,'behaviour','seedRandom');
dlmwrite(sprintf('%s.text',filename),behaviour,'delimiter',' ')

%return;

%%%
% Some functions...
%%%%%%%%%%

function xR = baRotateCoords(ori_deg,x,y)
% This function returns rotated coordinate axes.  Those input coordinate
% matrices need to be square. If not, caveat emptor.

  opRot = [cos(pi*ori_deg/180) -sin(pi*ori_deg/180); sin(pi*ori_deg/180) cos(pi*ori_deg/180)];
  aa = opRot * [x(:)'; y(:)'];
  xR = reshape(aa(1,:),size(x));

return;

function maskedGrating = baAnnulus(grating,x,y,innerRadius_stimwid,outerRadius_stimwid)
% Argument 'grating' should be in contrast space; as is the returned masked
% grating. The annulus is smoothed using a Gaussian kernel -- the constant
% below specified its standard deviation.

  SD_GAUSSIAN_STIMWID = 0.02;

  r = sqrt(x.^2 + y.^2);
  diskMask = double(r < outerRadius_stimwid & r > innerRadius_stimwid);
  kernelGauss = fspecial('gaussian',round(6*[SD_GAUSSIAN_STIMWID * size(x,2) SD_GAUSSIAN_STIMWID * size(x,2)]),SD_GAUSSIAN_STIMWID * size(x,2));
  diskMask = conv2(double(diskMask),kernelGauss,'same');
  maskedGrating = grating .* diskMask;

return;

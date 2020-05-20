%% Settings for image sequence experiment


imageTarget = 'T.png';
imageDistractor = 'DC.png';
imageDistractor2='T2.png';
rotateDistractor = 1;



% Length of display in seconds
DisplayTime=0.150;

% Set size 
setSize = [7];

% Number of trials per set size
nTrials = 10;

% Block trials by set size? 
% 1: each set size will be run in a separate block
% 0: set sizes will be randomly interleaved.
blockSetSize = 0;

% Response keys 
responseKeys = {'p','q'};

% Number of trials to show before a break (for no breaks, choose a number
% greater than the number of trials in your experiment)
breakAfterTrials = 100;

% Background color (0 to 255)
backgroundColor = 0;

% Text color
textColor = 255;

% How long to wait for subject response before the trial times out
trialTimeout = 5;

% How long to pause in between trials 
% if 0, the experiment will wait for the subject to press a key before every trial)
timeBetweenTrials = 1;
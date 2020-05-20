 function ConjunctionSearch(subID)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up the experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

settingsConjunctionSearch; 
rand('state', sum(100*clock)); % Initialize the random number generator

% Keyboard setup
KbName('UnifyKeyNames');
KbCheckList = [KbName('space'),KbName('ESCAPE')];
for i = 1:length(responseKeys)
    KbCheckList = [KbName(responseKeys{i}),KbCheckList];
end
RestrictKeysForKbCheck(KbCheckList);

% Screen setup
clear screen
whichScreen = max(Screen('Screens'));
[window1, rect] = Screen('Openwindow',whichScreen,backgroundColor,[],[],2);
slack = Screen('GetFlipInterval', window1)/2;

W=rect(RectRight);% screen width
H=rect(RectBottom); % screen height
Screen(window1,'FillRect',backgroundColor);
Screen('Flip', window1);

% Squeezing stimuli around fixation cross
cropFactorX=H/6;
cropFactorY=H/6;

% conditions
leftPresConCorr=0; leftPressConErr=0; leftAbsConCorr=0; leftAbsConErr=0;
rightPresConCorr=0; rightPressConErr=0; rightAbsConCorr=0; rightAbsConErr=0;

% Make beep
cf = 500;                        % carrier frequency (Hz)
sf = 22050;                      % sample frequency (Hz)
d = 0.7;                         % duration (s)
n = sf * d;                      % number of samples
s = (1:n) / sf;                  % time-dependent values
tone = sin(2 * pi * cf * s);     % sinusoidal modulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up stimuli lists and results file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the image files for the experiment
imageFolder = 'images';
imTarg = imread([imageFolder '/' imageTarget]);
imDist = imread([imageFolder '/' imageDistractor]);
imDist2 = imread([imageFolder '/' imageDistractor2]);

% Set up the trials
nSetSizes = length(setSize);
setSize = repmat(setSize,[nTrials 1]);
% Randomly assign half of trials to be target-present; rest will be
% target-absent
randomizedTrials = randperm(nTrials*nSetSizes);
targetPresent = zeros(size(setSize));
for i = 1:nSetSizes
    t = zeros(nTrials,1);
    t(1:round(nTrials/2)) = 1;
    t = t(randperm(nTrials));
    targetPresent(:,i) = t;
end


% Set up the output file
resultsFolder = 'results';
outputfile = fopen([resultsFolder '/resultConjunctionSearch_' num2str(subID) '.txt'],'a');
fprintf(outputfile, 'subID\t trial\t setSize\t targetPresent\t response\t loc\n');

% Squeezing stimuli around fixation cross
cropFactorX=H/6;
cropFactorY=H/6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start screen
Screen('DrawText',window1,'Press the space bar to begin', (W/2-300), (H/2), textColor);
Screen('Flip',window1);
% Wait for subject to press spacebar
while 1
    [keyIsDown,secs,keyCode] = KbCheck;
    if keyCode(KbName('space'))==1
        break
    end
end

% Run experimental trials
for n = 1:length(randomizedTrials)
    t = randomizedTrials(n);
    
    % making mask fram
    x=[];y=[];
    for i=1:60 %number of lines in mask frame
        y(1,i)= rand*H;
        x(1,i)= rand*W;  
    end
    
    xy=[x;y]; 
    
    % Make search display
    if rand <0.5
        img = makeSearchDisplayLeft(imTarg,imDist,imDist2,setSize(t),targetPresent(t),rotateDistractor,backgroundColor,H,W);
        loc=('L');
    else
        img = makeSearchDisplayRight(imTarg,imDist,imDist2,setSize(t),targetPresent(t),rotateDistractor,backgroundColor,H,W);
        loc=('R');
    end
    imageDisplay = Screen('MakeTexture', window1, img);
    % Calculate image position (center of the screen)
    imageSize = size(img);
    pos = [(W-imageSize(2))/2+cropFactorX (H-imageSize(1))/2+cropFactorY (W+imageSize(2))/2-cropFactorX (H+imageSize(1))/2-cropFactorY];

    % Screen priority
    Priority(MaxPriority(window1));
    Priority(2);
    
    % Show fixation cross
    fixationDuration = 0.5; % Length of fixation in seconds
    drawCross(window1,W,H);
    tFixation = Screen('Flip', window1);

    % Blank screen
    Screen(window1, 'FillRect', backgroundColor);
    Screen('Flip', window1, tFixation + fixationDuration - slack,0);
    
    %Display time
    DisplayTime = setSize(t)* 0.3; % allocating 30 ms per item
   
    % Show the search display
    Screen(window1, 'FillRect', backgroundColor);
    Screen('DrawTexture', window1, imageDisplay, [], pos);
    drawCross(window1, W,H);
    startTime = Screen('Flip', window1); % begining of stimulus 
    Screen('FillRect', window1, backgroundColor);
    Screen('Flip', window1,startTime+DisplayTime);%End of stimulus;
    Screen('FillRect', window1, backgroundColor);
    drawCross(window1, W,H);
    Screen('Flip', window1,startTime+DisplayTime+0.1); % 0.1s Delay after stimulus    
    Screen('DrawLines',window1,xy,7);%Masking Frame with 100ms delay
    drawCross(window1, W,H); 
    Screen('Flip', window1,startTime+DisplayTime+0.1);
    Screen(window1, 'FillRect', backgroundColor);
    drawCross(window1, W,H);
    Screen('Flip', window1,startTime+DisplayTime+0.1+0.075);
    
    
    
    % Get keypress response
    rt = 0;
    resp = 0;
    while GetSecs - startTime < trialTimeout
        [keyIsDown,secs,keyCode] = KbCheck;
        respTime = GetSecs;
        pressedKeys = find(keyCode);
                   
        % ESC key quits the experiment
        if keyCode(KbName('ESCAPE')) == 1
            clear all
            close all
            sca
            return;
        end
        
        % Check for response keys
        if ~isempty(pressedKeys)
            for i = 1:length(responseKeys);
                if KbName(responseKeys{i}) == pressedKeys(1)
                    resp = responseKeys{i};
                    rt = respTime - startTime;
                end
            end
        end
        
        % Exit loop once a response is recorded
        if rt > 0
            break;
        end

    end
    if resp == ['p'], answer = 1;
    elseif resp == ['q'], answer = 0;
    else answer = 2;
    end
    
    
        
    % Blank screen- making different delay fram after response
    d=[0.1 0.5 0.7]; %delay interval
    r=randperm(3);  
    delayAfterReponse=d(r(1));
    Screen(window1, 'FillRect', backgroundColor);
    drawCross(window1, W,H);
    Screen('Flip', window1, tFixation + fixationDuration - slack,0);

    % Save results to file
    fprintf(outputfile, '%s\t %d\t %d\t %d\t %d\t %c\n',...
        subID, n, setSize(t), targetPresent(t), answer, loc);
    
    % Count correct&error trials
  
   if loc==('L');
       kL=(targetPresent(t)*10+answer);
       switch kL
           case 11
               leftPresConCorr=leftPresConCorr+1;
           case 10 || 12
               leftPressConErr=leftPressConErr+1;
           case 0
               leftAbsConCorr=leftAbsConCorr+1;
           case 1
               leftAbsConErr=leftAbsConErr+1;
       end
   else
       KR=(targetPresent(t)*10+answer);
       switch  KR
           case 11
               rightPresConCorr=rightPresConCorr+1;
           case 10||12
               rightPressConErr=rightPressConErr+1;
           case 0
               rightAbsConCorr=rightAbsConCorr+1;
           case 1
               rightAbsConErr=rightAbsConErr+1;
       end
   end


%Clear textures
Screen(imageDisplay,'Close');
    
end

% analyse data

presConLeftPercent=leftPresConCorr/(leftPresConCorr+leftPressConErr)*100
absconleftPercent=leftAbsConCorr/(leftAbsConCorr+leftAbsConErr)*100
presConRightPercent=rightPresConCorr/(rightPresConCorr+rightPressConErr)*100
absconRightPercent=rightAbsConCorr/(rightAbsConCorr+rightAbsConErr)*100
result=[presConLeftPercent,presConRightPercent,absconleftPercent,absconRightPercent];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End the experiment dt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RestrictKeysForKbCheck([]);
fclose(outputfile);
Screen(window1,'Close');
close all
sca;
bar(result)
return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Draw a fixation cross (overlapping horizontal and vertical bar)
function drawCross(window,W,H)
    barLength = 16; % in pixels
    barWidth = 2; % in pixels
    Screen('FillRect', window, 255,[ (W-barLength)/2 (H-barWidth)/2 (W+barLength)/2 (H+barWidth)/2]);
    Screen('FillRect', window, 255,[ (W-barWidth)/2 (H-barLength)/2 (W+barWidth)/2 (H+barLength)/2]);
end
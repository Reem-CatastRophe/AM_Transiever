% Clear past run ( variables, figures, and command window )
clear FigureNum
clear plotFFT
clear
clc
close all

%% Transmitter Stage Overview
% This section of the code manages the signal processing and modulation steps
% for transmission. It encompasses the handling of audio signals, modulation,
% and the combination of modulated signals.

% **Parameters:**
% - `AudioFolder`: The directory containing stereo audio signals.
% - `CombinedSignal`: Variable to store the combined signal after modulation.
% - `CarrierFrequencyBase`: Base carrier frequency used for modulation.
% - `CarrierFrequencyIncrement`: Step between carrier frequencies for different signals.
% - `AudioFileNames`: Cell array specifying the names of audio files to be processed.

% - **Signal Processing:**
%   - `MaxSignalLength`: The maximum length among all audio signals to prevent aliasing.
%   - `MaxSamplingFrequency`: The maximum sampling frequency among all audio signals.
%   - `OriginalSignals`: Cell array to store the original stereo signals.
%   - `SamplingFrequencies`: Cell array to store the sampling frequencies of the signals.
%   - `MonoOriginalSignals`: Cell array to store mono versions of the original signals.
%   - `SampledOriginalSignals`: Cell array to store upsampled versions of the mono signals.
%   - `ModulatedSignals`: Cell array to store the modulated signals.

% **Signal Reading, Padding, and Modulation:**
% - The code iterates over each audio file, reading the stereo signal, converting
%   it to mono, and determining the maximum sampling frequency.
% - The mono signal is padded with zeros to match the maximum signal length,
%   and the spectrum of the original signal is plotted.
% - The mono signal is upsampled by a factor of 16, and a carrier signal is generated
%   at a specific frequency for modulation.
% - The upsampled signal is modulated with the carrier signal, and the modulated
%   signal is accumulated into the combined signal.
% - The spectrum of the modulated signal is plotted with x-axis limits focused on
%   the carrier frequency range.

% **Combined Signal Spectrum:**
% - After processing all audio files, the spectrum of the combined signal is plotted.

%% Transmitter Stage Code
% The detailed MATLAB code implementing the transmitter stage follows below:

% Define parameters for audio signals and modulation
% Folder containing stereo audio signals
AudioFolder = 'Signals';
% Variable to store the combined signal
CombinedSignal = 0;
% Base carrier frequency for modulation
CarrierFrequencyBase = 100e3;
% Step between carrier frequencies for different signals
CarrierFrequencyIncrement = 55e3;

% Specify audio file names
AudioFileNames = {
    fullfile(AudioFolder, 'Short_BBCArabic2.wav'),...
    fullfile(AudioFolder, 'Short_FM9090.wav'),...
    fullfile(AudioFolder, 'Short_QuranPalestine.wav'),...
    fullfile(AudioFolder, 'Short_RussianVoice.wav'),...
    fullfile(AudioFolder, 'Short_SkyNewsArabia.wav')
};

% Get maximum signal length to prevent aliasing
MaxSignalLength = max(cellfun(@(x) size(audioread(x), 1), AudioFileNames));

% Initialize variables for storing signals and frequencies
MaxSamplingFrequency = 0;
OriginalSignals = cell(1, length(AudioFileNames));
SamplingFrequencies = cell(1, length(AudioFileNames));
MonoOriginalSignals = cell(1, length(AudioFileNames));
SampledOriginalSignals = cell(1, length(AudioFileNames));
ModulatedSignals = cell(1, length(AudioFileNames));

% **Read, pad, and modulate signals**
figure('Name', 'Spectrum of Original Signal');
for i = 1:length(AudioFileNames)
    % Read stereo signal and store its sampling frequency
    [OriginalSignals{i}, SamplingFrequencies{i}] = audioread(AudioFileNames{i});

    % Convert to mono signal by summing channels
    MonoOriginalSignals{i} = sum(OriginalSignals{i}, 2);

    % Determine maximum sampling frequency
    MaxSamplingFrequency = max(MaxSamplingFrequency, SamplingFrequencies{i});

    % Pad signal with zeros to match longest signal length
    PaddingLength = MaxSignalLength - size(MonoOriginalSignals{i}, 1);
    MonoOriginalSignals{i} = [MonoOriginalSignals{i}; zeros(PaddingLength, size(MonoOriginalSignals{i}, 2))];

    % Plot spectrum of original signal (left subplot)
    subplot(length(AudioFileNames), 2, 2*i - 1);
    plotFFT(MonoOriginalSignals{i}, SamplingFrequencies{i}, ['Spectrum of Original Signal ' num2str(i)]);

    % Upsample mono signal by a factor of 16
    UpsamplingFactor = 16;
    SampledOriginalSignals{i} = interp(MonoOriginalSignals{i}, UpsamplingFactor);

    % Generate carrier signal at specific frequency
    CarrierFrequency = CarrierFrequencyBase + (i - 1) * CarrierFrequencyIncrement;
    t = (1:1:length(SampledOriginalSignals{i}))';
    FsCarrier = UpsamplingFactor * SamplingFrequencies{i};
    Carrier = cos(2*pi*CarrierFrequency*t/FsCarrier)/UpsamplingFactor;

    % Modulate upsampled signal with carrier signal
    ModulatedSignals{i} = Carrier .* SampledOriginalSignals{i};

    % Accumulate modulated signals into combined signal
    CombinedSignal = CombinedSignal + ModulatedSignals{i};

    % Plot spectrum of modulated signal (right subplot)
    subplot(length(AudioFileNames), 2, 2*i);
    plotFFT(ModulatedSignals{i}, FsCarrier, ['Spectrum of Modulated Signal ' num2str(i)]);

    % Set x-axis limits to focus on carrier frequency range
    xlim([-1.5 * CarrierFrequency , 1.5 * CarrierFrequency]); 
end

% Plot spectrum of combined signal
figure('Name', 'The spectrum of the output of the transmitter');
plotFFT(CombinedSignal, FsCarrier, 'The spectrum of the output of the transmitter');

%% RF Band Pass Filter Overview
% This section of the code filters the combined signal to isolate a specific
% chosen channel. It provides options to enable or disable the RF filter and
% prompts the user to choose a channel for listening.

% **Parameters:**
% - `RFFilterApplied`: Flag to enable (1) or disable (0) the RF filter.
% - `ChosenSignalIndex`: User-selected index for the desired channel.
% - `ChosenSignal`: Mono version of the chosen signal.
% - `FsChosenSignal`: Sampling frequency of the chosen signal.

% **Channel Selection:**
% - The user is prompted to choose a channel from a list of available options.

% **Bandpass Filter Design:**
% - The modulation frequency is calculated based on the chosen signal index.
% - The occupied bandwidth of the chosen signal is determined for designing
%   the bandpass filter.
% - Specifications for the bandpass filter are determined, including
%   stopband and passband frequencies, attenuations, and filter design parameters.

% **Filter Application:**
% - If the RF filter is enabled, a bandpass filter is designed using the
%   equiripple method.
% - The combined signal is filtered using the designed bandpass filter.
% - The spectrum of the filtered signal is plotted.

%% RF Band Pass Filter Code
% The detailed MATLAB code implementing the RF Band Pass Filter follows below:

% Enable or disable RF filter
RFFilterApplied = 1;  % Set to 1 to enable, 0 to disable

% **Prompt user for channel selection**
fprintf("Which Channel Would you Like To Listen To:\n");
fprintf("1. Short_BBCArabic2\n");
fprintf("2. Short_FM9090\n");
fprintf("3. Short_QuranPalestine\n");
fprintf("4. Short_RussianVoice\n");
fprintf("5. Short_SkyNewsArabia\n");
ChosenSignalIndex = input("Choose: ");

% **Retrieve chosen signal details**
ChosenSignal = MonoOriginalSignals{ChosenSignalIndex};
FsChosenSignal = SamplingFrequencies{ChosenSignalIndex};

% **Calculate parameters for bandpass filter design**

% Initialize variable to get the modulation frequency
NormalizedSignalNumber = 0.55 * ChosenSignalIndex + 0.45;
% Determine the bandwidth of the chosen signal to design BPF based on it
OccupiedBandwidth = 2 * obw(ChosenSignal, FsChosenSignal);
FrequencyModulationFactor = NormalizedSignalNumber * CarrierFrequencyBase;
HalfOccupiedBandwidth = OccupiedBandwidth / 2;


% **Design bandpass filter specifications**
StopbandAttenuation = 60;
StopbandStartFreq = FrequencyModulationFactor - OccupiedBandwidth;
PassbandStartFreq = FrequencyModulationFactor - HalfOccupiedBandwidth;
PassbandStopFreq = FrequencyModulationFactor + HalfOccupiedBandwidth;
StopbandStopFreq = FrequencyModulationFactor + OccupiedBandwidth;
StopbandAttenuation2 = 60;
PassbandAttenuation = 1;

% **Apply filter if enabled**
if RFFilterApplied ~= 0
    % Design bandpass filter using equiripple method
    BandpassFilter = fdesign.bandpass(StopbandStartFreq, PassbandStartFreq, PassbandStopFreq, StopbandStopFreq, ...
                                      StopbandAttenuation, PassbandAttenuation, StopbandAttenuation2, ...
                                      FsChosenSignal * UpsamplingFactor);
    BandpassFilter = design(BandpassFilter, 'equiripple');

    % Filter the combined signal
    FilteredCombinedSignal = filter(BandpassFilter, CombinedSignal);

    % Plot filtered signal spectrum
    figure('Name', 'The Output of the RF filter (before the mixer)');
    plotFFT(FilteredCombinedSignal, FsCarrier, ['The Output of the RF filter (before the mixer) for signal ' num2str(ChosenSignalIndex)])
end

%% Mixer to Demodulate Signal Overview
% This section of the code implements the mixing stage to demodulate the signal.
% It considers the creation of offset values, defines the intermediate frequency (IF),
% and generates the carrier frequency for mixing. The mixer demodulates the
% filtered or combined signal based on the user's choice, with an option to include
% an offset in the IF frequency.

% **Parameters:**
% - `Offset`: Array of offsets used in the code (e.g., [0, 0.1e3, 1e3]).
% - `IF`: Intermediate frequency at which the signal is centered after the IF stage.
% - `Fc`: Carrier frequency of the chosen signal.
% - `FCarrier`: Carrier frequency of the mixer.
% - `FsCarrier`: Sampling frequency of the carrier signal.

% **Signal Demodulation Process:**
% - The mixer creates a carrier signal at the IF frequency.
% - The chosen signal is demodulated by multiplying it with the carrier signal.
% - The demodulated signal is plotted, showcasing the output of the mixer.

%% Mixer to Demodulate Signal Code
% The detailed MATLAB code implementing the Mixer to Demodulate Signal follows below:

% Create array of offests that will be used in code --> Offset = [0,0.1e3,1e3];
Offset = 0;

% Signal is centered at IF frequency after IF stage 
IF = 27.5e3 + Offset; % which may suffer from frequency error that appear in offest

% Create carrier frequency of chosen signal
Fc = NormalizedSignalNumber * CarrierFrequencyBase;

% Carrier frequency of mixer 
FCarrier = Fc + IF;

% Sampling frequency of carrier
FsCarrier = FsChosenSignal * UpsamplingFactor; 

if RFFilterApplied ~= 0
    % Time vector for IFReceivedSignal
    TIF = (1:1:length(FilteredCombinedSignal))';
    
    % Generate carrier signal at IF frequency
    CarrierSignalIF = cos(2 * pi * FCarrier * TIF / FsCarrier);
    
    % The carrier signal is zero-padded to ensure its length matches that of FilteredCombinedSignal
    CarrierSignalIF(end + length(FilteredCombinedSignal) - length(CarrierSignalIF), 1) = 0;
    
    % Modulate signal
    IFReceivedSignal = FilteredCombinedSignal .* CarrierSignalIF;

    % Plot output of mixer
    figure('Name', 'The Output of the RF mixer');
    plotFFT(IFReceivedSignal, FsCarrier, ['The Output of the RF mixer for signal  ' num2str(ChosenSignalIndex) ' With offset ' num2str(Offset)])
else
    % Time vector for IFReceivedSignal
    TIF = (1:1:length(CombinedSignal))';
    
    % Generate carrier signal at IF frequency
    CarrierSignalIF = cos(2 * pi * FCarrier * TIF / FsCarrier);
    
    % The carrier signal is zero-padded to ensure its length matches that of CombinedSignal
    CarrierSignalIF(end + length(CombinedSignal) - length(CarrierSignalIF), 1) = 0;
    
    % Modulate combined signal
    IFReceivedSignal = CombinedSignal .* CarrierSignalIF;

    % Plot output of mixer without RF BPF
    figure('Name', 'Output of the RF mixer (no RF filter)');
    plotFFT(IFReceivedSignal, FsCarrier, ['The Output of the RF mixer (no RF filter) for signal  ' num2str(ChosenSignalIndex)])
end

%% IF Stage with BPF Overview
% This section implements the Intermediate Frequency (IF) stage with a Bandpass Filter (BPF)
% to extract the signal at the IF frequency. It defines the center frequency of the IF stage,
% sets BPF parameters, designs the bandpass filter, and filters the demodulated signal.

% **Parameters:**
% - `IF_BPF`: Center frequency of the IF stage with BPF, avoiding frequency errors.
% - `OccupiedBandwidth`: Adjusted occupied bandwidth to ensure filter validity.
% - `StopbandAttenuation`, `StopbandStartFreq`, `PassbandStartFreq`,
%   `PassbandStopFreq`, `StopbandStopFreq`, `StopbandAttenuation2`,
%   `PassbandAttenuation`: Parameters for bandpass filter design.

% **Bandpass Filter Design:**
% - The bandpass filter specifications are determined based on the IF parameters.
% - The bandpass filter is designed using the equiripple method.

% **Filter Application:**
% - The filtered signal at the IF frequency is obtained by applying the bandpass filter.
% - The output of the IF stage is plotted to showcase the signal at the IF frequency.

%% IF Stage with BPF Code
% The detailed MATLAB code implementing the IF Stage with BPF follows below:

% Center frequency of IF stage with BPF
IF_BPF = 27.5e3; % It doesn't suffer from frequency error

% To ensure that the space between the centered frequency of the filter and zero frequency is enough for the filter 'StopbandStartFreq'
if OccupiedBandwidth > IF_BPF
    OccupiedBandwidth = IF_BPF - 1;
end

% BPF parameters for IF
StopbandAttenuation = 60;
StopbandStartFreq = IF_BPF - OccupiedBandwidth;
PassbandStartFreq = IF_BPF - HalfOccupiedBandwidth;
PassbandStopFreq  = IF_BPF + HalfOccupiedBandwidth;
StopbandStopFreq  = IF_BPF + OccupiedBandwidth;
StopbandAttenuation2 = 60;
PassbandAttenuation = 1;

% Bandpass filter design for IF
BandpassFilter = fdesign.bandpass(StopbandStartFreq, PassbandStartFreq, PassbandStopFreq, StopbandStopFreq, ...
                                  StopbandAttenuation, PassbandAttenuation, StopbandAttenuation2, ...
                                  FsChosenSignal * UpsamplingFactor);
BandpassFilter = design(BandpassFilter, 'equiripple');

% Filter the signal
BPFFilteredSignal = filter(BandpassFilter, IFReceivedSignal);

if RFFilterApplied ~= 0
    % Plotting the output of IF stage
    figure('Name', 'The Output of the IF filter');
    plotFFT(BPFFilteredSignal, FsCarrier, ['The Output of the IF filter for signal  ' num2str(ChosenSignalIndex) ' With offset ' num2str(Offset)])
else
    % Plotting the output of IF stage without RF BPF
    figure('Name', 'Output of the IF filter (no RF filter)');
    plotFFT(BPFFilteredSignal, FsCarrier, ['The Output of the IF filter (no RF filter) for signal  ' num2str(ChosenSignalIndex)])
end

%% Baseband Detection 'Mixer' Overview
% This section implements the baseband detection 'mixer' to shift the signal to baseband.
% It defines the IF frequency used in the mixer for baseband detection, generates
% the carrier signal for mixing, and modulates the filtered signal to shift it to baseband.

% **Parameters:**
% - `IF_BB`: IF frequency used in the baseband detection mixer without frequency errors.
% - `FCarrier`: Carrier frequency for the baseband detection mixer.
% - `FsCarrier`: Sampling frequency of the carrier signal for the mixer.

% **Baseband Detection Process:**
% - The mixer uses the IF frequency to generate a carrier signal for modulation.
% - The filtered signal is modulated with the carrier signal to shift it to baseband.
% - The output of the baseband detection 'mixer' is plotted to showcase the signal at baseband.

%% Baseband Detection 'Mixer' Code
% The detailed MATLAB code implementing the Baseband Detection 'Mixer' follows below:

% IF frequency for mixer used in Baseband detection to shift signal to baseband
IF_BB = 27.5e3; % It doesn't suffer from frequency error

% IF carrier for mixer used in Baseband detection
FCarrier = IF_BB;
FsCarrier = FsChosenSignal * UpsamplingFactor; 

% Time vector for Carrier Signal
TIF = (1:1:length(BPFFilteredSignal))';

% Generate carrier Signal for IF mixer to shift signal to baseband
CarrierSignalIF = cos(2 * pi * FCarrier * TIF / FsCarrier);

% Ensure the length of CarrierSignalIF matches BPFFilteredSignal
CarrierSignalIF(end + length(BPFFilteredSignal) - length(CarrierSignalIF), 1) = 0;

% Modulate BPFFilteredSignal with Carrier IF Signal
IFReceivedSignal = BPFFilteredSignal .* CarrierSignalIF;

if RFFilterApplied ~= 0
    % Plot output of IF mixer before LPF to get message at baseband frequency
    figure('Name', 'The Output of the IF mixer (before the LPF)');
    plotFFT(IFReceivedSignal, FsCarrier, ['The Output of the IF mixer (before the LPF) for signal  ' num2str(ChosenSignalIndex) ' With offset ' num2str(Offset)])
else
    % Plot output of IF mixer before LPF to get message at baseband frequency without RF BPF
    figure('Name', 'Output of the IF mixer before the LPF (no RF filter)');
    plotFFT(IFReceivedSignal, FsCarrier, ['The Output of the IF mixer before the LPF (no RF filter) for signal  ' num2str(ChosenSignalIndex)])
end

%% Baseband Detection LPF Overview
% This section implements the Low-Pass Filter (LPF) for baseband detection
% to isolate the signal at the baseband. It defines LPF parameters, designs
% the LPF using designfilt, and filters the modulated signal to obtain
% the baseband signal.

% **Parameters:**
% - `LPF_cutoff_frequency`: Cutoff frequency of the LPF.
% - `LPF_order`: Order of the LPF.

% **LPF Design:**
% - The LPF is designed using the designfilt function with specified parameters.

% **Filter Application:**
% - The modulated signal is filtered using the designed LPF.
% - The output of the LPF for baseband detection is plotted.

%% Baseband Detection LPF Code
% The detailed MATLAB code implementing the Baseband Detection LPF follows below:

% LPF parameters
LPF_cutoff_frequency = HalfOccupiedBandwidth;  
LPF_order = 100;              

% Design the LPF using designfilt
LPF = designfilt('lowpassfir', 'CutoffFrequency', LPF_cutoff_frequency, 'FilterOrder', LPF_order, 'SampleRate', FsChosenSignal * UpsamplingFactor);

% Filter the signal using the designed LPF
LPFFilteredSignal = filter(LPF, IFReceivedSignal);

% Plotting signal after LPF of Baseband detection
if RFFilterApplied ~= 0
    figure('Name', 'The Output of the LPF');
    plotFFT(LPFFilteredSignal, FsCarrier, ['The Output of the LPF for signal  ' num2str(ChosenSignalIndex) ' at baseband With offset ' num2str(Offset)])
else
    figure('Name', 'Output of the LPF (no RF filter)');
    plotFFT(LPFFilteredSignal, FsCarrier, ['The Output of the LPF (no RF filter) for signal  ' num2str(ChosenSignalIndex) ' at baseband'])
end

%% Save Processed Signal Overview
% This section of the code handles the saving of the processed signal after
% the entire demodulation and filtering process. It specifies the folder,
% filename, and path to save the processed signal as a .wav file.

% **Folder and Filename Configuration:**
% - The code checks if the folder 'ReceivedSignals' exists; if not, it creates
%   the folder to store the saved signals.
% - The base filename is specified based on the presence of the RF filter and the offset value.

% **Save Process:**
% - The processed signal is saved in .wav format using the audiowrite function.
% - A message is displayed indicating the successful save.

%% Save Processed Signal Code
% The detailed MATLAB code implementing the Save Processed Signal follows below:

% Create a folder if it doesn't exist
%folderName = 'ReceivedSignals';
%if ~exist(folderName, 'dir')
    %mkdir(folderName);
%end

% Specify the base filename and path
%if RFFilterApplied ~= 0
    %if Offset == 0
        %baseFilename = 'Received_Without_Offset_Signal';
    %elseif Offset == 0.1e3
        %baseFilename = 'Received_With_Offset_0.1KHz_Signal';
    %else
        %baseFilename = 'Received_With_Offset_1KHz_Signal';
    %end
%else
    %if Offset == 0
        %baseFilename = 'Received_Without_Offset_Signal_Without_RF_Filter';
    %elseif Offset == 0.1e3
        %baseFilename = 'Received_With_Offset_0.1KHz_Signal_Without_RF_Filter';
    %else
        %baseFilename = 'Received_With_Offset_1KHz_Signal_Without_RF_Filter';
    %end
%end

% Save signal in .wav format
%filename = fullfile(folderName, [baseFilename '_' num2str(ChosenSignalIndex) '.wav']);

% Use audiowrite to save the audio signal
%audiowrite(filename, LPFFilteredSignal, FsCarrier);

% Display a message indicating the successful save and play audio
%disp(['Audio signal ', num2str(ChosenSignalIndex), ' saved as ', filename]);
play(LPFFilteredSignal)

%% FFT Plotting Function Overview
% The section defines a MATLAB function, `plotFFT`, which is used to plot the
% FFT (Fast Fourier Transform) of a signal. The function takes the signal, its
% sampling frequency, and a title as input, and it generates a plot with the
% FFT magnitude centered at zero frequency. The function also saves the
% generated figure with a unique filename.

% **Function Parameters:**
% - `signal`: The input signal for which the FFT is plotted.
% - `Fs`: The sampling frequency of the signal.
% - `titleText`: The title text for the plot.

% **FFT Plotting Process:**
% - The FFT of the signal is calculated and shifted to center zero frequency.
% - A frequency axis centered at zero is created based on the FFT result.
% - The FFT magnitude is plotted against the centered frequency axis.

% **Figure Saving:**
% - The function saves the generated figure with a unique filename in the
%   'FigurePlots' folder, creating the folder if it doesn't exist.

%% FFT Plotting Function Code
% The detailed MATLAB code defining the FFT Plotting Function follows below:

function plotFFT(signal, Fs, titleText)
    % plotFFT: Plots the FFT of a signal and saves the figure
    
    % Perform FFT and shift zero frequency to the center
    FFTResult = fftshift(fft(signal));
    
    % Create a frequency axis centered at zero
    FrequencyAxis = (-length(FFTResult)/2:1:length(FFTResult)/2 - 1)';
    Frequencies = FrequencyAxis * Fs / length(FFTResult);

    % Plot the FFT magnitude
    plot(Frequencies, abs(FFTResult), 'Color', '#139fff')
    
    % Set plot properties
    title(titleText);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    grid on;

    % Save the figure with a unique filename
    %persistent FigureNum;
    
    %if isempty(FigureNum)
        %FigureNum = 1;
    %else
        %FigureNum = FigureNum + 1;
    %end

    % Create a folder for saving figures if it doesn't exist
    %folderName = 'FigurePlots';
    %if ~exist(folderName, 'dir')
        %mkdir(folderName);
    %end

    % Define the filename based on the figure number and title
    %saveFilename = fullfile(folderName, [num2str(FigureNum),'- ', titleText, '.png']);
    
    % Save the figure only if the file doesn't already exist
    %if ~exist(saveFilename, 'file')
        %saveas(gcf, saveFilename);
    %end
end

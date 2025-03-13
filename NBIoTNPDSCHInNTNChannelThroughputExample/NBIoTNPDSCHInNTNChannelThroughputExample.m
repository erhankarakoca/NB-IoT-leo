%[text] # NB\-IoT NTN NPDSCH Throughput
%[text] This example demonstrates how to simulate a narrowband Internet of Things (NB\-IoT) narrowband physical downlink shared channel (NPDSCH) throughput in a non\-terrestrial network (NTN) channel. This example supports two narrowband NTN channels: the European Telecommunication Standards Institute (ETSI) Rician fading channel and the International Telecommunication Union Radiocommunication Sector (ITU\-R) P.681 land mobile\-satellite (LMS) channel.
%[text] ## Introduction
%[text] This example measures the NPDSCH throughput of an NB\-IoT link, as defined by the 3GPP NB\-IoT standards \[3\], \[4\], \[5\], \[6\], and \[7\].
%[text] The example models these features.
%[text] - Transport channel coding
%[text] - NPDSCH, narrowband reference signal (NRS), and synchronization signals (narrowband primary synchronization signal NPSS, narrowband secondary synchronization signal NSSS)
%[text] - ETSI Rician channel and ITU\-R P.681 LMS channel
%[text] - Single\-input\-single\-output (SISO) link
%[text] - Doppler pre\-compensation at the transmitter, and Doppler compensation at the receiver
%[text] - Optional power amplifier modeling \
%[text] The figure shows the implemented processing chain. For clarity, the figure does not show the NRS and synchronization signals.
%[text] ![BlockDiagram_png.png](text:image:556c)
%[text] For more information about the NB\-IoT NPDSCH transmitter and receiver processing units, see the [NB\-IoT NPDSCH Block Error Rate](docid:lte_ug#example-NPDSCHBlockErrorRateExample) example.
%[text] To reduce the total simulation time, you can use Parallel Computing Toolbox™ to execute the range of transmit power values of the transmit power loop in parallel.
%%
%[text] ## Configure Simulation Length, Transmitter, and Receiver
%[text] Set the length of the simulation in terms of the number of transport blocks. By default, the example uses 10 transport blocks, but you must use a large number of transport blocks to produce meaningful throughput results. Set the range of transmit power values to simulate. The transmitter power is defined as the power of the time\-domain waveform before performing Doppler pre\-compensation and includes the gain of the power amplifier. The receiver includes its noise figure and the antenna temperature. The noise figure models the receiver internal noise, and the antenna temperature models the input noise. This receiver specifies the noise per antenna element. In this example, you perform the simulation for different repetition values and compare the performance improvement with repetitions. The `iReps` variable is applicable only when NPDSCH does not carry the scheduling information block SIB1\-NB.
numTrBlks = 10;             % Number of simulated transport blocks
iReps = [0 5];              % Range of repetitions simulated
txPower = 30:5:45;          % Transmit power (dBm)
rxNoiseFigure = 6;          % Noise figure (dB)
rxAntennaTemperature = 290; % Antenna temperature (K)
%%
%[text] ## Power Amplifier Configuration
%[text] To configure a memoryless power amplifier nonlinearity, use `enablePA`. You can choose one of these power amplifier models, as defined in Annex A of TR 38\.803\.
%[text] - 2\.1 GHz Gallium Arsenide (GaAs)
%[text] - 2\.1 GHz Gallium Nitride (GaN) \
%[text] Alternatively, you can set `paModel` to `Custom` and use `paCharacteristics` to define the power amplifier characteristics in a matrix with three columns. The first column defines the input power in dBm. The second column defines the output power in dBm. The third column defines the output phase in degrees. When the `paCharacteristics` variable is set to empty and the `paModel` is set to `Custom`, this example uses a 2\.1 GHz laterally\-diffused metal\-oxide semiconductor (LDMOS) Doherty\-based amplifier.
%[text] The memoryless nonlinearity applied to the waveform follows this equation for power amplifiers (excluding a custom configuration).
%[text] $y\_P \\left(n\\right)=\\sum\_{k\\;\\in {\\;K}\_p } a\_k \\;x\\left(n\\right)\\;{\\left|x\\left(n\\right)\\right|}^{2k}${"editStyle":"visual"}
%[text] In this equation,
%[text] - $y\_P \\left(n\\right)${"editStyle":"visual"} is the output signal.
%[text] - $x\\left(n\\right)${"editStyle":"visual"} is the input signal.
%[text] - $K\_p${"editStyle":"visual"} is the set of polynomial degree(s).
%[text] - $a\_k${"editStyle":"visual"} is the polynomial coefficient. \
%[text] By default, the example sets `enablePA` to `false`.
enablePA = false;                 % true or false %[control:checkbox:057b]{"position":[12,17]}
paModel = "2.1GHz GaAs"; % "2.1GHz GaAs", "2.1GHz GaN", or "Custom" %[control:dropdown:089e]{"position":[11,24]}
paCharacteristics = [];         % Lookup table as empty or a matrix with columns: Pin (dBm) | Pout (dBm) | Phase (degrees)
%[text] When you set `enablePA` to `true`, use `scaleFactor` variable to modify the maximum input signal amplitude to excite the power amplifier nonlinearity. The `scaleFactor` controls the operating region of the power amplifier and applies the signal amplitude scaling at each transmit antenna. You can also use `scaleFactor` variable to set power backoff. For example, to provide a power backoff of 3dB to a signal passed through the power amplifier, set `scaleFactor` to \-3\. Ensure that the input signal is within the characterization range of the power amplifier model.
%[text] When both `scaleFactor` and `paCharacteristics` are set to empty and `paModel` is set to `Custom`, the example uses a default value of \-35 dB. In all other cases, when `scaleFactor` is set to empty, the example uses a default value of 0 dB.
scaleFactor = []; % Amplitude scaling, in dB

% If enablePA is set to true, visualize the power amplifier gain and phase
% characteristics
paModelImpl = paModel;
paInputScaleFactor = 0;                                 % in dB
if enablePA == 1
    % Set the power amplifier as applicable for the further processing
    if lower(paModel) == "custom"
        if isempty(paCharacteristics)
            tableLookup = getDefaultCustomPA;
            paInputScaleFactor = -35;
        else
            tableLookup = paCharacteristics;
        end
        % Use table look-up option of comm.MemorylessNonlinearity and provide
        % the power amplifier characteristics
        mnl = comm.MemorylessNonlinearity(Method="Lookup table", ...
            Table=tableLookup);
        plot(mnl)
        paModelImpl = mnl;
    else
        paMemorylessNonlinearity(paModel)
    end
end
% Update the power amplifier input scaling factor, based on scaleFactor
if ~isempty(scaleFactor)
    paInputScaleFactor = scaleFactor;
end
%%
%[text] ## Doppler Compensation Configuration
%[text] The example supports two Doppler compensation configurations: one at the transmitter end and the other at the receiver end. For compensation at the transmitter end, enable `txDopplerCompensator`. Setting the `txDopplerCompensator` variable to `true` pre\-compensates the transmitted waveform for Doppler caused by the satellite movement. This example assumes that the Doppler effect due to satellite movement is known when performing Doppler pre\-compensation. To compensate at the receiver end, enable `rxDopplerCompensator`. By setting the `rxDopplerCompensator` variable to true, the receiver performs Doppler compensation using the NPSS. It is important to note that Doppler estimates may be inaccurate at the receiver if the transmit power is low, resulting in an operating signal\-to\-noise ratio below \-10 dB.
txDopplerCompensator = true; % true or false %[control:checkbox:9af6]{"position":[24,28]}
rxDopplerCompensator = false; % true or false %[control:checkbox:7109]{"position":[24,29]}
%%
%[text] ## Set Up Higher Layer Parameters
%[text] To configure the eNodeB and NPDSCH parameters, set these higher layer parameters.
npdschDataType = "NotBCCH"; % "SIB1NB", "BCCHNotSIB1NB", or "NotBCCH" %[control:dropdown:407b]{"position":[18,27]}
iSF = 0;                               % Resource assignment field in DCI (DCI format N1 or N2)
schedulingInfoSIB1 = 0;                % Scheduling information field in MasterInformationBlock-NB (MIB-NB)
iMCS = 4;                              % Modulation and coding scheme field in DCI (DCI format N1 or N2)
%%
%[text] ## eNodeB and NPDSCH Configuration
%[text] Set these eNodeB parameters.
%[text] - NB\-IoT physical layer cell identity
%[text] - Operation mode \
enb = struct;
enb.NNCellID = 0;                         % NB-IoT physical layer cell identity
enb.OperationMode = "Standalone"; % "Standalone", "Guardband", "Inband-SamePCI", or "Inband-DifferentPCI" %[control:dropdown:9a1b]{"position":[21,33]}
%[text] Set the radio network temporary identifier in the NPDSCH structure.
npdsch = struct;
npdsch.RNTI = 1; % Radio network temporary identifier
%%
%[text] ## Propagation Channel Model Configuration
%[text] Create a channel model object for the simulation. The example supports the ETSI Rician channel and ITU\-R P.681 LMS channel. To obtain the NTN channel from these channels, apply an additional Doppler shift. Use this equation to calculate the Doppler shift due to satellite movement, as specified in TR 38\.811\.
%[text]  $f\_{d,\\textrm{sat}} =\\left(\\frac{\\nu\_{\\textrm{sat}} \\;}{c}\\right)\*\\left(\\frac{R}{R\+h}\\;\\cos \\left(\\alpha\_{\\textrm{model}} \\right)\\right)\*f\_c${"editStyle":"visual"}
%[text] In this equation,
%[text] - $\\nu\_{\\textrm{sat}}${"editStyle":"visual"} is the satellite speed.
%[text] - $c${"editStyle":"visual"} is the speed of light.
%[text] - $R${"editStyle":"visual"} is the Earth radius.
%[text] - $h${"editStyle":"visual"} is the satellite altitude.
%[text] - $\\alpha\_{\\textrm{model}}${"editStyle":"visual"} is the satellite elevation angle.
%[text] - $f\_c${"editStyle":"visual"} is the carrier frequency. \
%[text] By default, this example considers a satellite in low Earth orbit at an altitude of 600 km, with a carrier frequency of 2 GHz. The NB\-IoT user equipment (UE) is moving at a speed of 3 km/h. Also, this example assumes that satellites move in circular orbits.
channel = struct;
channel.NTNChannelType = "ETSI Rician"; % "ETSI Rician" or "ITU-R P.681" %[control:dropdown:34ce]{"position":[26,39]}
channel.CarrierFrequency = 2e9;                % Carrier frequency (in Hz)
channel.ElevationAngle = 50;                   % Elevation angle (in degrees)
channel.MobileSpeed = 3*1000/3600;             % UE speed (in m/s)
channel.MobileAltitude = 0;                    % Mobile altitude (in m)
channel.SatelliteAltitude = 600e3;             % Satellite altitude (in m)
channel.Seed = 73;                             % Random seed
channel.IncludeFreeSpacePathLoss = true;        % Include or exclude free space path loss %[control:checkbox:48e1]{"position":[36,40]}

% Set these fields based on the type of channel selected
if lower(channel.NTNChannelType) == "etsi rician"
    % For ETSI Rician channel, set KFactor
    channel.KFactor = 10;           % In dB
else
    % For ITU-R P.681, set Environment and AzimuthOrientation
    channel.Environment = "Urban";  % "Urban", "Suburban", "RuralWooded", or "Residential"
    channel.AzimuthOrientation = 0; % In degrees
end
%%
%[text] ## Channel Estimator Configuration
%[text] Configure a practical channel estimator by using the `cec` structure. By default, this example configures the channel with these specifications.
%[text] - Carrier frequency — 2 GHz
%[text] - Speed of NB\-IoT UE — 3 km/h \
%[text] This configuration results in a Doppler spread of 5\.5 Hz. Therefore, perform frequency averaging over pilot estimates with these settings.
%[text] - Time window — 1 resource element (RE)
%[text] - Frequency window — 25 REs, to ensure averaging over all subcarriers for the resource block \
% Configure channel estimator
cec.PilotAverage = "UserDefined";   % Type of pilot symbol averaging
cec.TimeWindow = 1;                 % Time window size in REs
cec.FreqWindow = 25;                % Frequency window size in REs
cec.InterpType = "Cubic";           % 2-D interpolation type
cec.InterpWindow = "Centered";      % Interpolation window type
cec.InterpWinSize = 3;              % Interpolation window size
cec.Reference = "NRS";              % Channel estimator reference signal
%%
%[text] ## Processing Loop
%[text] To determine the throughput at each repetition index and transmit power index, follow these steps.
%[text] 1. **Generate the transport block** — Get the transport block size depending on the configured higher layer parameters.
%[text] 2. **Generate the resource grid** —  Map the modulated bits, along with the NPSS, NSSS, and NRS, signals to the resource grid. The [`lteNDLSCH`](docid:lte_ref#mw_104bcf67-f870-4f4a-b31b-0af4d98ca362) function performs transport channel coding on the input transport block. The [`lteNPDSCH`](docid:lte_ref#mw_68108778-ba88-4606-b432-3f904efbb5f0) function then modulates the encoded data bits.
%[text] 3. **Generate the waveform** — Generate the NB\-IoT time\-domain OFDM waveform with half subcarrier shift using the [`lteSCFDMAModulate`](docid:lte_ref#bt0lmvu_1) function.
%[text] 4. **Apply power amplifier nonlinearities** — Apply the memoryless nonlinearities to the baseband OFDM signal.
%[text] 5. **Apply Doppler pre\-compensation** — Apply the Doppler shift due to satellite movement to the generated waveform to pre\-compensate the channel\-induced satellite Doppler shift.
%[text] 6. **Model and apply a noisy channel** — Pass the generated waveform through an ETSI Rician or ITU\-R P.681 LMS fading channel to get the faded waveform. Apply path loss and add thermal noise to the faded waveform.
%[text] 7. **Apply Doppler compensation** — Estimate the Doppler shift in the received waveform, and compensate the Doppler shift.
%[text] 8. **Perform synchronization and OFDM demodulation** — Perform timing synchronization by correlating the received waveform with the NPSS. The [`lteSCFDMADemodulate`](docid:lte_ref#bt1zks3) function then demodulates the synchronized signal.
%[text] 9. **Perform channel estimation** — Estimate the channel using NRS.
%[text] 10. **Decode the NPDSCH** — Decode the NPDSCH, with the estimated channel and noise variance, by using the [`lteNPDSCHDecode`](docid:lte_ref#mw_5cc19f14-aa08-49fe-b8b6-dfad0c9920c4) function.
%[text] 11. **Decode the transport block** — Decode the soft bits using the [`lteNDLSCHDecode`](docid:lte_ref#mw_99b0811b-04df-4037-bcd3-4a37a5b08535) function. The function decodes the codeword data and returns the block cyclic redundancy check (CRC) error. \
% Use the higher layer parameters, and check if the provided configuration
% is valid
numRep = numel(iReps);
npdschInfo = hNPDSCHInfo;
npdschInfo.NPDSCHDataType = npdschDataType;
npdschInfo.ISF = iSF;
npdschDataTypeLower = lower(npdschDataType);
if npdschDataTypeLower == "sib1nb"  % NPDSCH carrying SIB1-NB
    npdschInfo.SchedulingInfoSIB1 = schedulingInfoSIB1;
    % Store a copy of the information structure for all the repetitions
    npdschInfo = repmat(npdschInfo,numRep,1);
else % NPDSCH not carrying SIB1-NB
    npdschInfo.IMCS = iMCS;              % Modulation and coding scheme field in DCI (DCI format N1 or N2)
    % Store a copy of the information structure for all the repetitions
    npdschInfo = repmat(npdschInfo,numRep,1);
    for repIdx = 1:numRep
        npdschInfo(repIdx).IRep = iReps(repIdx); % Repetition number field in DCI (DCI format N1 or N2)
    end
end

% Initialize some parameters of enb
enb.NFrame = 0;
enb.NSubframe = 0;
enb.NBRefP = 1;
opMode = lower(enb.OperationMode);
inbandSamePCI = (opMode == "inband-samepci");
inbandDifferentPCI = (opMode == "inband-differentpci");
if inbandSamePCI
    enb.CellRefP = enb.NBRefP;     % Number of cell RS antenna ports (Fixed to 1 in this example)
    enb.NCellID = enb.NNCellID;
elseif inbandDifferentPCI
    enb.CellRefP = enb.NBRefP;     % Number of cell RS antenna ports (Fixed to 1 in this example)
    enb.NCellID = 1;
end
if ((npdschDataTypeLower == "bccnnotsib1nb") || (npdschDataType == "notbcch")) && ...
        (inbandSamePCI || inbandDifferentPCI)
    enb.ControlRegionSize = 3;     % The allowed values are 0...13
end

% Apply default window size according to TS 36.104 Table E.5.1-1a
if(~isfield(enb,"Windowing"))
    enb.Windowing = 6;
end

% Store enb structure with a name used for OFDM modulation and
% demodulation. The NB-IoT downlink waveform is a 1/2 subcarrier shift
% waveform. The lteSCFDMAModulate and lteSCFDMADemodulate functions use the
% NBULSubcarrierSpacing field to modulate and demodulate the NB-IoT
% downlink waveform, respectively.
enbOFDM = enb;
enbOFDM.NBULSubcarrierSpacing = "15kHz";

% Get the waveform information and set up the NTN channel
waveformInfo = lteSCFDMAInfo(enbOFDM);
ntnChannel = setupNTNChannel(channel,waveformInfo.SamplingRate);

% Compute the noise amplitude per receive antenna
kBoltz = physconst('Boltzmann');
NF = 10^(rxNoiseFigure/10);
T0 = 290;                                               % Noise temperature at the input (K)
Teq = rxAntennaTemperature + T0*(NF-1);                 % K
N0_ampl = sqrt(kBoltz*waveformInfo.SamplingRate*Teq/2.0);

% Compute path loss based on the elevation angle and satellite altitude
c = physconst("lightspeed");
d = slantRangeCircularOrbit( ...
    channel.ElevationAngle,channel.SatelliteAltitude,channel.MobileAltitude);
lambda = c/channel.CarrierFrequency;
pathLoss = fspl(d,lambda)*double(channel.IncludeFreeSpacePathLoss); % in dB

% Initialize throughput result
numTxPow = numel(txPower);
throughputPercent = zeros(numTxPow,numRep);

% Absolute subframe number at the starting point of the simulation
NSubframe = enb.NFrame*10+enb.NSubframe;

% Loop over repetitions
repVal = zeros(numRep,1);
for repIdx = 1:numRep %[output:group:271f66c2]
    % Add these fields to the npdsch structure
    npdsch.NSF = npdschInfo(repIdx).NSF;
    npdsch.NRep = npdschInfo(repIdx).NRep;
    npdsch.NPDSCHDataType = npdschDataType;
    repVal(repIdx) = npdsch.NRep;

    % Get the bit capacity and transport block length
    [~,info] = lteNPDSCHIndices(enb,npdsch);
    rmoutlen = info.G;                 % Bit length after rate matching (codeword length)
    trblklen = npdschInfo(repIdx).TBS; % Transport block size

    % The temporary variables 'enb_init', 'enbOFDM_init', and
    % 'channel_init' create the temporary variables 'enb', 'enbOFDM', and
    % 'ntnChannel' within the transmit power loop to create independent
    % simulation loops for the 'parfor' loop
    enb_init = enb;
    enbOFDM_init = enbOFDM;
    channel_init = ntnChannel;

    for txPowIdx = 1:numTxPow
    % parfor txPowIdx = 1:numTxPow
    % To enable the use of parallel computing for increased the speed,
    % comment out the 'for' statement and uncomment the 'parfor' statement.
    % This functionality requires the Parallel Computing Toolbox. If you do
    % not have Parallel Computing Toolbox, 'parfor' defaults to the normal
    % 'for' statement.

        % Reset the random number generator so that each transmit power
        % point experiences the same noise realization
        rng(0,"threefry");

        enb = enb_init;                        % Initialize eNodeB configuration
        enbOFDM = enbOFDM_init;                % Initialize eNodeB configuration related to OFDM waveform
        ntnChannel = channel_init;             % Initialize fading channel configuration
        txcw = [];                             % Initialize the transmitted codeword
        numBlkErrors = 0;                      % Number of transport blocks with errors
        estate = [];                           % Initialize NPDSCH encoder state
        dstate = [];                           % Initialize NPDSCH decoder state
        lastOffset = 0;                        % Initialize overall frame timing offset
        offset = 0;                            % Initialize frame timing offset
        subframeGrid = lteNBResourceGrid(enb); % Initialize the subframe grid
        foffsetRS = 0;                         % Initialize frequency offset using reference signal 

        N0 = N0_ampl;
        pl_dB = pathLoss;
        subframeIdx = NSubframe;
        numRxTrBlks = 0;
        reset(ntnChannel.BaseChannel);
        reset(ntnChannel.ChannelFilter);
        while (numRxTrBlks < numTrBlks)

            % Set current subframe and frame numbers  
            enb.NSubframe = mod(subframeIdx,10);
            enb.NFrame = floor((subframeIdx)/10);

            % Generate the NPSS symbols and indices
            npssSymbols = lteNPSS(enb);
            npssIndices = lteNPSSIndices(enb);
            % Map the symbols to the subframe grid
            subframeGrid(npssIndices) = npssSymbols;

            % Generate the NSSS symbols and indices
            nsssSymbols = lteNSSS(enb);
            nsssIndices = lteNSSSIndices(enb);
            % Map the symbols to the subframe grid
            subframeGrid(nsssIndices) = nsssSymbols;

            % Establish if either NPSS or NSSS is transmitted, and if so,
            % do not transmit NPDSCH in this subframe
            isDataSubframe = isempty(npssSymbols) && isempty(nsssSymbols);

            % Create a new transport block, and encode it when the
            % transmitted codeword is empty. The receiver sets the codeword
            % to empty to signal that all subframes in a bundle have been
            % received (it is also empty before the first transmission)
            if isempty(txcw)
                txTrBlk = randi([0 1],trblklen,1);
                txcw = lteNDLSCH(rmoutlen,txTrBlk);
            end

            if (isDataSubframe)
                % Generate NPDSCH symbols and indices for a subframe
                [txNpdschSymbols,estate] = lteNPDSCH(enb,npdsch,txcw,estate);
                npdschIndices = lteNPDSCHIndices(enb,npdsch);
                % Map the symbols to the subframe grid
                subframeGrid(npdschIndices) = txNpdschSymbols;
                % Generate the NRS symbols and indices
                nrsSymbols = lteNRS(enb);
                nrsIndices = lteNRSIndices(enb);
                % Map the symbols to the subframe grid 
                subframeGrid(nrsIndices) = nrsSymbols;
            end

            % Perform OFDM modulation to generate the time domain waveform.
            % Use NB-IoT SC-FDMA to get the 1/2 subcarrier shift on the
            % OFDM modulation.
            txWaveform = lteSCFDMAModulate(enbOFDM,subframeGrid);

            % Normalize the waveform with maximum waveform amplitude
            txWaveform0 = txWaveform./max(abs(txWaveform));

            % Apply power amplifier nonlinearities
            txWaveform1 = paMemorylessNonlinearity(paModelImpl,txWaveform0,...
                db2mag(paInputScaleFactor),enablePA);

            % Scale the waveform power based on the input transmit power
            wavePower = 10*log10(sum(var(txWaveform1)));
            powerScaling = (txPower(txPowIdx)-30)-wavePower;      % In dB
            txWaveform1 = db2mag(powerScaling)*txWaveform1;

            % Pad waveform with 25 samples. This covers the range of
            % delays expected from channel modeling (a combination of
            % implementation delay and channel delay spread)
            txWaveform1 = [txWaveform1; zeros(25,enb.NBRefP)]; %#ok<AGROW>

            % Apply Doppler pre-compensation
            txWaveform2 = compensateDopplerShift(enbOFDM,txWaveform1, ...
                ntnChannel.SatelliteDopplerShift,txDopplerCompensator);

            % Pass data through channel model
            rxWaveform = generateNTNChannel(ntnChannel,txWaveform2);

            % Apply path loss to the signal
            rxWaveform = rxWaveform*db2mag(-pl_dB);

            % Add thermal noise to the received time-domain waveform. Multiply
            % the noise variance with 2 as wgn function performs the scaling
            % within.
            noise = wgn(size(rxWaveform,1),size(rxWaveform,2),2*(N0^2),1,"linear","complex");
            rxWaveform = rxWaveform + noise;

            % Perform receiver Doppler compensation using reference signals
            if (enb.NSubframe == 5)
                % Use NPSS signal for estimating Doppler
                refInd = npssIndices;
                refSym = npssSymbols;
                foffsetRS = estimateDopplerShiftUsingRS(enbOFDM,rxWaveform,refInd, ...
                    refSym,rxDopplerCompensator);
            end
            rxWaveform1 = compensateDopplerShift(enbOFDM,rxWaveform,foffsetRS, ...
                rxDopplerCompensator);

            % In this example, the subframe offset calculation relies
            % on NPSS present in subframe 5, so we need to pad the
            % subframes before it so that the frame offset returned by
            % lteNBDLFrameOffset is the offset for subframe 5
            sfTsamples = waveformInfo.SamplingRate*1e-3;
            if (enb.NSubframe==5) 
                padding = zeros([sfTsamples*5,size(rxWaveform1,2)]);
                offset = lteNBDLFrameOffset(enb,[padding; rxWaveform1]);
                if (offset > 25) || (offset < 0)
                    offset = lastOffset;
                end
                lastOffset = offset;
            end

            % Synchronize the received waveform
            rxWaveform1 = rxWaveform1(1+offset:end,:);

            % Perform OFDM demodulation on the received data to recreate
            % the resource grid. Use NB-IoT SC-FDMA to get the 1/2
            % subcarrier shift on the OFDM demodulation.
            rxSubframe = lteSCFDMADemodulate(enbOFDM,rxWaveform1,0.55);

            % Channel estimation
            [estChannelGrid,noiseEst] = lteDLChannelEstimate( ...
                enb,cec,rxSubframe);

            % Data decoding
            if (isDataSubframe)
                % Get NPDSCH indices
                npdschIndices = lteNPDSCHIndices(enb,npdsch);

                % Get PDSCH resource elements from the received subframe.
                % Scale the received subframe by the PDSCH power factor
                % Rho. The PDSCH is scaled by this amount, while the
                % reference symbols used for channel estimation (used in
                % the PDSCH decoding stage) are not.
                [rxNpdschSymbols,npdschHest] = lteExtractResources(npdschIndices, ...
                    rxSubframe,estChannelGrid);

                % Decode NPDSCH
                [rxcw,dstate,symbols] = lteNPDSCHDecode( ...
                                     enb,npdsch,rxNpdschSymbols,npdschHest,noiseEst,dstate);

                % Decode the transport block when all the subframes in a bundle
                % have been received
                if dstate.EndOfTx
                   [trblkout,blkerr] = lteNDLSCHDecode(trblklen,rxcw);
                   numBlkErrors = numBlkErrors + blkerr;
                   numRxTrBlks = numRxTrBlks + 1;
                   % Re-initialize to enable the transmission of a new transport block
                   txcw = [];
                end
            end

            subframeIdx = subframeIdx + 1;

        end

        % Calculate the throughput percentage
        throughputPercent(txPowIdx,repIdx) = 100*(1-(numBlkErrors/numTrBlks));
        fprintf("Throughput(%%) for %d transport block(s) at transmit power %d dBm with %d repetition(s): %.4f \n", ... %[output:0221353b]
            numTrBlks,txPower(txPowIdx),npdsch.NRep,throughputPercent(txPowIdx,repIdx)) %[output:0221353b]

    end

end %[output:group:271f66c2]
%%
%[text] ## Results
%[text] Display the measured throughput, which is the percentage of the maximum possible throughput of the link given the available resources for data transmission.
% Set figure title
if npdschDataType == "SIB1NB"
    npdsch.NSF = 8;
end
figure; grid on; hold on; %[output:332ad7a1]
legendstr = repmat("",numRep,1);
for repIdx = 1:numRep
    plot(txPower,throughputPercent(:,repIdx),"-o") %[output:332ad7a1]
    legendstr(repIdx) = "NRep = " + repVal(repIdx);
end
hold off;xlabel('Input Transmit Power (dBm)'); ylabel('Throughput (%)'); %[output:332ad7a1]
title(npdsch.NPDSCHDataType + ": TBS=" + trblklen + ... %[output:332ad7a1]
    "; NSF=" + npdsch.NSF + "; " + enb_init.NBRefP + " NRS port") %[output:332ad7a1]
legend(legendstr,Location="southeast") %[output:332ad7a1]
%[text] This example figure shows the throughput results obtained by simulating 1000 transport blocks (`numTrBlks = 1000`, `txPower = 30:45`) for repetition indices 0 and 5\. The simulation setup includes the default higher layer, eNodeB, and NPDSCH configuration with an ETSI Rician channel. To obtain the lines corresponding to the `Tx Doppler Comp.`, you can set `txDopplerCompensator` to `true` and `rxDopplerCompensator` to `false`. For the lines corresponding to the `Rx Doppler Comp.`, you can set `txDopplerCompensator` to `false` and `rxDopplerCompensator`  to `true`.
%[text] ![23a_longRun.jpg](text:image:1d52)
%%
%[text] ## Further Exploration
%[text] This example shows the throughput simulation for NB\-IoT NPDSCH in an NTN channel. In addition to the default behavior, you can run the example for these cases.
%[text] - Analyze the throughput at each transmit power value and repetition index for a different satellite orbit by varying the satellite altitude and satellite speed.
%[text] - Observe the link performance without any Doppler compensation techniques by setting the `txDopplerCompensator` and `rxDopplerCompensator` variables to `false`.
%[text] - Observe the link performance with only Doppler compensation at only the receiver by setting `txDopplerCompensator` to `false` and `rxDopplerCompensator` to `true`.
%[text] - Check the throughput performance for different repetitions by changing the `iReps` value.
%[text] - Compare the throughput performance in an NTN and terrestrial network by using the [`lteFadingChannel`](docid:lte_ref#bt3f52s) channel as shown in [NB\-IoT NPDSCH Block Error Rate](docid:lte_ug#example-NPDSCHBlockErrorRateExample) example. \
%%
%[text] ## Supporting Files
%[text] This example uses this helper function.
%[text] - [hNPDSCHInfo](matlab:openExample('shared_satcom_lte/NBIoTNPDSCHInNTNChannelThroughputExample','supportingFile','hNPDSCHInfo.m')) — NB\-IoT NPDSCH information \
%%
%[text] ## References
%[text] 1. ETSI TS 101 376\-5\-5 V1\.3\.1 (2005\-02). GEO\-Mobile Radio Interface Specifications (Release 1); Part 5: Radio interface physical layer specifications; Sub\-part 5: Radio Transmission and Reception; GMR\-1 05\.005\.
%[text] 2. ITU\-R Recommendation P.681\-11 (08/2019). “Propagation data required for the design systems in the land mobile\-satellite service.” P Series; Radio\-wave propagation.
%[text] 3. 3GPP TS 36\.211, "Physical channels and modulation", *3rd Generation Partnership Project; Technical Specification Group Radio Access Network; Evolved Universal Terrestrial Radio Access (E\-UTRA)*.
%[text] 4. 3GPP TS 36\.213, "Physical layer procedures", *3rd Generation Partnership Project; Technical Specification Group Radio Access Network; Evolved Universal Terrestrial Radio Access (E\-UTRA)*.
%[text] 5. 3GPP TS 36\.321, "Medium Access Control (MAC) protocol specification", *3rd Generation Partnership Project; Technical Specification Group Radio Access Network; Evolved Universal Terrestrial Radio Access (E\-UTRA)*.
%[text] 6. 3GPP TS 36\.101, "User Equipment (UE) radio transmission and reception", *3rd Generation Partnership Project; Technical Specification Group Radio Access Network; Evolved Universal Terrestrial Radio Access (E\-UTRA)*.
%[text] 7. 3GPP TS 36\.104, "Base Station (BS) radio transmission and reception", *3rd Generation Partnership Project; Technical Specification Group Radio Access Network; Evolved Universal Terrestrial Radio Access (E\-UTRA)*.
%[text] 8. 3GPP TR 36\.763, "Study on Narrow\-Band Internet of Things (NB\-IoT) / enhanced Machine Type Communication (eMTC) support for Non\-Terrestrial Networks (NTN) (Release 17), *3rd Generation Partnership Project; Technical Specification Group Radio Access Network*.
%[text] 9. 3GPP TR 38\.803, "Study on new radio access technology: Radio Frequency (RF) and co\-existence aspects (Release 14)", *3rd Generation Partnership Project; Technical Specification Group Radio Access Network*.
%[text] 10. 3GPP TR 38\.811, "Study on New Radio (NR) to support non\-terrestrial networks (Release 15)", *3rd Generation Partnership Project; Technical Specification Group Radio Access Network*.
%[text] 11. O. Hammi, S. Carichner, B. Vassilakis, and F.M. Ghannouchi. “Power Amplifiers’ Model Assessment and Memory Effects Intensity Quantification Using Memoryless Post\-Compensation Technique.” IEEE Transactions on Microwave Theory and Techniques 56, no. 12 (December 2008): 3170–79\. \
%%
%[text] ## Local Functions
%[text] This example uses these local functions.
function chanOut = setupNTNChannel(channel,sampleRate)
% Setup NTN channel

    % Assign temporary variables for carrier frequency and maximum Doppler
    % shift due to mobile movement
    fc = double(channel.CarrierFrequency);
    c = physconst("LightSpeed");
    maxDoppler = (double(channel.MobileSpeed)*fc)/c;
    elevAngle = double(channel.ElevationAngle);
    h = double(channel.SatelliteAltitude);
    % Calculate the Doppler shift due to satellite movement
    maxDopplerSat = dopplerShiftCircularOrbit(elevAngle,h,channel.MobileAltitude,fc);
    % Check the maximum Doppler shift and sample rate
    if ((maxDoppler+maxDopplerSat) >= (sampleRate/10))
        error("satcom:setupNTNChannel:MaxDoppler", ...
            "The maximum Doppler shift (%d Hz) due to mobile and satellite " + ...
            "movement, must be less than %d Hz which is one-tenth of SampleRate.", ...
            (maxDoppler + maxDopplerSat),sampleRate/10)
    end

    chanOut = struct;
    chanTypeLower = lower(channel.NTNChannelType);
    if chanTypeLower == "etsi rician"
        channelName = "ETSI Rician";
        baseChannel = etsiRicianChannel;
        baseChannel.SampleRate = sampleRate;
        baseChannel.KFactor = channel.KFactor;
        baseChannel.MaximumDopplerShift = maxDoppler;
    elseif chanTypeLower == "itu-r p.681"
        channelName = "ITU-R P.681";
        baseChannel = p681LMSChannel;
        baseChannel.SampleRate = sampleRate;
        baseChannel.Environment = channel.Environment;
        baseChannel.CarrierFrequency = channel.CarrierFrequency;
        baseChannel.MobileSpeed = channel.MobileSpeed;
        baseChannel.ElevationAngle = channel.ElevationAngle;
        baseChannel.AzimuthOrientation = channel.AzimuthOrientation;
        baseChannel.FadingTechnique = "Sum of sinusoids";
        baseChannel.SatelliteDopplerShift = maxDopplerSat;
    end
    baseChannel.RandomStream = "mt19937ar with seed";
    baseChannel.Seed = channel.Seed;

    % Set the channel filter
    chanFilt = comm.ChannelFilter( ...
                SampleRate=sampleRate,PathDelays=0, ...
                NormalizeChannelOutputs=false);

    % Set the output structure
    chanOut.ChannelName = channelName;
    chanOut.CarrierFrequency = fc;
    chanOut.SatelliteAltitude = h;
    chanOut.ElevationAngle = elevAngle;
    chanOut.BaseChannel = baseChannel;
    chanOut.SatelliteDopplerShift = maxDopplerSat;
    chanOut.ChannelFilter = chanFilt;

end

function [out,sampleTimes] = generateNTNChannel(channel,in)
% Generate NTN channel

    if isprop(channel.BaseChannel,'SatelliteDopplerShift')
        [out,~,sampleTimes] = channel.BaseChannel(in);
    else
        % Get the channel information before channel processing
        prevInfo = info(channel.BaseChannel);
        numSamplesStart = prevInfo.NumSamplesProcessed;

        % Get the path gains of base channel
        [~,pathGainsBase] = channel.BaseChannel(in);

        % Get the channel information after channel processing
        postInfo = info(channel.BaseChannel);
        numSamplesEnd = postInfo.NumSamplesProcessed;

        % Get the channel sample times
        sampleTimes = (numSamplesStart:(numSamplesEnd-1)).'/channel.BaseChannel.SampleRate;

        % Apply satellite Doppler shift to the base channel path gains
        pathGains = pathGainsBase.*exp(1i*2*pi*channel.SatelliteDopplerShift*sampleTimes);

        % Perform channel filtering
        out = channel.ChannelFilter(in,pathGains);
    end

end

function out = compensateDopplerShift(enb,inWave,foffset,flag)
% Perform Doppler shift correction

    if flag
        % Correct frequency offset
        out = lteFrequencyCorrect(enb,inWave,foffset);
    else
        out = inWave;
    end

end

function out = estimateDopplerShiftUsingRS(enb,rxWave,refInd, ...
    refSym,flag)
% Estimate the Doppler shift using NPSS

    if flag
        % Set the Windowing field to 0, as this information is not known at
        % the receiver
        enb.Windowing = 0;
        ofdmInfo = lteSCFDMAInfo(enb);
        K = 12;                             % Number of subcarriers     
        L = 14;                             % Number of OFDM symbols in slot

        % Initialize temporary variables
        rxWave1 = [rxWave; zeros((mod(size(rxWave,1),2)),1)]; % Append zero, if required
        rxLen = size(rxWave1,1);

        % Generate reference waveform
        refGrid = complex(zeros([K L]));
        refGrid(refInd) = refSym;
        refWave = lteSCFDMAModulate(enb,refGrid);
        refWave = [refWave; zeros((rxLen-size(refWave,1)),1)];

        % Compute the correlation of received waveform with reference
        % waveform
        x_wave = rxWave1.*conj(refWave);

        % Compute FFT of the resultant waveform
        x_fft = fftshift(fft(x_wave));

        % FFT bin values
        fftBinValues = (-rxLen/2:(rxLen/2-1))*(ofdmInfo.SamplingRate/rxLen);

        % Use the FFT bin index corresponding to the maximum FFT value.
        % The FFT bin value corresponding to this bin index is the integer
        % frequency offset.
        [~,binIndex] = max(x_fft);
        out = fftBinValues(binIndex);
    else
        out = 0;
    end

end

function varargout = paMemorylessNonlinearity(paModel,varargin)
% Apply power amplifier nonlinearity (TR 38.803)
% out = paMemorylessNonlinearity(paModel,in,enable) returns the
% impaired output.
% paMemorylessNonlinearity(paModel) returns the plot with the gain and
% phase characteristics of the power amplifier

    if nargin == 1
        in_NoScale = randn(1e6,1,'like',1i);
        scaleFactor = 1/sqrt(2);
        enable = 1;
    else
        in_NoScale = varargin{1};
        scaleFactor = varargin{2};
        enable = varargin{3};
    end

    if enable
        in = scaleFactor*in_NoScale;
        if isa(paModel,"comm.MemorylessNonlinearity")
            % paModel is a comm.MemorylessNonlinearity System object
            out = paModel(in);
            paModelName = "";
        else
            absIn = abs(in);
            paModelName = paModel;
            switch lower(paModel)
                case "2.1ghz gaas"
                    % 2.1GHz GaAs
                    out = (-0.618347-0.785905i) * in + (2.0831-1.69506i) * in .* absIn.^(2) + ...
                        (-14.7229+16.8335i) * in .* absIn.^(2*2) + (61.6423-76.9171i) * in .* absIn.^(2*3) + ...
                        (-145.139+184.765i) * in .* absIn.^(2*4) + (190.61-239.371i)* in .* absIn.^(2*5) + ...
                        (-130.184+158.957i) * in .* absIn.^(2*6) + (36.0047-42.5192i) * in .* absIn.^(2*7);
                otherwise
                    % 2.1GHz GaN
                    out = (0.999952-0.00981788i) * in + (-0.0618171+0.118845i) * in .* absIn.^(2) + ...
                        (-1.69917-0.464933i) * in .* absIn.^(2*2) + (3.27962+0.829737i) * in .* absIn.^(2*3) + ...
                        (-1.80821-0.454331i) * in .* absIn.^(2*4);
            end
        end
    else
        out = in_NoScale;
    end

    if nargout > 0
        varargout{1} = out;
    end

    if nargin == 1 || (nargout == 0)
        % Gain Plot
        inpPower = 20*log10(absIn);
        gain = 20*log10(abs(out))-inpPower;
        figure
        subplot(211)
        plot(inpPower,gain,".")
        grid on
        ylim([-Inf 1])
        xlim([-30 0])
        xlabel("Normalized input power (dB)")
        ylabel("Gain (dB)")
        title("Gain Characteristics of PA Model " + paModelName)

        % Phase Plot
        phase = angle(out.*conj(in))*180/pi;
        subplot(212)
        plot(inpPower,phase,".")
        grid on
        xlim([-30 0])
        xlabel("Normalized input power (dB)")
        ylabel("Phase (deg)")
        title("Phase Characteristics of PA Model " + paModelName)
    end

end

function paChar = getDefaultCustomPA()
% The operating specifications for the LDMOS-based Doherty amplifier are:
% * A frequency of 2110 MHz
% * A peak power of 300 W
% * A small signal gain of 61 dB
% Each row in HAV08_Table specifies Pin (dBm), gain (dB), and phase shift
% (degrees) as derived from figure 4 of Hammi, Oualid, et al. "Power
% amplifiers' model assessment and memory effects intensity quantification
% using memoryless post-compensation technique." IEEE Transactions on
% Microwave Theory and Techniques 56.12 (2008): 3170-3179.

    HAV08_Table = ...
        [-35,60.53,0.01;
        -34,60.53,0.01;
        -33,60.53,0.08;
        -32,60.54,0.08;
        -31,60.55,0.1;
        -30,60.56,0.08;
        -29,60.57,0.14;
        -28,60.59,0.19;
        -27,60.6,0.23;
        -26,60.64,0.21;
        -25,60.69,0.28;
        -24,60.76,0.21;
        -23,60.85,0.12;
        -22,60.97,0.08;
        -21,61.12,-0.13;
        -20,61.31,-0.44;
        -19,61.52,-0.94;
        -18,61.76,-1.59;
        -17,62.01,-2.73;
        -16,62.25,-4.31;
        -15,62.47,-6.85;
        -14,62.56,-9.82;
        -13,62.47,-12.29;
        -12,62.31,-13.82;
        -11,62.2,-15.03;
        -10,62.15,-16.27;
        -9,62,-18.05;
        -8,61.53,-20.21;
        -7,60.93,-23.38;
        -6,60.2,-26.64;
        -5,59.38,-28.75];
    % Convert the second column of the HAV08_Table from gain to Pout for
    % use by the memoryless nonlinearity System object.
    paChar = HAV08_Table;
    paChar(:,2) = paChar(:,1) + paChar(:,2);

end
%[text] *Copyright 2022\-2023 The MathWorks, Inc.*

%[appendix]
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[text:image:556c]
%   data: {"imgAlign":"baseline","imgHeight":61,"imgWidth":932,"src":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAA5kAAAA9CAYAAAAnM61DAAANhUlEQVR42u3dTc9bRxnG8Wzytu6SjxWJj9F11qwjyIIlG4TEli0gseqOrgAJqWJFQSAiBcoDpGlK4nZCJ9y5OzNn7Mcvx\/bP0tU8tX3sY8\/f15lrXu\/c+fBnG7psrfmmfPCHP8If4Q9\/hD+6MJX\/fOfHn9KF6hxMRjnhD3+EP8If\/gh\/dDn8CZlMhskQ\/gh\/+CP8Ef7wR0ImMRnCH\/4If4Q\/\/BH+SMgkJsNk8Ic\/wh\/+8Ic\/wh8JmcRkmAzhj\/CHP8If4Q9\/JGQqZCZD+MMf\/vBH+MMf\/vBHQiYxGcIf\/gh\/hD\/8Ef5IyCQmw2QIf4Q\/\/BH+CH\/4EzJ9GUyGyRD+CH\/4I\/wR\/vBHQiYxGcIf\/gh\/hD\/8Ef5IyCQmw2Twhz\/CH\/7whz\/CHwmZxGSYDOGP8Ic\/wh\/hD38kZDIZJkP4wx\/+8Ef4wx\/+8EdCJjEZwh\/+CH+EP\/wR\/kjIJCbDZAh\/hD\/84Q9\/hD\/8CZnEZJgM4Y\/whz\/CH+EPfyRkEpMh\/OGP8Ef4wx\/hj4RMYjJMBn\/4I\/wR\/vBH+CMhs+rDj55v\/vXqzeaTz169d\/9Hf\/n83X0\/+M0\/N69ev3kPqnJMObb1muX5z1++7j5e9N1fPNv89T\/\/bQJb7i+Px+eXc8nvGd+nfo58q5+h9\/hP\/\/Dvb33ueCv\/PzrvfE7l\/fJ3OfN9MBlykbtM\/ma9Z58qr93yUfzxP8If\/8PfGviTP84nf9w6ZJZCLCfTK+TynPh468Ps+qFar59V3it\/6a1CjoVWC6Tc13q8wlvvK68dC60eU9+zvl783OXYeIyQSSpZ+Gt57Mh7VLLwR\/jjf\/zvGkOm\/LH+\/HGrkFlO4I83X773o1wq5NGPOH+o2BLRaoGYLeRyjvF5S4Vcjyufpfd4eax+jlYBxcdzgeb3EDJVRlSy8DdTycreUj01t6bWC8vHz16+89DqNdVXfvv8i7f3xwt19Ofc+lnPo\/pu1W38CX\/8j\/\/hj\/\/hT\/64zPxx65D5o9\/fvEvdsy0JuSt3lPDjl5DhmC3k8hr5vEaFHO\/rFXJ8jWpyPXhHrSeGy6pkqWThb5tKVvSF2KoZWyzr3\/XiEr2yHlP9Jl+QWhevcn+tjB26NwF\/xP\/wx\/\/wJ3+cf\/64dcgsJxKT8tKY6NFQhGweuVUhF+g2hRwLa2ZMdIVrppBbY57j5xyBnbvV8+22LWUucuQid7mVrNwyW734h7\/738U3+k71oeyb2cvLceX4+N5xCM+M7+KP\/xH++B\/+Dhky5Y\/154+9hMxRYm91V0cTiF9qNo88zji+37aFHLvKS+vHUnf1ktH1uqBbLQN6MlWyVLLwt69KVr6oxgtK9Yzqr\/HYWMmKvjKqZLUufvvwJfzxP\/6HP\/6HP\/njsvPH3kJm\/cLL+ONRIbcK6xgtCTHtl3O8bSFXc\/vex\/94+3r5HGILWw+I+qMQMlVGVLLwt8ucpNmW\/Op9rZb8+P+xkjXj0\/jjf4Q\/\/oe\/U4ZM+WO9+WNvITMuozvTktA68UOOic4rM9Vz2KWQW6s7xXPLqzn1VneKE84vNWQ+efJkc3Nzc5SLXKvLfy0T81v87nM1vGI2seW2fheFn6JTzR+5Jv5Osbri0pykylxrTlL1xZk5SfF9VLLOh79c1vm6GBdNybda1qPje++bl9RvDTs8JRfHWkWU\/\/E\/\/ndz8JApf6w3f+w1ZOYP0xoTnZccPsbqTvliGJf8nS3k1gV4dFHNBZbHTOfv4VJD5oMHDzZ3797dPH78uGk2+w6Zo32TTqVSjvWCt+95HK3KUvm71bp1jRe5Y\/J3qn3iRqsrxoaXWDErvlL9KHrraHXFvDqjStb6+cuVmtF1sXcNmj2+dX3NFSUhk\/\/xP\/53iJApf6wzf9wxXMdwnUPenj59+tZo7t27t7l\/\/\/63zObQIXNpSfJc6WntHxTH1PcuZn\/\/+sfY+vHHYQk58LaWUq\/vn1vklpZcrz2Zv\/rz5+99xvr\/ueW3dXGt+07ts\/J1TfytSaPK\/DlVkvB3+5BZK9TRe7YJmbPHL80hij7V2loie2z00\/JvHZURPbAVJrLH9vzuGkIm\/+N\/+KNT8idkMpmD3z744IOvQbvzVtlsDh0y8wToOCymVlTi2PQ69KBWYGqPYGy1ag3L6fWWxtcoF7Xy3NxSFodkxGXSt1lyPQ6XjRfY0ca+ve8DfypZ\/O8y+Ktl\/f1v\/s0LZcyGzJnjZxakyMMY89YQsWcg+1xuGGwNi4yvt+R31xAy+R\/\/wx8JmXTgQr6zOpUhFI8ePTr4nMzc45jnBtV5SfV5v\/7by80v\/\/Ti3QIBJSCW+1sbQMdwN9p7qb5GPce8qXOsVNVKXhwOMrPk+kzIzBXA0Wpq+CP+d\/78xd943sB7m5A5c\/xMD2H03aWFOEaLV1SvHL3ekt8dK2TyP+J\/+BMySUvWBfZkjoZxxT2MSoWjtNaX+8oS02Voanlu7D1sLSrUWip9NFa+NQw2V5xaIXNpyfXZkJlvtdX\/UJUuLanE\/07fk1k8I3rBLiFz6fhtejKrX7a2Img1EmZ\/yiGz9Xqn8Dv8Ef\/DHwmZTOaIY\/IfPnx4sjmZvcpNbgEvobL0YJZ\/S0WlvE7pgYwreLVeezQspxVsR+FxtiezteT6bMjsVawuNWQekz\/if2vjL3tMDYGt0Rcziz+Mjh\/NyawjOvIq2KOex9yTuUvIPLbf4Y\/4H\/5oBSGzNbSxt4rTqVYEbe2ZMxrew2TWt7psbxuRPAexVFrK4j0xdOb5kzHo1dfKlaalVv28X9dsyFxacn2XOZlxlbNLDZn75m\/N2+SsWa1NqQ9dwbe6Z3v+WWV4l5A5On5pddly\/1IorMe15lpuGzJP4Xf466+Meazf\/DX4Gv6Om00OsV7FUieFkLlnI1rzlyxk7n479j6Zo5DZWl22txx1q+cyDr\/KlaEWv73zqZWdn3\/6YqvhsqMl12dCZl5tsbVK7aWFzH3zt9Ztcs6pMsb\/Thsy8x5w24bM3vGjaQJ5Q\/JWKIz+XLypjCzpzaGcCZmn8Dv89RnJ6xrwNf631pDZusZfYiPJ1YbM1nYK8aKVNxEtw3DqhSSC0FrWfNttI37yyf8fe\/Hlm81nX7xunkO+mNbQUnWqno41jMkf3QyX2K2iaLjOafgbbZNTg37czmbUsDFq+c8Lloy2YSj+Vx+rPfTRi1rPq6\/RO7\/csBEf6+0fl8NFXL05PndpW6Gl9+d\/xP\/WzV+rUj5q8ORr+FtzyNxmD834WF4scrQgZW6wiI3Xx9567ipCZm51zHsXxm0p6mO1IPL8s7wVxey2ETPDZaOZxvePf7vIqWQJmZcbMut9re1s4vNH28TkCli9wOQegTzkO4fZOBww\/h0rTdk\/W+fX8snY4Bd73Ouw8iWfbg0r6w1hH70\/\/yP+t17+8noHMz2dfO0Z\/s6kJ3O0NVJrIbM6NWq0BV4rs5xy67mr6cnsVaqjGUQgZhY52WbbiKWQ2TLTFlQucipZKlmXEzLzLV8Y4nYK0dt22but1YIaF33KlZ84kqK3PU0epdE6v+xf0U9nhwdnn86VsTyfOc9X7r0\/\/yP+t\/6QOWpg52vHGbbN\/\/YzJ7PVe7jN1khLW+BlNmc4P+TWc1cRMmOB5KEHcVhr\/JG3egNaw1i33TZiKWTmW3m9tfQ+MRnC32F7MkfbM8SFL6IX5aH8vYtFvn+0euaoMhaf15oKkM+vtapobypCb7hZz6fra5Xtf3rbCi29P\/8j\/ne+PZl87bmQeQY9mS2WR1sj9eZuzmyBF+sPMZesYSumiw+ZueB6LUm9FrTckrXvkLlNz4SLnEqWStb1hMzZRcJ27cmcrYy1WmJLZax3frOVofgZZ316tsVfyCT+dzlzMqsPtOphfE3IXOtw2Tz1bZetkWa3wKvPq7+TXd9PyLxFyKzPXwqZsVDznMx9hczeubXMxEVOJUsl67pCZn7+aJ52a7uF2Crem7s0WxmLLfO9uUt5TnmvMhTfJ85dmvXp2blLQibxv\/Pkr7e6bA2FS3My+Rr+zmVOZuS6tdVdHdk4swVeZbW1BdQpt2K6+JCZV10sS5u3jKG1dcM+hsvW94\/d273VZfPqtEKmSpZK1vWGzNGKha1W9tYKckurMM5Uxsrfz79Z9XY0FaHnk\/G98jDg1vDg7NN1uNC2qzAKmcT\/zpO\/7BOjVTj5mpC55pCZG0V6WyP1VpGd3QKvtXjgGrZiOuuQSS5yygl\/+Duc1n4hwh+fwB\/++Br+yD6ZxGSYDOFPZQx\/+CP88TX88T8SMpkMkyH84Q9\/+CP84Q9\/+CMhk5gM4Q9\/hD\/CH\/4IfyRkEpNhMoQ\/wh\/+8Ic\/wh\/+hExiMkyG8Ef4wx\/hj\/CHPxIyickQ\/vBH+CP84Y\/wR0ImMRn84Q9\/hD\/CH\/4IfyRkEpNhMoQ\/wh\/+CH+EP\/yRkMlkmAzhD3\/4wx\/hD3\/4wx8JmcRkCH\/4I\/wR\/vBH+CMhk5gMk8Ef\/gh\/+MMf\/gh\/+BMyickwGcIf4Q9\/hD\/CH\/5IyCQmQ\/jDH+GP8Ic\/\/OGPhExiMsoJf\/gj\/BH+8Ef4IyGTmAyTIfwR\/vBH+CP84Y+ETCbDZAh\/+MMf\/gh\/ygl\/+CMhk5gM4Q9\/hD\/CH\/4If3QeIZMuW2s3GcIf\/gh\/hD\/8Ef7ocvQVymui22CiKRgAAAAASUVORK5CYII="}
%---
%[control:checkbox:057b]
%   data: {"defaultValue":false,"label":"enablePA","run":"Nothing"}
%---
%[control:dropdown:089e]
%   data: {"defaultValue":"\"2.1GHz GaAs\"","itemLabels":["2.1GHz GaAs","2.1GHz GaN","Custom"],"items":["\"2.1GHz GaAs\"","\"2.1GHz GaN\"","\"Custom\""],"label":"paModel","run":"Nothing"}
%---
%[control:checkbox:9af6]
%   data: {"defaultValue":false,"label":"txDopplerCompensator","run":"Nothing"}
%---
%[control:checkbox:7109]
%   data: {"defaultValue":false,"label":"rxDopplerCompensator","run":"Nothing"}
%---
%[control:dropdown:407b]
%   data: {"defaultValue":"\"NotBCCH\"","itemLabels":["NotBCCH","SIB1NB","BCCHNotSIB1NB"],"items":["\"NotBCCH\"","\"SIB1NB\"","\"BCCHNotSIB1NB\""],"label":"npdschDataType","run":"Nothing"}
%---
%[control:dropdown:9a1b]
%   data: {"defaultValue":"\"Standalone\"","itemLabels":["Standalone","Guardband","Inband-SamePCI","Inband-DifferentPCI"],"items":["\"Standalone\"","\"Guardband\"","\"Inband-SamePCI\"","\"Inband-DifferentPCI\""],"label":"Drop down","run":"Nothing"}
%---
%[control:dropdown:34ce]
%   data: {"defaultValue":"\"ETSI Rician\"","itemLabels":["ETSI Rician","ITU-R P.681"],"items":["\"ETSI Rician\"","\"ITU-R P.681\""],"label":"Drop down","run":"Nothing"}
%---
%[control:checkbox:48e1]
%   data: {"defaultValue":false,"label":"Check box","run":"Nothing"}
%---
%[text:image:1d52]
%   data: {"imgAlign":"baseline","imgHeight":389,"imgWidth":512,"src":"data:image\/jpeg;base64,\/9j\/4AAQSkZJRgABAgAAAQABAAD\/2wBDAAgGBgcGBQgHBwcJCQgKDBQNDAsLDBkSEw8UHRofHh0aHBwgJC4nICIsIxwcKDcpLDAxNDQ0Hyc5PTgyPC4zNDL\/2wBDAQkJCQwLDBgNDRgyIRwhMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjL\/wAARCALIA6kDASIAAhEBAxEB\/8QAHwAAAQUBAQEBAQEAAAAAAAAAAAECAwQFBgcICQoL\/8QAtRAAAgEDAwIEAwUFBAQAAAF9AQIDAAQRBRIhMUEGE1FhByJxFDKBkaEII0KxwRVS0fAkM2JyggkKFhcYGRolJicoKSo0NTY3ODk6Q0RFRkdISUpTVFVWV1hZWmNkZWZnaGlqc3R1dnd4eXqDhIWGh4iJipKTlJWWl5iZmqKjpKWmp6ipqrKztLW2t7i5usLDxMXGx8jJytLT1NXW19jZ2uHi4+Tl5ufo6erx8vP09fb3+Pn6\/8QAHwEAAwEBAQEBAQEBAQAAAAAAAAECAwQFBgcICQoL\/8QAtREAAgECBAQDBAcFBAQAAQJ3AAECAxEEBSExBhJBUQdhcRMiMoEIFEKRobHBCSMzUvAVYnLRChYkNOEl8RcYGRomJygpKjU2Nzg5OkNERUZHSElKU1RVVldYWVpjZGVmZ2hpanN0dXZ3eHl6goOEhYaHiImKkpOUlZaXmJmaoqOkpaanqKmqsrO0tba3uLm6wsPExcbHyMnK0tPU1dbX2Nna4uPk5ebn6Onq8vP09fb3+Pn6\/9oADAMBAAIRAxEAPwD3+iiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAoormvGPjbSfA9jb3mri4MVxL5SeRHuO7GeeR6UAdLRWLqvijTNJ8LSeJJJWm01IkmEkA3F0YjBUf8AAhV3StRh1fSLPU7bf9nu4Eni3jB2soYZHrg0AXaKKzIdf0u516bRILyOXUbeLzpoE5Ma5AG49AeRx1oA06K8xvPjt4Qsr+4s5Y9TMtvK0T7LcEblOD\/F6iul8IfEDw\/43Sf+x7lzNBgywTJsdQehx3H0oA6miuG8VfFbwz4S1Iabcvc3moDl7axiEjp3G7JABxzjOaveDviFoHjiOUaTPItzCN0trcJslQeuOQR9CaAOrorhvFXxW8M+EtSGm3L3N5qA5e2sYhI6dxuyQAcc4zmr3g74haB44jlGkzyLcwjdLa3CbJUHrjkEfQmgDq6K53w14y0vxXcapBpwuA+mz+RP5qBRuyw+XnkfKabb+NdKufHF14RjFx\/adrEJpCUHl7Sqtw2fR17UAdJRRXMeMfHWkeB7a1n1ZblkuXKJ9njDnIGTnkUAdPRXlI\/aD8FngR6r\/wCAy\/8AxVbmv\/Fbw74cstIu75L4x6rbC6txFCGOwgH5ueD8woA7qivN9K+N\/g3VdTgsRPd2kk7BEe6g2puPABIJx9TxWz4w+JGgeC7iC11Frie9nXdHaWke+Qr0zgkAd+\/NAHX0VxfhH4n+HvGV9JYWTXVrfxqWNreRCOQqOpGCQfpnNdpQAUVHNMlvBJNK22ONS7H0AGTXmK\/H7wQzAB9QyeP+Pb\/69AHqVFcv4s8faD4MNvHqs8jXNxzFbW8ZkkYeuOw+tWPCnjHRvGemve6POzrG+yWORNrxt6MKAOgoormvGPjbSfA9jb3mri4MVxL5SeRHuO7GeeR6UAdLRWLqvijTNJ8LSeJJJWm01IkmEkA3F0YjBUf8CFVNR8daFpHhG08TX9w8NhdxRyQKUzJJvXcqhR\/FigDpaK820j42+FNV1SCwlj1LTnuCFhkvrcJG5J45Vjj6nir\/AIr+K\/hzwbrQ0nVEvjc+Usv7iEMNpzjncPSgDuqK4Xwp8VvDvjLWf7K0yO+Fx5TS5nhCrgYzzuPrXdUAFFFczo3jrRNf8TajoWmzST3Ngu6aRVHlHkAhWzyQTj8DQB01FFFABRXMr450STxuPCUU0kup+WZH8tQY0wCSrNn72B09xXTUAFFcpa\/ELQbrxvceERJNHqkOeJEwkhADYU55ODn8DVy98W6bYeMNO8MTCf8AtDUImlhKplNqhicnPH3DQBv0V51r3xn8LeHNdu9HvU1A3Vq2yTyoAy5wDwd3vV7wr8VvC3jDU\/7N065njvSpZIbiLYZABk7eoJA7daAO3orhfFPxX8OeFNXOkzi9vtQUBpLexhEjR5GRuyQOnOK0fB\/j\/QfHEU50maVZ4MedbXCbJEB745BH0JoA6miiigAooooAKKKzIdf0u516bRILyOXUbeLzpoE5Ma5AG49AeRx1oA06KKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKr3d5bWEHnXU8cMe5U3O2BuY4A+pJpt7f2mm2j3V7cxW8Ef35JnCqPxNAFqiqOmatp+s2pudNvYLuEHaXhcMAfQ46VXh8S6Hcao2mQ6vZSXykgwLMpfI6jGevtQBrUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFeOftBwLc6H4egckJLqYRiOuCpFex1wnxM8F6j4ztNHi06a0iayvluZPtDMoKgYwNqnmgDyO71a88LeCPFvw21yQma0QS6ZM3SWIyKxUfh8w\/4EO1dLN4k1+ey+H\/gXw5erp1zqOjW9xPelAzJH5R4X8I3Pr05FdV8V\/hk3j6xtZtPkt7fVrVtqyzkhHiPVWKgng8jj19apav8MNYax8Kalomp21p4l0GxitC77jDMqrgjOM45ftyG7UAVLDV\/E\/gP4k6X4Z1vXH13S9XQ+RcTR7ZYpORjvnnA6\/xdq5\/wZoms\/wDC+dfi\/wCEjk8202y3Uv2cf6XHlD5ZGfl6gZ56V2Og+BPEupeNrXxZ43v7GW5sY9lnZ2Kny0PPzMSPcnvzjnjFA8E+KdL+Lt14n0e601tM1IxreJcFvMSMbN4UAY3fLwc9+fWgDg\/AGr+LtM1zxgvhnw3Fq6PqTGZpLpYvLIZ8cMRnPP5VtfDK6eT4weIrnxHbNpniW8gASxWLbF5fyMxDZO5vlU\/mcnPE+l\/D74l+GtW1m48P6xoEEGo3bTsJt7tjcxXrEccNW94P+Hmu2njWfxf4u1e2vtVMPkxJaqRGgxjPRf4eMY7k0AYPwCjTUD4m128QPq1xfbZZH5dVPzEe2WJ\/Ielek2vhbw5beL7rXba1iTW5UxO6SnO0gDlM45x1xzXD33w78VeHPFGoa34A1Wzt4tRbfdWF6p2bsk5Xg9ySOmMkVoeBvh5qml+JrzxZ4r1OHUNeuU8tfJH7uJeOnA5wAOAMDPXNAHPfAKNNQPibXbxA+rXF9tlkfl1U\/MR7ZYn8h6V6Ta+FvDlt4vutdtrWJNblTE7pKc7SAOUzjnHXHNcPffDvxV4c8UahrfgDVbO3i1Ft91YXqnZuyTleD3JI6YyRWh4G+HmqaX4mvPFnivU4dQ165Ty18kfu4l46cDnAA4AwM9c0AedeAtW8Zabr\/i9fC3h631WN9SYztLcLH5Z3PgDLLnPNaXgS61e9\/aH1a416wjsNTew\/fW0cgdUwkQX5gT1Xaeveu++HXgnUvCN\/4knv57SRNTvPtEIt2YlVy5+bKjB+YdM0yw8Dana\/GjU\/GLz2h067tVhSJXbzQQka8jbtxlD\/ABelAG34f0nxHY+INautX1lL3T7mXdYWyrg267mODx6FR+FdNXM+H7LxbbeINam13VLS60qaXOmQQqA8Kbm4b5Fydu3u3SumoA8Y+FH\/ACVn4i\/9fr\/+jpKg+MU97bfErwPNplot5exuzQW7OEEj71wu49K63wT4H1Pw1438V61eT2kltq9w0sCwuxdQZHb5gVABww6E1V+I\/gbxH4k8SaDrPh28062uNK3MDeM33twK4ARsjjvQB5\/8StU8T6xHpEXjjw1\/YugxXqtLdWpW5kGQRtDBuARn9OuMHq\/GPhrxE\/jm18ceBrixvbxrRVezmdcshHDLuI+UqR3B\/OodY+H3xK8YwQ6Z4n8S6QulCVZJVs423tj22Lnqe9aPiH4aa1ZeKbfxN4F1K2sr5LZLWa3uwSkqIoUc4P8ACqjBH8IOaAMnwx4rt774n2cfjLwf\/ZHit4ilteqzhZBtYY2k45BYBstnp2rL8Q6pNH4l1VB8YjYBbyYfZPssp8jDn93kddvT8K6nQ\/h\/4p1PxzZ+K\/G+p2UstguLS0slO1TzjJwOhJPfPHPGK7mfwZ4WuriW4uPDWjTTyuXkkksImZ2JySSVyST3oA868LNqWseDvFtnp3jP\/hK76a0WGAFGi+zs6yL1f16\/8ArnYJPG3wm8K6Tc6xpmh3GixyCKa1QbrhN7MxLP03dehYdK9gvfBulnw\/qWmaLBDoUl7GFNzpsKwOGHKk7MZx6ehI7153P8OfH\/AIi06w8N+Jtc019Bs5Vdp4N7XE6rkAEsvXBxk\/X5qAMS91DWtS\/aHvpPDsNlc3g0+MWcl8zCK3jaGNy528\/xsMDu9db4D1WfV18Y6KNM0\/RfF0JdLm7skxHNIwcJL3PDEnn+9nuRVnxV4A1qPxbZ+LPBV3Z2uow24tpba6B8qWMDaOgPbAx\/sjmneGfBmr+FtJ8U67qWr2Y8SapFJK10Plt7VgrFTlh90E5JK9B0PcA7LwtY6vpvhy1tNd1BdQ1KPf51yowHy7Few6KQPwrzX9oOBbnQ\/D0DkhJdTCMR1wVIruvh3qGrap4D0y+1yfz9QmEjPL5Yj8xfMbYwUBeCm0jgcVm\/EzwXqPjO00eLTprSJrK+W5k+0MygqBjA2qeaAPI7vVrzwt4I8W\/DbXJCZrRBLpkzdJYjIrFR+HzD\/gQ7V2mt+Drnxd8JPBY069gt9TsrS1mtknYBJT5S\/L9emOPbvW38V\/hk3j6xtZtPkt7fVrVtqyzkhHiPVWKgng8jj19ai8R\/DC61vwP4ZsIdQjtNe0GCJbe5QsY96oobnGcZRSDjPHSgDitf8V6zCLG1+K\/gWO406OcbL22ZlCMQRnKsVY4z8uV6dOKm8W6rqVl8frHUPD+mjVrptLUxQCTYJFZX+bP05rU1TwF8SfG0VtpfizW9Kh0mGVZJTZofMmI4zjaB6+g56GumHgK+i+Ldh4nt5bRNKtLAWghLt5owjKMDbjHI70Aafg7X\/FWsz3aeIvDH9jRxophfzxJ5hOcjj04rx3+15\/8AouJ\/8BJa+ka5\/wD4QTwf\/wBCpof\/AILof\/iaAOQ8M6RJ4w+Hl7pr+ObrVlmvdsmpWytE6oFQmEbv\/wBXzVz\/AML9HstA+M3i3StPQx2lraokasxY4ynUn3r2LTdK07R7ZrfTNPtbGBm3mO2hWJS2AMkKBzwOfauP8O+C9R0n4oeI\/E1xPaNZanGqQxo7GRSNv3gVAH3T0JoA2bTS\/EMXju\/1O41hJdAltwltp4X5opMJls4\/2X7\/AMVaus2Emq6NeWEV5NZyXETRrcQnDx5\/iX3rldP1PXpPjBq+myajDcaFFp6ypbRIpNtKfLAEjbchm\/eMF3H5Tn0ruaAPBtA8Maf4R\/aD03StOMrRrpjSPJM+55HZXyzH1PtXvNcHP4K1KT4yW3i8TWv9nxWRtzGXbzd21hnG3GPm\/vV3lAHzB4o0LU9U+KfjLVNGldNT0Vo7+FUHLbQuce4HPvjHeuj03xbbeNPjF4C1eDasj6fOlxEDnypRHNuX+o9iK9A8O+C9S0j4neI\/EtxNatZ6nGqwxozGRcbfvAqB\/D2JrndO+EcugfGCDxRp1zYxaN5kkn2VmZZUZ43UqgC7Su5sjkccduQDmodd1fw\/8cfGN1o\/h2bXJDGokgim8tkXCHd91s+mAO9W\/B0198Uviha+MXtrPTLbR08p7dJg87thtu4YB\/iPJA+7gV3fh7wTqWk\/FLxF4nnntGstSiCQxxsxkU\/J94FcD7p6E1n3fw61bTPijF4u8LXFlBb3PGpWdw7IJcn59u1T1+9z\/EM96AMPWvDPjPwx8QNW8UeCVsdVjvwPtVpMyl4z3ByVOMjIwc9scVb+HHiTStU8e6nHqPhZtC8XyQk3DFnxOny7vlb7p4U9Oeual1L4feLdB8Zal4i8C6pYxpqbeZd2V6DtLk5JGAc8knsRkjmr3gr4f65Z+L7rxh4u1K2vNYmi8qOO1UiOJcAdcDnAxjHryaALvxh1vUvD3w+uNQ0q7e1ulniUSoBkAtz1rifEt34\/8JeHtN8bzeKPthkaI3WmGALAiuMhR6\/3c9ec5rp\/j1\/yS26\/6+Yf\/Qq5a+8I+LfEHh\/QtF8QeJdGg8NKkMyzZ8u5lUL8qkHgsAcdfc5NAG9418X63qvijw54Q8L3g02XWLUXkt6yBnSIhjhfQ4Rj69ORVew1fxP4D+JOl+Gdb1x9d0vV0PkXE0e2WKTkY755wOv8Xatrxn4AvNY1LRfEXhPULex1jSohFAZcmKSLnCkjP95h0OQxqtoPgTxJqPja18WeN7+xlubGPZZ2dip8tDz8zEj3J78454xQBy3gS48c+OtT1It4ums7HS9RGVWEM8wLcpkYwAq+\/LVX8GaJrP8AwvnX4v8AhI5PNtNst1L9nH+lx5Q+WRn5eoGeeleifDTwXqXgy31mPUZrSU3181xH9ndmwpHRtyjmsweCfFOl\/F268T6Pdaa2makY1vEuC3mJGNm8KAMbvl4Oe\/PrQB6fSUVWur+0sUDXd1BbqehmkCA\/nQBaopkcqTRrJG6ujDKspyDT6ACiiigApKKpX+rWGmKrX15Bbhvu+bIF3fTPWgEm9i7S1Xtby3vYVmtZ45om6PGwYfmKnoC1txaKKKACkoqnf6pY6XEJL67ht0JwDK4XJ9s9aASb2LlFVbLUbPUoPOsrmK4izjdE4YZ9OKtUA01uLRRRQAlFFVL7UrPTIRNe3UVvHnAaVwoJ9B6mgEm9i3RVSx1Ky1KHzbK6huEHVonDY\/KrdANNbi0UUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAed\/EfRfMjt9XnvrmXyr61W2td22KHMihmwPvMeeT0zWl4mRb3x34SsJ1D2ubq6aNhlWkjRQhI9txNanizRbnXdIjtLV4kkW6gnJlJA2pIrHoDzgVH4k0S9v7rTNU0qWBNS02VmiW4z5ciOu10YjkZGOeeRQBzHiW4fQ\/F2uT6eDC9x4amuZPLGMzRMQkh9wGIzTdf0y1034L201rEkc9hbW93byqvzLMCjF8+p5z9a3bLw1falqOpan4j+zLNeWX9npbWrsyxQHJbLMBlmJ9OwrNPhfxLfaLa+GNRn046PCY0lu42fzriGNgVTYRhWO1QTuP40Ad8jb0VsEZGcHtWNL4glF\/d2lpomo3v2SRYpZYGgVQ5RXx+8lU\/dde3etusTQ\/+Qx4m\/wCwkn\/pJb0AaVvczTQJI9lPAzDJikZCy+x2sR+RqXzH\/wCeMn5r\/jUlFAEfmP8A88ZPzX\/GjzH\/AOeMn5r\/AI1JRQBH5j\/88ZPzX\/GjzH\/54yfmv+NSUUAR+Y\/\/ADxk\/Nf8aPMf\/njJ+a\/41JRQBH5j\/wDPGT81\/wAaPMf\/AJ4yfmv+NSUUAR+Y\/wDzxk\/Nf8aPMf8A54yfmv8AjUlFAEfmP\/zxk\/Nf8aPMf\/njJ+a\/41JRQBH5j\/8APGT81\/xo8x\/+eMn5r\/jUlFAEfmP\/AM8ZPzX\/ABo8x\/8AnjJ+a\/41JRQBH5j\/APPGT81\/xo8x\/wDnjJ+a\/wCNSUUAR+Y\/\/PGT81\/xo8x\/+eMn5r\/jUlFAEfmP\/wA8ZPzX\/GjzH\/54yfmv+NSUUAR+Y\/8Azxk\/Nf8AGjzH\/wCeMn5r\/jUlFAEfmP8A88ZPzX\/GjzH\/AOeMn5r\/AI1JRQBH5j\/88ZPzX\/Gm+cwYL5L5Iz\/D\/jU1Rn\/j5T\/cb+YoAPMf\/njJ+a\/40eY\/\/PGT81\/xqSigCPzH\/wCeMn5r\/jR5j\/8APGT81\/xqSjNAEXmP\/wA8X\/Nf8aTzjnHlPn0yv+NZXibWV0bRZ50dTcHCRLvAJLMFyM+mc\/hXNLp0clt5hnae7YcvHp7Omf8Aro2GJ9w4+g6VDkbQo8y5md35j\/8APF\/zX\/GjzG\/54yfmv+NcfpXi947Bk1GNVuIJGid5bmKMZHT7zA9COcf1q5\/wl3mH9wtjIvqly0uf++Iz\/Oqu+zE6Ek7HSeY\/\/PGT81\/xo8x\/+eL\/AJr\/AI1zX9uahMf3QkB9Bpdw2fozbBQbrXJhhIdUx\/ejit4\/0kcn9KNewvZ92jfurwWlrLcSxSCONSzYwTx6DNc1\/wAJFfzxtcRbxGDwIrJ5owPQyKfm9yowPeqXiCz1ufQL3db3Mn7slvNvFBIHJwsaYPApltva2jT7JDbfIAkTyXNxgY424wCPoaTT6\/mb06cVG+7+86vStXXVbTzo4wSrFW8uQMuevDdwQQfxq877kIeFipGCCVxj864Pw\/oN5enUJ9umQxm8ZVV7RpeiqDg+bgfMCCOfmDc1tDwgC2WksAfWPTowf\/Hi1O0l2\/r5Gc4U1J2ZvteRp97C\/V1H9aibVrNes0Y+sqf41kr4QjXpfyr7R2tsB\/6KqYeFoM5e+vGPqCif+gqKLMm1PuPvfE+l2dlPcNc27GKNn2CdMtgZwOa5i3vbjU4ftd3r1vbzSDckMcxPl+3yyKPwwceproL3wjZ3dhPb\/ab9XljZAxvJcKSOpUNg\/Q1z8Ei2duLbVdHu5r9V2nE80gkPqOuc+2foKHH+rG1Lks+Xc1dB8VpLbzW+qTRC6tpNhcOuJFwCGB4rXPiPTFGTcLj1yMfzrF8O+FoWhubzV9OgE9zLvSBhnykAAA+vUn6itoeGNBDBv7F0\/PqbZM\/yos+j\/r7yKjpcz0Im8WaOn3rtBjqcjj9aibxt4fTO7UYRj\/aH+NX10DR0xs0mxXHpboP6VMumWCfdsrde3EQ\/wos+5F6XZ\/eYc\/j\/AMOxW8ki38bsqlgoYZbA6Vz1rqnhvUIkudY1e4ur2TDsIp5EjiP91FUjgdMnk16BLZW00EkLQp5cilGAXGQeKwof+Eh0iNbKKxh1KCPCw3DXAifYBxvBByR6jrVLYuMofZTT9TK0Txfp1jJd2NxfTXMcbK8NwYWLurDo2ByRjGfTHpWz\/wAJpo\/AzdknsLSQ\/wDstWtF065tjcXmoSRvf3RBl8ofIirwqLnnAyevcmtfAodiJyhfb8TAHjDTCMrFqDfSxm\/+Jpw8V2R6WWrH6abN\/wDE1u4FFGhF49jBk8VW6RM66dq5wCRnT5Rn81rGs7nTr2zhu9S0W8v7qeNXkklsjIoz\/Cm7gLzgY69Tk812xAIwawv7P1LTYvI0+4tpLQOPLiuFIMQz0DL1X0GOPWqTQ00ZumawNP1C40+DTtUa0WJJIomhJaHJYFRk\/d+UY\/4EPYbH9vn\/AKA+qf8Afgf41Pp2nvaNLPcT\/aLufHmy7dowOiqvZRk8c9Tya0aTauDauY\/9uydRouqY\/wCuSf8AxVH9uTnldD1Mj\/djH83rYopXXYm6Oa1R4dbsTZal4Wvby2YhjFMkJUkdOslUfD8xXS4bn\/hG72aWeMFpNsH3SOEGXBCqMKBgdPXNdnWOdJvbWV\/7Kv4raB2LtBNb+aisTklcMpXJOcZI9qpNWsUmjJs9QmsNclt7bQr+OGeHzjbKYQFcNguB5mBnIz7j61r\/ANs3n\/Qv6l\/31B\/8cqew0xbOSWeWZ7i7mwJJ3ABIGcKAOAoyePfua0KTauDauZH9s3n\/AEL+p\/8AfUH\/AMco\/tm8\/wChf1P\/AL6g\/wDjla9FK6Juuxkf2zef9C\/qf\/fUH\/xysXR7uSeN9QuPD17c3U0j7pT5LbQGICLukyAvTHrk967HrWRJpd3BPLJpl8luJn3vDND5ke49WUAqQT164zzjk5pNbFJoxYNQuNP15Y7bQdSjiuoZJHtg0AG9WT51HmYH3\/m9crWudbv+3hvU\/wAXt\/8A47VnT9MNpNLdXNwbq9mAV5iu0BR0VV\/hXknuT3JrRpSkgcl2MQ63qHbw3qWf+ukH\/wAcpP7c1P8A6FnUf+\/sH\/xytyikLmXYwDrmrf8AQsX\/AP3+g\/8Ai6wdIvb6e5u9Sl8MXF3dvcSR+aZoT5So5URrluMY5x1OTzXe1h3OhXS38l5pWpvYvMd00bRCWJz\/AHtpIw3HUGi5pCcddEvvOe+3XuneJoJbTwzc28l5G6y26TQgTbcENgNgEZxn\/arc\/t3Wf+hWvf8AwIh\/+Lq3pui\/Y7l766upb2\/kQI08gChV67UUcKM89\/rWsKGwnUi3tf7zn\/7e1n\/oVb3\/AMCIf\/i6P7e1n\/oVb3\/wIh\/+LroaKRPPH+Vfic7\/AG7rP\/Qq3v8A4EQ\/\/F1haZf3txqt\/qUvhm5urtZ2gDGaH9wqgfIMtx6kjrmu+rEvdDmN9LfaZqL2FzMAJh5QkjkwMAlTj5u2QRx+FNMqNSOqsl95z0t9f2XiGzubXwzc20l0XjmhSWHFwApYEgN1GPvHscdwK2xrus\/9Cte\/+BEP\/wAXVnTdEe1vDf399LfXxTyxI6hEjUnJCIOBnAyeTxWxQ2Eqkeyf3nP\/ANvaz\/0Kt7\/4EQ\/\/ABdH9u6z\/wBCre\/+BEP\/AMXXQ0UieeP8q\/E53+3dZ\/6FW9\/8CIf\/AIqsG3vLu61zUdRuvDFzdXEDrDGrSwkWyhFbAy3UliSR6gdq7+uR1y3uJdXnbQZHi1IQgXbnHkFf4RID1bGcYGcdeMVUVfQaqwW6tf1Luo2sGn67pl\/aIsdxcz\/Z51XjzkKM3I7kbc5+vrXQ1yuj295FrML+IJHm1No3Fu6EfZwONwQcENjGdwzgccZrqqJJrczclLZ3sLRRRUiCiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKxND\/5DHib\/sJJ\/wCklvW3WJof\/IY8Tf8AYST\/ANJLegDbooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooASlpKzptYtllaGDfdTqcGOBd2D6E\/dX8SKTaW41FvYfq2pw6Rps17P8AcjxxnGSTgCuel1XW4iLyS2nWAIWJaFAAvGTt378f+Pf7Paq3i6\/1A6ZGrpDFMsiTrZR\/vp3CMCWHZduM5ww96qKUvFW8uTcwsEZ0VYfPbbkHPnSgxgHj7uFH96lvq9jqp00o3ep09j4o0+9sYrhHZnYDdHAjTFTnB+6D+fpzUFx4shicxCFEl7JcXCIT\/wABUs\/\/AI7VLwv4ctZ9FWfUYDObiWSVYpZWePYzkqdh+U5GDnaDzzXU21nbWcfl2tvFAn92NAo\/Sny92ZT9nGTSRzg1bWr1Q1paT7T3S22Y\/wCBTFf0Q\/Sj+ydcvD\/pE6RrnIM1w8hx6FIxGv6kc11VFFl2F7VrZWORvPBZuNOnijvmSd8OnlxJEnmKdykkDf1HdjWYloIo2t7jRWurxcjfc2Uk7MRkZ3tlD9d6\/Ra9Bop9LDVeXXUwfDOgLo9gRNHEbqWRpXZVHy56AEDsMe3pit6ilpt3MpScndhRRRSEIQCMHpWOfDdhysZuoYWJJhhuZEjOTk\/Kp45z061s0lNNrYabWxFbW8NpbxwW8SRQxrtREGAo9hUtFFIQtFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABUc33B\/vr\/6EKkqOb7g\/wB9f\/QhQBJRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRmgBKKy7jW4I5jb2sct7cg4MVuuQp\/2mOFX6E5qP7JqeoqDe3P2KE4Jt7RssfZpMZ\/75APvVcvV6Ec62WpZvdVtLBgkrs8xGVgiQySN9FHOPfpXOxXd59t1K2dk0zz5DMHn2tI\/yKNo5KZAUdycY4rft4LDSg0FlbjzWO5lTLOx9WY\/zY1Uu7cyx3EeoIq20mX2bQwJx0J\/Af41hUrRj7sfmyo05S1ZFc3BvNWsLMTQyyW9wJZJYuiYRxtI5wx54z610VZkMUDNBDYwJFbRNvJRNq9MAD1\/+tWnmnCV7pbLYfLbV7i0UUVoAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABWJof\/IY8Tf8AYST\/ANJLetusTQ\/+Qx4m\/wCwkn\/pJb0AbdFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQBBPNHbQPNM6pGgyWPaqDayyDzZNOvUtv+ezIvA9SoO4flTPEMsdva21xcHFrDcI8\/GflGcHHfDbT+FNP2u\/ga5muG0+y25CLtDle5djnb9ByPXsMpSd7I2hBcqbL1zqljawxzTXUapJjyyDnfnptA6\/hVf7bfXfFnZGJO0118ox6hB8x+h21B4csLa101Z47dUaVnIdl+cxliUyTz93b1pz6xJeTPb6PCtyyghrlziBG9Nw5Y+y\/iRVRUpq+wmoxdkr+o27s7aGF7rW79pYQMMjny4R7bB1z6MWpkTX18I4tOg\/szTgP9a8YErD0SM\/c+rDP+z3q1a6MiXC3l7M17ej7ssgwsfsidF+vX1JrUq1GMdtyXMpWGmWunI3kRkySYMszsWkkPqzHk\/wBO1VT4e0j7Yr\/2ba5OXx5QxuyDnHTPvWvUZ\/4+E\/3G\/mKd2TzMlooopCCiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAqOb7g\/31\/wDQhUlRzfcH++v\/AKEKAJKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigBKwibjVJ7qQ3slpYW7mMCLCtIV+8xc5woORgY6HJ5xWrd3cFjbtPcyrHEpAy3cngAepJ7Vx622pXd1LO9jM2mtcGVbQbCyuccuCy55ydvOCTnJ+7Sairsym7uyLcWt3NvqENpYtNqkFyxSKefEaI4UtgPj512qxyoOMYzzxqNo817k6reSTIR\/x7QkxRD64O5vxOPYVSS5uLvVLWS7t2Q25cxQLt3lyNoZsMQvyl8DPc\/Strybm45nk8pP8AnnEeT9W6\/lioliE9Ka16jhRf23p0ASWtki2ttEoKDCwQKBtH06D8aXyp7j\/XP5Sf884zyfq3+H51PFDHCm2NFVfQCpKjlctZM1ulsRxQpCu2NAo64Aon\/wCPeX\/dP8qkpk\/\/AB7y\/wC6f5VaSWiESUUUUwCiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAoorjfF8txdeIfDugfbZrOy1B52uJIJDG8nlqCsYYcjOTnHPFAHZUV51Nq\/8AwhviHVtLt7qe5s49Gk1OKC5maZoZEJG0MxLbWHOCe3vVS9srzQfA9t4wXVL+bWIkhu7nzLlzFOrld8ZjztC4Y4wBjAoA9QrE0P8A5DHib\/sJJ\/6SW9bKsGUMDkHkGuUstf0bSte8SQajq9hZzNqCOI7i5SNiv2W3GcMenB\/KgDraKrQX9ndwJPbXcE0LjKSRyqysPYipfPh\/56x\/99CgCSio\/Ph\/56x\/99Cjz4f+esf\/AH0KAJKKj8+H\/nrH\/wB9Cjz4f+esf\/fQoAkoqPz4f+esf\/fQo8+H\/nrH\/wB9CgCSio\/Ph\/56x\/8AfQo8+H\/nrH\/30KAJKKj8+H\/nrH\/30KPPh\/56x\/8AfQoAkoqPz4f+esf\/AH0KPPh\/56x\/99CgCSio\/Ph\/56x\/99Cjz4f+esf\/AH0KAJKKj8+H\/nrH\/wB9Cjz4f+esf\/fQoAkoqPz4f+esf\/fQo8+H\/nrH\/wB9CgB9Vb29hsIPNmJ5IVUUZZ2PQAdzUN\/qtvZooBE08hxFCjDc5\/oPUngVHaW6Cc3t7PFLdkYUBspCv91P6nqfpgCHLoi1FJXkZ2pLeMkM86Ibu4lWG0t3O5ICeS7D+JgoY\/hgdyaF9oem2CJb28t9dakoDxwo+8ZzwzRnEarkdSFHpzWn4huRMtvZ2Mkban5izQEt8sW08u\/ouMr6nOB7QWM2oWlu8MWmp9tmO6a8nuUMTvjG84O89Bhdo4wOK1jBJXNFOVtyTT7S41+0WfWZFMYdkaxhBWIMrFTvPWTkHg4X2PWuiRFjQIihVUYAAwAKpadDb6fYx2wuEcjLO5IBd2JZm\/FiT+NW\/Ph\/56p\/30KUnqZSd2S0VH58P\/PWP\/voUefD\/wA9Y\/8AvoUiSSoz\/wAfKf7jfzFHnw\/89Y\/++hUZnh+0J+9T7jfxD1FAFiio\/Ph\/56x\/99Cjz4f+esf\/AH0KAJKKj8+H\/nrH\/wB9Cjz4f+esf\/fQoAkoqPz4f+esf\/fQo8+H\/nrH\/wB9CgCSio\/Ph\/56x\/8AfQo8+H\/nrH\/30KAJKKj8+H\/nrH\/30KPPh\/56x\/8AfQoAkoqPz4f+esf\/AH0KPPh\/56x\/99CgCSio\/Ph\/56x\/99Cjz4f+esf\/AH0KAJKKj8+H\/nrH\/wB9Cjz4f+esf\/fQoAkoqPz4f+esf\/fQo8+H\/nrH\/wB9CgCSio\/Ph\/56x\/8AfQo8+H\/nrH\/30KAJKKj8+H\/nrH\/30KPPh\/56x\/8AfQoAkoqPz4f+esf\/AH0KPPh\/56x\/99CgCSio\/Ph\/56x\/99Cjz4f+esf\/AH0KAJKKj8+H\/nrH\/wB9Cjz4f+esf\/fQoAkqOb7g\/wB9f\/QhR58P\/PWP\/voVHNPDsH71Pvr\/ABD+8KALFFR+fD\/z1j\/76FHnw\/8APWP\/AL6FAElFR+fD\/wA9Y\/8AvoUefD\/z1j\/76FAElFR+fD\/z1j\/76FHnw\/8APWP\/AL6FAElFR+fD\/wA9Y\/8AvoUefD\/z1j\/76FAElFR+fD\/z1j\/76FHnw\/8APWP\/AL6FAElFR+fD\/wA9Y\/8AvoUefD\/z1j\/76FAElFR+fD\/z1j\/76FHnw\/8APWP\/AL6FAElFR+fD\/wA9Y\/8AvoUefD\/z1j\/76FAElFR+fD\/z1j\/76FHnw\/8APWP\/AL6FAElFR+fD\/wA9Y\/8AvoUfaIf+esf\/AH0KAH1S1DUodOjUuGklkO2KGMZeVvQD+p4HeoL\/AFmK3kW2tttxeyDKRBsBR\/edv4V\/n2yaq2qR2tw8zTLeapKoWSYnaiDrtH91fYZJ7560NxgryIu5O0SWGzbzhqWryI06n9zCvKQcY+UfxOc\/e\/AY7gMgdkuPNgtXZnUgctk5IJH3f889auRJErebNcJLN\/eJAC+wHarJnh\/56x\/99CsJxnVs3p5GsOWGi1KUSJLPD9ni8u3hJbO3aGJBGAPxPNaNR+fD\/wA9Y\/8AvoUefD\/z1T\/voVdOHIhN3JaKj8+H\/nrH\/wB9Cjz4f+esf\/fQrQRJUc\/\/AB7y\/wC4f5UefD\/z1j\/76FRzzxGCT96n3T\/EPSgCxRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFcr42igmtLJb\/AEOTVNN87NwbcMZrbj5ZEC\/N14O3nBrqqKAPMdN8Kwazf6vc2Gm3GnadLpUunwSXiuJZ5JfvSnfl8ABQC1JeXWo674JtvB\/9jajBqsiQ2l08luywwqhXfJ5n3WUhTjBPWvT6KAEChVAAwBwBWLof\/IY8Tf8AYST\/ANJLetusTQ\/+Qx4m\/wCwkn\/pJb0AbdFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUm4HIyOKACqF7qHkSpa28fn3kgysWcBV\/vMew\/n2zTLvUJHuDY2CiS5xmRz9yEHoW9T6L39hzU9lYx2UZVWaSVzuklc5aRvU\/4dB2qG29EWkoq7GWOni3ZriZ\/OvJP9ZMR\/wCOqP4VHp+eTzUWp6k9s6WllGs2oTAmONidqDu7nso\/XoKfqWotaBILeMT303EMOcZ9WY9lHc\/h1IpdO01bESSO\/nXc5DXE5GC57ADso6Advckk3GKihN31YmmaWmnRuzSNPdTHfcXDj5pW\/oB0A6AVo0UUNkt3CiiigAooooAKjP8Ax8p\/uN\/MVJUZ\/wCPlP8Acb+YoAkooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACo5vuD\/AH1\/9CFSVHN9wf76\/wDoQoAkooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKAEzRRWFNe6vf3ksWlrbwW0D+W89wpYyMPvBVBGAOmT1\/mrpbiZu1jX+quZZLPTynnR4E9xIMx24Pr\/ebHRR6jOO+NBrOq6pKdNWSD5sj7bakr5gHUJuBA7\/ADAn2GenQ2GkwWUca7FOz7igfKnuPU+rHk1POvs6v+v6sKzl5IhsdOEUbrF5iLK26aeQ\/vpz6k9vT2HAC4rTjijhj2RqFUdhUlHahR15nqytlZbC0UUVYBRRRQAUUUUAFRz\/APHvL\/uH+VSVHP8A8e8v+4f5UASUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFYmh\/wDIY8Tf9hJP\/SS3rbrE0P8A5DHib\/sJJ\/6SW9AG3RRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAZmqySE2lpFI0f2qbY7qcEKFZjg9iduPxrLuNLtJ7s22l20cM0RBmvk4eM+gPVnI9eMHnOcHQ105tYooi32x5V+z7eocdz\/sgZz7Z9agtdO1WC0WzFzAi5Je5UEyPnktg8BifqOenFYSV5WaudEHyxTv8A13LWhKg0a3KqFYjLkHO5+hOTyckdTT9S1E2axxQReffT5EEAON2OrMf4VGeW\/DkkAsu7uLR7OC2t4TLOw8u2tlPLke\/YDqSelP0\/T2tnkublxNezAebKBxgdFX0UZOB9T1JreCtFXMpNOTkJpemfY\/MuLiT7RfT4M05GM+iqP4VHYfickk1pUUUN3JbuFFFFAgooooAKKKKACoz\/AMfKf7jfzFSVGf8Aj5T\/AHG\/mKAJKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAqOb7g\/wB9f\/QhUlRzfcH++v8A6EKAJKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigBtcdf3yRzSaLdyfZ7PzWMkyZJdD82zjODzzntj1NdRdzuiiGDb9ok+4DyB6k+wogsYILUQbfMXqxk+Yse5PqaxleUvd6f1YdtNTI0me1v7uL+zbZk0+1QhZShQOx4woIzgAdfeuhpqqqgBQAB0Ap1XGPKK4tFFFWAUUUUAFFFFABRRRQAVHP\/AMe8v+4f5VJUc\/8Ax7y\/7h\/lQBJRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAViaH\/wAhjxN\/2Ek\/9JLetusTQ\/8AkMeJv+wkn\/pJb0AbdFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFACUySRIo2kdgqKMlicAD1p9ZN5nU73+z0J+zxYe6PZu6x\/j1PtgfxVMnZFRjd67FP7WU\/wCJvLA0txcN5FhbgYbaeR1+7uxuJ7KB6Uy+1PXNLRJZksbqaUkRWEAcSOcc4c5zjqSVA+lXddDwRWmoIoZLCYzyKeP3ex1Yj3AYn8Ky9M1mxmkn1N2FxqdwCsVpEuZYogflTH8JPDMTgZPJworSnC0e5pvrY0tBtVkhXVpphc3l3GrGUDCoh5EaA8qo\/Mnk+21VDRrSSx0m3t5iplVcybegYnJA9smr9KW5lLcWiiikIKKKKACiiigAooooAKjP\/Hyn+438xUlRn\/j5T\/cb+YoAkooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACo5vuD\/fX\/wBCFSVHN9wf76\/+hCgCSiiigAooooAKKKKACiiigAooooAKKKKACiiigBKrXt3FY2zzyk7V6Aclj2AHc1O7rGjO7BVUZJPYVg6eX16+GpyqRYQnFojf8tG7yEen938T6Gok3stw2NKxgkw11cr\/AKTL1Gf9WvZPw7+pq\/SUtOMeVWB6hRRRVAFFFFABRRRQAUUUUAFFFFABUc\/\/AB7y\/wC4f5VJUc\/\/AB7y\/wC4f5UASUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABXN+KtbvtMk0uw01YBeanc\/Z0muATHEApYsQCMnjgZrpK4Xx2sN9rPh3R9TlEOjXk0rTyH5d0qKDGm\/8Agyd3IwTjAoAtWHiW+07U9Z0zXpLe4Om2YvvtdrGUDRfNkMhY7WG315FZ7eJfE1hoVp4q1EWH9lTGOSeyjiYSwQyEBWEm75mG4EjaO9ZU+mRWuseJfD3hxjcWt9oks06eZ5rR3HKIDIcsdw\/hYn1qxr2rWWqfBy2tbOeOW7v7e3tIIFYbzKSilceowc\/SgD02sTQ\/+Qx4m\/7CSf8ApJb1souxFXJOBjJ71yllolrqOv8AiSaeW\/RhqCKBb6hPAuPstv8Awxuozz160AdbRVa3soraBIUedkQYBkuJJG\/FmYk\/ial8lfWT\/v43+NAElFR+SvrJ\/wB\/G\/xo8lfWT\/v43+NAElFR+SvrJ\/38b\/GjyV9ZP+\/jf40ASUVH5K+sn\/fxv8aPJX1k\/wC\/jf40ASUVH5K+sn\/fxv8AGjyV9ZP+\/jf40ASUVH5K+sn\/AH8b\/GjyV9ZP+\/jf40ASUVH5K+sn\/fxv8aa0caKWZnCgZJMjcfrQBX1K8NpABEFe5lYRwRkn5nP9ByT7A06wsxY2wj3F3JLSSHq7nkk\/546VS06H7fcNqcnmiNhstUZzwnduvVv5Ae9WNTuY9OsZLjbLK4wscSyHMjk4VRz3OKiK5nf7jR6e6vmVL8jVtVj0pCTbQbZ7wjoe6R\/iRuPsAP4q2wMVm6VphsrILPIZLqQmS4kRmAeQ9cc9BwB7AVe8lc9ZP+\/jf41o+yIfZEtFR+SvrJ\/38b\/GjyV9ZP8Av43+NIRJRUfkr6yf9\/G\/xo8lfWT\/AL+N\/jQBJRUfkr6yf9\/G\/wAaPJX1k\/7+N\/jQBJRUfkr6yf8Afxv8aPJX1k\/7+N\/jQBJRUfkr6yf9\/G\/xo8lfWT\/v43+NAElRn\/j5T\/cb+Yo8lfWT\/v43+NRmFPtCcyfcb\/lo3qPegCxRUfkr6yf9\/G\/xo8lfWT\/v43+NAElFR+SvrJ\/38b\/GjyV9ZP8Av43+NAElFR+SvrJ\/38b\/ABo8lfWT\/v43+NAElFR+SvrJ\/wB\/G\/xo8lfWT\/v43+NAElFR+SvrJ\/38b\/GjyV9ZP+\/jf40ASUVH5K+sn\/fxv8aPJX1k\/wC\/jf40ASUVH5K+sn\/fxv8AGjyV9ZP+\/jf40ASUVH5K+sn\/AH8b\/GjyV9ZP+\/jf40ASUVH5K+sn\/fxv8aPJX1k\/7+N\/jQBJRUfkr6yf9\/G\/xo8lfWT\/AL+N\/jQBJRUfkr6yf9\/G\/wAaPJX1k\/7+N\/jQBJRUfkr6yf8Afxv8aPJX1k\/7+N\/jQBJRUfkr6yf9\/G\/xo8lfWT\/v43+NAElFR+SvrJ\/38b\/GjyV9ZP8Av43+NAElRzfcH++v\/oQo8lfWT\/v43+NRzQpsHMn31\/5aN\/eHvQBYoqPyV9ZP+\/jf40eSvrJ\/38b\/ABoAkoqPyV9ZP+\/jf40eSvrJ\/wB\/G\/xoAkoqPyV9ZP8Av43+NHkr6yf9\/G\/xoAkoqPyV9ZP+\/jf40eSvrJ\/38b\/GgCSio\/JX1k\/7+N\/jR5K+sn\/fxv8AGgCSio\/JT1k\/7+N\/jR5K+sn\/AH8b\/GgB9Hao\/JX1f\/v43+NYmtXrRGOyst0l5M2xVMp64z69AOT6D3IqZO3qHmRahI+v6h\/ZNtIy2kWHvJUOMjPCA+px+XP93PQxxpDGscaKkagBVUYAHoKo6ZpMGmWa26NIzZLSSFyC7HqTz\/8AqGBV3yV9ZP8Av43+NNKxK7sloqPyV9ZP+\/jf40eSvrJ\/38b\/ABplElFR+SvrJ\/38b\/GjyV9ZP+\/jf40ASUVH5K+sn\/fxv8aPJX1k\/wC\/jf40ASUVH5K+sn\/fxv8AGjyV9ZP+\/jf40ASUVH5K+sn\/AH8b\/GjyV9ZP+\/jf40ASUVH5K+sn\/fxv8aPJX1k\/7+N\/jQBJUc\/\/AB7y\/wC4f5UeSvrJ\/wB\/G\/xqOeFRbycv90\/xt6fWgCxRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFVr2wtNStWtr62hubd\/vRTIHU\/gas0UAUtO0nTtHgMGm2NvaRE5KwRhAT6nHWoYvD+jW+ptqUOlWUd82SbhYFDknqd2M5NadFABWJof\/IY8Tf8AYST\/ANJLetusTQ\/+Qx4m\/wCwkn\/pJb0AbdFFFABRRRQAUUUUAFFFFABRRRQAUUUUAJWRf51K8GlpzAoD3Z\/2e0f1bv8A7I\/2hVzULw2dsWRPMmchIox\/G56D6dyewBpNOs\/sdqFd\/Mmcl5pSMb3PU\/0A7AAVEtXYuPurmLYAUAAYA6CsaDGr6y12Tm0sGaKAdnm6O\/8AwHlB7l\/aptYuZ44orOzOLy7by42xny1\/jkP+6P1KjvV2ytIbGzhtbdNkUShVX2FaLRC2VyxRRRSJCiiigAooooAKKKKACiiigAooooAKjP8Ax8p\/uN\/MVJUZ\/wCPlP8Acb+YoAkooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACo5vuD\/AH1\/9CFSVHN9wf76\/wDoQoAkooooAKKKKACiiigAooooAKKKKAG0MwUZJAHvS1zX9opPJI0enveuJNjO5CxoSeEXPVhxnA9cntUTk1sCt1NbUtQjsbdnLqpClix5CKOrH\/DucCqui2DrnULqNkuZlwkbnJhj67T\/ALRPLH146KKxtHlh1bWDG0MkUdsfMMJI2vIrEBgBkFAc7efvAnGRXZUQvvJWf5CdparYWiiirGFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFRz\/8AHvL\/ALh\/lUlRz\/8AHvL\/ALh\/lQBJRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAViaH\/yGPE3\/AGEk\/wDSS3rbrE0P\/kMeJv8AsJJ\/6SW9AG3RRRQAUUUUAFFFFABRRRQAUUUUAJSEgDPYVR1G5niMFva7BcTvtVnGQigZZiO\/09SKxdRtdQmn+wW+o3FwdnmXKyBAhXsmVUEFsEcds+1ZynbbU0hT5t3Y0rDOo3Z1N8GAApaD\/Z7v\/wAC7e31NajukcbO7BUUZLMcAD1qKxuI7qygniUrHIisqkdAR0rN1Q\/2nfR6MgJhKia9OOPLzxGf98g\/8BVvWrgtBS1lbsGjRteTy6zMGzcAJbIwx5cA5HHqx+Y\/8BH8NbdJilqm7kthRRRSEFFFFABRRRQAUUUUAFFFFABRRRQAVGf+PlP9xv5ipKjP\/Hyn+438xQBJRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFRzfcH++v\/AKEKkqOb7g\/31\/8AQhQBJRRRQAUUUUAFFFFABRRRQAUUUUAIelcPb\/aZ5Bosdw6Tq0vmHcCY03EM49NwIC+7t\/drqNY1FNL06S5bcW+6iqu4sx6ADufas\/TvDdq9mH1azgubuU75POUSBPRBnsBx7nJ71DjGbs+n9f18iJNr4RNPjT+34o7R1MFnatDKIhtijLFNiKPUBWJ7jI9RXQ96iggitoVhgiSKJRhURQoUewFS1ei2HFPqLRRRQUFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFRz\/8AHvL\/ALh\/lUlRz\/8AHvL\/ALh\/lQBJRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFY+t63Joxg8vR9T1Hzd2fsMSv5eMfe3MOuePoaANiiuW0fxtb6trF1pjaTqlhNaw+fO15GiJGvbcQ5wT1HsDUNv49tpvs1zJpWoQaTdyiG31KVVEbMThSV3blVj0YjuKAOvrE0P\/AJDHib\/sJJ\/6SW9bdYmh\/wDIY8Tf9hJP\/SS3oA26KKKACiiigAooooAKKKKACiiigDK1gPDDHqEZQPZ5dt5wGjx8wz245+oFUrCa\/NpJ5djKl\/csZJJJgFSMngZ5ycDAwPTtVzUcXt5b6aOUOJ5\/9xTwPxbH4A1q1ly80m0zbm5YJNGdJJBoOiAncYraNUQfxOeFVR7k4H1NLpFlLaWrPcsHvLhzLcMOm4\/wj\/ZUYUewqoANX1vcebPTn49JLjHP1CA\/99E91rbra3KrGbYtFFFIkKKKKACiiigAooooAKKKKACiiigAooooAKjP\/Hyn+438xUlRn\/j5T\/cb+YoAkooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACo5vuD\/fX\/wBCFSVHN9wf76\/+hCgCSiiigAooooAKKKKACiiigBKKWqd9K6osEJxPNlVP90d2\/D+eKmUuVXGldmWIDrOvi4kGbHT2xED0km7t9F\/n9K36htbaO0to7eIYjQYH+J96moirLUXW4tFFFUAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFRz\/APHvL\/uH+VSVHP8A8e8v+4f5UASUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABWfreqRaLol9qcwylrC0pGfvYHA\/E8VoVl6\/olt4i0W40q7eZLe427zCwVsKwbGSD1xg+xoA5ODSrjT\/hRrVzcEvq2o2FxeXcmOTI8ZO3\/gIwoHtUfinyv+FG\/u8bP7OtfKx65j24\/HFeguiPGY3UMjDaVI4I9K5S38A2MEltC+o6jNplrKJoNNllUwowOV\/h3MoPQEmgDrE3bF343Y5x61yllDrMmv+JDp1\/YW8P8AaCZW4snmYt9lt+dyypx04x+NdbWJof8AyGPE3\/YST\/0kt6ANK3S9W3QXM8EkwHzvHAyKT7KXOPzNS4m\/56J\/3wf8akooAjxN\/wA9E\/74P+NGJv8Anon\/AHwf8akooAjxN\/z0T\/vg\/wCNGJv+eif98H\/GpKKAI8Tf89E\/74P+NGJv+eif98H\/ABqSigCLbN\/z0T\/vg\/401vMVSzSRhRySUPH61MTWVrDG48jTEOGu2Ik46RLy\/wCfC\/8AAhUydlcqMbuxVtboWtncazeSJGLlgVHlEts6RqADkk9cermqmr+J7mytgn2Sa1lnby4J7iNQgJ6uQHJCqMscgdPetDV2S2v9Knn2raRysGYj5UcoQhPp1K\/VhVeKW2nl1HWdQx9kjD2sG\/keUOHIH+04x7hVrSEVFK+pd03zM1bGyNhZxW0MgKouNzKSWPck55JOSfrVnE399P8Avg\/41T0SOeHQ7CO63CdLeMSBzkhtozn3rQpPcze4zE3\/AD0T\/vg\/40Ym\/wCeif8AfB\/xqSikIjxN\/wA9E\/74P+NGJv8Anon\/AHwf8akooAjxN\/z0T\/vg\/wCNGJv+eif98H\/GpKKAI8Tf89E\/74P+NGJv+eif98H\/ABqSigCPE3\/PRP8Avg\/40Ym\/56J\/3wf8akooAjxN\/wA9E\/74P+NGJv8Anon\/AHwf8akooAjxN\/z0T\/vg\/wCNGJv+eif98H\/GpKKAI8Tf89E\/74P+NRkTfaE+dPuN\/AfUe9WKjP8Ax8p\/uN\/MUAGJv+eif98H\/GjE3\/PRP++D\/jUlFAEeJv8Anon\/AHwf8aMTf89E\/wC+D\/jUlFAEeJv+eif98H\/GjE3\/AD0T\/vg\/41JRQBHib\/non\/fB\/wAaMTf89E\/74P8AjUlFAEeJv+eif98H\/GjE3\/PRP++D\/jUlFAEeJv8Anon\/AHwf8aMTf89E\/wC+D\/jUlFAEeJv+eif98H\/GjE3\/AD0T\/vg\/41JRQBHib\/non\/fB\/wAaMTf89E\/74P8AjUlFAEeJv+eif98H\/GjE3\/PRP++D\/jUlFAEeJv8Anon\/AHwf8aMTf89E\/wC+D\/jUlFAEeJv+eif98H\/GjE3\/AD0T\/vg\/41JRQBHib\/non\/fB\/wAaMTf89E\/74P8AjUlFAEeJv+eif98H\/GjE3\/PRP++D\/jUlFAEeJv8Anon\/AHwf8aMTf89E\/wC+D\/jUlFAEeJv+eif98H\/Go5hNsHzx\/fX+A\/3h71YqOb7g\/wB9f\/QhQAYm\/wCeif8AfB\/xoxN\/z0T\/AL4P+NSUUAR4m\/56J\/3wf8aMTf8APRP++D\/jUlFAEeJv+eif98H\/ABoxN\/z0T\/vg\/wCNSUUAR4m\/56J\/3wf8aMTf89E\/74P+NSUUAVp5WtoJJ5ZokjjUu7FDgAck9aztIFzd79TnCxvcDEaMhykQ+734J6mqurO2s6tHo0JH2eJhJesD1\/iVP5E\/VfU10KqFUKoAA4AFQ1eVuwk+omJv76f98H\/GjE3\/AD0T\/vg\/41JRVjI8Tf8APRP++D\/jRib\/AJ6J\/wB8H\/GpKKAI8Tf89E\/74P8AjRib\/non\/fB\/xqSigCPE3\/PRP++D\/jRib\/non\/fB\/wAakooAjxN\/z0T\/AL4P+NGJv+eif98H\/GpKKAI8Tf8APRP++D\/jRib\/AJ6J\/wB8H\/GpKKAI8Tf89E\/74P8AjRib\/non\/fB\/xqSigCPE3\/PRP++D\/jRib\/non\/fB\/wAakooAjxN\/z0T\/AL4P+NGJv+eif98H\/GpKKAI8Tf8APRP++D\/jUc4l+zyZdPun+A+n1qxUc\/8Ax7y\/7h\/lQBJRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAViaH\/wAhjxN\/2Ek\/9JLetusTQ\/8AkMeJv+wkn\/pJb0AbdFFFABRRRQAUUUUAFFFFACVl6Z\/plzcakRlZD5UH\/XNT1\/E5P020\/V5nFqttA5W4um8mNh1XP3m\/Bcn64q7BDHbW8cEShI41CKo7ADAFRvL0LWkfUz\/EE\/laNOgRHlnxbxIyhgXc7RkHqBnJ9gai0\/wzpenJbrHA0htwPLMsjOFI\/iAJwD7gClk\/07xJFGATFp6GR\/QyuML+Sb8\/74rYrW7SsK7SsLRRRUkhRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAVGf+PlP9xv5ipKjP\/Hyn+438xQBJRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFRzfcH++v\/oQqSo5vuD\/fX\/0IUASUUUUAFFFFABRRRQAlUNVvzYWmYwj3MreXAjHAZz6+wGST2ANXiQASa522H9t6gL1iDblStuMf8sc8v\/wMjA\/2R7mk3yrmJersi9oVgtnZ7yzPJKS7OwwWJOdx9ySTjtnHatWiilFWWpWgtFFFUAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABUc\/\/HvL\/uH+VSVHP\/x7y\/7h\/lQBJRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRWfq2sWWiWYur6bYjOI41VSzyOeiIo5Zj6CuafxF4xv2\/4lXhAW8DjMdxql2sZH+9EuWH0zQB2tFcanizX9Pct4g8J3MFqDtN3YTLdKPVmjX51Xqc4NdTY3ttqNlFeWc8c9vKu6OSM5VhQBZooooAKKKKACsTQ\/+Qx4m\/7CSf8ApJb1t1iaH\/yGPE3\/AGEk\/wDSS3oA26KKKACiiigAooooASiiqWq3b2lizRAGeQiOEHu7HA\/DufYGk3ZXGld2K9pm91i4uzkw24+zw+hPWRvzAX\/gJ9a0Zpo7aCSeZwkUal3Y9AByTUdlarZWcVuhJWNcbieWPcn3PWs7XSLtrTSMZ+2SEyj\/AKYphnz7H5U\/4HRCPcptOWmxJoEMq6d9quFK3N45uZQRypYDC\/8AAVCr\/wABrVo6cUU27shu7FooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACoz\/AMfKf7jfzFSVGf8Aj5T\/AHG\/mKAJKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAqOb7g\/wB9f\/QhUlRzfcH++v8A6EKAJKKKKACiiigBKKDXLSv9r0n+2NQursRS4a2tbWZoshj+7HykFnbI6nGT04zTjG5Ep8pe1mYXcn9lrLtjKebeyDqkP936vgj6BvatCyiZIzI6bXk5K\/3B2X8B+uawdFtLq0vDp9\/dR3E0yrdyNzvJGFKn1UHZg\/pXU1EtZ2Wy\/MdP4bvdi0UUVRQUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAVHP8A8e8v+4f5VJUc\/wDx7y\/7h\/lQBJRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFAEbmMPGHZAxbEe48lsHp74z+tZfiLXBoOnxzLbNdXNxOlta26sF82V\/uruPCjg8n0pnipdBOiSSeI\/JFhGwfdJkFH7FCvzBuTjbzXEW93Fr1pPp+laivivTIyrvY3jNb30OCMPHKwXfg92we26gDq9M8Sak+tQ6VruiLps11E0lq8d2s6S7MF1OACpAYdsHB5rpERUGFUKCScAY5Jya4\/wAJL4bfUpZbOW+\/tlItksGqzyNdQoSG27ZGOF+7yvB45NdnQAVXu7210+1e6vLiK3t4xl5ZXCqv1JqxXn3jq7kbxh4W05dPbUQzT3C2hcKkkqKAhYnoFyzd\/oaAOx0vWtM1qF5tMv7e7jQ4ZoZA20+hx0qGHxLodxqjaZDq9lJfKSDAsyl8jqMZ6+1cVd6pJZ3PiWe80kaX4kTRZJkkt7gyRXES5w44HzK2ByM49qNf0y1034L201rEkc9hbW93byqvzLMCjF8+p5z9aAPS6xND\/wCQx4m\/7CSf+klvWyjb0VsEZGcHtXKWWo3Vnr3iSODRb++U6gjGS3eAKD9lt+D5kinP4Y560AdbRVa3uZpoEkeyngZhkxSMhZfY7WI\/I1L5j\/8APGT81\/xoAkoqPzH\/AOeMn5r\/AI0eY\/8Azxk\/Nf8AGgCSio\/Mf\/njJ+a\/40eY\/wDzxk\/Nf8aAH1lj\/TtbJODDY8D3lYf0U\/8Aj5qxfX32KzluGhc7RwuVyx6ADnqTgfjUelwSWenxxSRs8xy8rDb8zscsevqTUS1di46RbNCsfTM3urX+ok5jVvskAx0CE7z+L5H\/AAAVLrF\/LZaZNLDC32lh5cCsR80rfKg6+pH4ZqXT7Yadp8FnFC5SFAufl+b1PXqetabIWyL9FR+Y\/wDzxk\/Nf8aPMf8A54yfmv8AjSJJKKj8x\/8AnjJ+a\/40eY\/\/ADxk\/Nf8aAJKKj8x\/wDnjJ+a\/wCNHmP\/AM8ZPzX\/ABoAkoqPzH\/54yfmv+NHmP8A88ZPzX\/GgCSio\/Mf\/njJ+a\/40eY\/\/PGT81\/xoAkoqPzH\/wCeMn5r\/jR5j\/8APGT81\/xoAkoqPzH\/AOeMn5r\/AI0eY\/8Azxk\/Nf8AGgCSio\/Mf\/njJ+a\/40eY\/wDzxk\/Nf8aAJKKj8x\/+eMn5r\/jR5j\/88ZPzX\/GgCSoz\/wAfKf7jfzFHmP8A88ZPzX\/GozI\/2hP3Mn3G7r6j3oAsUVH5j\/8APGT81\/xo8x\/+eMn5r\/jQBJRUfmP\/AM8ZPzX\/ABo8x\/8AnjJ+a\/40ASUVH5j\/APPGT81\/xo8x\/wDnjJ+a\/wCNAElFR+Y\/\/PGT81\/xo8x\/+eMn5r\/jQBJRUfmP\/wA8ZPzX\/GjzH\/54yfmv+NAElFR+Y\/8Azxk\/Nf8AGjzH\/wCeMn5r\/jQBJRUfmP8A88ZPzX\/GjzH\/AOeMn5r\/AI0ASUVH5j\/88ZPzX\/GjzH\/54yfmv+NAElFR+Y\/\/ADxk\/Nf8aPMf\/njJ+a\/40ASUVH5j\/wDPGT81\/wAaPMf\/AJ4yfmv+NAElFR+Y\/wDzxk\/Nf8aPMf8A54yfmv8AjQBJRUfmP\/zxk\/Nf8aPMf\/njJ+a\/40ASUVH5j\/8APGT81\/xo8x\/+eMn5r\/jQBJRUfmP\/AM8ZPzX\/ABo8x\/8AnjJ+a\/40ASVHN9wf76\/+hCjzH\/54yfmv+NRzSPsH7mT7y91\/vD3oAsUVH5j\/APPGT81\/xo8x\/wDnjJ+a\/wCNAElFR+Y\/\/PGT81\/xo8x\/+eMn5r\/jQA+uIie7S4itYoXn03TwHiaFlDozZ8tTuODhG9\/vL3FdVf6gmn2M93NDJ5cSFj059hz1PSsfTIp7SOO3ntGLqfOn8k7wZWHAP+70HsFrKrNwWuz0+YlHmlp0\/In0ezuYtRnudQI+1yx\/KFbIC7j7fe+7ntwMVvVRh86S5NxJbuny7UXKkgZySee\/H5Va8x\/+eL\/mv+NFL4bFS3JaKj8x\/wDnjJ+a\/wCNHmP\/AM8ZPzX\/ABrURJRUfmP\/AM8ZPzX\/ABo8x\/8AnjJ+a\/40ASUVH5j\/APPGT81\/xo8x\/wDnjJ+a\/wCNAElFR+Y\/\/PGT81\/xo8x\/+eMn5r\/jQBJRUfmP\/wA8ZPzX\/GjzH\/54yfmv+NAElFR+Y\/8Azxk\/Nf8AGjzH\/wCeMn5r\/jQBJRUfmP8A88ZPzX\/GjzH\/AOeMn5r\/AI0ASUVH5j\/88ZPzX\/GjzH\/54yfmv+NAElFR+Y\/\/ADxk\/Nf8aPMf\/njJ+a\/40ASUVH5j\/wDPGT81\/wAaPMf\/AJ4yfmv+NAElRz\/8e8v+4f5UeY\/\/ADxk\/Nf8ajnkb7PJ+5f7p7r6fWgCxRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFAHP8AiuzW8ttPKajbWN7DepLZvcqGjkmCsAhUkZypbocjqOlQaNpGuvrSax4hutPa5ht3toYNPiZUCsyszMzncT8g46D8a3b\/AE+01Szls7+2iubaQYeKVdwNcq\/gS7sj\/wASHxXrGnR4CpbyuLqGNf7qrJyPzoA6ua1tXuIryeGFprcMY5nUbogR82G7ZHWltLy2v7VLqzniuIHzsliYMrYODgjryK5aLwBDcSrL4g1nU9bOQxt7mXZbbh0PlLgduhyK62ONIYkijRUjRQqqowFA6ACgCSud8SaHe393pmq6VNDHqemyO0S3GfLlRxtdGxyMgDnnGK6KigDj7Tw3f6vqt7qniRbVGnsG06K1tGLqkTHLlmYDLH6cCs8+F\/Et9otr4Y1GfTjo8JjSW7jZ\/OuIY2BVNhGFY7VBO4\/jXoFFABWJof8AyGPE3\/YST\/0kt626xND\/AOQx4m\/7CSf+klvQBt0UUUAFFFFACUUVDdXMdpay3EzbY41LsT2Ao2BLoULg\/btahtRzDaATy\/75yEH\/AKE3thfWtWs\/SLeSK0M067bm5bzpQeoJ6L+AAX8KuyyJDE8sjBY0UszHoAKmC69y572RlT5vvEVvb4\/c2KfaZPeRtyoPwG8\/981sVk6BG5sWvZlKzX0huWDDBUHhFPuECj8DWt0q5diX2FooopCCiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKjP\/Hyn+438xUlRn\/j5T\/cb+YoAkooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACo5vuD\/fX\/0IVJUc33B\/vr\/6EKAJKKKKACiiorieK1tpZ5nCRRIXdz0UAZJoE2Y+qSpc6lFakFre0AurgDuwP7pPxYbv+AD1rUtImigzJjzXO+Qj+8f84\/CsfR7eWb99coyzTv8Aa7hW6qTxHH\/wFQPxGe9dBUP3p+S0+fUI6Ru92LRRRVjCiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAqOf\/j3l\/wBw\/wAqkqOf\/j3l\/wBw\/wAqAJKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAo6rqH9l6bNem0ursRAEw2ke+RucfKuRmuLsfif51xfJL4Y8QusNxsjEGnMzKuxTiT5uHyTx6ba9CrkdQ0DxRJql1caT4nt9PtZ3DiD+zI3IO0AktnLHjqfagBdC8Vajr3iR4E0LULHSY7Qu02oWrQuZt4AUc4I25P4V1tcR4el8T2PjObSvEOtR30MlkZ7TyrRY1fDqrEkchlyPl5BD5zxiu3oAKKKKACiiigArE0P\/kMeJv8AsJJ\/6SW9bdYmh\/8AIY8Tf9hJP\/SS3oA26KKKACiiigCne3qWaJlXklkbZFEg+Z2xnHoOAeTxWFqGo3ck8Ftf2aQWyss9wYZfN2xqeN\/yjA3YzjPCt2rT1KRbTUbO9mOLaNZI3bshbaQx9vlIz71RhvIRY3d8zLNPfsfIhBGXXG1FH1Hzf8CNYTk22rnRTikr2v8A5nRBgQGByD0NZGvFrmO20tACb2XZKPSEcyfgR8v\/AAMVoWMDW1hbwO25o41Qse5Axms3TwL7Xb7UCCUg\/wBDgPrjmQj\/AIFhf+2ddEO5itG2bQAAwOKWiikSFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFRn\/AI+U\/wBxv5ipKjP\/AB8p\/uN\/MUASUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABUc33B\/vr\/AOhCpKjm+4P99f8A0IUASUUUUAJ3rE1qZZ7iDTjzFj7TdY5\/dqeF\/wCBNge4DVsuwRSzEAAZJPasDSka\/uZL6Uc3LibBGMRLkQr\/ADf2LUNuMbrfoS1zNRNizhaOEtJ\/rpDvk+p7fgMD8KtUlLUxjyqxTd2FFFFUAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFRz\/8AHvL\/ALh\/lUlRz\/8AHvL\/ALh\/lQBJRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABXC+JkvLTVHnm+IkWiwTcw2ssEPygAA4LHLc\/zruq8wtbrwpa+NvEUHiuXTJNUe5DwTXgVkW22LsQFuEYc5HBPB5oA3vCWh2f22XxCfET+IL6WL7L9rDp5caAhtiqnC9if6ZNdjXl2l+IPBunfEK5n0jUdPsbA2Hl3myVIoJZt4KbV43MF35YcfMO5Nel21xDd20dxbypLBKoeORGyrKehBoAmooooAKKKKACsTQ\/+Qx4m\/7CSf8ApJb1t1iaH\/yGPE3\/AGEk\/wDSS3oA26KKKACiiigDL1p3NktrGxWW7cQIVOCM\/eIPqFDH8Ks22n2dmc21rDEcYJSMKSKqxn7Zr0smQ0dknlrj\/no2C35Ls\/76NalQlduRcm0lEpatetp+lzXKIHlACxJ\/fkYhUX8WIFO0yyGnabb2gcuYkAZz1du7H3JyfxqldkX\/AIgtbQHMVmv2mYdt5+WMH\/x9v+ArWxWj0ViXorC0UUUhBRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABUZ\/4+U\/3G\/mKkqM\/wDHyn+438xQBJRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFRzfcH++v\/oQqSo5vuD\/AH1\/9CFAElFFFAGJ4glEkUOnc7brcZyP4YFwZPzyqf8AA60bKExQbnUCSQ7mA7eg\/AYH4Vjaeo1TVJtQIyjkLHn\/AJ4oTt\/76fc30C10VS9ZW6L8\/wCtBR2v3FoooqhhRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFRz\/8AHvL\/ALh\/lUlRz\/8AHvL\/ALh\/lQBJRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABXm96+peI9f1QaD4b8OTQWM\/2ae81WMs00wALABVz8uQOevGK7TXE1t7FRoMtlHd+YNxvFZk2YOfu85zivOQfFGi69fTnxZ4PtJ7lle6s5JWVDJtA37WO5WKhehGcCgDo9C8N6ib1x4g8P+EltPLOw2FsS+\/IxncMYxu\/Suyhhit4UhgjWOJBtREXCqPQCuP8Na\/qU+q+Tq\/iLwtdROhWKLTpv3pkyMdW5GN36V2tABXD\/EbUvssGiWP228tFvL9RK9kziYxKrFlXZ82SdoruK5\/xFqd\/pNzp13DpbX1iHdbwwxmSeEFflZFHUZzu74oA5zSNU0nS9O12\/sNU1u5vdPs3llstXml3JhSwOx8dcYyKpXtleaD4HtvGC6pfzaxEkN3c+ZcuYp1crvjMedoXDHGAMYFW5tNn8aeIdR1CG0ubKwbRpdNSW7haJp3kOchTzsX1PfpVW8utR13wTbeD\/wCxtRg1WRIbS6eS3ZYYVQrvk8z7rKQpxgnrQB6arBlDA5B5BrlLLX9G0rXvEkGo6vYWczagjiO4uUjYr9ltxnDHpwfyrrAoVQAMAcAVi6H\/AMhjxN\/2Ek\/9JLegDSgv7O7gSe2u4JoXGUkjlVlYexFS+fD\/AM9Y\/wDvoVJRQBF58P8Az1j\/AO+hUNzf21rbSzySrsjUscNk8VarL1M\/abuz08AlZH86XnoiYP6tsH0zUydkVBXeoukBbfTo\/OkjFxJmWb5h99uT+WcfhV1rmBFLNNGABySw4qasTXZo7qCLSopQZL6XyW2tyIx80nTp8oK\/VhVQj0D4pXDQJEktJdRldRLfyefgnBVMARjB6fIFz7k1r+fD\/wA9U\/76FPUBVAAwBwAKdTbuxN6kfnw\/89Y\/++hR58P\/AD1j\/wC+hUlFIRH58P8Az1j\/AO+hR58P\/PWP\/voVJRQBH58P\/PWP\/voUefD\/AM9Y\/wDvoVJRQBH58P8Az1j\/AO+hR58P\/PWP\/voVJRQBH58P\/PWP\/voUefD\/AM9Y\/wDvoVJRQBH58P8Az1j\/AO+hR58P\/PWP\/voVJRQBH58P\/PWP\/voUefD\/AM9Y\/wDvoVJRQBH58P8Az1j\/AO+hR58P\/PWP\/voVJRQBH58P\/PWP\/voUefD\/AM9Y\/wDvoVJRQBH58P8Az1j\/AO+hR58P\/PWP\/voVJRQBH58P\/PWP\/voVGZ4ftCfvU+438Q9RVioz\/wAfKf7jfzFAB58P\/PWP\/voUefD\/AM9Y\/wDvoVJRQBH58P8Az1j\/AO+hR58P\/PWP\/voVJRQBH58P\/PWP\/voUefD\/AM9Y\/wDvoVJRQBH58P8Az1j\/AO+hR58P\/PWP\/voVJRQBH58P\/PWP\/voUefD\/AM9Y\/wDvoVJRQBH58P8Az1j\/AO+hR58P\/PWP\/voVJRQBH58P\/PWP\/voUefD\/AM9Y\/wDvoVJRQBH58P8Az1j\/AO+hR58P\/PWP\/voVJRQBH58P\/PWP\/voUefD\/AM9Y\/wDvoVJRQBH58P8Az1j\/AO+hR58P\/PWP\/voVJRQBH58P\/PWP\/voUefD\/AM9Y\/wDvoVJRQBH58P8Az1j\/AO+hR58P\/PWP\/voVJRQBH58P\/PWP\/voUefD\/AM9Y\/wDvoVJRQBH58P8Az1j\/AO+hR58P\/PWP\/voVJRQBH58P\/PWP\/voVHNPDsH71Pvr\/ABD+8KsVHN9wf76\/+hCgBPPh\/wCesf8A30KyPEF8iWAtYp9k123lB0YZRcZd\/wAFDY98DvW0a523P9sa7LcA5t4CYY\/TCt8x\/wCBOuPpF70m+VXRMlf3e5p6asFraKoMcZIHyBh8gxgL+AAFXPPh\/wCesf8A30Kko70oqysUxnnw\/wDPWP8A76FHnw\/89Y\/++hUlFUBH58P\/AD1j\/wC+hR58P\/PWP\/voVJRQBH58P\/PWP\/voUefD\/wA9Y\/8AvoVJRQBH58P\/AD1j\/wC+hR58P\/PWP\/voVJRQBH58P\/PWP\/voUefD\/wA9Y\/8AvoVJRQBH58P\/AD1j\/wC+hR58P\/PWP\/voVJRQBH58P\/PWP\/voUefD\/wA9Y\/8AvoVJRQBH58P\/AD1j\/wC+hR58P\/PWP\/voVJRQBH58P\/PWP\/voUefD\/wA9Y\/8AvoVJRQBH58P\/AD1j\/wC+hR58P\/PWP\/voVJRQBH58P\/PWP\/voUefD\/wA9Y\/8AvoVJRQBH58P\/AD1j\/wC+hR58P\/PWP\/voVJRQBH58P\/PWP\/voVHPPEYJP3qfdP8Q9KsVHP\/x7y\/7h\/lQBJRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABXmFpD4Vk8a+JH8WQ6WmqfaVECXwRUa22LsZd3DMcNuPJ7cV6fXDeJtH8S6tqTiPRvCl\/Yxn\/R21JJGlUYGe2Bz6UAcudK8DXPizW\/LttMbw7HpiveXUQXy4LgvhRG4+6Sp6L1IHevQPA899ceB9Fl1LJuntULMSSWGPlJJ7lcE+5rn\/AAw2p\/21deFte0Pw\/Z2S232qK3s4SUm\/eD5wD8uARznDZ2nGOa9BoAKKKKACiiigArE0P\/kMeJv+wkn\/AKSW9bdYmh\/8hjxN\/wBhJP8A0kt6ANuiiigBD0rL0zN1eXl+w+V38iH\/AHEJBP4sW\/DFTarcyW2nSPDjzmIjiz\/fY7R+pz+FT2dqllZw20edkSBFJ6kAd6h6yt2LWkb9zN1ZDe6pY6Y7MtrMks0wVipk2bAEyOx35Prtx0JrJ\/si2uL3Ur2wihs\/scf2e1kgUJiRfmdjjqM7Vwf7reta3iRFTTHvUdo7q1Be2dME+YflC4PB3Ehce\/Y81UtND1L+zU0++vYBbHJnFvGwecsSzZYn5QxJzgd+ordOyGnZG1pty15plpdOmxpoUkK\/3SQDirdMRVRQqgBQMADoKfWZmFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABUZ\/4+U\/3G\/mKkqM\/8fKf7jfzFAElFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAVHN\/qx\/vL\/6EKkqOb\/Vj\/eX\/wBCFAGfrt69jpjmFgtzKfKhJ6Bj\/EfZRlj7Kado1ktjp0cYVgdo4Y5IAGAD74HPvms1idY8TMvBtbH5B7uQCx\/9BX\/vsV0dS9ZW7fmSu4tFFFUUFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFRz\/APHvL\/uH+VSVHP8A8e8v+4f5UASUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAVyerr43Gozvpl54fi08keULtJfMHHOSDjrmusrzi08P6T4y8V+IW8SRG5u7G7FvbWTzsqwW+xSrhVI+\/knJ9KANvQNA1oa82v+INVt7q6+zG2gt7OMrBChYMSMnLE7V\/+vxjrK4LQNNsPDXxCn0TQ5nSxm043NzY72dYJQ6qrgsTgsrNx\/sg+ld7QAUUUUAFFFFABWJof\/IY8Tf8AYST\/ANJLetusTQ\/+Qx4m\/wCwkn\/pJb0AbVLRTHdY0Z3YKqjJJOABQBmTD7Zr8EXPl2aec3oXbKr+QDn8RWrWZoil7Nr1wRJeOZzkYIU\/dH4KFH4Vp9KmG1+5c97djG1DF\/rljp+3McH+mT8cccRg\/VssP+udbNY2gn7XHcasQP8ATpN8f\/XJRtT8wC3\/AAM1s1cuxL7C0UUUhBRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAVGf8Aj5T\/AHG\/mKkqM\/8AHyn+438xQBJRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAMkdY0Z2ICqMk+grjrnUL270k6q91c28UhX7PDFsRQS3yAkhmZuRnHyjnrjNdfKgljdGAKspBBHWuIjuL6a8t7O6tJJoNNYD\/R1\/1rHgZBxjaPl9DuJz6Q5taImSua3hQm0WbTp\/muVAnZ24Z9xOSw9d2T\/wIV0tYulG6udSu725iWIFVijjBDFQMkgkd8nnBI6DtW1RB3Q0rC0UUVYwooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAqOf\/AI95f9w\/yqSo5\/8Aj3l\/3D\/KgCSiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACvKdUm8I6t4w1Ww8bf2dDeWsoFlKkrRFoCoIDyKw+b\/ZbGM8ZBr1avN719R8R+INUGg+G\/Dc0FjP8AZp7zVYyzTTBQWACrn5cgc9eMUAbvhBPBlj51h4VuNPaRh5sq29x5sjAcZYkliBu\/DPvXV1yHhbSddsNTkl1PSfDNpCYSqyaVE6yltw4OR93AP4gV19ABRRUU80VtBJPPIscUal3djgKo5JNAEtFcp4X8aQ+IrHV9Qe3a0s7CdlV5M5aMIG3sO3BziqUfjbU47C01y+0SO30C6dFWYXO6aJHICSOm3G05HRiRmgDuKxND\/wCQx4m\/7CSf+klvW3WJof8AyGPE3\/YST\/0kt6ANqsvWW82KHT1yWvJPLbHaMcv\/AOOgj6sK1Ky7X\/S9aubnkx26\/Z4+P4urkf8Ajo\/4CaiXbuXDfm7GmOAB6VleIZJDpws4JDHcX0gto3A5XIJZvwQMfwrWNY8JF\/4lnlBDRafH5K\/9dXAZ\/wAl2f8AfRrSJK3uasUaQxJFGoVEUKqjoAKkoopCCiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKjP8Ax8p\/uN\/MVJUZ\/wCPlP8Acb+YoAkooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKAKOqXw07T5bjbucYWNM43uThR+JIrO0zSoo9MV7hFlmmcO7suCxZsk\/jyce9R3n\/ABOPEEdoObaz5k9DIV5\/JTj\/ALaf7Nb033B\/vr\/6EKiUVJ2ey\/MUW90ORFRAiKFUcAAYAp9FFUlYYUUUUwCiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACo5\/8Aj3l\/3D\/KpKjn\/wCPeX\/cP8qAJKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAxPEt1qNtYQDS73TLS7luFiRtRJCSZB+RcHJYnGB9a5O18PfEOx1a7v7bVNAj+2MHntxHL5TSBQu\/HUNgDOCM45q38UZZoNN0Ca3tzcTx67atHAHC+Yw3YXJ4GTxmpf+Er8Z\/wDRPZv\/AAaQ\/wCFAGroMfi5L6Q6\/caPJaeWdgso5FffkYzu4xjd+ldFXN6Drev6lfSQ6r4Wk0qARFlna9Sbc2QNuFGRwSc+1dJQAVxHjee9lv8AT9OOj6lfaO37+8+wxBzKVPyxHLLhc8t64A7mu3ooA8h0nUJ9S0r4h2kWlahA9x9qlUzRKoQ+SF8psMcP7enetnxNPDcfBBPJYP8AaLG2jhVerMxQBR75\/lXb2OlWWmzXktpD5b3kxnnO9jvfAGeTxwB0rKtvA3h6z1GO8hsSGik82KIzOYYpDn5ljLbFPPYcdqAOhQEIoY7mA5PrXKWWiWuo6\/4kmnlv0YagigW+oTwLj7Lb\/wAMbqM89etdbWJof\/IY8Tf9hJP\/AEkt6ALM6RaVpLtGZ3EKfIJLh3Zj2BZiSSTxyadplh9k0+GKRmMuN0rK5AZzyx69ySaivv8AS9Us7IfcQ\/aZh7KfkH\/fXP8AwCtTNQtZNlvSNu5TvpYbGxnupTJshQuQHOTgdBz1qvotg9rpUK3G4XUmZpwJG4kc7mHXoCcfhUerZvNQsNM\/gkf7TNx\/BGQQPxcp9QGrYrTZC2QzyV9ZP+\/jf40eSvrJ\/wB\/G\/xqSikSR+SvrJ\/38b\/GjyV9ZP8Av43+NSUUAR+SvrJ\/38b\/ABo8lfWT\/v43+NSUUAR+SvrJ\/wB\/G\/xo8lfWT\/v43+NSUUAR+SvrJ\/38b\/GjyV9ZP+\/jf41JRQBH5K+sn\/fxv8aPJX1k\/wC\/jf41JRQBH5K+sn\/fxv8AGjyV9ZP+\/jf41JRQBH5K+sn\/AH8b\/GjyV9ZP+\/jf41JRQBH5K+sn\/fxv8aPJX1k\/7+N\/jUlFAEfkr6yf9\/G\/xo8lfWT\/AL+N\/jUlFAEfkr6yf9\/G\/wAaPJX1k\/7+N\/jUlFAEfkr6yf8Afxv8aPJX1k\/7+N\/jUlFAEfkr6yf9\/G\/xqMwp9oTmT7jf8tG9R71YqM\/8fKf7jfzFAB5K+sn\/AH8b\/GjyV9ZP+\/jf41JRQBH5K+sn\/fxv8aPJX1k\/7+N\/jUlFAEfkr6yf9\/G\/xo8lfWT\/AL+N\/jUlFAEfkr6yf9\/G\/wAaPJX1k\/7+N\/jUlFAEfkr6yf8Afxv8aPJX1k\/7+N\/jUlFAEfkr6yf9\/G\/xo8lfWT\/v43+NSUUAR+SvrJ\/38b\/GjyV9ZP8Av43+NSUUAR+SvrJ\/38b\/ABo8lfWT\/v43+NSUUAR+SvrJ\/wB\/G\/xo8lfWT\/v43+NSUUAR+SvrJ\/38b\/GjyV9ZP+\/jf41JRQBH5K+sn\/fxv8aPJX1k\/wC\/jf41JRQBH5K+sn\/fxv8AGjyV9ZP+\/jf41JRQBH5K+sn\/AH8b\/GjyV9ZP+\/jf41JRQBF5K+sn\/fxv8ap6pOmn6fJcBZJJBhY4xIcu7HCr17kitCudvT\/a2uLbKcwWhwcd5WXn\/vlD+cg9KL2V30JfZdSxoGmi2sA7u0kkhLNIGI3knJb8SS30I9K05oU8scyfeX\/lo3qPepUQIAo4AGABTZvuD\/fX\/wBCFTFWWpXoL5K+sn\/fxv8AGjyV9ZP+\/jf41JRVAR+SvrJ\/38b\/ABo8lfWT\/v43+NSUUAR+SvrJ\/wB\/G\/xo8lfWT\/v43+NSUUAR+SvrJ\/38b\/GjyV9ZP+\/jf41JRQBH5K+sn\/fxv8aPJX1k\/wC\/jf41JRQBH5K+sn\/fxv8AGjyV9ZP+\/jf41JRQBH5K+sn\/AH8b\/GjyV9ZP+\/jf41JRQBH5K+sn\/fxv8aPJX1k\/7+N\/jUlFAEfkr6yf9\/G\/xo8lfWT\/AL+N\/jUlFAEfkr6yf9\/G\/wAaPJX1k\/7+N\/jUlFAEfkr6yf8Afxv8aPJX1k\/7+N\/jUlFAEfkr6yf9\/G\/xo8lfWT\/v43+NSUUAR+SvrJ\/38b\/GjyV9ZP8Av43+NSUUAR+SvrJ\/38b\/ABo8lfWT\/v43+NSUUAR+SvrJ\/wB\/G\/xqOeFRbycv90\/xt6fWrFRz\/wDHvL\/uH+VAElFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAZWs6Fa66LAXUkyfYryO9j8ogZdM4DZB+Xn\/AOvWrRRQAUUUUAFFFFABRRRQAlYmiHGr+Jv+wkn\/AKSW9bdcokzpceJoomxNcapHBHjqC1rbgkfQZb\/gNJuyuOKu7GvpK+e91qJ63Mm1P+uaZC\/geW\/4FWp3qK3hjtreOCJQscahVUdgBgVS1y6kttMkW3JF1ORBBjs7HAP0H3j7A0QjokNvmloQaOBd3l\/qhJIlk8iHnjy4yVyPq2857gitmobS1jsrOG1hGIoUWNAfQDAqaqbE3qLRRRSEFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAVGf8Aj5T\/AHG\/mKkqM\/8AHyn+438xQBJRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQBR1S9\/s\/T5bgJ5jjCxxg4LuThV\/EkVU0GyNtbFnbzJCTmTGN7E5dvxbP4Yqteu+pa0IYj+6sztBHedl5P\/AEJP1celb0aLFGsaDCqAAPQVMtZKPYmOt5ElRzfcH++v\/oQqSo5vuD\/AH1\/9CFUUSUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABUc\/\/HvL\/uH+VSVHP\/x7y\/7h\/lQBJRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQBn3t7LFPFaWyK11MCV3n5UUYyxx9QMdyfqRylj9vs\/E+t3l0izW1rcq8jou3cWt4gXUZP3VUDHuTznFdNqAktb+DUEieWJUMUyoNzAEghgO+COR7+2KwLO6F7P4ks7OOWRru+Co3lkJGptoAzEkYGDuO3qfSsJ6vc6KeiukdkrBlBByDWO\/+neJo0wDBp8e8\/8AXZxgfkm7\/v4K02eKztGd2CRQpkseiqBVDQIZV077VcKVub1zcyqRypYDC\/8AAVCr\/wABrojormO12a1FFFIkKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKjP\/AB8p\/uN\/MVJUZ\/4+U\/3G\/mKAJKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooASqWqXo0\/T5LgJ5kgwscYPLuThV\/EkCrtc7eTPf61sjIMdiQi9w1y4\/9kQlv+Bj0ouknJ9CZX2W7LWiWX2eH538x13AyY\/1khOZH\/Fv0ArYqOGJYIUiT7qgAVJURTtd7srRaIWo5vuD\/AH1\/9CFSVHN9wf76\/wDoQqwJKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAqOf8A495f9w\/yqSo5\/wDj3l\/3D\/KgCSiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKzdY1yw0K1SfUJigkcRxoiM7yOeiqqgljWlVG80+C4mgvDbRS3loHa1aQ42My4PPbPSgCDSPEGn64s4s5JBLbsFngniaKWInpuVgCM1QtvHOg3d7FbRXMu2eUwwXDW8iwTSD+FZCNrH8ea428udQtNQ8ZNrCR2+tT6I8tsLRi0BgRWGQxwxcM3OQPatHxPBDb\/A9BEAggsLaSFl\/hZShDD3zQB6LWHoY\/wCJx4m\/7CSf+klvW0hJRSw2sRyPSuUsodZk17xIdOv7C3h\/tBMrcWTzMW+y2\/O4Spx04x+NAGnrwF59k0kDIvZP3o9IU+Z8+x+VP+B1s1z+jLd3eoXl\/cTQSPF\/oaPHCVVthJkYKWJGXO3GT\/qxW5ib++n\/AHwf8ab7DfYloqPE3\/PRP++D\/jRib\/non\/fB\/wAaQiSio8Tf89E\/74P+NGJv+eif98H\/ABoAkoqPE3\/PRP8Avg\/40Ym\/56J\/3wf8aAJKKjxN\/wA9E\/74P+NGJv8Anon\/AHwf8aAJKKjxN\/z0T\/vg\/wCNGJv+eif98H\/GgCSio8Tf89E\/74P+NGJv+eif98H\/ABoAkoqPE3\/PRP8Avg\/40Ym\/56J\/3wf8aAJKKjxN\/wA9E\/74P+NGJv8Anon\/AHwf8aAJKKjxN\/z0T\/vg\/wCNGJv+eif98H\/GgCSio8Tf89E\/74P+NGJv+eif98H\/ABoAkoqPE3\/PRP8Avg\/40Ym\/56J\/3wf8aAJKKjxN\/wA9E\/74P+NGJv8Anon\/AHwf8aAJKKjxN\/z0T\/vg\/wCNGJv+eif98H\/GgCSoz\/x8p\/uN\/MUYm\/56J\/3wf8ajIm+0J86fcb+A+o96ALFFR4m\/56J\/3wf8aMTf89E\/74P+NAElFR4m\/wCeif8AfB\/xoxN\/z0T\/AL4P+NAElFR4m\/56J\/3wf8aMTf8APRP++D\/jQBJRUeJv+eif98H\/ABoxN\/z0T\/vg\/wCNAElFR4m\/56J\/3wf8aMTf89E\/74P+NAElFR4m\/wCeif8AfB\/xoxN\/z0T\/AL4P+NAElFR4m\/56J\/3wf8aMTf8APRP++D\/jQBJRUeJv+eif98H\/ABoxN\/z0T\/vg\/wCNAElFR4m\/56J\/3wf8aMTf89E\/74P+NAElFR4m\/wCeif8AfB\/xoxN\/z0T\/AL4P+NAElFR4m\/56J\/3wf8aMTf8APRP++D\/jQBJRUeJv+eif98H\/ABoxN\/z0j\/74P+NAFXVb8abp0tzs8yQYWKMHBkdjhVH1JAqno1ibaNUd\/MaLd5kmP9ZMx3SN+ZwPTpVW7llvNZYqyNFpxAUbeHuXGAOv8KnP\/Awe1bVvBJBAkSyKQo5JQ5J7nrUz1ah83+hMNby+SLVFR4m\/56J\/3wf8aMTf89E\/74P+NUUSVHN9wf76\/wDoQoxN\/wA9E\/74P+NRzCbYPnj++v8AAf7w96ALFFR4m\/56J\/3wf8aMTf8APRP++D\/jQBJRUeJv+eif98H\/ABoxN\/z0T\/vg\/wCNAElFR4m\/56J\/3wf8aMTf89E\/74P+NAElFR4m\/wCeif8AfB\/xoxN\/z0T\/AL4P+NAElFR4m\/56J\/3wf8aMTf8APRP++D\/jQBJRUeJv+eif98H\/ABoxN\/z0T\/vg\/wCNAElFR4m\/56J\/3wf8aMTf89E\/74P+NAElFR4m\/wCeif8AfB\/xoxN\/z0T\/AL4P+NAElFR4m\/56J\/3wf8aMTf8APRP++D\/jQBJRUeJv+eif98H\/ABoxN\/z0T\/vg\/wCNAElFR4m\/56J\/3wf8aMTf89E\/74P+NAElFR4m\/wCeif8AfB\/xoxN\/z0T\/AL4P+NAElFR4m\/56J\/3wf8aMTf8APRP++D\/jQBJRUeJv+eif98H\/ABoxN\/z0T\/vg\/wCNAElRz\/8AHvL\/ALh\/lRib\/non\/fB\/xqOcS\/Z5Mun3T\/AfT60AWKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigArB8QeH5dWnsL6yvjY6lYOzQTeX5iEOMOrLkZBHvxW9RQBzGneFZTqN3qmvXy6le3Nt9jwkPlRRQk5ZFXJPJ6kms+PwPqD2dro15r32nQLV0ZLb7MBLIiHKRvJu5UEL0UZxXb0UAJXLRXjWD+LJ4wDMNRjSJT\/FI1rbKg\/FiB+Nbeo6lHp0UZaOWaaVtkUMK5eRuuB2HHUnAHrXFWOp+b4z1K3vraW3hi1CO5lyQypI1tCkQYj\/dY56Z21UYt6lRi3qdxplkNP023tN5cxIFZz1du7H3Jyfxq5SUtTuS9QooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKjP\/Hyn+438xUlRn\/j5T\/cb+YoAkooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKAI5ZUgheWRgqIpZmPYDrXOXuo6wNLfU43t7ZGA+zWrwmSSRmwEDHcACxI4HT1rd1C3S7065t5G2pLEyM3oCCK5G31mS+EU12rRGwz5KiF3S4lClWYHHY5A56k+gp86gr2uZTTk7Xtc1fD8UiZtLny\/tNmS8xSTcJZJCSZPUfxcHpyOwroqwfD0czT3d5cRtDJKkaCJ1w20bm3sOxLO3HbArerOL5ve7miVkl2FoooqxhUc33B\/vr\/AOhCpKjm+4P99f8A0IUASUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABUc\/\/AB7y\/wC4f5VJUc\/\/AB7y\/wC4f5UASUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAYl\/Ilr4isLm4IWBoZIFkY8LIzIQPbIUj8AO9ZEs0Ua+Lo2CvPcXot4Y+A0rNaQbVH4sfoOe1dbLDFcQtDNGkkbgqyOuQw9CDXOeG9KsbTW\/EUkFnBG8d8kSMqAFU+zQNtB7DLE496pNFJo6OFGjgjR3LuqgMx\/iPrUtFFSSFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABUZ\/4+U\/3G\/mKkqM\/8fKf7jfzFAElFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABSUtVb68i0+zluZtxVB91RksegAHck4AHqaFqJuyuylq1zJLLHpdqy\/aLgEyMefKi\/iYj3+6Pc+xqGzlgRYY2Tybe2XbFGqHaxHG4ccj0\/P0pLG0nkMwuT\/pdzh7tgeI1\/hiU+wP6k\/xVtqoVQqgAAYAHasaic5e49F+L\/wCAVDRXlu\/wK1srPPLcsrIHAVVbg7Rnkjt1P6Vbooq4R5VYG7i0UUVYgqOb7g\/31\/8AQhUlRzfcH++v\/oQoAkooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACo5\/+PeX\/cP8qkqOf\/j3l\/3D\/KgCSiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKazKilmICgZJPQUAOopkciSoHjdXU9GU5FIJo2lMQkQyAZKBuR+FAElYmh\/8hjxN\/2Ek\/8ASS3rbrE0P\/kMeJv+wkn\/AKSW9AG3RRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFRn\/j5T\/cb+YqSoz\/x8p\/uN\/MUASUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAJXOyzNquqLJGFeC1kKWysvEkw4aQ\/7KZIHvn2qxrt8yhNOt5ClxcA73Q\/NFH0LD\/aJ+Vfc57GrmnWSWcCqsYQ7QqoOiKOij6fqeaUm0uVbv8AAm3M9dkWIIVgiCKSe7MerHuTU9JS0JJKyKbuFFFFMAooooAKjm+4P99f\/QhUlRzfcH++v\/oQoAkooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACo5\/wDj3l\/3D\/KpKjn\/AOPeX\/cP8qAJKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigArD1\/w5F4huNPW8mLafbSNJNZFcpcnGE3c9FPOOQa3K43x4niS7gtLDQ7Cee0mY\/b5Le5jhk8sfwKzkY3c5I9PegDBjuLfR9b8V3\/hiFIdKsdJfzxAMQNfLll2r0yF4bHtUep+H7HQvhpbeIrOFU1y0ihvTf8A\/LWWRipfe3VlbcwweOa6LSbW8vdFuvDlx4UfQ9MezkhSQXcUo+YbcbVOc\/MTk+nvWVPp3ibV\/C9r4PvNH8hFEUF3qX2hDE0KEZZFB3bmCjggYzQB6KjB0Vl6MMiuUstRurPXvEkcGi398p1BGMlu8AUH7Lb8HzJFOfwxz1rre1Ymh\/8AIY8Tf9hJP\/SS3oA0re5mmgSR7KeBmGTFIyFl9jtYj8jUvmP\/AM8ZPzX\/ABqSigCPzH\/54yfmv+NHmP8A88ZPzX\/GpKKAI\/Mf\/njJ+a\/40eY\/\/PGT81\/xqSigCPzH\/wCeMn5r\/jR5j\/8APGT81\/xqSigCPzH\/AOeMn5r\/AI0eY\/8Azxk\/Nf8AGpKKAI\/Mf\/njJ+a\/40eY\/wDzxk\/Nf8akooAj8x\/+eMn5r\/jR5j\/88ZPzX\/GpKKAI\/Mf\/AJ4yfmv+NHmP\/wA8ZPzX\/GpKKAI\/Mf8A54yfmv8AjR5j\/wDPGT81\/wAakooAj8x\/+eMn5r\/jR5j\/APPGT81\/xqSigCPzH\/54yfmv+NHmP\/zxk\/Nf8akooAj8x\/8AnjJ+a\/40eY\/\/ADxk\/Nf8akooAj8x\/wDnjJ+a\/wCNHmP\/AM8ZPzX\/ABqSigCPzH\/54yfmv+NHmP8A88ZPzX\/GpKKAI\/Mf\/njJ+a\/41GZH+0J+5k+43dfUe9WKjP8Ax8p\/uN\/MUAHmP\/zxk\/Nf8aPMf\/njJ+a\/41JRQBH5j\/8APGT81\/xo8x\/+eMn5r\/jUlFAEfmP\/AM8ZPzX\/ABo8x\/8AnjJ+a\/41JRQBH5j\/APPGT81\/xo8x\/wDnjJ+a\/wCNSUUAR+Y\/\/PGT81\/xo8x\/+eMn5r\/jUlFAEfmP\/wA8ZPzX\/GjzH\/54yfmv+NSUUAR+Y\/8Azxk\/Nf8AGjzH\/wCeMn5r\/jUlFAEfmP8A88ZPzX\/GjzH\/AOeMn5r\/AI1JRQBH5j\/88ZPzX\/GjzH\/54yfmv+NSUUAReY3\/ADwf81\/xqtfaglhZS3U0UgSMZwNpLHoABnkk4A+tXa5p0PiLWxkZ0uwfjI\/103IP4L0+pPsaTdhO+yJdFtrh2k1K9idridt+AVIUfwgeyg4H1J\/ipx8a+GlJVtd0sEcEG+h4\/wDH63cYwB0rF8G\/8iP4f\/7Btv8A+ilpRXfce2gz\/hNvDP8A0HtK\/wDA+H\/4uj\/hNvDP\/Qe0r\/wPh\/8Ai636KoDA\/wCE28M\/9B7Sv\/A+H\/4uj\/hNvDP\/AEHtK\/8AA+H\/AOLreqmNV083ZtBfWxuc48rzV359MZzQkxNpbmb\/AMJt4Z\/6D2lf+B8P\/wAXR\/wm3hn\/AKD2lf8AgfD\/APF1v0UDMD\/hNvDP\/Qe0r\/wPh\/8Ai6jl8a+Gigxruln5l\/5f4fUf7ddHRQBgf8Jt4Z\/6D2lf+B8P\/wAXR\/wm3hn\/AKD2lf8AgfD\/APF1v0UAYH\/CbeGf+g9pX\/gfD\/8AF0f8Jt4Z\/wCg9pX\/AIHw\/wDxdb9FAGB\/wm3hn\/oPaV\/4Hw\/\/ABdH\/CbeGf8AoPaV\/wCB8P8A8XW\/RQBgf8Jt4Z\/6D2lf+B8P\/wAXR\/wm3hn\/AKD2lf8AgfD\/APF1v0UAYH\/CbeGf+g9pX\/gfD\/8AF0f8Jt4Z\/wCg9pX\/AIHw\/wDxdb9FAGB\/wm3hn\/oPaV\/4Hw\/\/ABdH\/CbeGf8AoPaV\/wCB8P8A8XW\/RQBgf8Jt4Z\/6D2lf+B8P\/wAXR\/wm3hn\/AKD2lf8AgfD\/APF1v0UAYH\/CbeGf+g9pX\/gfD\/8AF0f8Jt4Z\/wCg9pX\/AIHw\/wDxdb9FAGB\/wm3hn\/oPaV\/4Hw\/\/ABdH\/CbeGf8AoPaV\/wCB8P8A8XW\/RQBgf8Jt4Z\/6D2lf+B8P\/wAXR\/wm3hn\/AKD2lf8AgfD\/APF1v0UAYH\/CbeGf+g9pX\/gfD\/8AF0f8Jt4Z\/wCg9pX\/AIHw\/wDxdb9FAGB\/wm3hn\/oPaV\/4Hw\/\/ABdH\/CbeGf8AoPaV\/wCB8P8A8XW\/RQBgf8Jt4Z\/6D2lf+B8P\/wAXR\/wm3hn\/AKD2lf8AgfD\/APF1v0UAYH\/CbeGf+g9pX\/gfD\/8AF0f8Jt4Z\/wCg9pX\/AIHw\/wDxdb9FAGB\/wm3hn\/oPaV\/4Hw\/\/ABdXbfVbTVLCWfT54rqABl82CVJFzjplWNaVYGlf8hDxV\/1\/r\/6R29AG\/RRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAViaH\/yGPE3\/YST\/wBJLetusTQ\/+Qx4m\/7CSf8ApJb0AbdFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAVGf+PlP9xv5ipKjP\/Hyn+438xQBJRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAJRUc0qwQvK\/CopY\/QVzVzfa+bNtR8y0s7fGYrd0Z5H\/ug+hPoM1Ep2A2dSkllAsbZ9s8w+Z1ODEndvr2Hv9Ks2trDZ20dvAgSNBtUCqOizyXEU73MTRXgkKzI2Plx0AxnjH8zWrU0\/e95jemgtc94VkeLwBoUiQyTMumW5EcZUM37peBuIH5muhrD8G\/8iN4f\/wCwbb\/+ilrURL\/a1\/8A9C7qP\/fy3\/8AjtH9rX\/\/AELuo\/8Afy3\/APjta9FPmXYjlfc5bWtX1E6f5S6PqNt50scRlEkG4BnAIGJDgkHAPYnNSq22w+xDwnd\/ZcEeXm2IP\/kWtu7tYr21ktp13xSLtYZxWcdM1TyxAutv5PQsYFM2P9\/OM++2tFJWtt95lKEr31f3FHRtX1D7C0TaLqM4hlkiWTzYCSFYgAkyckYwT7dTWj\/a1\/8A9C7qX\/fy3\/8AjtXrK0hsLSK1gXbFGuACcn6k9z71Y\/lUykm7pFxhJJJsyf7Wv\/8AoXdR\/wC\/lv8A\/Hatxn+0bGRLuylgSVWR4ZipJU8H7jEc\/WrlFS35FpNdTyO\/+HnhOL4maNpqaLCtnPp9zLJEHfDMrIFPXtk10+p+IbH4dvZ2M+ni28Om3cW9xCWcpMuW8or\/ALQztOeTxWneaFdXHjzS9dR4Ra2llPbupJ3lnZCCBjGPlPeqXiXwjP4v1ZINVnVNBt4WMUEEjCWS5YECRuMAICdoBPzHJ6YpFG7od1qF9o9td6nZLZXcy72tlfcYgeik4HzYxn3rz7w14L8PeJNY8XXmr6Yl1cx69NEkjOwKqEjIAwR3Y\/nXe+HoNXttEgttcnguL6EFGuICcTKPuuQQMMRjI5Gc81ylnovjjQtU1xtJTw7NZ6jqMl6hu55lkXcqrghUx\/D60AUdN8QP4Ik8V6TI1zqNrpU1qdNikl3St9pGEg3t2DjgnOFNdBbeI9csNe07TPEenWMK6nvW1nsbhnVJFXd5bhlXqoOGHp0rOb4fXd1oGrG91KJ\/EOpXUV694kREUckRUxIqk52LtxzzyauxaN4i1rxDpOoeIE060t9KZ5Y4bKZ5TPMyFNxLKu1QGJA5PrQBl2\/j7X38MDxXPoVpHokbHz0W5Y3Hlh9rSqu3bgcnaTkgdq1JvEuv33iXVdI0LS7CVNPSGQ3V1csqSeYm4KoVTz156D8a5Pwxo\/ibxF8NLbQ9+nR6PeGRZLzzH+0JD5rbkEe3aW6gNu6ds1rWv\/CQW3xC8VjQYNNmiCWcZivJXi8tvJ+VgVVsj1XjoMGgBq+IIfEHizwLqoja1DR6ms8LtkwyIqq6k+zKefxpg+Kc39mf8JB9n0oaHv8A9R9v\/wBO8ndt83y8Y\/2tmc471o6b4BuNPuvDTPdQ3CWC3zX7nKtNJc4JKADpu3dxxisuz8A61p+mRaDbWPhs28TBY9ZlgD3Ih3ZwYmjKtJjjJbHegDoR4m1q+8Z3mjaXpdpJZ2DW7XV3NcMuY5F3fIoU5bqfTj3pPiHe3kGl6XptjcyWsmr6nDp73MRw8Ub7mcqezbUIH1rS0nRJ9P8AFGvak7Q\/Z7\/7P5KITuXy49p3DGB7YpfFnh3\/AISTRhbR3Btby3nS6s7kLu8mZDlWx3HUEehNAGPcfDLRIEgn0BDo2qW8iyR38O53bDAssgLfvAwyDuPep7jxHrmoa7qWneHNOspk0wqlzPfXDRh5WUN5aBVbopGWPc9Kp3umeOfEFsml6nNpOmWTOpurvTp5WmlUEErGCq+XnGCSTj3qw+i+ItD1\/Vr7w+mm3dtqrpNJDezPEYJggQsCqtvUhQSODmgCp\/wsC8v4PD66Po6SXmrPcwtDdT7Bayw\/fDsFOQCG6DnA9arr418VyWeteXoOmG50JmF8WvXEcwCb8Q\/JnJX+9jHHXnF3R\/BF3pN14al+1wznT3vJr2UgqZZZ+SUXGMbie44Aq1B4Xvo08ZK0tv8A8TuRmtsM3yAwLH8\/HHI7Z4oAhm8Xatfa1p+maFpltIb3S01Lz7uZkSFWbGGCglj04HqfStfwxr8mu2t4l3araahYXT2l3Asm9Q6gEMrYGVZWUjgdfauOS117SfHelWulpYXNzaeGY4Z4biVo0k2ybcq4UkEEd15GeldFouh63ounXtysmnz6zqV+Lu73F1hRTtUohA3Hai8Ejk9aANnxL\/yKur\/9eU3\/AKAa8\/8AhhqV5ommaR4c1aZpIdQ06O+0m5f+MNGrSwH\/AGkLZH+yfavRtXtJNQ0W+s4iokuLeSJS+cAspAz7c1zdx4LkvPh3pOhSXCQ6rplrb\/ZryIkiG5iQKHXIyVyCOnKk0AeZx\/2Y\/gz4YLrVrNd6cftPnQxxPKz\/ALttvyp8x+bB49K29Hi0K78Y6T\/wgOmX9lLZ3X\/E2kaOWGFYChykiv1cnbt44610Wh+BtT0y28DRTz2jHQPP+1bHY798bKNny88sOuK1tc8NXsvijS\/EehyQQ38B8i8SZmVLq1PVSQD8ynleOtAGLdfEa6ebVLjTIdIaw02aSFo7vUPKuLpo\/v8AlrggDOQu772O1X38Z6jqmr6dYeHNPtp0v9LXUluruVkSJGbGGVQST04Hv6VlP4G1bTZdUttL0zw5e297cSXFveagn760Mh3EFdjeYASdvzD3pZrfWtM+I1jBowsbuW38PKk0dx\/o6TATYyuxSEOecbcckcUAbdh4p1jUtJuTaaJbSaxY3rWd7aSX3lRqQu7ekmwllIKkfKOp9Kv6TqXiW5vlj1PQLKytipJmh1Pz2B7Db5S\/zrO0vwNazaber4otrLUrvUL439wmzdFHJt2KqbucKoxk8nJrT0vwX4a0W9W90zRLK0uVBUSwxBWAPXmgDlr+11eH4u+G59R1MTw3Av8A7PaRR7I4Y1RdpPOWchuT+QxXo9c9quhXN74x8PaxE8It9NS6WZWJ3sZVULt4x\/Cc5IroaACsDSv+Qh4q\/wCv9f8A0jt636wNK\/5CHir\/AK\/1\/wDSO3oA36KKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACsTQ\/+Qx4m\/wCwkn\/pJb1t1iaH\/wAhjxN\/2Ek\/9JLegDbooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACoz\/wAfKf7jfzFSVGf+PlP9xv5igCSiiigAooooAKKKKACiiigAooooAKKKKACiiigBKKKy9S1eOzTamXlLbAFGSW67VH8TcdPTk4FS5WDzE1+8gttJuFlkCtJGUUYzkngcfUge\/TvXNW+s20mofbfEF01tLAc29i6FdpxjcR\/EfTsOvPGN7TNIladdR1T57vO6OIkEQZGOv8TYzz0GSFA5zt4Gc4qfZ68z3Ju2tNjL0aWa7E9\/LC8CXDDyopF2uEHALDsT6fStWilqkktiteoVieDf+RG8P\/8AYNt\/\/RS1t1z\/AITkaLwDoTrE8rLpluRGhGW\/dL0yQP1qgN+isz+073\/oB33\/AH3B\/wDHKP7Tvf8AoB33\/fcH\/wAcqOdFcjNPtVFtZ01bn7O19bibds2GQZ3en19qytY1O8OnFH0m9hWSSONn8yEHazAEAiTIJzjPbOasi4lW0+zDw5dCDbt8sNBtx6Y8ypdTWyLVOyu\/zRtZornNI1K8SzMQ0m9lSKR40bzIc7QxABJk5I6fhWh\/ad7\/ANAO+\/77g\/8AjlONRNXFKm4uxqUmazP7Tvf+gFff9\/IP\/jlWV\/4mFlJHdWkkKSqyPFKVJKng\/dJH61SkmQ4tFuivI7\/4eeE4viZo2mposK2c+n3MskQd8MysgU9e2TXT6n4hsfh29nYz6eLbw6bdxb3EJZyky5byiv8AtDO055PFUI7Wis3Q7rUL7R7a71OyWyu5l3tbK+4xA9FJwPmxjPvXn3hrwX4e8Sax4uvNX0xLq5j16aJJGdgVUJGQBgjux\/OgD1OivMNN8QP4Ik8V6TI1zqNrpU1qdNikl3St9pGEg3t2DjgnOFNdBbeI9csNe07TPEenWMK6nvW1nsbhnVJFXd5bhlXqoOGHp0oA6qCCK2iEUEaRRr0RFAA+gFCQxJLJKkSLJJje4UAtjpk9689t\/H2vv4YHiufQrSPRI2PnotyxuPLD7WlVdu3A5O0nJA7VqTeJdfvvEuq6RoWl2EqaekMhurq5ZUk8xNwVQqnnrz0H40AdnRXmi+IIfEHizwLqoja1DR6ms8LtkwyIqq6k+zKefxpg+Kc39mf8JB9n0oaHv\/1H2\/8A07yd23zfLxj\/AGtmc470AenUVx48Ta1feM7zRtL0u0ks7Brdrq7muGXMci7vkUKct1Ppx70nxDvbyDS9L02xuZLWTV9Th097mI4eKN9zOVPZtqED60AdjRXE3Hwy0SBIJ9AQ6NqlvIskd\/Dud2wwLLIC37wMMg7j3qe48R65qGu6lp3hzTrKZNMKpcz31w0YeVlDeWgVW6KRlj3PSgDr6K4H\/hYF5fweH10fR0kvNWe5haG6n2C1lh++HYKcgEN0HOB61XXxr4rks9a8vQdMNzoTML4teuI5gE34h+TOSv8AexjjrzgA9B8mLz\/tHlJ523Z5m0btuc4z6ZqWuKm8Xatfa1p+maFpltIb3S01Lz7uZkSFWbGGCglj04HqfStfwxr8mu2t4l3araahYXT2l3Asm9Q6gEMrYGVZWUjgdfagDeorL8S\/8irq\/wD15Tf+gGvP\/hhqV5ommaR4c1aZpIdQ06O+0m5f+MNGrSwH\/aQtkf7J9qAPVKK8Bj\/sx\/BnwwXWrWa704\/afOhjieVn\/dtt+VPmPzYPHpW3o8WhXfjHSf8AhAdMv7KWzuv+JtI0csMKwFDlJFfq5O3bxx1oA9jqLyYvP+0eUnnbdnmbRu25zjPpmvP7r4jXTzapcaZDpDWGmzSQtHd6h5VxdNH9\/wAtcEAZyF3fex2q+\/jPUdU1fTrDw5p9tOl\/pa6kt1dysiRIzYwyqCSenA9\/SgDtqK5Cw8U6xqWk3JtNEtpNYsb1rO9tJL7yo1IXdvSTYSykFSPlHU+lX9J1LxLc3yx6noFlZWxUkzQ6n57A9ht8pf50AdBRXkXibSLTws+nPGdQbXpr+GWXxJMGEUatL8yyvnaFK5QJ93leldL4n8O+H0uNT8ReLrn7VYrGi20MxYLagLgiMBuXducgbugFAHcVgaV\/yEPFX\/X+v\/pHb1F4Ci1SHwPpSaz5324RksJzmRU3EoHP94JtB981LpX\/ACEPFX\/X+v8A6R29AG\/RRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFZutX15YWHmWGnS6hduwSOFGCjJ7sx+6ox1oA0qK4rwv4hv8A\/hBtR1jWnWe5s5btpRGMLiNm+VfbjArMnvfEej+FLXxhd61JcMRFcXenGKMQiGQrlUwNwZQ3Xcc4oA9IrE0P\/kMeJv8AsJJ\/6SW9bdclZa\/o2la94kg1HV7CzmbUEcR3FykbFfstuM4Y9OD+VAHW0VWgv7O7gSe2u4JoXGUkjlVlYexFS+fD\/wA9Y\/8AvoUASUVH58P\/AD1j\/wC+hR58P\/PWP\/voUASUVH58P\/PWP\/voUefD\/wA9Y\/8AvoUASUVH58P\/AD1j\/wC+hR58P\/PWP\/voUASUVH58P\/PWP\/voUefD\/wA9Y\/8AvoUASUVH58P\/AD1j\/wC+hR58P\/PWP\/voUASUVH58P\/PWP\/voUefD\/wA9Y\/8AvoUASUVH58P\/AD1j\/wC+hR58P\/PWP\/voUASUVH58P\/PWP\/voUefD\/wA9Y\/8AvoUASUVH58P\/AD1j\/wC+hR58P\/PWP\/voUASUVH58P\/PWP\/voUefD\/wA9Y\/8AvoUASUVH58P\/AD1j\/wC+hR58P\/PWP\/voUASUVH58P\/PWP\/voUefD\/wA9Y\/8AvoUASUVH58P\/AD1j\/wC+hR58P\/PWP\/voUASVGf8Aj5T\/AHG\/mKPPh\/56x\/8AfQqMzw\/aE\/ep9xv4h6igCxRUfnw\/89Y\/++hR58P\/AD1j\/wC+hQBJRUfnw\/8APWP\/AL6FHnw\/89Y\/++hQBJRUfnw\/89Y\/++hR58P\/AD1j\/wC+hQBJRUfnw\/8APWP\/AL6FHnw\/89Y\/++hQBJRUfnw\/89Y\/++hR58P\/AD1j\/wC+hQBJRUfnw\/8APWP\/AL6FH2iH\/nrH\/wB9CgB9ISFGSQAO5qCW8t4U3vKnpgHJJ9AKxr7UWllW3QRy3TgNHal\/ljH9+Ujt147479RDk2+WOrB2SuybUNUYKkVujyNKD5MSHEk+OuM\/dX1c+vHaptM0k2ri6u2SS9KlcoMJEpOdiDsOmT1YjJ9A7TrOCyLzS3S3F5NjzrhyAWx2A\/hUc4H8zk1f8+H\/AJ6p\/wB9CqUVHzfcnWWr+4loqP7RD\/z1j\/76FHnw\/wDPWP8A76FMokoqPz4f+esf\/fQo8+H\/AJ6x\/wDfQoAkrD8G\/wDIjeH\/APsG2\/8A6KWtjz4f+esf\/fQrH8G\/8iN4f\/7Btv8A+iloA3KKKKAIbi3jurd4JkDxuCrKe4rP\/sy68vyP7VuBBwOg8zHpv6\/j1961arfb7bdjz0POM54z9aznyr4nYuLkth9tbxWlvHBCgSNBhQOwqaiitEQ2LRRRQBzt5oV1cePNL11HhFraWU9u6kneWdkIIGMY+U96peJfCM\/i\/Vkg1WdU0G3hYxQQSMJZLlgQJG4wAgJ2gE\/McnpiuvooAyfD0Gr22iQW2uTwXF9CCjXEBOJlH3XIIGGIxkcjOea5Sz0XxxoWqa42kp4dms9R1GS9Q3c8yyLuVVwQqY\/h9a9BooA8\/b4fXd1oGrG91KJ\/EOpXUV694kREUckRUxIqk52LtxzzyauxaN4i1rxDpOoeIE060t9KZ5Y4bKZ5TPMyFNxLKu1QGJA5PrXZ0UAeQeGNH8TeIvhpbaHv06PR7wyLJeeY\/wBoSHzW3II9u0t1Abd07ZrWtf8AhILb4heKxoMGmzRBLOMxXkrxeW3k\/KwKq2R6rx0GDXokEEVtEIoI0ijXoiKAB9AKEhiSWSVIkWSTG9woBbHTJ70AcPpvgG40+68NM91DcJYLfNfucq00lzgkoAOm7d3HGKy7PwDrWn6ZFoNtY+GzbxMFj1mWAPciHdnBiaMq0mOMlsd69RooAwNJ0SfT\/FGvak7Q\/Z7\/AOz+SiE7l8uPadwxge2KXxZ4d\/4STRhbR3Btby3nS6s7kLu8mZDlWx3HUEehNb1FAHCXumeOfEFsml6nNpOmWTOpurvTp5WmlUEErGCq+XnGCSTj3qw+i+ItD1\/Vr7w+mm3dtqrpNJDezPEYJggQsCqtvUhQSODmuzooA4bR\/BF3pN14al+1wznT3vJr2UgqZZZ+SUXGMbie44Aq1B4Xvo08ZK0tv\/xO5Ga2wzfIDAsfz8ccjtniuvooA8wS117SfHelWulpYXNzaeGY4Z4biVo0k2ybcq4UkEEd15GeldFouh63ounXtysmnz6zqV+Lu73F1hRTtUohA3Hai8Ejk9a6jyYvP+0eUnnbdnmbRu25zjPpmpaAKOr2kmoaLfWcRUSXFvJEpfOAWUgZ9ua5u48FyXnw70nQpLhIdV0y1t\/s15ESRDcxIFDrkZK5BHTlSa7KigDz7Q\/A2p6ZbeBop57RjoHn\/atjsd++NlGz5eeWHXFa2ueGr2XxRpfiPQ5IIb+A+ReJMzKl1anqpIB+ZTyvHWurooA80fwNq2my6pbaXpnhy9t724kuLe81BP31oZDuIK7G8wAk7fmHvSzW+taZ8RrGDRhY3ctv4eVJo7j\/AEdJgJsZXYpCHPONuOSOK9KqLyYvP+0eUnnbdnmbRu25zjPpmgDk9L8DWs2m3q+KLay1K71C+N\/cJs3RRybdiqm7nCqMZPJya0tM8FeGtFvBeaZollZ3IUqJoYgrAHrzW\/RQB55f+GvGGraBL4V1G70+fT5W2S6s0jG4eHdnb5Wzb5mON27HfGabrPhvxjeeNP7Xhj0C9sLQKNNt764mAt2A+aQqqEGQn+LJwOlei0UAUdKbUzpsJ1hLRL\/5vNW0ZmiHzHG0sAemO3XNZ2lf8hDxV\/1\/r\/6R29b9c7pciLqPilWdQTfrwT\/06W9AHRUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQB594csJNU+HPiHT4sebc3GoQpn+8zuB+tZep+ILPXvhzbeGbNy2vXUUNk1iVPmQupUOXH8KqFY5Neq0zYocuFG4jBOOaAFRAiKq9FGBWNof\/IY8Tf8AYST\/ANJLetusTQ\/+Qx4m\/wCwkn\/pJb0AbdFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAVGf+PlP9xv5ipKjP\/Hyn+438xQBJRRRQAUUUUAFFFFABRRRQAUlFMkkSJC7sFUdSTSbtqwI7uVoLaSRAC4Hyg9M9qxrjU9Jh8xZLvzrlGCyMJGAVvdhwuP0o8Qak0Gjzzlkt4SuA0owz5PRRkHpn3qC0jvZtMFnpdoLC2OMXMiBGZO5VOSG92H4VhKLm+Zr3e4OfL7qevYek93NMbazYTXgwDeMN0UKEA5H95vbvjJwMVsWGnw6fEyRF3d23yyyHLyN6sazPDFulkmoWUDK1vb3OxGUdTsRn57\/MW\/l2rfxW0I+zjyr+vUm\/O+Zi0UUVRQUUUUAFFFFABWD4R3jwJoJQBm\/s23wGOB\/ql71vVh+Df+RG8P8A\/YNt\/wD0UtIDQ8y\/\/wCfa3\/7\/n\/4ijzL\/wD59rf\/AL\/n\/wCIq5RWfJL+Z\/gO\/kZV7JfG22tBCqsyqxExPBIz\/DUwN4Itgs7bZjG3zzjH\/fFXHRZUZHAKsMEHvVb7E+Nn2ufyv7uRnHpnGf1zWUqc1K6bf3FKSasVrGS8FttSCFkVmVSZj0BIx93t0\/CrPmX\/APz7W\/8A3\/P\/AMRVlEWNAiAKqjAA7U6rhTlGKXMxOSbvYqebf\/8APrb\/APf8\/wDxFTxGVkBlRVb0Vtw\/PAqSitIxaeruK58++HJfh5b+G1n8R6LdzX\/2qdZLlbOdlYmdwgDr8ucbRwetd3oV\/rnhP4fLf3mnXMtpDeySG3uZC1zbWBc7Seu5kXkqe30xW94X8KnTfA\/\/AAj2sJb3Su1wJlQlkZJZXfHIB6MPxrIg8J+KE8KDwrLqls1gLgwm9Er\/AGk2H\/PPG3HmY+TOcBferEdB4a8RSeJnvby2tgujpJ5VndFjm6x99wuOEzwD3wTXM+IdA0vxJ8XrKy1ezW7tk0GWVY3YgBxOgB4Pox\/Ouh8L+H7zwzPfafDLC+gF\/N0+Iuxltt334uRgpnlecjJFUdf0TxN\/wm1t4h8PjSZNmnPZSR6hLIn3pFfI2Kf7ooAzRp9t4E8faHaaS8lvpGsRXMdxZvKzxxPFH5iyIGJ28AqccdKk\/wCE71waAPFh0S1Hhz\/WbftLfa\/Izjzdu3b0+bbnOO9X7DwzrGpeIY9b8Vz2Dy21vJb2lnYB\/KiEmBI5ZsMzEAL0AArK\/wCEQ8Unwl\/whRn03+x9n2f+0vMf7R9mz93ytu3ft+Xdux3xQBqy+JfEF\/4l1XSNC0ywlj09IZDdXdyyJJ5ibgoCqTnrz0H41V\/4T69v7fw8NH0hJLvV2uIniuZ9i2skPD7iAcgEN0HOB61TtB4gtfiD4rTQINOmhCWcZivZni8s+T8rBlVtw9V47YNU7nQtU8N6r4IsNNntrnUke\/nme4BSO4d13yDKglRljg4OMDg0AL4z8QXt74N8SaPqtrDbapp8lk7\/AGeQtHNHJOm10yAR91gQehHXmtvVPG90PEN\/pOkjSEOnBBcS6nfeQHkYbgkYAPQEZY8AnGKpal4I1rXNM1+5vpbGLV9VNqkcUcjNDbwwSBwu\/blmPzEnb1xU2o+D7+28Tanqmm6VoWqw6mUkeLVMq1vKqhcqwR9yEAErxz0oAD8QrzUrLw62gaVFPda0biMR3NxtS3kh+\/uZQcqMNyOuBjrXTatql3pPgy+1a4hiF7aafJcyRIxZPMSMsVB4JGRjNY1h4S1C1vPC881xZyNpi3RuzDCIFdpVwPLRVxgHjnBwM8musvLSC\/sp7O5QPBPG0Uin+JWGCPyNAHBaF8PND1jw1Z3+uwPqOr31ulxPqEsreaHcbv3bZ+QLnChcdK0ZdS1bw9Bofhe2kTWtduIpD9ru2MKCKPGZJMbiT8yrgfePORVWw0zx5oGlRaHp0miXlpbp5NpfXckiSxxjhfMjVSGIHHBGcc4pT4M1bSF0G+0m\/TUNU0yKWGc6jIyi8SUhn+YBihDDK8EY4oAbe+Pr7SNJ17+09KiGr6MsMjQ28xaK4jlbaroSAf7wwR1HvSv4p8Wwa\/Boc2haWL2\/t3ubRlvnMcKoRuEp2ZJG4fdHJP41BfeCtb1mw8QXd\/LYxatqq28MUMTs0NvDE+4Lv25ZjliTt9K6O70O6uPHema4rwi1tbKe3dCTvLOyEEDGMfKe9AGC\/j3VItBjmOixS6qNZ\/seW2juP3bSc\/Mrlfun5TyOMn0rX0nxBqY8Sv4f12ztYLx7U3dtNaSs8cqBgrL8wBDKWX6g54rlPEmi6ppFnb+RPai7vfFyXlsTuZAGB2q\/APbnHrxXS6douuXXiZ\/EWsrYQXEFm1pZWltK0qLuYMzu5VTklVGAOB70AddXjOi6je+GPGHiTxBNMz6DPr0lhfoelqdsZjn\/AN3c5VvYr6V65p5vjp8H9pC3W92DzhbFjHu77S3OPrWBofhY2sHia11Rbe5tdY1Ka5EakkGGREXa2QOflPT25oA4Hxa5\/sz4skN3s8YP\/TCOqWrxeBZbKSw8M+H9Vg8TTxN\/Z3kW09vIJB0fc2AFU4JPpW7B8MNZtfCni\/Rv7Strl9UMUdjNM7ZWKNVVBL8vUKoHGc4rtfF\/hkeJdC+zQzfZtRt2FxYXY6wTryrZ9Ox9jQBkan4vvtN1Gx8PRNpZ1VbFLm9utQuvJgT+HC4GWZmDHHGAKr\/8LHuJdBiubPS4bnUv7WXSZbWO6BjMhGQySY5Qjac46E+lLqHhPWJtZtfELadoWo38lilrqFldM3lFlJIkikKMR94jBXpWXr2jaxpOgaKm\/S4tUn8RwTRpbWwS2i4bbH8oVmXjlj83J9qAOntPEetRaxc6Hqmm2X9pmya8sjbXLeTcBTtKFmXKMGK84PDZ7Yp0es+MmkRZPCumohYBmGtZIHrjyeahtPDmpav4gm1fxPBYKn2F9PhsLaRpk8uQgyM7Mq5LbQMY6VZj+HPgyKVJI\/DWmq6EMrCAZBFAGB8XLbV5PDU91FqYttMgNvvtoo\/3lxI06qQ754QAj5QOT1OK9HrnfG2h3PiPwpd6VZvEk8zwsrTEhfklRznAJ6Ke1dFQAVgaV\/yEPFX\/AF\/r\/wCkdvW\/WBpX\/IQ8Vf8AX+v\/AKR29AG\/RRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAViaH\/yGPE3\/AGEk\/wDSS3rbrE0P\/kMeJv8AsJJ\/6SW9AG3RRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFRn\/j5T\/cb+YqSoz\/x8p\/uN\/MUASUUUUAFFFFACUVXvLuKxtJbqYkRxrk7Rkn2A7mss3Wslg4\/s6JjyLSR2L\/i44B+in6mnbuS5WNyiuej8Vw3Fp5tnaT3EoJR1TARHBwV3nhjkfw5PQ45qW3tNWvpfO1G6FvBnK2ttkf8AfT\/eP4Y+lTL3dOoRkparYn1DXrWxmW2TdcXbfdghUu31OOg96px2mt384luXiso+qgYllX6D7in3+etyG3it02QxKi9woxmpalLrLcbTfXQzbbRLC1dpvJEtwwIaec+ZIwPUbjyB7Dj2qFdCaH91a6nfW1r2t0ZCqj0UspYD2B47YrYorXmZPs49iC1tYrO3SCBNsa9BnPuSSepJ5JPWp6KWpK20QUUUUDCiiigAooooAKw\/Bv8AyI3h\/wD7Btv\/AOilrcrnvCjTL4A0I28cckw0y32pI5RSfKXqwBx+RoA6DNFZH2nxD\/0CtN\/8GMn\/AMZo+0+If+gVpv8A4MZP\/jNOzHymsWAGScCskeJtILsBeDYpwZvLbyh\/20xt\/WsvX7nWzpZS40\/T0tpJY0mK6hJyjOAQT5Qwpzgn0Jq+H1ww+V\/Y+leVt27Pt74x6Y8jpVKOmpXLpqbgIIyKK5Pw\/ca2umeXbadp720c0scRbUHAChyAB+5OQMYB9q1ftPiHP\/IK03\/wYv8A\/GaTi0xONnY16KyPtPiH\/oFaZ\/4MZP8A4zV60e7eAG8ghhlz92GYyLj6lV\/lSaE0WqK+fPDkvw8t\/Daz+I9Fu5r\/AO1TrJcrZzsrEzuEAdflzjaOD1ru9Cv9c8J\/D5b+8065ltIb2SQ29zIWubawLnaT13Mi8lT2+mKQj0iiue8NeIpPEz3t5bWwXR0k8qzuixzdY++4XHCZ4B74JrmfEOgaX4k+L1lZavZrd2yaDLKsbsQA4nQA8H0Y\/nQB6PRXnI0+28CePtDtNJeS30jWIrmO4s3lZ44nij8xZEDE7eAVOOOlSf8ACd64NAHiw6Jajw5\/rNv2lvtfkZx5u3bt6fNtznHegDvUhiSWSVIkWSTG9woBbHTJ70NBE8scrxI0kedjlQSueuD2rkZfEviC\/wDEuq6RoWmWEsenpDIbq7uWRJPMTcFAVSc9eeg\/Gqv\/AAn17f2\/h4aPpCSXertcRPFcz7FtZIeH3EA5AIboOcD1oA72ivK\/GfiC9vfBviTR9VtYbbVNPksnf7PIWjmjknTa6ZAI+6wIPQjrzW3qnje6HiG\/0nSRpCHTgguJdTvvIDyMNwSMAHoCMseATjFAHc0V5+fiFealZeHW0DSop7rWjcRiO5uNqW8kP39zKDlRhuR1wMda6bVtUu9J8GX2rXEMQvbTT5LmSJGLJ5iRlioPBIyMZoA2qK870L4eaHrHhqzv9dgfUdXvrdLifUJZW80O43fu2z8gXOFC46Voy6lq3h6DQ\/C9tImta7cRSH7XdsYUEUeMySY3En5lXA+8ecigDs6K4G98fX2kaTr39p6VENX0ZYZGht5i0VxHK21XQkA\/3hgjqPelfxT4tg1+DQ5tC0sXt\/bvc2jLfOY4VQjcJTsySNw+6OSfxoA7iWGKbZ5sSSbHDpvUHaw6Ee9S1wD+PdUi0GOY6LFLqo1n+x5baO4\/dtJz8yuV+6flPI4yfStfSfEGpjxK\/h\/XbO1gvHtTd201pKzxyoGCsvzAEMpZfqDnigDqKKK8Z0XUb3wx4w8SeIJpmfQZ9eksL9D0tTtjMc\/+7ucq3sV9KAPZqK8Y8Wuf7M+LJDd7PGD\/ANMI6pavF4FlspLDwz4f1WDxNPE39neRbT28gkHR9zYAVTgk+lAHulRSwxTbPNiSTY4dN6g7WHQj3ri9T8X32m6jY+Hom0s6qtilze3WoXXkwJ\/DhcDLMzBjjjAFV\/8AhY9xLoMVzZ6XDc6l\/ay6TLax3QMZkIyGSTHKEbTnHQn0oA9BorkrTxHrUWsXOh6pptl\/aZsmvLI21y3k3AU7ShZlyjBivODw2e2KdHrPjJpEWTwrpqIWAZhrWSB648nmgDq6K4LxnoGmob\/XdXsNT10tGsVpZWqMxtcKxLRhT8pY8l+o4FM0jQz4q8G+F5tZ1lb+yt7YyXkaMTHePtwvmNkZCEHIPU9elAHoFYGlf8hDxV\/1\/r\/6R29YPw7SL7f4gl0fzF8MNPGumgsTGWVSJmiz\/wAsy23GOMhsVvaV\/wAhDxV\/1\/r\/AOkdvQBv0UUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFYmh\/8AIY8Tf9hJP\/SS3rbrE0P\/AJDHib\/sJJ\/6SW9AG3RRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAJUZ\/4+E\/3G\/mKg1G9FhZPPsMj5CpGvV2JwBWDdHX7VRfXGp2MO1SxhdQsQXjILHnPvnFK6vYTdjqcUVzsPiuG4tomitJ3ncfMikBFPs7YDD3XP07Uqz69fZwI7OM\/wDPJC7gf7z7QD\/wBv8AE5o9xX7I3ndY0Lu4VQMkk4ArMbxDYyFlsjLfyA422kZkGfQv9xfxYVXTQPNcSXZFxIDndcsZ8H1CnCL+C1pLYJtAkllkAAAUttUfguBS5\/5U3+CHyye7sYOt3GpXdgyNDb2oEkciRNOGmdldWUcfKvIGeWqgk2iQELc3b3V5Hz5N2QChPXEKcuefRvr3rsRZ26wvEsKBHGGUDGaiFvdqNguxs7Fo8v8AnnH6VnJzcveWnl\/wQ5FbR6+Zn+HbRY4ZppLbypTKQFbGVXAAHHTgDj8OcVu1HFEsMYjTgDPXqT3NSVVOPLGxTYtFFFaCCiiigAooooAKKKKACiiigAooooAKw\/Bv\/IjeH\/8AsG2\/\/opa3K5208MXNhZQWdp4m1eK3t41iijCWp2IowBkwknAHegDoqKxP7C1H\/oa9Y\/79Wn\/AMYo\/sLUf+hr1j\/v1af\/ABigDWlijniaKVFeNwVZWGQR6VlL4eCRiBdU1EWo4EHmrjb6b9u\/H\/AqT+w9R\/6GvWP+\/Vp\/8YpP7D1H\/oa9Y\/79Wn\/ximm0NNo14IIraCOCCNY4o1CIijAUDoBUtYf9h6j\/ANDXrH\/fq0\/+MUf2HqP\/AENesf8Afq0\/+MUhG5RWJ\/YWo\/8AQ16x\/wB+rT\/4xTJNE1FFBHirV\/vAf6q07n\/rhQBT8L+FDpvgf\/hHtYS3uldrgTKhLIySyu+OQD0YfjWRB4T8UJ4UHhWXVLZrAXBhN6JX+0mw\/wCeeNuPMx8mc4C+9dN\/YWo\/9DXrH\/fq0\/8AjFH9haj\/ANDXrH\/fq0\/+MUAVPC\/h+88Mz32nwywvoBfzdPiLsZbbd9+LkYKZ5XnIyRVHX9E8Tf8ACbW3iHw+NJk2ac9lJHqEsifekV8jYp\/uitn+wtR\/6GvWP+\/Vp\/8AGKP7C1H\/AKGvWP8Av1af\/GKAMew8M6xqXiGPW\/Fc9g8ttbyW9pZ2AfyohJgSOWbDMxAC9AAKyv8AhEPFJ8Jf8IUZ9N\/sfZ9n\/tLzH+0fZs\/d8rbt37fl3bsd8V1v9haj\/wBDXrH\/AH6tP\/jFH9haj\/0Nesf9+rT\/AOMUAcjaDxBa\/EHxWmgQadNCEs4zFezPF5Z8n5WDKrbh6rx2wap3Ohap4b1XwRYabPbXOpI9\/PM9wCkdw7rvkGVBKjLHBwcYHBrtU8OXiSySp4m1VZJMb3EFmC2OmT5HNDeHLx5Y5X8Taq0kedjmCzJXPXB8jigDmdS8Ea1rmma\/c30tjFq+qm1SOKORmht4YJA4XftyzH5iTt64qbUfB9\/beJtT1TTdK0LVYdTKSPFqmVa3lVQuVYI+5CACV456V0v9haj\/ANDXrH\/fq0\/+MUf2FqP\/AENesf8Afq0\/+MUAY9h4S1C1vPC881xZyNpi3RuzDCIFdpVwPLRVxgHjnBwM8musvLSC\/sp7O5QPBPG0Uin+JWGCPyNZf9haj\/0Nesf9+rT\/AOMUf2FqP\/Q16x\/36tP\/AIxQBzthpnjzQNKi0PTpNEvLS3TybS+u5JEljjHC+ZGqkMQOOCM45xSnwZq2kLoN9pN+moappkUsM51GRlF4kpDP8wDFCGGV4IxxXQ\/2FqP\/AENesf8Afq0\/+MUf2FqP\/Q16x\/36tP8A4xQBy194K1vWbDxBd38tjFq2qrbwxQxOzQ28MT7gu\/blmOWJO30ro7vQ7q48d6ZrivCLW1sp7d0JO8s7IQQMYx8p71N\/YWo\/9DXrH\/fq0\/8AjFH9haj\/ANDXrH\/fq0\/+MUAcP4k0XVNIs7fyJ7UXd74uS8tidzIAwO1X4B7c49eK6XTtF1y68TP4i1lbCC4gs2tLK0tpWlRdzBmd3KqckqowBwPer0vhy8m2eb4m1WTY4dN8FmdrDoR+461L\/YWo\/wDQ16x\/36tP\/jFAGjp5vjp8H9pC3W92DzhbFjHu77S3OPrWBofhY2sHia11Rbe5tdY1Ka5EakkGGREXa2QOflPT25q9\/YWo\/wDQ16x\/36tP\/jFH9haj\/wBDXrH\/AH6tP\/jFAHBwfDDWbXwp4v0b+0ra5fVDFHYzTO2VijVVQS\/L1CqBxnOK7Xxf4ZHiXQvs0M32bUbdhcWF2OsE68q2fTsfY1P\/AGFqP\/Q16x\/36tP\/AIxR\/YWo\/wDQ16x\/36tP\/jFAHNah4T1ibWbXxC2naFqN\/JYpa6hZXTN5RZSSJIpCjEfeIwV6Vl69o2saToGipv0uLVJ\/EcE0aW1sEtouG2x\/KFZl45Y\/Nyfau5\/sLUf+hr1j\/v1af\/GKil8OXk2zzfE2qybHDpvgsztYdCP3HWgDOtPDmpav4gm1fxPBYKn2F9PhsLaRpk8uQgyM7Mq5LbQMY6VZj+HPgyKVJI\/DWmq6EMrCAZBFXf7C1H\/oa9Y\/79Wn\/wAYo\/sLUf8Aoa9Y\/wC\/Vp\/8YoAoX9v4usdavbnSZLC\/srtU8uC+neL7G4GCV2o25D1I4Oa5nUfAXiWPwvo3h\/TLrTbmxheSfU47uWWAXbu5fZ+7VsR7mb5c88V2v9haj\/0Nesf9+rT\/AOMUf2FqP\/Q16x\/36tP\/AIxQAzw3H4hhhkh1u10a3ijVFtU0ySRlAGcghlXAHy4x70mlf8hDxV\/1\/r\/6R29Sf2FqP\/Q16x\/36tP\/AIxU1jpQ0u1vibu4u5rqQzSzXGwMzbFQcIqqPlRe1AGrRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFAENxcRWlrNczuEhhRpHY\/wqBkmuNj8banHYWmuX2iR2+gXToqzC53TRI5ASR0242nI6MSM1veLoZLjwXrkMQLSPYTqoHUkxtxXKeJp4bj4IJ5LB\/tFjbRwqvVmYoAo98\/yoA7y5uDCUREMkr52pnGfUk9h\/jWJpFw8GqeIHkjBX+0E85lP3G+ywdPVcbefrxWrLm3uYZ5GLIIzG7AdDwc+w4rLvdMsLia4ntp7j7TdEMyQXTqrttVAzKpxjCqCT2FctScldp6roaRSdjoqKx\/+Edt+97qefa\/lH\/s1H\/CO2\/8Az+6p\/wCDCb\/4qtVKf8pNl3Niisj\/AIRy3\/5\/dT\/8GE3\/AMVR\/wAI5b\/8\/up\/+DCb\/wCKo5p\/yisu5oTziEIFUvI5wqA9f\/rVEbi5hXzJ4k8sfeMbklfwI5rLuNCggaORrrUTEMh2+2yll987s49fw9KZLpdh5ZEWoahNIw+SNdQlbd+G7pWE6s7vpY0UVZHRAgjIIxS1ip4cgCKGvNRBA5230oGfYBuKf\/wjtt\/z+6p\/4MJv\/iq2jObWsSLLua9FZH\/COW\/\/AD+6n\/4MJv8A4qj\/AIRy3\/5\/dT\/8GE3\/AMVT5p\/yisu5oXE4gQHaXdjtRB1Y1C091CvmSwoyDlhGxJA9sjmsy50KGHZILrUmjB+f\/TJWYDHUfNn8qbJpmnrHlNQ1GVj91E1GUlj6ferCVSd3fSxooqy6nQI6yIrKcqRkEd6dWJF4chESCS81EOFAYJeyqucdgG4FP\/4R23\/5\/dT\/APBhN\/8AFVtGc2k3EhpX3Niisj\/hHLf\/AJ\/dT\/8ABhN\/8VR\/wjlv\/wA\/up\/+DCb\/AOKp80\/5RWXc0p5hBHuIJYnCqvVj6VCZrtFLvboVAyVjky36gZrLudAhjCSLc6k6q2WH22VmAwRkfN157ds0g0\/TY1Ei6hqEmfuoL+Vix9AN1YSqTu09C1FW7m6kiyRq6EFWGQafWZa6REtqizNOWx8w+0PgH0HNRy+GtLmJ3xznPXF1KP5NW0ZVGk2l95LUU9GaxYKMsQB71XbUbJDhry3Uj1kFY7+B9AkbLW1xn1F7OP5PTB4F0IdI9QH01S5\/+OVd59l9\/wDwBWXct6hrlmkAS21G0EzsF3eYreWD1bH+PGcVmPLo0X7yHWL43CjPmpJLLn6jBBHtip5fBem+RIttLqMUpXCu2pXL7fwMh\/PrVD\/hHLkr5KWt\/FIflM\/9v3ZQe+3zM\/hWUue7uOyNa38Qs1sn\/EuvriUACQww4Xd3xuIp51TWJQTbaE3\/AG8XKx\/yDVHH4SsEjRWudVLADJTVbmME9ztWQAZ61J\/wimnf8\/Gsf+Dm7\/8AjtaQU+VczJaV9BGbxNMTsTTbYdg++U\/oVpp03X5zmXXhB3229qn\/ALNmn\/8ACK6d\/wA\/Osf+Dm7\/APjtH\/CK6d\/z8ax\/4Obv\/wCO1ST6sVkU73w5dm0kcapf3M6jKo05jUn\/AID3\/lWAsmjw67ZldMvZiIJi6TQsZRMGi2KGPHQydDiur\/4RXTv+fnWP\/Bzd\/wDx2oz4P0kzrMZNVMqqUV\/7Xu9wU4JGfM6HaPyFR7PVtPcLIv6RbPb2KmeFI55CXkROi5OcZ9ulaFYn\/CK6d\/z8ax\/4Obv\/AOO0f8Ipp3\/PzrH\/AIObv\/47VQgopJFN3NyisT\/hFNO\/5+dY\/wDBzd\/\/AB2j\/hFNO\/5+dY\/8HN3\/APHasRt0Vif8Ipp3\/PzrH\/g5u\/8A47R\/wimnf8\/Osf8Ag5u\/\/jtAG3RWJ\/wimnf8\/Osf+Dm7\/wDjtH\/CKad\/z86x\/wCDm7\/+O0AbdFYn\/CKad\/z86x\/4Obv\/AOO0f8Ipp3\/PzrH\/AIObv\/47QBt0Vif8Ipp3\/PzrH\/g5u\/8A47R\/wimnf8\/Osf8Ag5u\/\/jtAG3RWJ\/wimnf8\/Osf+Dm7\/wDjtH\/CKad\/z86x\/wCDm7\/+O0AbdFYn\/CKad\/z86x\/4Obv\/AOO0f8Ipp3\/PzrH\/AIObv\/47QBt0Vif8Ipp3\/PzrH\/g5u\/8A47R\/wimnf8\/Osf8Ag5u\/\/jtAG3RWJ\/wimnf8\/Osf+Dm7\/wDjtH\/CKad\/z86x\/wCDm7\/+O0AbdFYn\/CKad\/z86x\/4Obv\/AOO0f8Ipp3\/PzrH\/AIObv\/47QBt0Vif8Ipp3\/PzrH\/g5u\/8A47R\/wimnf8\/Osf8Ag5u\/\/jtAG3RWJ\/wimnf8\/Osf+Dm7\/wDjtH\/CKad\/z86x\/wCDm7\/+O0AbdFYn\/CKad\/z86x\/4Obv\/AOO0f8Ipp3\/PzrH\/AIObv\/47QBt1la\/dy2WnQyw7Q7XtpEcjPyvcRo36Mah\/4RTTv+fnWP8Awc3f\/wAdqKbwfpU6BJpNVkUMrgPq92QGVgyn\/W9QQCD6igDoKKxP+EU07\/n51j\/wc3f\/AMdo\/wCEU07\/AJ+dY\/8ABzd\/\/HaANuisT\/hFNO\/5+dY\/8HN3\/wDHaP8AhFNO\/wCfnWP\/AAc3f\/x2gDborE\/4RTTv+fnWP\/Bzd\/8Ax2j\/AIRTTv8An51j\/wAHN3\/8doA26KxP+EU07\/n51j\/wc3f\/AMdo\/wCEU07\/AJ+dY\/8ABzd\/\/HaANuisT\/hFNO\/5+dY\/8HN3\/wDHaP8AhFNO\/wCfnWP\/AAc3f\/x2gDborE\/4RTTv+fnWP\/Bzd\/8Ax2j\/AIRTTv8An51j\/wAHN3\/8doA26KxP+EU07\/n51j\/wc3f\/AMdo\/wCEU07\/AJ+dY\/8ABzd\/\/HaANuisT\/hFNO\/5+dY\/8HN3\/wDHaP8AhFNO\/wCfnWP\/AAc3f\/x2gDborE\/4RTTv+fnWP\/Bzd\/8Ax2j\/AIRTTv8An51j\/wAHN3\/8doA26KxP+EU07\/n51j\/wc3f\/AMdo\/wCEU07\/AJ+dY\/8ABzd\/\/HaANuisT\/hFNO\/5+dY\/8HN3\/wDHaP8AhFNO\/wCfnWP\/AAc3f\/x2gDborE\/4RTTv+fnWP\/Bzd\/8Ax2j\/AIRTTv8An51j\/wAHN3\/8doA26KxP+EU07\/n51j\/wc3f\/AMdo\/wCEU07\/AJ+dY\/8ABzd\/\/HaANuisT\/hFNO\/5+dY\/8HN3\/wDHaP8AhFNO\/wCfnWP\/AAc3f\/x2gDbrK8S3cun+FtXvYCBNb2U0sZIyNyoSP1FQ\/wDCKad\/z86x\/wCDm7\/+O1HJ4Q0qeJ4ZpNUlikUq8cmrXTK6nqGUyYIPoaAN+iiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACuctvA3h6z1GO8hsSGik82KIzOYYpDn5ljLbFPPYcdq6OigBOoxSBFXO1QM9cCvLvE8A0DTpbm+8R3j+LbiQyWSQXDqhJfCRrDnbsxgHcPXmt3UvtOv+NLTw\/dzz29nb6aL66jtZmjM0jPsCFlIbaMMeD6UrIDtqK5HwnNcWeva94dluprmDT3hltZLiQvII5VJ2FjydpU4zngiuupgFFFFACU1UVSdqgZ64HWn0UrIAooopgFFFFACYpojRWLBVBPUgU+ilZAFFFFMAooooASmhFDFgoDHqcc0+ik0gCiiimAUUUUAFFFFABRRRQAUUUUAFFcb4l0+1W+udV8Sa5Na6NHGkdrDBcSQBJOdzMUILN\/d\/lXPXWu6jD8H0ubm\/uY57q5W3trwsRMYWm+VyV53eWD70Aep0VwPhZ9F\/4SEQwax4kN8sZdbTVpplEidNwSQDdXfUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRWVr+tJoemG48tp7mRhDa26n5ppW+6g\/qewBPauL8K6pqml6H40vdVuvtl5p9zNI2Sdm5Yg21R2XPAoA9Jory+9srzQfA9t4wXVL+bWIkhu7nzLlzFOrld8ZjztC4Y4wBjAr05WDKGByDyDQA6iiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAoqK4uIbS2lubiRY4YlLyOxwFUckmvO9G1XVdU+Jtpd3MksOn3mmyy2dmSRtiDqFdx\/eblvYECgD0mivJdR1WzvvF3iF9R1XxHDa2UscMceltP5cKqg3u+wED5s\/lXpGhG0Oh2Zsb2W9tWjDRXMspkeQHuWPJoA0qKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooA871PXr7UfDV\/oeseGrybWpUkgWK3tGa3lJyEkWQ5VV5U5JyOaWKzv8Awlq+japewXV\/ANGTTbyS1jaV45EO4OVHzFTlhn6Zr0OigDkfCdvc3Wu694huLSa0i1B4YraK4TZJ5cSkbmXqu4seD6CuuoooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKAOWv\/El5pGu3VvqWlXM2luiNZ3NlbPOS2MOkgXODnpxjBrC0aDW9C8OXeoW2hlorjWGvF0xlzLBatjOxQeJON23tn1r0aigDg3mk8WeMtBvbKwvrez0ozSzXN3bNDuLptEaBgCfftxXeUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFAHM694YvNX1qy1Oz1ybT5rSNkjVbeOUAt95gHyAccZ9K5jw74W1qVfGdnqGo3LRXzzQL51osazs8agTggfhgfLXptFAHmF5dajrvgm28H\/2NqMGqyJDaXTyW7LDCqFd8nmfdZSFOME9a9OChVAAwBwBS0UAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFAGJ4o0FvEmjnT1vns1aVHZ0jV9wU52lW4Izjg+lciPDviSH4laZcTa7eXUMdk2+6NhGqbfMXMJKrtG71+96V6TRQBx58U3ekX2pWmqaDes4nZrN9PtHlS6jI+XLDID9jnFW\/Aej3Wh+ELOzvUEdwWkleJTkRb3LbPwz+ea6WigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACsvWdcs9DFibsS\/6beR2UWxc\/vHzjPtx1rUriPiR9zwp\/2Mdn\/7PQB29FcBcHWta+IuraLFrd1Y6Zb2sEpFsFEm5s8KxB256n6CsuXXdb0vwv41sn1SW4udFdBa3rqvmbXAIDcYJHr70Aep0V5vq8fiHw54YXxS3iG7u7yBY5ruykCfZpEYjciKFyuM8NnPHvW3Lq103xJ0yxiuHFhPpclw0WOC28YP5GgDq2YIpZjgAZJNVdN1K01ewivrGUTW0uTHIAQGwSO\/0rk767v73xxrOjDULiC0XR0mQQ7QyOXOWGQewxR8KLR4PAGnTNeXEwmjLLHIV2w4ZhhMAHH1JoA1ta8Ww6JeNby6PrV1tQOZbOxaWMD\/AHhVHS\/iFY6y9r9j0bXnhuXVEufsDeTycbi+cBR3NdRf\/wDIPuf+uTfyrmfhf\/yTbRP+uTf+htQB11Fec+A49e15E1rUPEN0YYLuaJLRFUJIisw+c45Of0UetZ+p6oqSXv8Aavj6e31dXk8my0jEkUIB+VGURlmPTOcUAerUV5Y3iTxBq2ieB5LXUBaXeqyPFcyrECDhSC23pngsB0z7Vr2h1Pw74\/07SH1q91Ow1K2mcre7WeORMHKsqjg+lAHeUV43B4gup2lXVvF9\/oXicStts7qNUslAY7RymGXbj5t1dNrdzrl5450jRbTVmsre5055bp7ZVbow5j3A4PYH0NAHfVT1TUI9K02e+limljhXcyQJvdvoveuC0yy1688T6z4Zk8T6gNOsBFMlyuz7U3mLkIZNuMAhu2enSoU8Q6zD8NPFEkl+8moaPez2Ud4VAdwjLhj2zhqAPQIdXsptSGnJN\/pn2cXJhKkMsZOAT6c1clkWKF5W+6iljj2rzbSrCeb4tSzvqt9uOkw3DKGTDAv\/AKs\/L9z9fevRL\/8A5B9z\/wBcm\/lQByNt8S9PvLVLq20LxFNbuMrLHprMjD1DA4rodC1\/TvEdh9t0ycyRBzG6spV43HVWU8g1574C8WX9h4I0q0g8KaxeLHGQs8KL5b\/M3QlqEPiDQLS6uvJi07V\/FOsRxQxEiQWiEH527M2M\/p9KAPV6K8\/1RNW8E3Ol3669fanZXN5HaXdvfbWPz8B4yqjbg\/w1JP4lvNAv\/GMN\/cmRbO3W\/sN4H3GUjYPYSLj8aAO8oryu08Ta9D4BvbC6umPiZNRXTI5GA3F5GVlb8FZv++al1nxItx4vvdDvPEd5pOn6ZFEpe1TM11Ky7iS+1toUY4xyTQB6fRXkp8YainhXxZBbarLeNpkccljqbQ7HkRz0YEAFlPGcc12\/hbTNUtYDfaprdxfzXkSO0LqqxwseSEA7c4\/CgDT1jVrbRNOe8uQ7KGVEjiGXldjhUUd2JNXl5UEgg+hrjtcc3\/xK8OaY+TBawzagyHozgbEP\/AcsaTW7rUda8bQ+GbLUp9OtIbL7bdz22BK+X2qikg7fXOKAN5vEFmut3WkkS\/aba1F3J8vy7CSOD68Vz1v8TdNurRbuHRPEMlqw3CdNOZkI9dwOKydNsL7TfiHr1ve38t9\/xJFME8oAkMe9sBsAAkHcM454qt4H1fxdbeANNi0zwtBdwrC3kzvqKpv+Zv4CvH50Aej6Tq1jrumQ6jp1wJ7WYZVwCOnBBB5BzV+vG7TVLvSvhyv9kXUltrP9tiK+jliC+XO7\/Mm3J+Tp359ug29UtfEeieJtFsLTxRe3A1kyxXDXUcbeUVUOXiUKAvG7A5H1oA9JrO0vV7fVvtQhWSOa1na3nhlADow+hPBBDA9wa5jQpNR0f4g3XhyfVbvUrOTTFv4nvCrSRt5nlldwAyD1qRnOnfF1ETiHVdLJkX+9LE\/Df98tigDtKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigArnPFug3WvLootZIU+w6tBeyeaxGY03ZC4B+bn\/wCvXR0UAc9Y6HdW3jjVtbd4jbXltDFGgJ3gpnORjGOfWuT8W6DdaZ4c8e6jNJC0OpCKSEISWUKAp3cfyzXptV72yttRs5bS8hSe3lG143GVYe9AHCyaF4q8Q6FaaHqdxpyaSyxfaLqBn8+4jXBC7cYUnAycn29K1fEOg6sfEGm6\/oJszc2kL20ltdFlSSJsHhlB2kEV1aIsaKigKqjAA7CnUAcZo3hvW08W32u6zcWT\/bLFbYw2u7ERDfdG4cjHO71J4xVrwNo2s+HtCXSNUaxkhtfltZbZn3OpLE7wwwDyOma6migCG5jM1rNEuMujKM+4rH8G6Nc+HvCWn6TdvE89tGVdoiSpyxPGQD39K3qKAOc8HaDdeHvDzaddyRNKbiaXdAxIw7lh1A55rntC8MeKtC0OTw7atpCWpaQf2nlzMVck5Me3BcZ7tjp1r0SigDgNL8E6pZ2Hg+CaezL6JNK85R2IdWDBdvy9eRnOK3NR0G6vfGujazHJEttZQTxSqWIclwANvGO3rXR0UAee3vh\/xjdaLc+H7l9H1K1lDxx6jeM\/nIjZ5ZNuC4zwQewrUs\/CdxYeJNEvIp45LPTdK+wHeSJHYbcNjGMYX1rrqKAOf03RLmz8Ya7q8jxG3v47ZYlUnepjVg27jH8XHJrCfwVqTeEvFmlCe18\/V9RuLuBtzbURypAb5cg\/KegNd7RQBx\/\/AAj2s2fjOz1mwewe2axisryOZnDgK2S0eBg8djiuquYzNazRLjLoyjPuKmooAwfBujXPh7wlp+k3bxPPbRlXaIkqcsTxkA9\/Sk8WeH5PEGmwLbXItr+zuUu7SZhuVZU6bh3U5IrfooA4ifQ\/EviS+04eIP7NtdPsbhLporKR5GuJU+7ksBtTPOOTT\/GHgy48Ra5pN7bTwxQxkRagjk5mgEiSBRgH+JO+OtdpRQBxM\/gq4m+JUfiH7RF\/ZoRZpLck72uVRo1bGMYCt1z1HSpr7Qtb03xRd694eazmF\/HGl5Z3jMgZkGFdGUHBxxgiuwooA4e+8M+I9Y8La9a6nqVs9\/qYUQwIWFtaquPlBxuOe5xXY2kTQWcELkbkjVTjpkCp6KAON8TRnTfGfhzXzxb7n0+5b+6Jf9WT6DeMf8CFWNc0LVB4itvEegyWxvo7c2lxb3RZY54t24fMoO1ge+K6G9srbUbKW0u4Vmt5V2vG3IIqdRtUDngY5JNAHF6Z4c19vEupa1q9xYl73TxapDbs+ISGOFyRyMc59SeKz9A0b4i+HtDtdJtW8LPBbKVV5WuCx5J5wAO9ejUUAeeH4f3\/APYXktfwXGqXOrx6leTuCiHa2SqAA9B0\/pXR6xotzqHibw9qUTxLDpsk7zKxO5g8ZQbePX1xXQUUAc7\/AGFdf8LEPiHzIfsn9k\/YdmT5nmebvzjGNuPf8Kz7SNtW+KN5frza6TZCyDdmnkO9sfRdoP1rsqq2dlbWELRWsIiR5GlbH8TscsSe5JNAFqiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigD\/2Q=="}
%---
%[output:0221353b]
%   data: {"dataType":"text","outputData":{"text":"Throughput(%) for 10 transport block(s) at transmit power 30 dBm with 1 repetition(s): 0.0000 \nThroughput(%) for 10 transport block(s) at transmit power 35 dBm with 1 repetition(s): 0.0000 \nThroughput(%) for 10 transport block(s) at transmit power 40 dBm with 1 repetition(s): 80.0000 \nThroughput(%) for 10 transport block(s) at transmit power 45 dBm with 1 repetition(s): 100.0000 \nThroughput(%) for 10 transport block(s) at transmit power 30 dBm with 32 repetition(s): 40.0000 \nThroughput(%) for 10 transport block(s) at transmit power 35 dBm with 32 repetition(s): 100.0000 \nThroughput(%) for 10 transport block(s) at transmit power 40 dBm with 32 repetition(s): 100.0000 \nThroughput(%) for 10 transport block(s) at transmit power 45 dBm with 32 repetition(s): 100.0000 \n","truncated":false}}
%---
%[output:332ad7a1]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAjAAAAFRCAYAAABqsZcNAAAAAXNSR0IArs4c6QAAIABJREFUeF7tfQ1wVtd55huIBLYk5BDbCooNBWRENwbbaA3xji3baT3NNLJmrAFtS5yadKsEorJdYxx2MgwsDHWXGtu7QxU70bYDW4dtZVp3ZJl6yW5tCw8ODlIVOt2uvKQqJkaV7UoWkrF+Atl5Lz4fR5f7fffnu7\/vee4MY0v33PPznOe+99H7nvOeT82fP\/8XhAsIAAEgAASAABAAAhlC4FMQMBmaLXQVCAABIAAEgAAQsBCAgAERgAAQAAJAAAgAgcwhAAGTuSlDh4EAEAACQAAIAAEIGHAACAABIAAEgAAQyBwCEDCZmzJ0GAgAASAABIAAEICAAQeAABAAAkAACACBzCEAAZO5KUOHgQAQAAJAAAgAAQgYcCBWBA4cOEC1tbU0ODhIjz32GA0MDFjtr1mzhvbs2UOlpaXU1tZGHR0dgfr15JNP0t133+347Pj4OG3fvp1OnDgx4\/7mzZtp\/fr1M3536NAh2r9\/\/4zfLV68mJ566ilasGBB7vf2caj233jjDXr88cdz5fRnnep26rBTe\/Zy7777LlVWVlJ5ebnjmO394ELNzc3U2tpqYc1XPlwKTYC9Di47NTU1Y+50XO333CZX8YHHZR+DqteOvRNeOtZBuFGon6qPXMaJV25j1HnvVofX94br0efW3gcn7oXBBy9jRRkgEDYCEDBhI4r6CiKgDDEX0o1pUAHDH60nnniCXnzxRUv0FPpIcZv2j57eH3vH+\/v7acOGDdavnT7Yqrz+cU6bgLHjnG8cfkWMk+hzwkHH1E8buoCxP+ckYAqJPcWzsAWM4o6fcdk55rUOr++Nm4CxC8mw+BC22ePx9vT0XPVHRNjtoL5sIwABk+35y1zvdUOsG\/4gAsbpmXwCwqms+hDajbrqo\/r9j3\/845znJZ+oUb8PU8DYJ9ep7kJeAFVe77P9d7pQ8OoZ4n7lG6fuVVDeNB0\/r23o\/eI69TE4CRj1IeY5U94QNY9uc+P3JbKLpSACxm8dXt8bJWD+5V\/+ZYaHk3+v6tA9WmHxwS+G+coH8VSG1TbqyR4CEDDZm7NM91gZ0eHhYZo\/f34uPJBPwNg9JMr42j9w6iP3\/vvvWyEke9jB6UOv6rZ\/VJURffXVV62\/ANXH0emjwB\/T+++\/P\/ex8CNgCokAp0kOKmCcPlj5PrpOuOp94bq++93vXhVK00WgwovDg8qD5UWM6WV0PnDdemixkIDhcvlEkhvebiE7FkK7du3KjZ0x5PZ00eQm7vi+3k6+OvJ5atzem0ICppCg9SPC1HvzyiuvWO+aCl\/acXcLTal6WODedtttFpYfffQRlZWV5YbPHtOmpqZM2zx0PjoEIGCiwxY1OyBgN37qw3TmzJkZa2D0v9rt1fCH5Nlnn7XK62s\/+PdKwOQDXxlZ9RH57Gc\/67rmRn0wdS9AvvrdQlj8nD2k4bROxa+AybcGxh4ys39U8q1FKrSm5vDhw1dhz\/1VH8HVq1dba4pOnTpFS5YssebI7uVyExO64GR8vvzlL+fCfw0NDVb9et\/t4sNpzU2hueE56erqukqY6fOgCxgWt0NDQ9Z6k6ACplAd+QSMEg353pt8AkbHR+ebGx+ceFgo7Kq47RRi5Lr0ebHXw\/c+\/vhja02XuiBg8BkphAAEDPgRKwK614Mb5g+RLkiUYa6qqrLuOf1lz8\/xQl+76PGyBkbV5\/SsmygJW8D4Bb6QByaf4OA27P12+rh4GZvqrz1c89577+U+\/PxxfOedd65aFG3\/eLmN3e4x27Rpk7X4mz+QijdexJfuWXATMPZF2259dApbuT1jv++1Dq\/vjRIwaoG2vT0nT4tfPjgtKNZ\/98wzz9Cjjz5qLXbXeWUP6TnVgxCSXwaZXR4Cxuz5j330uiHW\/+K1\/2V55513OoaC9Di+8gTo4QUvIRz+yDo9G6aACWMXkr0\/fkNI6sOULzzgFJ7L511RfcnnLdLDE7z4ksWn3q794+VGPLuAWbRoUc7bYffIqJ1sqk77B1kJnW9961uOnFLPeQkh6SExr+Kj0Fi91uH1vSkkYNxEar5wrb3\/TqFXfRx\/+qd\/Sl\/72tesx\/TdWfaxKlGqcwoCxu3NwH0dAQgY8CFWBOzGz+7CVh6SsAUMD9K+BsBpUSOXy7cGxh4q4LLFrIHxC7xfAeN1YbQuLpxCc3o\/+WPDF699yLc4+MiRI5bY0NfA+AnDcf1Oa5bsHhS7B8aOp13ArV27VoSA0ddlKU+Lem+UgFHrtfhntfXf6zoXN7EZtoDR185AwPi1CmaXh4Axe\/5jH72T8dM\/TMoQewkhqXUy+joWrx4YztGSbxeSqsNtF5K+4NVtp0sYhtmvgLF\/wPVQj1N+FLe\/0BVZlOjknzmUZ99l9NZbb1lrZPS\/wPOJxXwEdBIwdg+JLmCcFqja\/+J3EzB+Xwav3pOwPTAq1OX03tgFDItIHbd864b88MEp9KP6wvX7DSHZF\/\/mW1zvd35QXj4CEDDy5zhVI3QyTroQcBIN9gGoD639g+ZlEW++LdNOIEWdB8ZtIau9T2GsgQlrHYjTQk4nQaGPwUu+HFU+3\/ZwPTykt1coT4+buAz6guQTMH7m1asI8vreOAkY\/p2Oj9OuNCcM8u3mCnsRbz4Bw33CIt6g7DTjOQgYM+Y5NaPM99eVW04WNQD7Ggz9g8zG7ty5c1RXV+c43nzZYJPKxOvnQ8cDCiJgnLwqXsdbiDSFvCHqOScPgcqw7Db2QvltnDwA3Ga+rfVq3Ypbm35fkqQFDPfX\/t7kEzA6fwplTFYYFMrXo\/D\/q7\/6KwtzlZk66DZq+3M6Pzns9cADD\/idGpQ3BAEIGEMmGsMEAmlDgD+EvO1dP3IhbX0M2h\/pY1M7wvzu3AqKJ54DAk4IiBAw9oWUaqC6q9Ou8vW\/Dr3m4QCFgAAQCAcBfv9uuOGGvInuwmklmVokj40RxRqVZHiFVq9GIPMCJt\/Bbvz7xsZGaxsfJ9ZS\/88H+bHrt6Wlhdrb2y1E1P8HPUAQxAICTgisW3YN3bWgJHfrhbcn6M3BKYAFBDwjkEYOQcB4nr7EC6aRP2GCkmkBwy8Sx+J53cOcOXNmnP3B9\/ji+Ld9Wyz\/hbR06dJcecnu3jDJgrq8I8CG46byWfRM70e5hx5dVUY3V8ymLa+f914RShqLADhk7NSHMnAT+JNpAcNeFo7B2gWJEiw\/\/elPrfi6\/Wdd3CiXqBI7oTAHlRiNwF0LSumLC0pmiBcFCIuYHw1OwxNjNEPcBw8OuWOEEvkRMIU\/mRYwavryCRh1GJ8SKWrBoN3j4hazVqvsVXt8oBouIJAPgafvraQtr49at28o\/TmxMdFDSfwzQkngTyEE3Djidh\/omo2AEz90z+\/T984T4QmGgPlke2q+BYUsXngdzapVq3JvxNGjR4n\/4QICTghsqh6iZ89V0Q0l03TvdWP0wA0T9O5Hv8gVvaF0mt6furI2BigCATsCOkc4225JaQl9NH4lHAkOgTOFELhu3jw6N+vGXJGSiWFq\/ItBujQ6aP2uo+Ez1Nw1knkQRQuYMEJILFw42yhnFuU8I3yxB2ZkJPuTn3n2pnQAT6yZQ985MUkcLmpc9Cn64\/dvphd+\/I80OTlJN5XPpp13\/IJajl1Kae\/RrTQgoDjEfeHTmZcvX059fX3gUBomJ8V9KFm4iuasfJC+M\/8kbb\/0GzR9poeu\/fA0faHkHB0\/fvwyfypmE3tgIGBSMpH2EBJ3Sw8TOS3i1T0uhRbxKgHT1NSUEzApGTa6kVIE2H37aF2ZFTr649Oz6YdTi+nkyZOW8cAamJROWsq6pa9hYAFz6623gkMpm6O0dGdW5QKau7LBEi6zKxfQxdFBuu30Abrj\/f9prcOTzB+RHhgmVljbqCFg0vKaZqcfvPqf\/8L52dhF+vJfT1sfn6H+v6UHF37KGoS+Myk7o0JP40aAOXR27CL9yU8\/DQ7FDX7K22PRUrKojuauaLD+y6KFvS2Tp7po+p0eq\/cm8EesgFFeGM4YyVfQRHYQMCl\/k1PWPXbPvvkb19MLb39MnPflvsUV9IXPf4b4IMU3B6et3+MCAl4RYE8MOOQVLfnlWLhce883LI8Li5ZLH56jib\/rsoSL0yWdPyIETJS0hYCJEl1ZdSvxwjuMVHzZ7r6VNWKMJg4EwKE4UE5vGypEVLKwLudtmTz1Ek2c6sotyi3Ue8n8gYBx4S0ETHpf7DT1jMXLC1\/5jNWlu\/7sg1zXJBuPNOEvuS\/gkOTZdR6blxCRV1Qk8wcCBgLG63uAcnkQUKv6by6fTeteHrHWvqhLsvEAIeJBAByKB+c0tKK8LRwm4ovDRBeOfT9viMhLnyXzBwIGAsbLO4AyAcQLPyLZeIAU8SAADsWDc1KtFBsicuu3ZP5AwEDAuPEf9wsgwNuim5ddY2W1dMquK9l4gBjxIAAOxYNznK2oEBGva1ELcu27iMLqj2T+QMBAwIT1nhhXD4uXLXXl9HTPeN6t0ZKNh3ETntCAwaGEgI+g2ShCRG7dlMwfCBgIGDf+474DAirXSyHxghASqBMGApI\/QGHgk\/Y6nBbk8i6i6TO9uZwtUY5BMn8gYCBgonx3RNbNuRX4LBE38QIBI3L6Yx+U5A9Q7GDG1CCLltmV1TRnZUPkISK3IUnmDwQMBIwb\/3FfQ0BPVKef7poPJMnGA8SIBwFwKB6cw2jFKa0\/e1suHGsPo\/pAdUjmDwQMBEygl8LEh\/LleimEhWTjYSIHkhgzOJQE6t7bDDNni\/dWvZeUzB8IGAgY72+CwSV18WLP9QIBYzAxYhi65A9QDPBF1oQ9rT\/vIuJziPKl9Y+sIy4VS+YPBAwETFLvVWbaLZSozm0Qko2H29hxPxwEwKFwcAyjlnwhIq9p\/cPog986JPMHAgYCxu\/7YFx5t1wv8MAYR4lYByz5AxQrkAEbS3uIyG1YkvkDAQMB48Z\/o+8r8dLx9sd5c71AwBhNkcgHL\/kDFDl4RTSgh4i4Gg4RFTr5uYimIn1UMn8gYCBgIn15sly511wvEDBZnuX0913yByht6Eed1j+J8UrmDwQMBEwS71Tq2wxDvPAgJRuP1E+ikA6CQ9FOZNZDRG7oSOYPBAwEjBv\/jbvvN9cLPDDGUSTWAUv+AMUKpK2xJNL6JzFeyfwRK2AWL15M+\/bto+rqaosz3d3dtG3bthx\/9u7dS\/X19Y73dJKtWrWK2traqKmpiQYHB5PgH9qMEYEguV4gYGKcIAObkvwBins6k07rH\/d4pXuBxQqYAwcOUEVFBW3dutXiDIuZ1157jfbv30\/Nzc3U0tJC7e2XsyOq\/+\/o6LiKXxAwSbxyybQZtniRbjySmSXzWoWAKW7O4zz5ubieRvO0ZP6IFDBr1qyhPXv2UGdnpyVY+Nq8eTPV1dXRhg0biL0vNTU1lrgZGBggFjtDQ0MzPDSKShAw0bxUaau1mFwv8MCkbTZl9UfyByjKmUpjWv8ox5uvbsn8MUrA3HfffZZo2blzpzXXLGb4YgGj\/6wTAQImiVcu3jZZvKy7ZS41L7uG+HyjNwenQuuAZOMRGkioqCAC4JB3gpgYInJDRzJ\/RAoYtf6FJ1YPIc2bN4+2b99OmzZtmuFxYY9MVVVVTtA4CZjW1tbcGpjh4WE3zuB+hhD41hc+TV9bed1l8XIuPPHCEHAYc8WKFdTb20sTExMZQgVdTQsC4JD7TFyqujV38nPJxLCV0t\/K23Kqy\/1h4SWYP7W1tXTy5EmanJwUNVqRAoZnSIWRysvLaWpqin7yk5\/QL\/\/yLwcWMPqsHz16lPgfruwjsPaGYeJ\/h9+fb\/0L+yorK6MlS5ZQf3+\/xUNcQMAvAuCQM2LTc+fTaNVqOv+51cT\/z8Jl3j+\/RdefecUvxKLLM3\/435EjRyBgsjrT+rqXICEkXlOjdiGxB2ZkZCSrUKDfnyCwbtlc2nk70dM94\/QnP\/10JLiw+3b58uXU19cnznhEAhgqvQoBcOgKJFaIaGEdzV3ZQBeuq6GLo4OWp4UPULx29DTY44AA82fp0qV0\/PhxcTZIrAfGPo\/6Ohd7yAiLeM1778PM9VIIPcnxZ\/NYk8yIwSEi+8nPlz48l8m0\/kkwSDJ\/xAoYXZTwDqS1a9da+Vx4qzS2USfxGqWnTbVd+uz4RWruitaTJtl4pGdGZffEVA5JTOufBFMl80esgNHXwDBpDh06lNtSzT8jkV0Sr1LybUaR6wUemOTnVXIPJH+A7PMmPa1\/EjyVzB+xAiYsomAbdVhIJl9PVLleIGCSn1vJPZD8AVLzZkpa\/yR4Kpk\/EDAujIKASeKVC79NFi9bVpXRXQtKad3LI\/SzsYvhN+JQo2TjEQuAaETsgaAIEcVDbsk2CAIGAiaetyjhVh5dVRZJojq3YUk2Hm5jx\/1wEJDEIdPT+ofDCH+1SOKPfeQQMBAw\/t6GDJZm8bKlrtzaLv1M70exjkCy8YgVSIMbk8AhhIiSI7AE\/uRDDwIGAia5NyuGltctu4aevndeIuKFhyfZeMQwfWgiwxxCWv900FeyDYKAgYBJx1sWQS94vUtHw2cSEy8QMBFMqoFVZukDxKJldmV1Lq2\/SjTHqf052Ryu+BHIEn\/8ogMBAwHjlzOZKB9Xojo3MCQbD7ex4344CGSBQzj5OZy5jqKWLPAn6LghYCBggnIntc\/FneulEBCSjUdqCSCsY2nlEHK2ZINoaeVPGOhBwEDAhMGj1NShi5c4t0vnA0Cy8UjNpAvvSNo4ZE\/rz2cRIUSUXhKmjT9hIgUBAwETJp8SrSuJRHVuA5ZsPNzGjvvhIJAGDuULEU2c6qJLo4PhDBS1RIJAGvgTycCICAIGAiYqbsVeb1K5XhBCin2qjWowqQ8QQkQyaJYUf+JADwIGAiYOnkXehhIvHW9\/HHuuFwiYyKfX6Abi\/gDpISIGnkNEE3\/XhV1EGWVh3PyJEyYIGAiYOPkWSVtJ53qBgIlkWlHpJwjE8QFCWn+5dIuDP0mhBwEDAZMU90JpN83ihQco2XiEMoGoxBWBqDiEtP6u0IsoEBV\/0gAOBAwETBp4GKgPacn1Ag9MoOnDQx4RCPsDhLT+HoEXUixs\/qQJFggYCJg08dFzX9KU6wUCxvO0oWAABML4ACGtfwDghTwSBn\/SCgUEDARMWrmZt19ZES8IIWWOWqnscNAPEEJEqZzO2DsVlD+xdzRAgxAwEDABaJPcI2nM9QIPTHJ8MKFlvx8gpPU3gRXex+iXP95rTr6kaAFz4MABqq2ttVDu7++nDRs25BDfu3cv1dfXWz93d3fTtm3bHGdj1apV1NbWRk1NTTQ4iIRNSVKWxcu6W+ZS87JraMvr5+nNwakku+OpbcnGwxMAKFQUAnNWNlDFsn9DlfMqaei996ytzJz11n4hRFQUzKIflmyDxAoYFig1NTW0detWuvHGG2nPnj3U29trCZXm5mZqaWmh9vZ2i7jq\/zs6Oq4iMgRMet7tNCaqc0NHsvFwGzvuF4cAi5fZlQuo5FQH3XrrrXTy5Emavfq3aFZlNY137bIqL1lYh5Ofi4NZ\/NOSbZBYAcPeF76U10X\/WRc3AwMDxPeGhoYcvTAQMOl4v1m8bKkrp6d7xlOVqM4NHcnGw23suB8cARYmJYtW0YVj7Vdtxb\/2nhZLxJQsqrMEzsXRQZo89ZJVFhcQsCMg2QaJFTBOHpjOzk7av3+\/JVjyiRv75CsB09ramgshDQ8P4y2JEYHGRbPo6Xvn0Xf\/\/ueWgMnSVVFRQStWrLC8fxMTE1nqOvqaIAIVD+6ksZcue1mYQ0tuu4v6z5cSLXvAEi4lE8NWKInPIuJMubiAQD4EmD+8lII9eJOTk6KAEitgeJY2b95M69evp6mpKWsdiwoR2T0uLHaqqqpmrJFRs6wEjD7rR48eJf6HK3oEbiiZpv23nKHXP6ygZ89VRd9gyC2UlZXRkiVLrDVYzENcQMALAv9cu54+13\/IKvqpis\/RO7f\/Lk1PTVnCZd7QW3S+ajXd\/JM\/8lIVyhiOANsg\/nfkyBEImKxwgUUKK09eA8PXvn376PTp01aYKIiA4TU0ahEve2BGRkayAkVm+3lT+Wz60weuodMfXKCWY5cyOQ523y5fvpz6+vrEGY9MTkhGOl3ywDaa\/uFesnYUrf2vlnAZ69plnfzMv7t0z38g6nLeeJCRIaKbMSHANmjp0qV0\/PhxcTZIpAdmzZo11qJdFTJS3pjGxkbavn07bdq0yaKO0\/oYO6ewBiamt8zWTJZyvRRCSHL8ORlmmNGqWgPDa114FxJ7W37yxlHrA8RrYKbP9DruRjIDHYzSDwKSbZCRAoa3ROshIyzi9fM6RF82a7leIGCi54SJLZQ37KS5KxtoVvd\/oVt+8U\/0t6fPES37VQsKLNg1kRHBxgwBEwy3RJ9yCiGNjY1ZXhdso050ago2zuJly6oyumtBKa17eYR+NnYxvZ310DPJxsPD8FEkIAIcJqp8+HuXn+7\/Ic1fuNzKA8MLdzkXDC4g4BUByTZIpAeGJ3bx4sXWupfq6urLNgCJ7LzyPdFyWcz1Ag9MopQR2TiHieasfJBGn\/8mVdCFXB4YabtIRE5eygYFAZOyCYmzO1gDEx\/aWc31AgETH0dMaInXv1Q+\/Jy1aJe9LZI\/QCbMZ9JjlMwfsR6YsEgDARMWkoXrWbfsGivXS9YS1bmhI9l4uI0d9\/0joEJHlz48R6M\/2GhVAA75xxFPXEFAMn8gYFyYDgETvSng9S4dDZ8RJ17w8YmeO9Ja0ENHvGUaHJI2w\/GPBwImfsxT0yIETLRTwYt23\/yN6+mFtz+2DmiUdkk2HtLmKunxqNDRhWPfn7HLCBxKemay3b5k\/sADAw9MYm+nlFwvhQCUbDwSI47Ahjl0VNGw0xqZCh2pYYJDAic8xiFJ5g8EDARMjK\/SlaZ08SJhu3Q+ECUbj0SII7RRp9ARBIzQyY55WJJtEAQMBEzMrxORpER1buBJNh5uY8d9bwiw92V+ayfZQ0cQMN7wQ6nCCEi2QRAwEDCxv\/\/Scr0ghBQ7hcQ0qEJHs66rppG2RsdxSf4AiZnIFA9EMn8gYCBgYn31lHjpePtjeqb3o1jbTqIxycYjCTyltalCR+Mv7cp7thE4JG3W4x2PZP5AwEDAxPY2Sc31Ag9MbBQS1ZBb6AghJFHTndhgIGASgz75hrGNOpw5MFG8MHKSjUc4zDCzFi+hIwgYM7kR9qgl2yB4YOCBCft9uao+6ble4IGJnELiGvASOoKAETftiQwIAiYR2NPRKDwwxc2DCbleIGCK44hpT3sNHUHAmMaMaMYLARMNrpmoFQIm+DSZLl4QQgrOHalP+gkdQcBIZUG844KAiRfvVLUGARNsOkzK9QIPTDCOmPiUn9ARBIyJDAl\/zBAw4WOamRohYPxPFYuXdbfMpeZl11jnG705OOW\/EiFPSDYeQqYotmH4DR1BwMQ2NaIbkmyDsIjXhboQMP7fbZMS1bmhI9l4uI0d968gECR0BAEDBoWBgGQbJFLArFmzhvbs2UPl5eUz5n98fJy2b99OJ06coL1791J9fb11v7u7m7Zt2+bIFQgYf68Qi5ctdeX0dM+4EYnq3NCRbDzcxo77VxCYs7KBrr3nG1QoYV0+vMAhMKkYBCTzR6SAsU\/24sWLad++fXT69GlLqDQ3N1NLSwu1t7dbRdX\/d3R0XMUTCBjvr46puV4KISTZeHhnhtklg4aO4IExmzdhjV6yDTJCwLC3paamhrZu3UoDAwOW90X\/+cCBAzQ0NOTohYGA8fYamZzrBQLGG0dMLFVM6AgCxkTGhD9mCJjwMY2tRhVO6uzspP3791vtsmDha8OGDY4\/652DgHGfKrVd+uz4RWruGnF\/wKASko2HQdMYeKjFhI4gYALDjgc1BCTbIPEeGLu3RQkY3ePCZaqqqnKCxknAtLa20uDgoHVreHgYL8gnCNxQ+nN64SufoTlz59Bd\/+MD4GJDoKKiglasWEG9vb00MTEBfAxCYHrufJrf2kn09g9p7KVdgUcODgWGDg8SEfOntraWTp48SZOTk6IwES1g7Gtf1MzZQ0ZeBIw+60ePHiX+Z\/p1Q8k0bap+j24onabd\/\/R5en+6xHRIrhp\/WVkZLVmyhPr7+2lqytzt5CYS4+xtv0ssYpac2F3U8MGhouAz\/mHmD\/87cuQIBEyW2MDhox07dtDBgwdJX6AbJITEu5p0D8zIiNmhkpvKZ9O3vvBpWjX\/En3tf12AeMnzYrD7dvny5dTX1yfOeGTJFsTd1zkrGmj6tn9L0z\/8zzR76O+Lah4cKgo+4x9m\/ixdupSOHz8uzgaJ9sBs3ryZGhsbc1unFZPtHhcs4vX\/jiPXizfMJMefvSFgXim162jiVBeNdwUPHSnkwCHzOBTmiCXzR7SAyRcawjbq4l4P5Hrxjp9k4+EdBbNKVn71OZp1XTWNtDWGMnBwKBQYja1EMn9ECxh7qEhnMBLZBXufkevFH26SjYc\/JMwozbuOKhp20ujzG2n6nZ5QBg0OhQKjsZVI5k9iAubJJ5+ku+++25FUb7zxBj3++OOpIBy2UV+ZhrsWlFJHw2eQZdcHMyUbDx8wGFE07NARQkhG0CbyQUq2QbEKGD3FfyGRosQN79poa2ubsQA38tm2NQABcxkQJKoLxjzJxiMYInKfCjt0BAEjlytxjkyyDYpNwLB44TT+HLrhs4i8XOqZpqYmL8UjKQMBc1m8cK4Xvu76M+R68UM0ycbDDw7Sy0YROoKAkc6aeMYn2QbFJmDimarwWzFdwOjiZd3LI\/SzsYvhgyy4RsnGQ\/C0+RpaVKEjCBhf04DCeRCQbIMgYFxob7KAYfHy9L3z6Oby2QTxEsw+SjYewRCR91RUoSMIGHlcSWJEkm1Q4gJGXxfDk5uGdS86yUwWMMj1Ury5kWw8ikcn+zVEGTqCgMk+P9IwAsk2KHEBw1ude3p6cgctco6WRx55hHbv3u15rUyUJDFVwCjhF0SLAAAgAElEQVTx0vH2x\/RM70dRQiy6bsnGQ\/TEeRhc1KEjCBgPk4AirghItkGxCRinRbx8VtFTTz1Fr776KgSMKw3jK4BcL+FhLdl4hIdSNmuKOnQEAZNNXqSt15JtUGwChidVhYvGxsboscceo4GBgdzvysvLrXlHCClZ+kO8hIu\/ZOMRLlLZqi2O0BEETLY4kdbeSrZBsQoYNcEcJmptbbUEzIYNG9I671a\/TAohIddL+FSUbDzCRysbNXLoqPLh79H0mZ5QzjpyGzU45IYQ7hdCQDJ\/EhEwCmw+bHH9+vWUpsy7diKYImCQ6yUaIyjZeESDWPprLW\/YSSWL6kI768htxOCQG0K4DwGTIAdU5t1Dhw7l1sIk2J0ZTZsgYCBeomMbPj7RYZtEzSp0NNa1iyZPdcXSBXAoFpjFNiKZP7F6YFToqLS01CLL+Pg4bd++PbfbiHck8cLepI8P0JksXcAg10u0dkuy8YgWufTVHnfoSCEADqWPC1nqkWT+xCZgeAHvjh076ODBg7mzjTiEVFdXN2MdDAuYnTt30rPPPott1BG\/JSxe1t0yl5qXXUNbXj9Pbw5ORdyiedVLNh6mzaYKHY0+\/026NDoY2\/DBodigFtmQZP6kTsCkjUGSPTBIVBc92yQbj+jRS08LJQvrqPLh5yjO0BE8MOmZ\/yz3RLINik3AMAHcQkhpJIlUAcPiZUtdOT3dM45EdREST7LxiBC2VFWdVOgIAiZVNMhsZyTboFgFTBYZIFHAINdLfEyUbDziQzHZlpIKHUHAJDvvUlqXbIMgYFxYKk3AINdLvGZJsvGIF8lkWksydAQBk8ycS2tVsg2KTcA4HSXgRhT1TFNTk1tRx\/t79+6l+vp66965c+do69atVvI8vvR73d3dtG3bNsc6JAkYtV367PhFau4aCYQpHvKHgGTj4Q+J7JVWoaNLH56j0R9sTGwA4FBi0ItoWDJ\/YhMwzAT95OlCOV9UXphijhXgHU5r167NbcnmLdp8ceZfXovT0tJC7e3t1u\/U\/3d0dFxFWCkCBrlekrFFko1HMojG12rSoSN4YOKba8ktSbZBsQoYnSQsKGprax15E0ZmXq5\/aGjI0bPC3peampqcR6ZQWQkCBrlekjNPko1HcqhG33IaQkcQMNHPswktSLZBiQmYKImjPD2dnZ2OmX11bwz3w\/6z3jclYPjspsHBy7kfhoeHo+x+qHXfUPpz2rKqjO5bXGGFjc6OXQy1flRWGIGKigpasWIF9fb20sTEBODKAALTc+dbZx2VTAwT53xJ+gKHkp6BbLfP\/GFnwcmTJ2lycjLbg7H1XqyA4aR5b731Fn3pS18izvyrr4Gxe1zYI1NVVeV4sKQSMDpuR48eJf6XhWvtDcN073Xn6dl3q+j\/XLgmC10W1ceysjJasmQJ9ff3Wyet40o\/Ah8s+jKd\/9xquvknf2SJmKQvcCjpGch2+8wf\/nfkyBEImCxMpfLAnD9\/3goT8bVv3z4aGxuzREoQAbNnz54ZHpiRkfQvguVcL\/+u5qKV6+VPfvrpLEyduD6y+3b58uXU19cnzniImywiKlm4iujBP6TpH\/5norf\/VyqGCA6lYhoy2wnmz9KlS+n48ePibJBYDwwLDj2ExIt6GxsbrbOXNm3aZJGRxQxfXkJIvBNKhZCywGTkeknHLEmOP6cD4fB6kZZdR\/YRgUPhzbGJNUnmj0gBw+cpscfltddey62BYQFz3333WR6ZjRs3zggZSVvEe9eCUupo+Ayy7KbAWkk2HimAN9QuXHtPC81Z+aC17iXOs47cBgEOuSGE+4UQkMyfRARMoUW2uqfkxIkTgZmp7zRSIaTTp09bu5Ikb6NGorrAlInkQcnGIxLAEqpU7Tq6cOz7dOHY5fQKabnAobTMRDb7IZk\/YgUMU01PVseLKFXIyH5PSiI75HpJn4GRbDzSh3awHnHoqKJhp\/Vwkgnr8vUeHAo2r3jqMgKS+ROrgFEJ6tyIZRcbbuWjvJ+VPDC6eFn38gj9DNulo6SF57olGw\/PIKS8YFpDRwo2cCjlBEp59yTzJ1YBo+bZLU9LmviQBQGDRHVpYszMvkg2HulF3XvP2Psyv7WT0hg6goDxPo8omR8ByTYoEQGTJbJlQcDwdunmZdfQltfP05uDyDWSJn5JNh5pwjlIX1ToaNZ11TTS1hikilieAYdigVlsI5L5k4iA0c9EcmLN+Pi4td25mEW8YbEx7QJGiZeOtz+mZ3o\/CmvYqCckBCQbj5AgSqwaFToaf2kXTb\/Tk1g\/3BoGh9wQwv1CCEjmTyICJh\/YLGw4g+7BgwfJ6WDFJGiaZgGDXC9JMMJfm5KNhz8k0lU6C6EjhJDSxZms9kayDUqVgGGC8Dbq+++\/nx577DEaGBhInDNpFTAQL4lTw1MHJBsPTwCksFBWQkcQMCkkTwa7JNkGpVLAqIy5CCE5vy3I9ZIdKyLZeGRnFmb2NCuhIwiYrDIsXf2WbINSJ2A4K+68efPggcnzDiDXS7qMg1tvJBsPt7Gn8X6WQkcQMGlkUPb6JNkGJSJgCi3i5RN729rasAbG4T2BeIHxyB4C6elx1kJHEDDp4U6WewIBk+XZK7LvaVkDg1wvRU5kQo9LNh4JQRq42ayFjiBgAk81HtQQkGyDEvHAZIldaRAwLF7W3TIXuV6yRJxP+irZeGRpOrIYOoKAyRLD0ttXyTYoUQHjdLTAoUOHcidIp4ESaRAwSFSXBiYE64Nk4xEMkfifymroCAImfq5IbFGyDUpMwLB4uf3222ckrFNrY\/r6+ujxxx9PBZeSFjAsXrbUldPTPeNIVJcKRvjrhGTj4Q+J5ErPWdlA197zDUp7wrp8CIFDyXFHQsuS+ZOIgCl0FhLngcE26suvDXK9ZN98SDYeWZidLIeO4IHJAsPS30fJNggCxoV\/SXlgkOsl\/YbBSw8lGw8v40+yTNZDRxAwSbJHTtuSbVAiAoapwZ6WtWvXztgyjRDS5ZdGbZc+O36RmrtG5LxJBo5EsvFI+3RmPXQEAZN2hmWjf5JtUCICxu0wR50WfLDjAw88kBhT4vbAINdLYlMdScOSjUckgIVUqQodTZzqovGuXSHVmkw14FAyuEtpVTJ\/EhEwcRDDSSR1d3fTtm3brOb37t1L9fX11v\/rv7f3LU4Bg1wvcTAj3jYkG494kfTXWuVXn6NZ11XTSFujvwdTWBocSuGkZKhLkvkjVsA0NzfTI488Qrt37yb7mUp8r6Wlhdrb2y0aqv93OgE7LgHD4mXLqjK6a0EprXt5hH42djFDrwi6mg8BycYjrbMuJXSk8AWH0sq0bPRLMn\/EChheY3PffffR1q1brzrVmr0vNTU1uXt8\/tLQ0FDOO6PTMi4Bg1wv2TAGfnsp2Xj4xSKO8pJCRxAwcTBGfhuSbVBiAoa9IK2trVRaWnoVg3jdy\/bt26\/ynPihGouUqqoq2rBhw1WPsWDhS92z\/xy3gEGuFz8zm62yko1HGmdCUugIAiaNDMtenyTboEQEzOLFi+mpp56i8+fPOwqMMCjCoqS2tjZX1blz5\/J6XAqJHeWBYbE1ODho1Tc8PBxGF606GhfNoqfvnUff\/fufW8nqcMlCoKKiglasWEG9vb00MTEha3BpG03tA1TRsJMmDv8eTZ\/pSVvvAvcHHAoMHR4kIuYPfwtPnjxJk5OTojBJRMAUSmQXBroskPbt20djY2OWQLL\/bA8ZeREwer+OHj1K\/K\/Y619d+zHt+KV36fD7861\/uOQhUFZWRkuWLKH+\/n7ik9ZxRYPA9Nz59I9rdlDlP79Fn+s\/FE0jCdUKDiUEvJBmmT\/878iRIxAwYcyp8sC8+uqrsZ17pGf43bRpkzUMPyGkPXv2zPDAjIwUl5\/lpvLZ9Mqvl9ALb39Mu\/o+FQasqCOFCLD7dvny5cTHY0j76ydNcM9Z+1+JRcyl57+Wpm6F0hdwKBQYja2E+bN06VI6fvy4OBuUiAeGmRT3kQH6ot6NGzfOWB8T9yJe5Hoxx5ZIjj+nZRZ51xGHjkaf30jT78gJHSl8waG0MC2b\/ZDMn9gEjN\/kdcUs4rWHqNTPvA6B88AkuY1aFy\/YLp1Ng+Cn15KNhx8coiorcdeRHStwKCr2mFGvZP7EJmDipopdMPEaBH1HUhKJ7JCoLm4WJN+eZOORPLpEEncdQcCkgVly+iDZBokVMGHRL8w8MMj1EtasZKceycYj6VmQHjpCCClphsloX7INSkTAeA0n2b0mSdApLAGjxEvH2x\/TM70fJTEUtJkAApKNRwJw5po0IXQEAZMkw+S0LdkGJSJgmBpPPvkk3X777TMS1umnUX\/3u9+NPFeMF4qGIWDWLbvGyvXCeV4gXrygLqeMZOOR5CyZEDqCgEmSYXLalmyDEhEwhfLA6LuTVq9eTY2NjZk+jRriRY4hCDISycYjCB5hPGNK6AgCJgy2oA7JNggCxoXfxXhgeNHum79xvZXrZcvr5\/EmGYiAZOORxHRy6Kjy4e9ZmXbHu3Yl0YXY2wSHYodcVIOS+ZOIgPESQnr88cetMBMn4GlqakqMUEEFDHK9JDZlqWpYsvFIAujyhp1UsqiORtoak2g+kTbBoURgF9OoZP4kJmCYHRwuWr9+\/QyiHDp0yMrOy+KFQ0htbW3U0dGRGJmCCBiIl8SmK3UNSzYecYOtQkdjXbto8lRX3M0n1h44lBj0IhqWzJ9EBUwW2OFXwCDXSxZmNb4+SjYe8aFIZGLoSOELDsXJNHltSeYPBIwLX\/0IGBYv626ZS83LrrHWvLw5iMP75JkDfyOSbDz8IVFcaRU6Gn3+m3Rp9PKp8KZc4JApMx3NOCXzJxEB45YHZnx8fMb26mim1VutfgQMEtV5w9SkUpKNR1zzWLKwjioffo5MCx3BAxMXw2S3I9kGJSJg8tGFhc2OHTvo4MGDia570fvnVcCweNlSV45cL7Jtge\/RSTYevsEI8IDJoSMImACEwSNXISDZBqVKwDDyvLD3\/vvvp8cee4wGBgYSp6MXAYNcL4lPU2o7INl4xAG6yaEjCJg4GCa\/Dck2KJUChpPXFXMadZiUdBMwyPUSJtry6pJsPKKeLdNDRxAwUTPMjPol26DUCZgDBw7QvHnzMuGBUdulz45fpOauETPeBozSFwKSjYcvIHwWVqGjSx+eo9EfbPT5tKzi4JCs+Yx7NJL5k4iAKbSId2pqKvHcLzrB8nlgkOsl7tcwm+1JNh5RzghCR1fQBYeiZJr8uiXzJxEBkyXKOAkY5HrJ0gwm21fJxiMqZBE6moksOBQV08yoVzJ\/IGBcOGwXMCxetqwqo7sWlNK6l0foZ2MXzXgLMMpACEg2HoEAcXkIoaOrAQKHomCaOXVK5k+iAoaPC7j77rtnMOmNN94gPgcpzEs\/4frEiRNW1Xv37qX6+nrr\/7u7u2nbtm2OTdoFDHK9hDkz8uuSbDyimL1r72mhOSsfJBMT1uXDExyKgmnm1CmZP4kJGBYvt99++4zdRmptTF9fX2giRtXJdFU7m5qbm6mlpYXa29stFqv\/t5+5xNujm1YtpDtWraIjL79slbW2TPeM0zO9H5nzBmCkgRGQbDwCg5LnQYSOnIEBh8Jmmln1SeZPIgJGiYrOzk7r4Eb9cvKWFEM39rTceuutVFpamhMw\/LuamhraunWrlWuGdz4NDQ3N8MKwULmpfBa9TrXWouL\/8R\/X0+PLLlgho7v+7INiuoRnDUJAsvEIcxoROsqPJjgUJtPMq0syf0QLGPa0PPTQQ8QenV\/91V\/NCRgWLHxt2LDB+q\/9Z17f8sUFJZaXhUNI\/237Jrrhr7dbnhe+fjQ4jXOOzLMDgUYs2XgEAiTPQwgdQcCEySfUdQUByTYoEQHD0MYRQmJh0tPTY82knhzP7nFhj0xVVVVO0Dx97zzrMEa+7rzlJvrL+6bpjbHrcr\/bvuIibXl9FO8IEHBFoKKiglasWEG9vb00MTHhWt7EApeqbrXOOrr41n+nC8e+byIEBccMDoESxSDA\/KmtraWTJ0\/S5ORkMVWl7tnEBIwSMVEt4uVQVF1dnSVK7GEpPwKmsf5O+sOGGhq559\/nJm\/2D75Fj\/7ttambTHQofQiUlZXRkiVLqL+\/nzjHEa6ZCEzPnU\/\/XLve+uXNP\/kjwOOAADgEWhSDAPOH\/x05cgQCphgg43p28eLF9MQTT9CLL75oHQrpJGC4L\/lCSLoHRu1C2rNnDw0ODtJN5bPp4XmnqeXYpbiGg3YyjAC7b5cvX26FMaX99RPGtHDoaGJRPU0e\/j0qmRgOo0pxdYBD4qY01gExf5YuXUrHjx8XZ4MS8cCwwHjqqafo1VdfvWoRbxgzy2tfWltbrYW7+qWy\/LJnRg8Z2T0y9jUwvIi3qanJEjC8jRprYMKYJTPqkBx\/LnYGeeHu\/NZOK2x04djlHYG4rkYAHAIrikFAMn8SETCFdiEVM1H5nrV7YLxso2YvzNmxi7ldSP\/+t9bR3RWXzzvCFuooZklmnZKNRzEzxuKlomEnzbqumkbaGoupSvyz4JD4KY50gJL5k4iA4dniRbzs1nrsscesrcxRXkET2bEn5qFVC6nx3jutPDBvDk7TC29\/HGVXUbcwBCQbj2KmSu06Gn9pF02\/c3mhPS5nBMAhMKMYBCTzJxEBU+gwR56o8fHxGQnuipm8Yp\/Nd5hjsfXieTMQkGw8gs4gQkf+kAOH\/OGF0jMRkMyfRARMlggGAZOl2UpfXyUbjyBoI3TkHzVwyD9meOIKApL5AwHjwnQIGJiCYhCQbDyC4ILQkX\/UwCH\/mOEJCJhIOZBvpxBCSJHCjspjRgAfnyuAI3QUjHzgUDDc8NRlBCTzJxEPjNpGff78+VwulrSSDR6YtM5MNvol2Xj4mQGEjvygNbMsOBQcOzwJARM6B+LeRl3MACBgikEPz+Ljc5kDCB0FfxfAoeDY4UkImNA5EHUiuzA7DAETJprm1YWPDxFCR8XxHhwqDj\/Tn5bMn0RCSEwop9wsaSQaBEwaZyU7fZJsPLzMAkJHXlAqXMZ0DhWPoNk1SOZPbALGLfeLTjHkgTH7hZM0esnGw8s8zVnZQNfe8w1CwjovaDmXMZ1DwZHDk4yAZP7EJmCySiV4YLI6c+not2Tj4YYwQkduCHm7bzKHvCGEUoUQkMyf2ARMlhbu6mSAgIFxKAYBycajEC4IHRXDmpnPmsqh8BA0uybJ\/IGAceE2BIzZL3+xo5dsPAphg9BRscy58rypHAoPQbNrkswfCBgIGLPf7ohHL9l45INOhY4mTnXReNeuiBGWX72JHJI\/q\/GNUDJ\/IGAgYOJ7kwxsSbLxyDedlV99jmZdV00jbY0Gznj4QzaRQ+GjaG6NkvkTu4ApLy93ZRJ2IblChAIZQUCy8XCaAoSOwiemaRwKH0Gza5TMn9gFTGdnJ+3fvz8zjMIamMxMVSo7Ktl42AFH6CgaCprEoWgQNLtWyfyBgEEIyey3O+LRSzYedugQOoqGTCZxKBoEza5VMn8gYCBgzH67Ix69ZOOhQ8eho4qGnTT6\/EaafqcnYlTNqt4UDpk1q\/GNVjJ\/xAoYe+bf7u5u2rZtW441e\/fupfr6eutn+z2dWgghxfeiSWxJsvFQ84XQUbTMNYFD0SJodu2S+RObgImTQnxY5L59++j06dOWaGlubqbW1lY6fPiwtf6Gf25paaH29narW+r\/Ozo6ruomBEycMyevLcnGQ80WQkfR8tYEDkWLoNm1S+aPSAFjp6td0LD3paamhrZu3UoDAwN04MABGhoamuGhUXVAwJj98hc7esnGg7FB6KhYhrg\/L51D7gigRDEISOaPEQJG97iwl4UFC18bNmyw\/mv\/2SmExB6cwcFB69bw8HAxfMKzBiFQUVFBK1asoN7eXpqYmBA18um582l+ayfR2z+ksZeQsC6qyZXMoagwQ71XEGD+1NbW0smTJ2lyclIUNKIFjPK8VFdXU39\/\/wzBontc2CNTVVWVu+8kYPTfHT16lPgfLiDghkBZWRktWbLE4t\/U1JRb8UzdP3vb7xKLmCUndmeq31nrrGQOZW0usthf5g\/\/O3LkCARMFidQCZmxsTFLpNhDRl4EzJ49e2Z4YEZGRrIIBfocMwLsvl2+fDn19fWJMh5zVjTQpfr\/QBOHf49mD\/19zKia1ZxUDpk1i8mNlvmzdOlSOn78uCgbxIiK9sDolNm8eTM1NjbS9u3badOmTdYtPyGkpqamnIBJjopoOWsISIw\/Y9dRvCyUyKF4ETS7Ncn8MUbA6At3N27cOCNkhEW8Zr\/gUY5eovHArqMoGXN13RI5FC+CZrcmmT8iBYx915HKCcMLKdW2amyjNvuljmv00owHdh3FxZwr7UjjUPwImt2iZP6IFDBMV3siO30RL99HIjuzX+q4Ri\/JeHDoqPLh79H0mR4a78KuI3AoLgTQTjEISLJBdhzECphiJlx\/FnlgwkLSzHokGY\/yhp1UsqiORtoazZzMhEYtiUMJQWh0s5L5AwHjQm0IGKPf\/aIHL8V4qNDRWNcumjzVVTQuqMA7AlI45H3EKBkmApL5AwEDARPmu4K6bAhIMB4IHSVLawkcShZBs1uXzB8IGAgYs9\/uiEcvwXio0NHo89+kS6OXs1Hjig8BCRyKDy20ZEdAMn8gYCBg8MZHiEDWjUfJwjqqfPg5QugoQpK4VJ11DiWHHFpmBCTzBwIGAgZveYQIZNl4IHQUITF8VJ1lDvkYJopGhIBk\/kDAQMBE9Nqg2qz\/9YPQUTo4LPkDlA6EZfdCMn8gYCBgZL+9CY8uq8YDoaOEiaM1n1UOpQdBs3simT8QMBAwZr\/dEY8+i8ZDhY4ufXiORn+wMWKEUL0bAlnkkNuYcD8+BCTzBwIGAia+N8nAlrJoPBA6ShdRs8ihdCFodm8k8wcCBgLG7Lc74tFnzXggdBQxIQJUnzUOBRgiHokQAcn8gYCBgInw1UHVWTIeCB2lk69Z4lA6ETS7V5L5AwEDAWP22x3x6LNkPK69p4XmrHyQkLAuYlL4rD5LHPI5NBSPAQHJ\/IGAgYCJ4RUyt4msGA8VOrpw7Pt04Vi7uROWwpFnhUMphA5dQiI7szmAwxzNnv9iR5+Fjw9CR8XOcrTPZ4FD0SKA2otBQDJ\/4IGBB6aYdwPPuiCQBeOB0FG6aZwFDqUbQbN7J5k\/EDAQMGa\/3RGPPu3Gg70v81s7CaGjiIlQRPVp51ARQ8OjMSAgmT9iBUxzczO1trZSaWmpRZH+\/n7asGFDji579+6l+vp66+fu7m7atm2bI5UQQorhDRPcRJqNB4uXioadFvpIWJdeEqaZQ+lFDT1TCEjmj0gBs2bNGtqzZw91dnbS\/v37Sf3c29trCRUWNy0tLdTefnmxovr\/jo6Oq1gPAQNDUAwCaTYeCB0VM7PxPZtmDsWHAloKioBk\/ogUME4TfeDAAevX7IVh70tNTQ1t3bqVBgYGiO8NDQ05emEgYIK+NniOEUir8UDoKDv8TCuHsoOg2T2VzB8jBYwuZpja9p91ukPAmP3yFzv6NBoPFTqadV01jbQ1FjtEPB8xAmnkUMRDRvUhIiCZP0YIGLUe5vDhw1ZIye5xYY9MVVXVjDUyij9KwPB6msHBQevXw8PDIdILVUlGoKKiglasWEEcvpyYmEjFUGev\/i2qvHMtjXXtoukzPanoEzqRH4E0cgjzlR0EmD+1tbV08uRJmpyczE7HPfRUvIBR61\/efffdnEAJImB0LI8ePUr8DxcQcEOgrKyMlixZYi0in5qacise+f3pufPpH9fsoM\/+0yt0\/ZlXIm8PDRSPQNo4VPyIUEOcCDB\/+N+RI0cgYOIEvti2nMQL1xkkhMSLgnUPzMjISLHdw\/MGIMDu2+XLl1NfX1\/ixoNDRyUP\/EdiEXPp+a8ZgL6MIaaJQzIQNWsUzJ+lS5fS8ePHE7dBYSMv1gNj33mkA2cPGWERb9i0Qn0KgTTFn9Wuo\/GXdtH0OwgdZYWlaeJQVjBDP68gIJk\/IgXM4sWLad++fTQ2Nua4rgXbqPF6x4VAWowHdh3FNePht5MWDoU\/MtQYBwKS+SNSwNiT2CmSnDt3Lrd1Gons4nh10EYajAd2HWWbh2ngULYRNLv3kvkjUsCESVdsow4TTfPqSoPxQOgo27xLA4eyjaDZvZfMHwgYF25DwJj98hc7+qSNB0JHxc5g8s8nzaHkEUAPikFAMn8gYCBgink38KwLAkkaD4SOZNAzSQ7JQNDsUUjmDwQMBIzZb3fEo0\/SeMxZ2UDX3vMNwq6jiCc54uqT5FDEQ0P1MSAgmT8QMBAwMbxC5jaRlPFQoaOJU1003rXL3AkQMPKkOCQAOgwhxeexhTE5EDAQMGHwCHXkQSCJjw9CR7LomASHZCFo9mgk8wcCBgLG7Lc74tEnYTwQOop4UmOuPi4OLViwIOaRobliEFCZ4d3qiIs\/bv2I4j4EDARMFLxCnZ8gELfxQOhIHvXi4BCLl+3btxPvusSVDQT4gFj9iJt8vY6DP0khBgEDAZMU94xoN27jUfnV52jWddU00tZoBL4mDDIODql0EV4+iCZgnvYx3nHHHfQ7v\/M71NTUlDujDwIm7bOWQP+QByYB0AU1GcfHR8GF0JEg4mhDiYNDsHPZ4o6f+YqDP0mhBw8MPDBJcc+IduMyHggdyaVTHBzy80GUi3R2RuZnvuLgT1LIQcBAwCTFPSPajct4IHQkl05xcMjPB9EP0uwVLFlYl3tk8lRXJCehr1mzxloP8u677844wJfPvKuqqrJ+l++MvP7+fsdDf\/2Ms5iy3K9HHnmEdu\/eTSdOnPBUlZ\/5ioM\/njodQSEIGAiYCGiFKhUCcRgP\/khUNOyk0ec3RvJxwGwmi0AcHPLzQfSKBvNyduUCunCsPfcIn8s1q7I69NxESsCUlpZSW1sbdXR0WG3aBUxLSwu1t7fn7qvnOjs7af\/+\/V6HFlo5JaqmpqasRdQQMFT0Xo0AABVSSURBVP6ghYCBgPHHGJT2hUDUHx+EjnxNRyYLR80hBiVsAcNel5JFq2aIFwU+i5jpM72him0WIjt27KDz588Ti5itW7fSwMCAq4DhPh04cICGhoZo27ZttHjxYtq3bx9VV1db3T106JAlbJSX5OzZs3TbbbcRCw5dKAUh1ubNm2nt2rX0N3\/zN7R69Wp4YAKACAEDAROANnjEKwJRf3wQOvI6E9ktFzWHnAQMC+NiLvYIjhXIAO12363tS6ODM4ooAfMXf\/EX9JWvfIVee+01S3j49cCwmOGLQ05KYLBQ4au1tdUSRXyP662pqckJJdUZfmb9+vVXdV8JIadxIYTkNtv570PAQMAEZw+edEUgyo8PQkeu8IsoECWHFEB2Dwx7SfgcrbReF459f4Z3RwmYgwcPWmteGhsbrZAMbzN2WwPT3d1teV\/s4STljWExxB4aPfykt6fCVUGxgoAJihwRBAwETHD24ElXBKL6+CB05Aq9mAJRcUgHyC5gmF+zKy+HUYJcLID0tS96HbOuW0BzVzTkve+lvYuj50j3wtgFhQoLcV26gNFFCHtLlNDhtSdKwJSXl8\/oAgucnp6eGQttIWC8zFL0ZSBgIGCiZ5nBLUT18UHoyBxSRcWhQgKmWHSTWgPDHhj2iLBXg8XKBx98QJOTk7ldSPZFvCx0KioqrFDQjTfeaK2jUXXoGKj61AJg+8+qLEJIxTLH3\/PiBYzuBtRXmXMMs76+3kJLuRCdoAt7cZu\/6UHprCMQxccHoaOss8Jf\/6PgkL0HUdi58oaddGn0XM7Twl6duSsbrKbzeWf8IXOltJNHRNl4tU3aSXTYw0b6Ghi1Q+jw4cNWCInXwPzoRz+ywk351sAE6T9CSEFQu\/yMaAGjryjXF1HpRGYQ7Ko8yr9Mgk8VnswiAmF\/fBA6yiILiutz2ByK8w815YnhrdN8Tb\/TQ5wLJuzLScAo+z82NpbXA8P9YDHCAo7XzLz33nszdiGpP26VyJiYmLB2KI2Pj\/va9lxovBAwwdkgVsAo9cwuxHnz5pG+z9+unvVtdHH8ZRJ8uvBk1hAI++OD0FHWGFB8f8PmUJwCpvjRp6OGICIjyp778ZjFwZ8ox1qobrEC5td+7descX\/44YdWhkZdwOhuQi5j\/9nJA8PuQ3V8+fDwcFLzhXYzhgDH11esWEF8ciz\/9VbUVfuAlbBu4vDv0fSZnqKqwsPZQSBUDuUZ9he+8AUrr4mXwwGzg1x4PU2rgPnN3\/xN18McmT+1tbV08uRJaz2QpEusgFGT5JRp0e5x0XMF5PPA6L8\/evQo8T9cQMANgbKyMlqyZAlxHJ6TXwW9pufOp7O3\/S5d++Fp+lz\/oaDV4LkMIhAWhwoNfenSpbRp0yYImIzwQ3lgfv\/3f59GRkYK9pr5w\/+OHDkCAZOR+c11MywBox8zzx4YN9JkDSf0NxoE2H27fPly6uvrK8p4lDywjS5W3UqXnv9aNB1FralFICwOFRrgHXfcQX\/wB38AAZNaFszsmBIwX\/\/61109MMwfFqjHjx8vygalERpjPTA8GZxRkS8vISS4VtNI3\/T3KYz4s9p1xJlNo1gAmX4Uze5hGBxyQ9DPmgq3unA\/egT8zFcc\/Il+xM4tGClg7CEjLOJNin7y2y3WePCuo8qHv2eteRkvkJpdPpLmjrBYDnlBzs8H0Ut9cZfJ4mnU+unY9l1N\/E3idSt8Oe148jNfcfAn7vlW7RkpYLCNOim6mdduscaDc2mULKqj0ee\/OSPzqHlImjviYjnkBTk\/H0Qv9aky65ZdQ3ctKMk98sLbE\/TmYPC1YPnaztpp1LzF+4knnqAXX3zRSryn74xtaGigurq6ghECP\/MVB3\/8cCLMskYKGAYQiezCpBHqyodAMcaDc2hUPvycdSgeQkfmcqwYDnlFzc8H0WudLF5uKp9Fz\/R+lHvk0VVldHPFbNry+nmv1Xgql8XTqPWBFdrlZD\/ygJ\/zM19x8MfTJEVQSLyAKRYzP0Qpti08Lw+BoMYDoSN5XAg6oqAc8tNe2HburgWl9MUFJTPEi+oPi5gfDU6H6onJ8mnU6g9qp9Ot893zM19x8McP18IsCwHjgqYfooQ5MahLBgJBjQdCRzLmP4xRBOWQn7btdu6mitl+Hr+q7NP3zivoZXG779b4z8YuziiS1dOo9QMk9WzxanD5zlzy812Kgz9u8xXVfQgYCJiouIV6iSiI8UDoCNTREQjCIb8I2j+I7CXZUjfzVGa\/dUZZ\/ume8RnenSyfRs04OaX70M9i0s\/x4\/IQMJfZBQEDAROlnTG+br8fHxU6uvThORr9wUbj8QMAwUSwX9ycPDA3lwf3wjxaV0bP9FxZ+6L3h70765bNzXvfS9\/Pjl8k3QtjFzBZOo2ax6vObTp9+rR1WCSve1m7dq2VHZkX+dovCBgIGC\/viS+l66lCFDIKAb8CBqEjo+jhabB+OeSpUlshPx9EL\/UntQbm4MGDuQ9+mk+jzie42tvbLXgLHTAMD8wVBsIDAw+MF3uEMgER8PPxQegoIMjCH\/PDoaBQhC1guB+8zuXs2MVcqMfyvNwy1+qivjMpaJ\/157J4GrWeB4bHotbA6Dtk1RjtuWD8zFcc\/AljDoPUAQEDAROEN3jGIwJejQdCRx4BNbCYVw4VA42fD6KfdpQnhrdO8\/Xm4DS98PbHfqpIRdm0HuboJUN8HPxJapIgYCBgkuKeEe16NR7X3tNCc1Y+iIR1RrDC3yC9cshfrTNLRyVgiulTmp6FgEnTbFzpCwQMBEw6mSmkV14+Pip0dOHY9+nCscsxcFxAQCHghUPFogUBUyyC8T7vZ77i4E+8o4eA8Yy3H6J4rhQFjUHAzXggdGQMFQIP1I1DgSvWHlR2rrW11fV04zDaQx3FIbBgwQJrhxJCSPPn\/6I4KGU\/DQEje36jHp3bxweho6hnIPv1u3EojBHyB3H79u3Wrktc2UCgt7eXWHC6XXHwx60PUd1HCMkFWQiYqKhnRr2FjAd7X+a3dhJCR2ZwIego4\/oAsYjhf7iygcDg4KAnb1lc\/EkCNQgYCJgkeGdMm\/mMB4uXioadFg5IWGcMHQINVPIHKBAgeMgXApL5AwEDAePrZUBhfwjkMx4IHfnD0eTSkj9AJs9rXGOXzB8IGAiYuN4jI9txMh4IHRlJhcCDlvwBCgwKHvSMgGT+QMBAwHh+EVDQPwJ246FCR7Ouq6aRtkb\/FeIJ4xCQ\/AEybjITGLBk\/kDAQMAk8EqZ06TdeKjQ0fhLu2j6nR5zgMBIAyMg+QMUGBQ86BkByfwxVsDo5010d3dbJ4DarzkrG2jhmgZateoOevnlIzR5qgsfHc+vDQoq\/lRWzqP\/+3\/7Le7wwl3sOgI3\/CDwS7\/0S\/T1r3+dnnvuOU+7TvzUjbLyEZDMHyMFjDpqvdDJn\/zxmV25gJZ\/1JNLGDRa00CzKqtpvGuXfNZjhEUhYOfP2q+30oVf+U8Wpz544s6i6sbDZiGAVA5mzXfYo5XMHyMFDHtfampqaOvWrTQwMEAHDhygoaGhnBeGU7uXLFplpXW3Tz6HAKbP9MITE\/ZbJqg+J\/48\/Mxf0\/AN\/5oufXjO4hXCR4ImPOKhSP4ARQwdqie66hsmCRQjBQwLFr42bNhg\/df+c3nDzpyXxSnFtn5fEhkwlnAQ0PnBicF2PtlG\/+5\/f5oqT79kiRfwJxycTalFpY1Hmn9TZjzccfo5diDclqOvzVgBo3tc2CNTVVWVEzT2D1D9jr+aMRP\/++ws+pWbL0U\/O2ghkwj83b98ilZ89soJHfwzX3\/8Kxet\/37n+Gx64t9c\/n9cQAAIAIGoEfB67EDU\/Qi7fggYIiokYBjwyq8+l8Odt7\/yxaEAXEDACQHmiM6P66fP0cSpLqq6loi3Uf+\/6gdp6T98H+ABASAABGJBwOuxA7F0JsRGjBUwjGG+EJK+hsGONdbAhMg+oVWBP0InFsMCAkAgVQgYKWDsHhf7Il6eIQ4jXRq9vOCSL\/7Lee7KBuv\/1e9SNZPoTKoQAH9SNR3oDBAAAgIRMFLAeNlGzXOt\/pLmrdN88c4RzgWDCwh4QQD88YISygABIAAEgiFgpIBhqLwksgsGKZ4CAkAACAABIAAEokbAWAETNbCoHwgAASAABIAAEIgOAQiY6LBFzUAACAABIAAEgEBECEDAfALs4sWLad++fVRdfXm9S39\/f26XEv+8Zs0a2rNnD5WXl9P4+Dht376dTpw4EdG0oNosIuDGoc2bN9P69etzQ5uamrKOqejo6MjicNHnCBHQ1+kpfsAGRQi4wKqdOCTNBkHAfEJcPRuvMhSc\/Ecd8qjft2fuFch9DCkAAm4csu9+C9AEHjEAASWEr7\/++hkCFzbIgMkPaYj5OCTNBkHA5CGM08eos7OT9u\/fT6xsH3nkEdq9eze8MCG9cBKrcTqyQs8ALXHMGFPxCPBfyY2NjVZFfOAse2DUH1WwQcXja0INThzicTulDMkyHhAwDrPnZCxaWlpyxsTJNZdlEqDv4SNg55D6i+i1116zRDAuIOCEAPPm29\/+NjFPWMQoAWO3ObBB4E8+BPJxSKINgoCxsYAVam1tLZ07dy53WrXd48IE2bFjBx08eBDrF2BHrkLAiUP6+gX1QHd3dy5ECRiBACPALn6+enp6yP5Hk+71hQ0CX\/IhkI9DEm0QBEweFvBHqKKiwhIxd95554yQEYwHjIcXBOwc4tOEDx8+nAtD6j97qQ9lZCPAfyg99NBD9J3vfMeyORAwsuc7itG5cUiaDYKAycOiQtl64b6N4tWTV6cbT7AYXN6cFzMi5gN7XtQ6u0JhazduFdMPPJtdBApxyGlUWbdBEDB5uKoWQfF2ab70kBEW8Wb3BY+z5zqHnLbcS1tQFye20tpycu+rMR46dIjeeust2CBpkx7yeNw45LT2Lus2CALmExLpSlQtdhobG3M8sTrrqjXk9wbVeeCQXczwz2vXrkUeGLDHEQEnDwu2UYMsfhCwc0iiDYKA+YQR9iRk+iJeLoIkUn5eHTPLunHInkSK\/7LGjiQzueI2aiSyc0MI990Q8JLILus2CALGjQW4DwSAABAAAkAACKQOAQiY1E0JOgQEgAAQAAJAAAi4IQAB44YQ7gMBIAAEgAAQAAKpQwACJnVTgg4BASAABIAAEAACbghAwLghhPtAAAgAASAABIBA6hCAgEndlKBDQAAIAAEgAASAgBsCEDBuCOE+EAACQAAIAAEgkDoEIGBSNyXoUFoQSCpd+2\/\/9m\/Tq6++SgMDAzOgsOeRseM0NTWVicR4fNhcTU1N7rBUzrG0aNEix4NRC2UXHR8fJ86U7ZTlOC0c4jm77777cmO190tPTpdvfu05qfyOjfGuqqrKJeX0+zzKA4G0IgABk9aZQb8SRyAJAWP\/uOcDQSXNO336dKZPtFYCpbOz0zGpX6H7+mGZdrGXOHk+SX7pdmq9XcA0NjbOEGVqnnk8fLBskHGqOl577TUkTkwDMdCH0BCAgAkNSlQkDQEImOhntBgBk8T8+EHEi+fDTcBwe25nannpk5snyEsdKAME0oYABEzaZgT9SQ0C9g+k+vmNN96gL33pS1RaWmr1VaXjVvf\/4R\/+ge68807rnh7WcfpY67+76aabqL6+\/qrnnADJ54FRB42ePXuWbrvtNlLhh507d1JtbW2uKr1fbuNSH9H169fnnu\/v759xTticOXMsPKqrq60y3d3d1snKra2t1u\/19pSX6eDBg9bHuby83HpGr1M1VEjgOM2Pak+vz6kOuyiw42kP5\/B4tm3bZnWL+79w4cLcePV7bv1mwaLmgeeGcZmcnLSwzCdU7L\/nOvzgzX2CFyY1ZgUdCREBCJgQwURVshDI94H84IMPcu58\/pitWrXKcvvzOg7+gPJHSa3N0O8zOnv27CE9XGL\/uBYbQuI+cx841MAfRfXB1dec2MMSLLb4mULjamlpofb2dmudir3P6qOshJz6+OtrN\/Rwz8aNG3NrYG688carMNFZlE\/A2MfQ0NBALLBUH9Rz7777roWD\/QBW\/pnraGtrs8akz3VdXV1uTnl9jaqrt7fXEjE8Ryw0C50j4yRGdC5wvaoeJdycnnESHn7xVmEnHEIryz5hNEQQMGABEMiDQD4Bc\/jw4dxaAr0MV8NCQL+v\/2X\/l3\/5l7EJGL0PTsPThZISMPnGxQtAC52cbV+L4iQ69PaCCBjlpdHHYhdIfE+JNuU1UmtKVq9eTer\/+d63v\/1tev\/992l0dDQnSljkPfPMM\/Too4+Sfb2ILi6amppmCBw3fFlAOGFiP\/W+0CJt3cvjF28lYLyKYxgEIJAVBCBgsjJT6GfsCOQLISlPBHfILmB0T4XqsPrL99lnn41FwDj1QX3Q9TCQEgAsYOzP6OP68Y9\/TPv27cuFh+yhHvtf9lEImHyLfPXwiH1BM\/dDLaI9c+ZM7v\/5mYceeoj6+vqIhQ0vjuUQ29DQECmR6SSY1K4nFjC6R8uLgMm3XsfLGhi7x8cv3krAhLGWJvaXEA0CgQIIQMCAHkDAhwem0Ieeq0mjgFEfQKd1KPzxdhMwHGJRl+4lUELG7wc1iAcmiIBxEmEsctTFYoU9MX\/+539OjzzyCPGaHCcvmp0eXjwZ9jLFCBj7+hy\/eEPAwMRJRQACRurMYlxFIxDEA+M3hKTWrKjwjZePYyGvg73P+Rb76mEIPwJGFzIqJLNp0ybr1yp8E7cHhtt2Wt9h9zgwtrfeeitNTExYooWFmXru85\/\/vLVuiS9ep6TWu3jxrngp4zWEZN9G7TTXQQWMV24V\/eKgAiAQEwIQMDEBjWayh0BQAZNvEe97771nhWL4Ujk9nBZkOn3E7OgV2oWke4Hs6yy4HuWR8RpCcloDo39E\/X5QdQ8M94cxyZfPxm2btS6oCi3i5XJKLNrnhxfk6mExxueLX\/zijKSA+RYh58vLonaD7d69O5doL8giXu630y4kP4IRi3izZ3vQY28IQMB4wwmlDEQgiIBh8aBvo7Zni7Vnln3llVfo7rvvzu1M0u8X2uXiVcDwtNnbZOHy1ltv0a\/\/+q9bH2m+CoXG2FOhRI+iQaEFtH48MPxxtQsqXRR4FTC6QFHb2+1rdZzEnN0DpsZXaLxePBn5ti3r9Tpto9bXKKm+2DnkVzDmW0Rs4CuNIQtDAAJG2IRiOMkhkPbEaskhY2bLXhLZxYUMEtnFhTTaiRMBCJg40UZbohGAgBE9vb4Hp++C0hdC+66oyAeQxK5IAPF4ahGAgEnt1KBjWUMAAiZrMxZ9f9Pg+UiTJyh6xNGCSQj8f\/04W15rCIuGAAAAAElFTkSuQmCC","height":337,"width":560}}
%---

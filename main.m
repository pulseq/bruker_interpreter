%% Pulseq interpeter for Bruker systems
% This script reads a Pulseq file and converts it to a PPG file necessary
% for Bruker MRI scanners. Additional RF shapes are also written in the
% appropriate format. Currently, ParaVision variables must be set manually
% to execute the sequence correctly.
%
% TODO:
%  - Arbitrary gradient shapes
%  - Different readout events (e.g. freq offset for RF spoiling)
%  - Program a Paravison Method (PVM): This would read the file, write the
%    PPG files, WAVE files and set the low-level variables. It would remove
%    the dependency on MATLAB.
%
% Kelvin Layton 2015
% email: pulseq.mr@uniklinik-freiburg.de 

%% Read the Pulseq file
seq = mr.Sequence();
seq.read('gre.seq');

%% Write PPG file
blocks=seq.blockEvents;
numBlocks=size(seq.blockEvents,1);


gamma         = 2*pi*42.576e6;         	% rad/s/T
max_T_m       = 32e-3;                	% T/m
max_rad_ms_mm = gamma*max_T_m*1e-6;     % rad/ms/mm (JEMRIS units)

PREEMP_rise_time = 0.114;
PREEMP_grad_cal_const = 630e-3*gamma*1e-2/(2*pi);
PREEMP_grad_cal_const = 313000;
Pulseq_to_Trim = 100/(PREEMP_grad_cal_const)*1e-2;

    
Pulseq_to_Rabi = 1;

fprintf('Writing PPG file...\n');

fid=fopen('arbitrary.ppg','w');
fprintf(fid,'; PPG file\n');
fprintf(fid,';\n');
fprintf(fid,'; Created from a Pulseq file\n');
fprintf(fid,';\n\n');

fprintf(fid,';-------include files-------\n');
%fprintf(fid,'#include <Avance.incl>\n');
%fprintf(fid,'#include <DBX.include>\n\n');
fprintf(fid,'#include <MRI.include>\n\n');

fprintf(fid,';-------definitions-------\n');
fprintf(fid,'preset off\n\n');
%fprintf(fid,'INIT_DEVICES\n\n');

fprintf(fid,';-------sequence-------\n');
fprintf(fid,'start,');


for iB=1:numBlocks
    if mod(iB, round(numBlocks/10))==0 || iB==numBlocks
        fprintf('Processing %d blocks (%d%%)\n', numBlocks, round(100*iB/numBlocks))
    end
        
    block = seq.getBlock(iB);
    isGrad = cellfun(@(x)~isempty(x),{block.gx,block.gy,block.gz});
    isRF=~isempty(block.rf);
    isADC = ~isempty(block.adc);
    
    amp=zeros(3,1);
    gradNames={'gx','gy','gz'};
    for iC=1:3
        grad = block.(gradNames{iC});
        if isGrad(iC) && strcmp(grad.type,'trap')
            % Trapezoid gradient
            tRampUp(iC) = grad.riseTime*1e6;
            tFlat(iC) = grad.flatTime*1e6;
            tRampDown(iC) = grad.fallTime*1e6;
            amp(iC) = grad.amplitude*Pulseq_to_Trim;
        end
    end

    % Delay
    if ~isempty(block.delay)
        delay = round(block.delay.delay*1e6);   % us
        fprintf(fid,'\t%.0fu',delay);
    end
    
    % Gradients on
    if any(isGrad)
        ramp=max(tRampUp);
        assert(norm(ramp-PREEMP_rise_time*1e3)<1e-10,'Invalid ramp time!');
        fprintf(fid,'\t%.0fu grad{(%.3f)|(%.3f)|(%.3f)}',ramp,amp(1),amp(2),amp(3));
    end
    
    
    % RF
    if isRF
        % Gate pulse (not needed?)
        %fprintf(fid,'\n\t10u gatepulse 1\n');
    end
    if isRF
        rf=block.rf;
        rfId = blocks(iB,2);
        freqOffset = rf.freqOffset;
        txPhaseOffset = rf.phaseOffset;

        fprintf(fid,'\n\t(p%d:sp%d ph%d):f1',rfId-1,rfId-1,rfId-1);
    end
    fprintf(fid,'\n');
    
    % Gradient flat-top
    if any(isGrad) && ~isRF && ~isADC
        flat=max(tFlat);
        fprintf(fid,'\t%.0fu \t\t\t; flat top\n',flat);
    end
    
    % ADC
    if isADC
        adc = block.adc;
        delay = ceil(adc.delay*1e6/10)*10;
        rxPhaseOffset = adc.phaseOffset*180/pi;

        fprintf(fid,'\tADC_INIT(ph0, ph1)\n');
        fprintf(fid,'\t%.0fu ADC_START\n',round(1e6*adc.dwell*adc.numSamples));
    end
    
    % Gradients off
    if any(isGrad)
        ramp=max(tRampDown);
        fprintf(fid,'\t%.0fu groff\n',ramp);
    elseif isADC
        fprintf(fid,'\t5m\n');  % No gradients but ADC, add delay
    end
    
    if isADC
        fprintf(fid,'\t10u ADC_END\n');
    end
end

fprintf(fid,'SETUP_GOTO(start)\n');
fprintf(fid,'exit\n\n');

fprintf(fid,';-------phase definitions-------\n');
fprintf(fid,'ph0 = %d\t\t\t; excitation phase\n',txPhaseOffset);
fprintf(fid,'ph1 = %d\t\t\t; reference phase\n',rxPhaseOffset);

fprintf('Done\n');

%% Write RF pulse files
fprintf('Writing pulse files...\n');
fprintf('MUST SET THE FOLLOWING PARAMETERS MANUALLY:\n');
fprintf('=== GENERAL ===\n');
fprintf('PULPROG = arbitrary.ppg\n');

keys=seq.rfLibrary.keys;
for iPulse=1:length(keys)
    data = seq.rfLibrary.data(keys(iPulse)).array;
    amp = data(1);
    magShape = data(2);
    phaseShape = data(3);
    freqOffset = data(4);
    phaseOffset = data(5);
    
    shapeData = seq.shapeLibrary.data(magShape).array;
    compressed.num_samples = shapeData(1);
    compressed.data = shapeData(2:end);
    mag = mr.decompressShape(compressed)*100;   % Pecentage of maximum
    
    shapeData = seq.shapeLibrary.data(phaseShape).array;
    compressed.num_samples = shapeData(1);
    compressed.data = shapeData(2:end);
    phase = mr.decompressShape(compressed)*360; % 0-360 degrees
    
    integralRatio=sum(mag)/100/length(mag);
    
    filename=sprintf('pulse%d.exc',iPulse-1);
    fid=fopen(filename,'w');

    fprintf(fid,'##TITLE=%s\n',filename);
    fprintf(fid,'##JCAMP-DX= 5 BRUKER JCAMP library\n');
    fprintf(fid,'##DATA TYPE= Shape Data\n');
    fprintf(fid,'##ORIGIN= MATLAB\n');
    fprintf(fid,'##OWNER= <klayton>\n');
    fprintf(fid,'##DATE= 20/10/14\n');
    fprintf(fid,'##TIME= 14:29:34\n');
    fprintf(fid,'##MINX= %.6e\n',min(mag));
    fprintf(fid,'##MAXX= %.6e\n',max(mag));
    fprintf(fid,'##MINY= %.6e\n',min(phase));
    fprintf(fid,'##MAXY= %.6e\n',max(phase));
    fprintf(fid,'##$SHAPE_EXMODE= Excitation\n');
    fprintf(fid,'##$SHAPE_TOTROT= 9.000000e+01\n');
    fprintf(fid,'##$SHAPE_BWFAC= 2.025000e+01\n');
    fprintf(fid,'##$SHAPE_INTEGFAC= %.6e\n',integralRatio);
    fprintf(fid,'##$SHAPE_REPHFAC= 50\n');
    fprintf(fid,'##$SHAPE_TYPE= conventional\n');
    fprintf(fid,'##$SHAPE_MODE= 0\n');
    fprintf(fid,'##NPOINTS= %d\n',length(mag));
    fprintf(fid,'##XYPOINTS= (XY..XY)\n');
    
    fprintf(fid,'%.6e, %.6e\n',[mag,phase].');
    
    fprintf(fid,'##END=\n\n');
    
    % Must set reference attenuation manually for each scan to ensure
    % correct flip angles
    PVM_RefAttCh1 = 17.092;   % dB 

    flipRef=pi/2;           % rad
    tauRef=1e-3;            % ms
    refB1=flipRef/tauRef;   % rad/s
    B1=2*pi*amp;            % rad/s
    
    atten = PVM_RefAttCh1-20*log10(B1/refB1);
    
    fprintf('=== PULSE #%d ===\n',iPulse-1);
    fprintf('Set TPPQ[%d].name to %s\n',iPulse-1,filename);
    fprintf('Set TPPQ[%d].offset to 0 Hz\n',iPulse-1);
    fprintf('Set TPPQ[%d].power to PVM_RefAttCh1-20*log10(B1/refB1)\n',iPulse-1);
    fprintf('Set TPPQ[%d].power to PVM_RefAttCh1-%.3f\n',iPulse-1,20*log10(B1/refB1));
    fprintf('Set TPPQ[%d].power to %.3f dB (PVM_RefAttCh1=%.2f dB)\n',iPulse-1,atten,PVM_RefAttCh1);
    fprintf('Set P[%d] to %d us\n',iPulse-1,length(mag)*1);

end

%% Print ADC parameters
fprintf('=== ADC ===\n');
fprintf('ACQ_O1B_list[0] = %.2f Hz\n', adc.freqOffset);
fprintf('PREEMP_grad_cal_const=%.0f Hz/cm\n',PREEMP_grad_cal_const);
fprintf('SW_h = %.2f Hz\n', 1/(adc.dwell));
fprintf('ACQ_size[0] = %d (samples = %d)\n',2*adc.numSamples,adc.numSamples );
fprintf('\nFinished\n');


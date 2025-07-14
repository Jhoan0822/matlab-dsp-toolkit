classdef GetSyntheticSignals < handle
    properties
        audioData
        time
        fs %Frecuencia de muestreo 
        f
    end
    methods
        function Chirp = GenChirp(obj)
            dlg_title = 'Chirp Parameters';
            prompt = { 'Amplitude:','Sampling Frequency(Hz):', 'Initial Instant(s):','Final Instant(s):','Start Frequency(Hz):','Final Frequency(Hz):','Method:'};
            num_lines = 1;
            defaultans = {'1','100','0','10', '0.5','5','linear'};

            answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
            if size(answer,1) == 0
                Chirp = 0;
                return;
            end

            A = str2double(answer{1});
            Fs = str2double(answer{2});
            ti = str2double(answer{3});
            tf = str2double(answer{4});
            fi = str2double(answer{5});
            ff = str2double(answer{6});
            method = answer{7};

            t = ti:1/Fs:tf-1/Fs;

            x = A*chirp(t,fi,tf,ff,method);
            
            obj.audioData = x;
            obj.time = t;
            obj.fs = Fs;
            Chirp = 1;
            obj.f = fi;
        end
        function Sinc = GenSinc(obj)
            dlg_title = 'Sinc Parameters';
            prompt = { 'Amplitude:','Sampling Frequency(Hz):', 'Initial Instant(s):','Final Instant(s):','Frequency(Hz):'};
            num_lines = 1;
            defaultans = {'1','100','-10','10','5'};

            answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
            if size(answer,1) == 0
                Sinc = 0;
                return;
            end

            A = str2double(answer{1});
            Fs = str2double(answer{2});
            ti = str2double(answer{3});
            tf = str2double(answer{4});
            freq = str2double(answer{5});

            t = ti:1/Fs:tf-1/Fs;

            x = A*sinc(freq.*t);
            
            obj.audioData = x;
            obj.time = t;
            obj.fs = Fs;
            Sinc = 1;
            obj.f = freq;
        end
        function Sin = GenSin(obj)
            dlg_title = 'Sin Parameters';
            prompt = { 'Amplitude:','Sampling Frequency(Hz):', 'Initial Instant(s):','Final Instant(s):','Frequency(Hz):'};
            num_lines = 1;
            defaultans = {'1','100','0','10','5'};

            answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
            if size(answer,1) == 0
                Sin = 0;
                return;
            end

            A = str2double(answer{1});
            Fs = str2double(answer{2});
            ti = str2double(answer{3});
            tf = str2double(answer{4});
            freq = str2double(answer{5});

            t = ti:1/Fs:tf-1/Fs;

            x = A*sin(2*pi()*freq.*t);
            
            obj.audioData = x;
            obj.time = t;
            obj.fs = Fs;
            Sin = 1;
            obj.f = freq;

        end
        function Cos = GenCos(obj)
            dlg_title = 'Cos Parameters';
            prompt = { 'Amplitude:','Sampling Frequency(Hz):', 'Initial Instant(s):','Final Instant(s):','Frequency(Hz):'};
            num_lines = 1;
            defaultans = {'1','100','0','10','5'};

            answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
            if size(answer,1) == 0
                Cos = 0;
                return;
            end

            A = str2double(answer{1});
            Fs = str2double(answer{2});
            ti = str2double(answer{3});
            tf = str2double(answer{4});
            freq = str2double(answer{5});

            t = ti:1/Fs:tf-1/Fs;

            x = A*cos(2*pi()*freq.*t);
            
            obj.audioData = x;
            obj.time = t;
            obj.fs = Fs;
            Cos = 1;
            obj.f = freq;

        end
        function Saw = GenSaw(obj)
            dlg_title = 'Sawtooth Parameters';
            prompt = { 'Amplitude:','Sampling Frequency(Hz):', 'Initial Instant(s):','Final Instant(s):','Frequency(Hz):'};
            num_lines = 1;
            defaultans = {'1','100','0','10','5'};

            answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
            if size(answer,1) == 0
                Saw = 0;
                return;
            end

            A = str2double(answer{1});
            Fs = str2double(answer{2});
            ti = str2double(answer{3});
            tf = str2double(answer{4});
            freq = str2double(answer{5});

            t = ti:1/Fs:tf-1/Fs;

            x = A*sawtooth(2*freq.*t);
            
            obj.audioData = x;
            obj.time = t;
            obj.fs = Fs;
            Saw = 1;
            obj.f = freq;

        end
        function Ramp = GenRamp(obj)
            dlg_title = 'Ramp Parameters';
            prompt = { 'Slope:','Sampling Frequency(Hz):', 'Initial Instant(s):','Final Instant(s):'};
            num_lines = 1;
            defaultans = {'1','100','-10','10'};

            answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
            if size(answer,1) == 0
                Ramp = 0;
                return;
            end

            A = str2double(answer{1});
            Fs = str2double(answer{2});
            ti = str2double(answer{3});
            tf = str2double(answer{4});

            t = ti:1/Fs:tf-1/Fs;

            x =A*t.*(t>=0);
            
            obj.audioData = x;
            obj.time = t;
            obj.fs = Fs;
            Ramp = 1;
            obj.f = 0;

        end
        function Trig = GenTrig(obj)
            dlg_title = 'Trig Parameters';
            prompt = { 'Amplitude:','Sampling Frequency(Hz):', 'Initial Instant(s):','Final Instant(s):','Frequency(Hz)'};
            num_lines = 1;
            defaultans = {'1','100','0','10','2'};

            answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
            if size(answer,1) == 0
                Trig = 0;
                return;
            end

            A = str2double(answer{1});
            Fs = str2double(answer{2});
            ti = str2double(answer{3});
            tf = str2double(answer{4});
            freq = str2double(answer{5});
            t = ti:1/Fs:tf-1/Fs;

            x = A * sawtooth(2*pi*freq*t,1/2);
            
            obj.audioData = x;
            obj.time = t;
            obj.fs = Fs;
            Trig = 1;
            obj.f = freq;

        end
        function Step = GenStep(obj)
            dlg_title = 'Trig Parameters';
            prompt = { 'Amplitude:','Sampling Frequency(Hz):', 'Initial Instant(s):','Final Instant(s):'};
            num_lines = 1;
            defaultans = {'1','100','-10','10'};

            answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
            if size(answer,1) == 0
                Step = 0;
                return;
            end

            A = str2double(answer{1});
            Fs = str2double(answer{2});
            ti = str2double(answer{3});
            tf = str2double(answer{4});
            t = ti:1/Fs:tf-1/Fs;

            x = A*(t>=0);
            
            obj.audioData = x;
            obj.time = t;
            obj.fs = Fs;
            Step = 1;
            obj.f = 0;

        end
    end
end
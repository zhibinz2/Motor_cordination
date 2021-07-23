function Sig=rampdamp(Sig,tc,fs)

[nch, dim] = min(size(Sig));
if nch >2 || length(Sig)<fs*0.5e-3
    error('Wrong stimulus size!');
end

if dim ==2
    Sig = Sig.';
end

if tc~=0
    ramp=sin(2*pi*(1/(4*tc))*(1:round(fs*tc))/fs).^2;
    damp=fliplr(ramp);
    for ich = 1:nch
        Sig(ich,1:length(ramp))=Sig(ich,1:length(ramp)).*ramp;
        Sig(ich,end-length(ramp)+1:end)=Sig(ich,end-length(ramp)+1:end).*damp;
    end
end

%eof
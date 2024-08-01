from scipy import signal
from matplotlib import pyplot as plt 
from numpy import log10, abs
#
#standard make filter with b,a coefficients for teaching. 
def makefilter(sr,fp,fs,gp=3,gs=20):
	""" 	Wrapper function around scipy filter functions.  
	Makes it convenient by providing frequency parameters in terms of 
	frequencies in Hz.   
	INPUT: 	sr - sampling rate in Hz. 
		fp - pass frequency in Hz
		fs - stop frequency in Hz
		gp - pass band ripple in dB, default 3 dB
		gs - stop band attenuation in dB, default 20 dB
		doPlot - make a plot of filter gain versus frequency, default 'no'
	OUTPUT: b,a filter coefficients. 
			w,h for making bodeplot 
	Automatically detects the type of filter.  if fp < fs the filter
	is low pass but if fp > fs the filter is highpass.  
	It is recommended that you DO NOT use this filter, but instead use 
	makefiltersss to improve stability in high filter order scenarios """
#
#set up filter parameters

	fn = sr/2
	wp = fp/fn
	ws = fs/fn
#get the filter order

	n,wn = signal.buttord(wp,ws,gp,gs);                                                            
#design the filter

#lowpass 
	if fp < fs:
		b,a = signal.butter(n,wn,btype='lowpass')
#highpass
	if fs < fp:
		b,a = signal.butter(n,wn,btype='highpass')
#get filter respons function	
	w,h = signal.freqz(b,a,fs=sr)
	return b,a,w,h
#
#sos version of makefilter.  Use this one expecially for high filter orders. 
def makefiltersos(sr,fp,fs,gp=3,gs=20):
	""" 	Wrapper function around scipy filter functions.  
	Makes it convenient by providing frequency parameters in terms of 
	frequencies in Hz.   
	INPUT: 	sr - sampling rate in Hz. 
		fp - pass frequency in Hz
		fs - stop frequency in Hz
		gp - pass band ripple in dB, default 3 dB
		gs - stop band attenuation in dB, default 20 dB
		doPlot - make a plot of filter gain versus frequency, default 'no'
	OUTPUT: sos filter coefficients. 
			w,h for making bode plot 
	Automatically detects the type of filter.  if fp < fs the filter
	is low pass but if fp > fs the filter is highpass.  """
#
#set up filter parameters

	fn = sr/2
	wp = fp/fn
	ws = fs/fn
#get the filter order

	n,wn = signal.buttord(wp,ws,gp,gs);                                                            
#design the filter

#lowpass 
	if fp < fs:
		sos = signal.butter(n,wn,btype='lowpass',output='sos')
#highpass
	if fs < fp:
		sos = signal.butter(n,wn,btype='highpass',output='sos')
#get filter respons function	
	w,h = signal.sosfreqz(sos,fs=sr)
	return sos,w,h


def narrowfilter(fp,fs,gpass=0.2,gstop=20,samplingrate=200):
	""" 	Wrapper function around scipy elliptical filter functions
			to make narrow band filters.   
			Makes it convenient by providing frequency parameters in terms of 
			frequencies in Hz.
			To be useful, very high sampling rates should be avoided and signal should be
			downsampled to sampling rates 200 Hz or less.  Dont forget to lowpass filter 
			using makefilter before you downsample. 
			Note that the passband must be inside the stopband, or this wont work   
		INPUT: 	
			fp - passband frequencies in Hz, [lo,hi]
			fs - stopband frequencies in Hz, [lo,hi]
			gp - pass band ripple in dB, default 0.2 dB
			gs - stop band attenuation in dB, default 20 dB
		    sr - sampling rate in Hz (recommended <= 200 Hz). 
		OUTPUT: esos elliptical filter coefficients. 
				w,h for making bode plot """

	n,wn = signal.ellipord(fp,fs,gpass,gstop,fs=samplingrate)
	esos = signal.ellip(n,gpass,gstop,wn,btype='bandpass',output='sos',fs=samplingrate)
	w,h = signal.sosfreqz(esos,fs=samplingrate)
	return esos,w,h

#make plot of filter response
def bodeplot(w,h):
	plt.plot(w,20*log10(abs(h)))
	plt.xlabel('Frequency(Hz)')
	plt.ylabel('Gain (dB)')
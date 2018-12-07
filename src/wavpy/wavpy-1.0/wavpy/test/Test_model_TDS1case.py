import matplotlib.pyplot as plt
import numpy as np
import sys
import os
HOME = os.getenv("HOME")
wavpyDIR = HOME + "/wavpy-1.0/wavpy/"
sys.path.insert(1,wavpyDIR)
import wavpy
import time 

print "0) Setting parameters"

#Declare waveform-model object
modelWAV = wavpy.ZaVoModel_GNSSR()

#Set waveform-model aspects
modelWAV.sampling_rate = 4091750.0 # 244 nanoseconds per lag (rate in samples/sec)
modelWAV.wav_length = 64
modelWAV.exponent_wav_model_length = 8
modelWAV.gnss_signal.weight_PY = 0.0
modelWAV.gnss_signal.weight_M = 0.0
modelWAV.gnss_signal.weight_IM = 0.0

#Receiver aspects
BW = 5000000.0 #Baseband bandwidth in Hz
F = 3.0 #Noise figure in dB
T_ant = 200.0 #Antenna temperature in K
Gain = 0.0 #Directivity is given, for this particular case, in the antenna pattern. No additional gain is applied.
Isotropic_ant = 0 #Non-isotropic antenna
modelWAV.receiver_Down.set_receiver_params(Gain, T_ant, F, BW, Isotropic_ant)
#Approximation of UK-TDS antenna pattern
angles_deg = [0.0, 20.0, 40.0, 45.0, 50.0, 65.0, 80.0, 90.0]
pattern_dB = [12.8, 8.26, -6.3, -15.0, -7.44, -3.3, -5.7, -15.0]
#Computes the whole antenna pattern by means of spline interpolation of the segments above a minimum value (here set as -15.0 dB)
modelWAV.receiver_Down.set_antenna_pattern_interp(angles_deg, pattern_dB, -15.0)
#If you have the whole antenna pattern as a 2D array antpattern[phi, theta] with phi=0...180 and theta=0..359 with 1 degree resolution, then:
#modelWAV.receiver_Down.set_antenna_whole_pattern(antpattern)

#Geometry aspects
#Distances in km and velocities in km/s
modelWAV.geometry.set_ECEFpos_Receiver([-3147.5815062814811, 3434.4546494885967, -5238.6507503218688])
modelWAV.geometry.set_ECEFvel_Receiver([5.1817276946851334, -2.7036802379824494, -4.8853405725848206])
modelWAV.geometry.set_ECEFpos_Transmitter([-14849.239465341682, 15882.877420247279, -14561.494119085249])
modelWAV.geometry.set_ECEFvel_Transmitter([0.398198421670646, -1.8727892375129452, -2.4768106814572866])
modelWAV.geometry.compute_specular_point(0)
# If the user is not interested in particular locations of the transmitter and receiver, but only the altitude of the receiver and the elevation angle of the observation, then the geometry can be set as: 
# modelWAV.geometry.set_geometry_from_ElevHeightsSpec(elevation_deg, heightReceiver, heightTransmitter, lonSpecular, latSepcular, azimTransmitter, heightSpecular)
# and the tangencial velocities:
# modelWAV.geometry.set_tangEarthVel_Receiver(velocity, specAzim_deg)
# modelWAV.geometry.set_tangEarthVel_Transmitter(velocity, specAzim_deg)
#
# If we need the specular point information stored in a variable, then do:
# SPp = modelWAV.geometry.get_ECEFpos_Specular()
# in general, if you want to check at any moment all the geometric parameters of your object, simply use:
# modelWAV.geometry.dump_parameters()

#Surface aspects
# Several parameters are set by default, such as Temperature and Salinity (permitivity of the water), some Gram-Charlier parameters, ... if you want to check the values at which your object is set, simply use:
#modelWAV.surface.dump_surf()
# There are several ways for the user to enter the roughness conditions.
# 1) by typing the mss or mss in up-wind and cross-wind directions:
# In two directions:
# modelWAV.surface.mss_1 = myvalue_of_mss_upwind
# modelWAV.surface.mss_2 = myvalue_of_mss_crosswind
# or omnidirectional:
# modelWAV.surface.mss = myvalue_of_total_mss
# 2) by having entered the wind speed and using Katzberg'06 relationship:
modelWAV.surface.wind_U10_speed = 6.0 #In m/s
modelWAV.surface.compute_mss_from_wind()
# 3) extracted from a spectrum. 
# The spectra can be either given by the user or computed by wavpy using Elfouhaily's 97.
# See Test_Elfouhaily_spectrum.py for examples on how to set different threshold values, how to run Elfouhaily's or enter by hand the spectrum, etc.

#Compute waveform model
print "1) Computing Waveform + DDM"
modelWAV.delta_doppler = 500.0
modelWAV.ddm_half_dopplers = 10
start_time = time.time()
modelWAV.compute_waveform(0, 1, 1, 0, 0)
print("=> Computation Time for waveform and DDM: %s seconds" % (time.time() - start_time))

#Get waveform and associated range
range_wav = modelWAV.waveform_POW.get_range_waveform(modelWAV.wav_length)
range_spec = range_wav - modelWAV.geometry.geometric_delay*1000.0
wav_pow = modelWAV.waveform_POW.get_waveform(modelWAV.wav_length)

#First sample has the noise_level value
noise_level = wav_pow[0]

#Plot SNR waveform
plt.figure(1)
plt.plot(range_spec/1000.0, 10.0*np.log10(wav_pow/noise_level), '.-')
plt.grid()
plt.title("Example of UK-TDS waveform model (elevation = %3.2f deg)" % (modelWAV.geometry.elevation))
plt.xlabel("Range from specular [km]")
plt.ylabel("SNR [dB]")

#Plot DDM
print "2) Plotting DDM"
DDM = np.zeros([(modelWAV.ddm_half_dopplers*2 + 1), modelWAV.wav_length])
  
for freq_ind in range(0, (modelWAV.ddm_half_dopplers*2 + 1)):
  DDM[freq_ind, :] = modelWAV.get_DDM_doppler_slice((modelWAV.ddm_half_dopplers - freq_ind), modelWAV.wav_length)
  
plt.figure(2)

plt.imshow(DDM[:,:], extent=[range_spec[0]/1000.0, range_spec[-1]/1000.0, -modelWAV.delta_doppler*float(modelWAV.ddm_half_dopplers), modelWAV.delta_doppler*float(modelWAV.ddm_half_dopplers)], vmin=noise_level, vmax=max(DDM[modelWAV.ddm_half_dopplers,:]), cmap='jet', aspect='auto')
plt.title("Example of UK-TDS DDM model (elevation = %3.2f deg)" % (modelWAV.geometry.elevation))
plt.xlabel("Range from specular [km]")
plt.ylabel("Doppler [Hz]")
plt.colorbar()

file_out = open("DDM_norm_Test_model_TDS1case.txt","w")
file_out.write("Doppler [Hz] | Range from specular [m] | DDM [norm]\n")
ddm_max = max(DDM[modelWAV.ddm_half_dopplers,:])
for freq_ind in range(0, (modelWAV.ddm_half_dopplers*2 + 1)):
  for del_ind in range(modelWAV.wav_length):
    doppler_freq = modelWAV.delta_doppler*float(freq_ind) - modelWAV.delta_doppler*float(modelWAV.ddm_half_dopplers)
    file_out.write("%f %f %f\n" % (doppler_freq, range_spec[del_ind], DDM[freq_ind, del_ind]/ddm_max))
file_out.close()

print "3) Normalized DDM stored at file DDM_norm_Test_model_TDS1case.txt for comparison against DDM_norm_Test_model_TDS1case_REF.txt"


#Compute ddm covariance model
print "4) Computing Waveform with COVARIANCE model"
start_time = time.time()
modelWAV.compute_waveform(0, 1, 1, 1, 0)
print("=> Computation Time for waveform: %s seconds" % (time.time() - start_time))

#Get waveform and associated range
range_wav = modelWAV.waveform_POW.get_range_waveform(modelWAV.wav_length)
range_spec = range_wav - modelWAV.geometry.geometric_delay*1000.0
wav_mean = modelWAV.waveform_POW.get_waveform(modelWAV.wav_length)

#Get waveform covariance
cov_wav = np.zeros([modelWAV.wav_length, modelWAV.wav_length])
  
for cov_ind in range(0, modelWAV.wav_length):
  cov_wav[modelWAV.wav_length - 1 - cov_ind, :] = modelWAV.get_cov_slice(cov_ind, modelWAV.wav_length)

plt.figure(3)
plt.imshow(cov_wav[:,:], extent=[range_spec[0]/1000.0, range_spec[-1]/1000.0, range_spec[0]/1000.0, range_spec[-1]/1000.0], vmin=min(cov_wav.min(axis=0)), vmax=max(cov_wav.max(axis=0)), cmap='jet', aspect='auto')
plt.title("Waveform covariance model")
plt.xlabel("Range from specular [km]")
plt.ylabel("Range from specular [km]")
plt.colorbar()

print "5) Computing noisy waveforms using stored covariance model and mean waveform"
plt.figure(4)
num_wavs = 1000
noisy_wavs = np.zeros([num_wavs, modelWAV.wav_length])
start_time = time.time()
for i in range(num_wavs):
  seed_random_generator = i
  noisy_wavs[i,:] = modelWAV.get_noisy_waveform(modelWAV.wav_length, i)
  if(i < 10):# Plot the 10 first noisy waveforms
    plt.plot(range_spec/1000.0, 10.0*np.log10(noisy_wavs[i,:]/wav_mean[0]), '--')

print("=> Computation Time for 1000 noisy waveform: %s seconds" % (time.time() - start_time))

plt.plot(range_spec/1000.0, 10.0*np.log10(wav_mean/wav_mean[0]), 'k.-', label="Mean waveform")
plt.grid()
plt.title("Example of 10 simulated noisy waveforms")
plt.xlabel("Range from specular [km]")
plt.ylabel("SNR [dB]")
plt.legend(loc=1)

plt.figure(5)
cov_sim = np.cov(noisy_wavs, rowvar=False)
plt.imshow(np.flipud(cov_sim[:,:]), extent=[range_spec[0]/1000.0, range_spec[-1]/1000.0, range_spec[0]/1000.0, range_spec[-1]/1000.0], vmin=min(cov_sim.min(axis=0)), vmax=max(cov_sim.max(axis=0)), cmap='jet', aspect='auto')
plt.title("Covariance from 1000 simulated noisy wavs.")
plt.xlabel("Range from specular [km]")
plt.ylabel("Range from specular [km]")
plt.colorbar()

print "6) Computing Waveform + DDM with COVARIANCE model"
start_time = time.time()
modelWAV.compute_waveform(0, 1, 1, 1, 1)
print("=> Computation Time for waveform and DDM: %s seconds" % (time.time() - start_time))

#Get DDM
DDM_mean = np.zeros([(modelWAV.ddm_half_dopplers*2 + 1), modelWAV.wav_length])
for freq_ind in range(0, (modelWAV.ddm_half_dopplers*2 + 1)):
  DDM_mean[freq_ind, :] = modelWAV.get_DDM_doppler_slice((modelWAV.ddm_half_dopplers - freq_ind), modelWAV.wav_length)
  
#Plot DDM  
plt.figure(6)
plt.imshow(DDM_mean[:,:], extent=[range_spec[0]/1000.0, range_spec[-1]/1000.0, -modelWAV.delta_doppler*float(modelWAV.ddm_half_dopplers), modelWAV.delta_doppler*float(modelWAV.ddm_half_dopplers)], vmin=min(DDM_mean[modelWAV.ddm_half_dopplers,:]), vmax=max(DDM_mean[modelWAV.ddm_half_dopplers,:]), cmap='jet', aspect='auto')
plt.title("Mean DDM from covariance model")
plt.xlabel("Range from specular [km]")
plt.ylabel("Doppler [Hz]")
plt.colorbar()

print "7) Computing noisy DDM using stored covariance model and mean DDM"
seed_random_generator = 1
start_time = time.time()
DDM_noise = modelWAV.get_noisy_DDM(modelWAV.wav_length*(modelWAV.ddm_half_dopplers*2 + 1), seed_random_generator)
print("=> Computation Time for 1 noisy DDM: %s seconds" % (time.time() - start_time))

DDM_noise_plot = np.flipud(DDM_noise.reshape((modelWAV.ddm_half_dopplers*2 + 1), modelWAV.wav_length))

plt.figure(7)
plt.imshow(DDM_noise_plot[:,:], extent=[range_spec[0]/1000.0, range_spec[-1]/1000.0, -modelWAV.delta_doppler*float(modelWAV.ddm_half_dopplers), modelWAV.delta_doppler*float(modelWAV.ddm_half_dopplers)], vmin=min(DDM_mean[modelWAV.ddm_half_dopplers,:]), vmax=max(DDM_mean[modelWAV.ddm_half_dopplers,:]), cmap='jet', aspect='auto')
plt.title("Example of noisy DDM from covariance model")
plt.xlabel("Range from specular [km]")
plt.ylabel("Doppler [Hz]")
plt.colorbar()

plt.show()

exit()

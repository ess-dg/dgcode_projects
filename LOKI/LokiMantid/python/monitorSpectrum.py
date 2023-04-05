import numpy as np
import G4Units.Units as Units
import math as m

def neutron_angstrom_to_meters_per_second(l_aangstrom):
  if l_aangstrom==0:
    return m.inf
  h_Planck = 6.62606896e-34 *Units.joule*Units.s
  c_light = 2.99792458e+8 *Units.m/Units.s
  neutron_mass_c2 = 939.56536 *Units.MeV
  return h_Planck*Units.s * c_light * c_light / (Units.m*Units.angstrom * neutron_mass_c2 * l_aangstrom)

def createTofSpectrumWithUniformPart(totalIntensity, tMin, tMax, spectrumStart, spectrumEnd, verbose=1):
  binNumberWithIntensity = tMax-tMin
  binValue = round(totalIntensity / binNumberWithIntensity)
  binError = round(np.sqrt(binValue)) #rounding probably not needed

  spectrumBinNumber = spectrumEnd - spectrumStart
  intensity = np.zeros(spectrumBinNumber)
  intensity[np.arange(tMin, tMax)] = binValue
  error = np.zeros(spectrumBinNumber)
  error[np.arange(tMin, tMax)] = binError
  tof = np.arange(spectrumStart, spectrumEnd)

  sumIntensity = np.sum(intensity)
  if verbose:
    print(f'    Created uniform TOF spectrum with:')
    print(f'        sum intensity: {sumIntensity} (intended intensity: {totalIntensity})')
    print(f'        number of bins: {binNumberWithIntensity} (tmin: {tMin}, tmax: {tMax})')
    print(f'        intensity per bin: {binValue}')

  if abs(totalIntensity - sumIntensity) > binNumberWithIntensity: #probably not a rounding issue
    print(f'WARNING: Unexpected difference in intended intensity: {totalIntensity} and sum intensity in the spectrum: {np.sum(intensity)}')

  return tof, intensity, error

def createFloodSourceTofSpectrum(totalIntensity, lambdaMin, lambdaMax, sourceToMonitorDistance, spectrumStart=0, spectrumEnd=100000, verbose=1):
  velocityMin = neutron_angstrom_to_meters_per_second(lambdaMax)
  velocityMax = neutron_angstrom_to_meters_per_second(lambdaMin)
  tMin = round((sourceToMonitorDistance / velocityMax)*1e6) #ms
  tMax = round((sourceToMonitorDistance / velocityMin)*1e6) #ms

  return createTofSpectrumWithUniformPart(totalIntensity, tMin, tMax, spectrumStart, spectrumEnd, verbose)


def extractMcStasMonitorData(filename, verbose):
    if verbose:
       print(f"    Reading McStas monitor file: {filename}")
    with open(filename) as monitorFile:
        data = monitorFile.readlines()

        headerLines = 0
        for line in data:
            if line.startswith('#'):
                headerLines +=1
        header = data[headerLines-1].split(" ")
        monitorData = np.array([x.split(" ") for x in data[headerLines:-1]])

        if header[2] == 't':
            tof = monitorData[:, 0].astype(np.float32) #values in McStas monitor file are already in microseconds so no conversion factor
            wavelength = None
        else:
            wavelength = monitorData[:, 0].astype(np.float32)
            tof = None
        intensity = monitorData[:, 1].astype(np.float32)
        error = monitorData[:, 2].astype(np.float32)

        return tof, wavelength, intensity, error
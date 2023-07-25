
from Core import FindData
from tempfile import NamedTemporaryFile

class IdfCreator:
  '''Creates instrument definition file (IDF) for mantidpython processing'''

  def __init__(self, params):
    self.bankIds = params.get('aiming_bank_id')
    self.nominalSourceSampleDistance = params.get('nominal_source_sample_distance_meters')
    self.preSampleMonitorDistance = params.get('source_monitor_distance_meters')
    self.rearDetDistance = params.get('rear_detector_distance_m')
    self.pixelNumber = params.get('analysis_straw_pixel_number')
    self.xmlDetector = self.readXmlFile("LOKI_Definition_template_detectors.xml")
    self.xmlMonitor = self.readXmlFile("LOKI_Definition_template_monitors.xml")

    self.uncommentBanks()
    self.setDistances() #sample, monitor, rear detector bank distance
    self.setPixelNumber()

  # def __del__(self):

  def getDetectorIDF(self): #TODO getDetectorIDF and getMonitorIDF are way too similar to be 2 functions
    if not hasattr(self, 'tempDetectorIdfFile'):
      self.tempDetectorIdfFile = NamedTemporaryFile(prefix='LOKI_Definition_detector_', suffix='.xml')
      with open(self.tempDetectorIdfFile.name, 'w') as f:
        f.write(self.xmlDetector)
    return self.tempDetectorIdfFile.name
  def getMonitorIDF(self):
    if not hasattr(self, 'tempMonitorIdfFile'):
      self.tempMonitorIdfFile = NamedTemporaryFile(prefix='LOKI_Definition_monitor_', suffix='.xml')
      with open(self.tempMonitorIdfFile.name, 'w') as f:
        f.write(self.xmlMonitor)
    return self.tempMonitorIdfFile.name

  def readXmlFile(self, filename):
    filepath = FindData("LokiMantid", filename)
    with open(filepath, 'r') as file:
      return file.read()

  def uncommentBanks(self):
    def isLineToRemove(line):
      return any([(f'<!--BANK{id}_START' in line) or (f'BANK{id}_END-->' in line) for id in self.bankIds])

    newXmlLines = [line for line in self.xmlDetector.splitlines() if not isLineToRemove(line)]
    self.xmlDetector = '\n'.join(newXmlLines)

  def setDistances(self):
    ## Monitor IDF ##
    self.xmlMonitor = self.xmlMonitor.replace('<PLACEHOLDER_NOMINAL_SOURCE_SAMPLE_DISTANCE>', str(self.nominalSourceSampleDistance))
    self.xmlMonitor = self.xmlMonitor.replace('<PLACEHOLDER_PRESAMPLE_MONITOR_DISTANCE>', str(self.preSampleMonitorDistance))
    beamstopMonitorDistance = self.nominalSourceSampleDistance + self.rearDetDistance - 0.05 #hardcoded -5cm
    self.xmlMonitor = self.xmlMonitor.replace('<PLACEHOLDER_BEAMSTOP_MONITOR_DISTANCE>', str(beamstopMonitorDistance))

    ## Detector IDF ##
    self.xmlDetector = self.xmlDetector.replace('<PLACEHOLDER_NOMINAL_SOURCE_SAMPLE_DISTANCE>', str(self.nominalSourceSampleDistance))

    rearDetFrontToPackCentreOffset = 0.048877 #hardcoded offset to position the detector front
    rearDetPositionZ = self.nominalSourceSampleDistance + self.rearDetDistance + rearDetFrontToPackCentreOffset
    self.xmlDetector = self.xmlDetector.replace('<PLACEHOLDER_REAR_DETECTOR_POSITION_Z>', str(rearDetPositionZ))
    #note (wrong distance) <location x="0" y="0.008651" z="28.648877"/> <!-- +0.048877 detector front to pack centre -->

  def setPixelNumber(self):
    self.xmlDetector = self.xmlDetector.replace('<PLACEHOLDER_PIXEL_NUMBER>', str(self.pixelNumber))
    self.xmlDetector = self.xmlDetector.replace('<PLACEHOLDER_PIXEL_LENGTH_500MM>', "{:.9f}".format(0.5/self.pixelNumber))
    self.xmlDetector = self.xmlDetector.replace('<PLACEHOLDER_PIXEL_LENGTH_1000MM>', "{:.9f}".format(1.0/self.pixelNumber))
    self.xmlDetector = self.xmlDetector.replace('<PLACEHOLDER_PIXEL_LENGTH_1200MM>', "{:.9f}".format(1.2/self.pixelNumber))

    numberOfPacksInBank = [28, 8, 6, 8, 6, 14, 16, 10, 16]
    numberOfPixelsInBank = [numberOfPacksInBank[i] * 8 * 7 * self.pixelNumber for i in range(9)] #8 tube/pack, 7 straw/tube, pixelNuber pixel/straw

    def firstPixelOfBank(bankId):
      offset = 1 #pixel ids start from 1 in the IDF file
      for id in range(bankId):
        offset += numberOfPixelsInBank[id]
      return offset

    def lastPixelOfBank(bankId):
      return firstPixelOfBank(bankId) + numberOfPixelsInBank[bankId] - 1

    for bankId in range(9):
      self.xmlDetector = self.xmlDetector.replace(f'<PLACEHOLDER_BANK{bankId}_PIXEL_ID_START>', str(firstPixelOfBank(bankId)))
      self.xmlDetector = self.xmlDetector.replace(f'<PLACEHOLDER_BANK{bankId}_PIXEL_ID_END>', str(lastPixelOfBank(bankId)))

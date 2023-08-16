
from Core import FindData
from tempfile import NamedTemporaryFile
import LokiMantid.mantidApi as api

class workspaceCreator:
  '''Creates Mantid workspaces for mantidpython processing'''

  def __init__(self, params, singleNexus, idfCreation):
    self.bankIds = params.get('aiming_bank_id')
    self.nominalSourceSampleDistance = params.get('nominal_source_sample_distance_meters')
    self.preSampleMonitorDistance = params.get('source_monitor_distance_meters')
    self.rearDetDistance = params.get('rear_detector_distance_m')
    self.pixelNumber = params.get('analysis_straw_pixel_number')
    self.singleNexus = singleNexus
    #read idf templates
    self.xmlMonitor = self.readXmlFile("LOKI_Definition_template_monitors.xml")
    self.xmlDetectorsTemplate = self.readXmlFile("LOKI_Definition_template_detectors.xml")
    #edit idf content based on input params using the placeholders
    self.setDistances() #sample, monitor, rear detector bank distance
    self.firstPixelOfBank = [None]*9 #filled up in setPixelNumbers()
    self.lastPixelOfBank = [None]*9 #filled up in setPixelNumbers(), used later for getting bankId(workspace id) for a certain pixelId in getDetectorWorkspaceKey()
    self.setPixelNumbers()

    if self.singleNexus: #create one xml, containing all detector banks
      print("    Creating a single IDF/workspace for all banks")
      self.xmlDetectors = {'All': self.getXmlForBanks(self.bankIds)}
    else:
      self.xmlDetectors = {id: self.getXmlForBanks(id) for id in self.bankIds}
    self.tempDetectorIdfFiles = {id: self.getTempDetectorIdfFile(id) for id in self.xmlDetectors}
    if not idfCreation: #no need for long workspace creation if only the idf files are required
      self.detectors = {id: api.createDetectorWorkspace(self.tempDetectorIdfFiles[id].name, id) for id in self.tempDetectorIdfFiles}

  def getDetectorWorkspaceKey(self, pixelId):
    #note: conversion(e.g. offset) is already applied to pixelId (also referred to as detectorId)
    if self.singleNexus:
      return 'All'
    else:
      for id, (firstPixel, lastPixel) in enumerate(zip(self.firstPixelOfBank,self.lastPixelOfBank)):
        if firstPixel <= pixelId <= lastPixel:
          return str(id)
    return 'outOfValidPixelRanges'

  def getTempDetectorIdfFile(self, id):
    self.tempFile = NamedTemporaryFile(prefix=f'LOKI_Definition_detector_bank{id}', suffix='.xml')
    with open(self.tempFile.name, 'w') as f:
      f.write(self.xmlDetectors[id]) #TODO we might want to create the idf file here?
    return self.tempFile

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

  def getXmlForBanks(self, bankIds):
    commentLineToRemove = lambda line : any([(f'<!--BANK{id}_START' in line) or (f'BANK{id}_END-->' in line) for id in bankIds])
    newXmlLines = [line for line in self.xmlDetectorsTemplate.splitlines() if not commentLineToRemove(line)]
    return '\n'.join(newXmlLines)

  def setDistances(self):
    ## Monitor IDF ##
    self.xmlMonitor = self.xmlMonitor.replace('<PLACEHOLDER_NOMINAL_SOURCE_SAMPLE_DISTANCE>', "{:.6f}".format(self.nominalSourceSampleDistance))
    self.xmlMonitor = self.xmlMonitor.replace('<PLACEHOLDER_PRESAMPLE_MONITOR_DISTANCE>', "{:.6f}".format(self.preSampleMonitorDistance))
    beamstopMonitorDistance = self.nominalSourceSampleDistance + self.rearDetDistance - 0.05 #hardcoded -5cm
    self.xmlMonitor = self.xmlMonitor.replace('<PLACEHOLDER_BEAMSTOP_MONITOR_DISTANCE>', "{:.6f}".format(beamstopMonitorDistance))

    ## Detector IDF ##
    self.xmlDetectorsTemplate = self.xmlDetectorsTemplate.replace('<PLACEHOLDER_NOMINAL_SOURCE_SAMPLE_DISTANCE>', "{:.6f}".format(self.nominalSourceSampleDistance))

    rearDetFrontToPackCentreOffset = 0.048877 #hardcoded offset to position the detector front
    rearDetPositionZ = self.nominalSourceSampleDistance + self.rearDetDistance + rearDetFrontToPackCentreOffset
    self.xmlDetectorsTemplate = self.xmlDetectorsTemplate.replace('<PLACEHOLDER_REAR_DETECTOR_POSITION_Z>', "{:.6f}".format(rearDetPositionZ))

  def setPixelNumbers(self):
    self.xmlDetectorsTemplate = self.xmlDetectorsTemplate.replace('<PLACEHOLDER_PIXEL_NUMBER>', str(self.pixelNumber))
    self.xmlDetectorsTemplate = self.xmlDetectorsTemplate.replace('<PLACEHOLDER_PIXEL_LENGTH_500MM>', "{:.9f}".format(0.5/self.pixelNumber))
    self.xmlDetectorsTemplate = self.xmlDetectorsTemplate.replace('<PLACEHOLDER_PIXEL_LENGTH_1000MM>', "{:.9f}".format(1.0/self.pixelNumber))
    self.xmlDetectorsTemplate = self.xmlDetectorsTemplate.replace('<PLACEHOLDER_PIXEL_LENGTH_1200MM>', "{:.9f}".format(1.2/self.pixelNumber))

    numberOfPacksInBank = [28, 8, 6, 8, 6, 14, 16, 10, 16]
    numberOfPixelsInBank = [numberOfPacksInBank[i] * 8 * 7 * self.pixelNumber for i in range(9)] #8 tube/pack, 7 straw/tube, pixelNuber pixel/straw

    idfPixelOffset = 1 #pixel ids start from 1 in the IDF file

    def getFirstPixelOfBank(bankId):
      if self.firstPixelOfBank[bankId] is None:
        if bankId == 0:
          self.firstPixelOfBank[0] = idfPixelOffset
        else:
          self.firstPixelOfBank[bankId] = getFirstPixelOfBank(bankId-1) + numberOfPixelsInBank[bankId-1]
      return self.firstPixelOfBank[bankId]

    def getLastPixelOfBank(bankId):
      if self.lastPixelOfBank[bankId] is None:
        if bankId == 0:
          self.lastPixelOfBank[0] = getFirstPixelOfBank(0) + numberOfPixelsInBank[0] - idfPixelOffset
        else:
          self.lastPixelOfBank[bankId] = getLastPixelOfBank(bankId-1) + numberOfPixelsInBank[bankId]
      return self.lastPixelOfBank[bankId]


    for bankId in range(9):
      self.xmlDetectorsTemplate = self.xmlDetectorsTemplate.replace(f'<PLACEHOLDER_BANK{bankId}_PIXEL_ID_START>', str(getFirstPixelOfBank(bankId)))
      self.xmlDetectorsTemplate = self.xmlDetectorsTemplate.replace(f'<PLACEHOLDER_BANK{bankId}_PIXEL_ID_END>', str(getLastPixelOfBank(bankId)))

  def saveIdfFiles(self, savename, printOnly=False): #note: saving the idf files is not needed for producing the nxs files
    from shutil import copyfile
    if not printOnly:
      filename = f'{savename}_monitor.xml'
      copyfile(self.tempMonitorIdfFile.name, f'{savename}_monitor.xml')
      print(f'    Saved monitor idf file: {filename}')
      for key in self.tempDetectorIdfFiles:
        filename = f'{savename}_bank{key}_detector.xml'
        copyfile(self.tempDetectorIdfFiles[key].name, filename)
        print(f'    Saved detector idf file: {filename}')

    else: #print output for testing
      f = open(self.tempMonitorIdfFile.name, "r")
      print(f.read())
      for key in self.tempDetectorIdfFiles:
        f = open(self.tempDetectorIdfFiles[key].name, "r")
        print(f.read())

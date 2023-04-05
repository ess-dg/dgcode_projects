
import mantid.simpleapi as api
from mantid.kernel import DateAndTime
import numpy as np
import sys
import MCPL
from datetime import datetime

def createMonitorWorkspace(monitorTOF, monitorIntensity, monitorError, instrumentDefinitionFile_monitor):
  mcStasMonitorWS = api.CreateWorkspace(OutputWorkspace="mcStasMonitorWS", DataX=np.array(monitorTOF), DataY=np.array(monitorIntensity), DataE=np.array(monitorError), NSpec=len(monitorTOF), UnitX='TOF', YUnitLabel='Counts')
  
  api.LoadInstrument(Workspace=mcStasMonitorWS, Filename=instrumentDefinitionFile_monitor, RewriteSpectraMap=True)
  return mcStasMonitorWS

def rebinMonitorWorkspace(workspace, wsBinParams, workspaceName):
  return api.Rebin(InputWorkspace=workspace, OutputWorkspace=workspaceName, Params=wsBinParams)

def rebinDetectorWorkspace(workspace, wsBinParams, workspaceName):
  tofmin = int(workspace.getTofMin())
  tofmax = int(workspace.getTofMax())
  #width = tofmax - tofmin
  #SortEvents(InputWorkspace=Geant4DataWS, SortBy='Pulse Time')
  #Geant4DataWS = Rebin(InputWorkspace=Geant4DataWS, OutputWorkspace=Geant4DataWS,
  #             Params=str(tofmin) + "," + str(width) + "," + str(tofmax), PreserveEvents=True)
  if (tofmin < wsBinParams[0]):
      print(f"Lowest detection event TOF ({tofmin}) is lower than the workspace TOF binning limit ({wsBinParams[0]}).", file=sys.stderr)
  if (tofmax > wsBinParams[2]):
      print(f"Highest detection event TOF ({tofmax}) is higher than the workspace TOF binning limit ({wsBinParams[2]}).", file=sys.stderr)

  return api.Rebin(InputWorkspace=workspace, OutputWorkspace=workspaceName, Params=wsBinParams, PreserveEvents=False)

def addProperties(workspace, params):
  run = workspace.run()
  run.addProperty("TimeUnit", "Micro Seconds", True) #required
  mcplMetadata = params.getMcplMetadata()
  for (k,v) in mcplMetadata.items():
    run.addProperty(k, v, True)
  
def saveNexus(workspace, filename):
  api.SaveNexus(workspace, filename)
  print(f'    Saved nexus file: {filename}')

def createDetectorWorkspace(instrumentDefinitionFile_detector):
  Geant4DataWS = api.LoadEmptyInstrument(instrumentDefinitionFile_detector, MakeEventWorkspace=True)
  return Geant4DataWS

def addMcplDetectionEventsToWorkspace(Geant4DataWS, mcplFile, idConverter=(lambda id:id), idFilter=(lambda _:True), verbose=False):
  
  numHists = Geant4DataWS.getNumberHistograms()
  detToIndex = {}
  for i in range(numHists):
      eventList = Geant4DataWS.getSpectrum(i)
      id = eventList.getDetectorIDs()
      detToIndex[id[0]] = i
      eventList.clear(False)

  pulsetime = datetime.now()
  dateTime = DateAndTime(pulsetime.isoformat(sep="T"))
  numHistograms = Geant4DataWS.getNumberHistograms()
  allEventList = [ Geant4DataWS.getSpectrum(idx) for idx in range(numHistograms)]

  readBlockLength = 100000000
  if verbose:
    print(f'    Loading detection events from {mcplFile}')
  
  tof = np.array([])
  detids = np.array([])
  with MCPL.MCPLFile(mcplFile, blocklength=readBlockLength) as myfile:
    for p in myfile.particle_blocks:
      detids = np.append(detids, p.userflags.astype(int))
      tof = np.append(tof, p.time * 1000.0)  # convert to microseconds
    countAddEventError = 0
    countFilteredOutEvents = 0 #e.g., Simulation is done for the full LOKI rear bank geometry instead of the 'reduced' geometry used for the experiment
    for time,detId in zip(tof, detids):
      try:
        if idFilter(detId):
          id = idConverter(detId) + 1 # +1 for 1-based spectrum indexing
          allEventList[detToIndex[id]].addEventQuickly(time, dateTime)
        else:
          countFilteredOutEvents += 1
      except:
        countAddEventError += 1
    if countAddEventError:
       print(f'    ERROR: Number of addEventQuickly errors in file {mcplFile} is: {countAddEventError}', file=sys.stderr)
  return countFilteredOutEvents


import mantid.simpleapi as api
from mantid.kernel import DateAndTime
import numpy as np
import sys
import mcpl
from datetime import datetime
from math import inf as infinity

colWarning = '\033[93m'
colEnd = '\033[0m'

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
      print(f"    {colWarning}WARNING: Lowest detection event TOF ({tofmin}) is lower than the workspace TOF binning limit ({wsBinParams[0]}).{colEnd}", file=sys.stderr)
  if (tofmax > wsBinParams[2]):
      print(f"    {colWarning}WARNING: Highest detection event TOF ({tofmax}) is higher than the workspace TOF binning limit ({wsBinParams[2]}).{colEnd}", file=sys.stderr)

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

def createDetectorWorkspace(instrumentDefinitionFile_detector, id=''):
  return api.LoadEmptyInstrument(instrumentDefinitionFile_detector, OutputWorkspace=f'Geant4DataWS_bank{id}', MakeEventWorkspace=True)

def addMcplDetectionEventsToWorkspaces(workspaces, filename, idConverter=(lambda id:id), idFilter=(lambda _:True)):
  pulsetime = datetime.now() #dummy pulse time
  dateTime = DateAndTime(pulsetime.isoformat(sep="T"))

  spectraNumber = {}
  allEventLists = {}
  detectorIdToWorkspaceIndex = {}
  for key,ws in workspaces.detectors.items():
    spectraNumber[key] = ws.getNumberHistograms()
    allEventLists[key] = [ws.getSpectrum(wsIndex) for wsIndex in range(spectraNumber[key])]
    detectorIdToWorkspaceIndex[key] = {ws.getDetector(wsIndex).getID():wsIndex for wsIndex in range(spectraNumber[key])}
  readBlockLength = 100000000
  tof = np.array([])
  detids = np.array([])
  print(f'    Loading detection events from {filename}')
  with mcpl.MCPLFile(filename, blocklength=readBlockLength) as myfile:
    oldFile = myfile.opt_userflags #This could break if new detection files are created with userflags for some reason
    for p in myfile.particle_blocks:
      detids = np.append(detids, (p.ekin.astype(int) if not oldFile else p.userflags.astype(int)))
      tof = np.append(tof, p.time * 1000.0)  # convert to microseconds
    countAddEventError = 0
    addEventErrorIdMin = infinity
    addEventErrorIdMax = -infinity
    countFilteredOutEvents = 0 #e.g., Simulation is done for the full LOKI rear bank geometry instead of the 'reduced' geometry used for the experiment
    for time,detId in zip(tof, detids):
      try:
        if idFilter(detId):
          id = idConverter(detId)
          wsKey = workspaces.getDetectorWorkspaceKey(id)
          allEventLists[wsKey][detectorIdToWorkspaceIndex[wsKey][id]].addEventQuickly(time, dateTime)
        else:
          countFilteredOutEvents += 1
      except:
        countAddEventError += 1
        if detId < addEventErrorIdMin:
          addEventErrorIdMin = detId
        if detId > addEventErrorIdMax:
          addEventErrorIdMax = detId
    if countAddEventError:
       print(f'    {colWarning}WARNING: Number of addEventQuickly errors in file {filename} is: {countAddEventError}. Possibly wrong IDF for the simulation (simulation pixel index range out of the pixel range defined in the IDF file). Minimum of ids causing error: {int(addEventErrorIdMin)}, maximum of ids causing error: {int(addEventErrorIdMax)} {colEnd}', file=sys.stderr)
  return countFilteredOutEvents

#Legacy version needed for Larmor2020 and Larmor2022 processing
def addMcplDetectionEventsToWorkspace(workspace, filename, idConverter=(lambda id:id), idFilter=(lambda _:True), verbose=False):
  pulsetime = datetime.now() #dummy pulse time
  dateTime = DateAndTime(pulsetime.isoformat(sep="T"))

  spectraNumber = workspace.getNumberHistograms()
  allEventLists = [workspace.getSpectrum(wsIndex) for wsIndex in range(spectraNumber)]
  detectorIdToWorkspaceIndex = {workspace.getDetector(wsIndex).getID():wsIndex for wsIndex in range(spectraNumber)}

  readBlockLength = 100000000
  tof = np.array([])
  detids = np.array([])
  if verbose:
    print(f'    Loading detection events from {filename}')
  with mcpl.MCPLFile(filename, blocklength=readBlockLength) as myfile:
    oldFile = myfile.opt_userflags #This could break if new detection files are created with userflags for some reason
    for p in myfile.particle_blocks:
      detids = np.append(detids, (p.ekin.astype(int) if not oldFile else p.userflags.astype(int)))
      tof = np.append(tof, p.time * 1000.0)  # convert to microseconds
    countAddEventError = 0
    countFilteredOutEvents = 0 #e.g., Simulation is done for the full LOKI rear bank geometry instead of the 'reduced' geometry used for the experiment
    for time,detId in zip(tof, detids):
      try:
        if idFilter(detId):
          id = idConverter(detId)
          allEventLists[detectorIdToWorkspaceIndex[id]].addEventQuickly(time, dateTime)
        else:
          countFilteredOutEvents += 1
      except:
        countAddEventError += 1
    if countAddEventError:
       print(f'    {colWarning}WARNING: Number of addEventQuickly errors in file {filename} is: {countAddEventError}. Possibly wrong IDF for the simulation (simulation pixel index range out of the pixel range defined in the IDF file){colEnd}', file=sys.stderr)
  return countFilteredOutEvents


def copyMonitorSpectraToDetectorWorkspace():
  gw = api.mtd['Geant4DataWS2D']
  mw = api.mtd['mcStasMonitorWS_rebin']
  for i in range(3):
    gw.setY(i, mw.readY(i))

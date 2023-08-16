
from mantid.kernel import *
from mantid.simpleapi import *
from mantid.api import *

class LoadLokiDetectionEvents(PythonAlgorithm):
  def PyInit(self):
    self.declareProperty(FileProperty(name="filename", defaultValue="", action=FileAction.Load, extensions=["mcpl", "mcpl.gz"]),
                                      doc="Input detection MCPL file.")
    self.declareProperty('savename', defaultValue="",
                         doc="Output nexus filename base (must be an absolute path), automatically extended by '_detector.nxs' and '_monitor.nxs'") #TODO the . might cause an isse
    self.declareProperty('rear_detector_distance_m', defaultValue='',
                         doc="Rear detector distance.")
    self.declareProperty('analysis_straw_pixel_number', defaultValue='',
                         doc="Number of pixels per detector straw.")
    self.declareProperty('detectorWorkspaceTofStart', defaultValue=11000,
                         doc="TOF start for the detector workspace spectra binning. (default: %(default)d microSeconds)")
    self.declareProperty('detectorWorkspaceTofBinWidth', defaultValue=1000,
                         doc="TOF bin width for the detector workspace spectra binning. (default: %(default)d microSeconds)")
    self.declareProperty('detectorWorkspaceTofEnd', defaultValue=113000,
                         doc="TOF end for the detector workspace spectra binning. (default: %(default)d microSeconds)")
    self.declareProperty('monitorWorkspaceTofStart', defaultValue=9000,
                         doc="TOF start for the monitor workspace spectra binning. (default: %(default)d microSeconds)")
    self.declareProperty('monitorWorkspaceTofBinWidth', defaultValue=1000,
                         doc="TOF bin width for the monitor workspace spectra binning. (default: %(default)d microSeconds)")
    self.declareProperty('monitorWorkspaceTofEnd', defaultValue=100000,
                         doc="TOF end for the monitor workspace spectra binning. (default: %(default)d microSeconds)")
    self.declareProperty('detectorOnly', defaultValue=False,
                         doc="Create only the detector spectrum file.")
    self.declareProperty('monitorOnly', defaultValue=False,
                         doc="Create only the monitor spectrum file.")
    self.declareProperty('noOutput', defaultValue=False,
                         doc="Developer option, no output is generated.")
    self.declareProperty('verbose', defaultValue=False,
                         doc="Enable more verbosity.")
    self.declareProperty('showMetadata', defaultValue=False,
                         doc="Only output metadata stored in the detection MCPL file (without producing nxs files).")
    self.declareProperty('idfCreation', defaultValue=False,
                         doc="Save the idf files created based on inputs.")
    self.declareProperty('singleNexus', defaultValue=False,
                         doc="Save a single nexus file (or idf) for all detector banks together, instead of one for each.")
    self.declareProperty('mcstasDir', defaultValue='',
                         doc="Name of the McStas directory containing the monitor spectra.")
    self.declareProperty('mcstas_monitors', defaultValue='',
                         doc="Name of the McStas monitor files (this option overrides the default ones)") #TODO
    self.declareProperty('neutronNumber', defaultValue='',
                         doc="The total number of neutrons used to create the (probably merged) detection MCPL input file.")
    self.declareProperty('neutron_wavelength_min_aangstrom', defaultValue='',
                         doc="Minimum neutron wavelength used for the flood source sampling.")
    self.declareProperty('neutron_wavelength_max_aangstrom', defaultValue='',
                         doc="Maximum neutron wavelength used for the flood source sampling.")
    self.declareProperty('sampling_cone_opening_deg', defaultValue='',
                         doc="Maximum opening angle of the cone used for the direction sampling.")
    self.declareProperty('sampling_cone_opening_min_deg', defaultValue='',
                         doc="Minimum opening angle of the cone used for the direction sampling.")
    self.declareProperty('source_monitor_distance_meters', defaultValue='',
                         doc="Source to preSample monitor distance, the preSample monitor TOF spectrum is generated for this distance")
    self.declareProperty('aiming_bank_id', defaultValue='',
                         doc="Id(s) of the banks to include.")
    self.declareProperty('nominal_source_sample_distance_meters', defaultValue='',
                         doc="Nominal source(moderator) to sample position distance [m].")

    # self.declareProperty(EventWorkspaceProperty("OutputWorkspace", self.OUTPUT_NAME, Direction.Output),
    #                       "Output event workspace")

  def category(self):
      return 'ESS'

  def PyExec(self):
    self.load()

  def _getParser(self):
    import argparse
    parser = argparse.ArgumentParser(description = 'Process detection MCPL files from loki simulations to produce nexus files (one for the detectors, one for the monitors) readable with Mantid.')
    parser.add_argument('filename', help = 'Input detection MCPL file. Assumed to be a path relative to the directory defined by the G4PROC_MCPL_BASEDIR_LOKI environment variable, unless a full path is given.')
    parser.add_argument('-s','--savename', help = 'Output nexus filename base. Automatically extended by "_detector.nxs" and "_monitor.nxs". Assumed to be a path relative to the directory defined by the G4PROC_SAVEDIR_LOKI environment variable, unless a full path is given.')

    parser.add_argument('--rear_detector_distance_m', help = 'Rear detector distance.')
    parser.add_argument('--analysis_straw_pixel_number', help = 'Number of pixels per detector straw.')
    parser.add_argument('--detectorWorkspaceTofStart', default = 11000, type=int, help = 'TOF start for the detector workspace spectra binning. (default: %(default)d microSeconds)')
    parser.add_argument('--detectorWorkspaceTofBinWidth', default = 1000, type=int, help = 'TOF bin width for the detector workspace spectra binning. (default: %(default)d microSeconds)')
    parser.add_argument('--detectorWorkspaceTofEnd', default = 113000, type=int, help = 'TOF end for the detector workspace spectra binning. (default: %(default)d microSeconds)')
    parser.add_argument('--monitorWorkspaceTofStart', default = 9000, type=int, help = 'TOF start for the monitor workspace spectra binning. (default: %(default)d microSeconds)')
    parser.add_argument('--monitorWorkspaceTofBinWidth', default = 1000, type=int, help = 'TOF bin width for the monitor workspace spectra binning. (default: %(default)d microSeconds)')
    parser.add_argument('--monitorWorkspaceTofEnd', default = 100000, type=int, help = 'TOF end for the monitor workspace spectra binning. (default: %(default)d microSeconds)')
    parser.add_argument('-det', '--detectorOnly', action = 'store_true', help = 'Create only the detector spectrum file.')
    parser.add_argument('-mon', '--monitorOnly', action = 'store_true', help = 'Create only the monitor spectrum file.')
    parser.add_argument('-no', '--noOutput', action = 'store_true', help = 'Developer option. No output is generated.')
    parser.add_argument('-v', '--verbose', action = 'store_true', help = 'Enable more verbosity.')
    parser.add_argument('-show', '--showMetadata', action = 'store_true', help = 'Only output metadata stored in the detection MCPL file (without producing nxs files).')
    parser.add_argument('-i', '--idfCreation', action = 'store_true', help = 'Save the idf files created based on inputs.')
    parser.add_argument('--singleNexus', action = 'store_true', default=False, help = 'Save a single nexus file (or idf) for all detector banks together, instead of one for each.')

    mcstasParamGroup = parser.add_argument_group('McStas parameters', 'Supplement McStas simulation parameters')
    mcstasParamGroup.add_argument('--mcstasDir', help = 'Name of the McStas directory containing the monitor spectra. Assumed to be a path relative to the directory defined by the G4PROC_MCSTAS_BASEDIR_LOKI environment variable, unless a full path is given.')
    mcstasParamGroup.add_argument('--mcstas_monitors', nargs='*', help = 'Name of the McStas monitor files (this option overrides the default ones).')

    floodPramGroup = parser.add_argument_group('Flood source parameters', 'Supplement flood source simulation parameters missing from the detection MCPL files (legacy files), or override parameters stored in the files.')
    floodPramGroup.add_argument('-n', '--neutronNumber', type=float, help = 'The total number of neutrons used to create the (probably merged) detection MCPL input file.')
    floodPramGroup.add_argument('--neutron_wavelength_min_aangstrom', help = 'Minimum neutron wavelength used for the flood source sampling.')
    floodPramGroup.add_argument('--neutron_wavelength_max_aangstrom', help = 'Maximum neutron wavelength used for the flood source sampling.')
    floodPramGroup.add_argument('--sampling_cone_opening_deg', help = 'Maximum opening angle of the cone used for the direction sampling.')
    floodPramGroup.add_argument('--sampling_cone_opening_min_deg', help = 'Minimum opening angle of the cone used for the direction sampling.')
    floodPramGroup.add_argument('--source_monitor_distance_meters', help = 'Source to preSample monitor distance. The preSample monitor TOF spectrum is generated for this distance.')
    floodPramGroup.add_argument('--aiming_bank_id', help = 'Id(s) of the banks to include.')
    floodPramGroup.add_argument('--nominal_source_sample_distance_meters', help = 'Nominal source(moderator) to sample position distance. [m]')
    return parser

  def _getInputsForArgParser(self):
    """Transforms the input of the Mantid Agorithm GUI into a format suitable for the command line argument parser"""
    inps = [self.getProperty("filename").value] #filename is the only positional argument
    propPairs = [[f"--{p.name}", str(p.value)] for p in self.getProperties() if p.name!='filename' and p.value!='' and p.value is not False]
    inps.extend([pair[0] for pair in propPairs if pair[1]=='True']) #only the name for "store_true" arguments
    inps.extend([item for pair in propPairs if pair[1]!='True' for item in pair]) # flat list of name,value,name,value
    # print(f"{inps=}")
    return inps

  def _resolveFilename(self, filename, baseEnvVar):
    """Ensure that the filename is an absolute path"""
    import os
    from pathlib import Path
    base = os.environ.get(baseEnvVar, None)
    if os.path.isabs(filename):
      return Path(filename)
    elif base:
      return Path(base) / filename
    else:
      raise Exception(f"    The provided filename ({filename}) is not an absolute path, and the {baseEnvVar} env var is not set.")

  def load(self, cli=False):
    parser = self._getParser()
    if cli:
      print("    Assuming dgcode script execution")
      args = parser.parse_args()
    else: #if not called from the loki dgcode script, get input arguments from Mantid algorithm interface
      print("    Assuming Mantid GUI execution")
      args = parser.parse_args(self._getInputsForArgParser())

      import dgbuild.cfg as cfg
      import os
      os.environ["ESS_DATA_DIR"]=str(cfg.dirs.datadir) #This is set by the bootstrap.sh script for dgcode cli execution
    mcplFile = self._resolveFilename(args.filename, 'G4PROC_MCPL_BASEDIR_LOKI')
    if not (args.noOutput or args.showMetadata):
      saveFileBase = self._resolveFilename(args.savename, 'G4PROC_SAVEDIR_LOKI')
    if args.mcstasDir:
      mcStasFolder = self._resolveFilename(args.mcstasDir, 'G4PROC_MCSTAS_BASEDIR_LOKI')
    #TODO check if files exists and stop with error, or overide?

    import LokiMantid.mantidApi as api
    from LokiMantid.monitorSpectrum import createFloodSourceTofSpectrum, extractMcStasMonitorData
    import numpy as np
    from glob import glob
    import math as m

    isFloodSourceSimulation = False
    if not args.savename and not (args.noOutput or args.showMetadata):
      parser.error("Save filename (-s, --savename) is required, unless no output (-no, --noOutput) or metadata (-show, --showMetadata) option is selected")
      raise Exception #parser error doesn't stop Mantid algorithm execution
    if not args.showMetadata and (not args.neutronNumber and args.mcstasDir is None):
      parser.error("Either --neutronNumber (for Flood source simulation) or --mcstasDir (for McStas+Geant4 simulation) is required, unless the --showMetadata option is used to examine the metadata stored in the detection MCPL file!")
      raise Exception #parser error doesn't stop Mantid algorithm execution
    if args.neutronNumber:
      if args.mcstasDir:
        parser.error("Both --neutronNumber (for Flood source simulation) and --mcstasDir (for McStas+Geant4 simulation) options cannot be provided!")
        raise Exception #parser error doesn't stop Mantid algorithm execution
      if not args.neutronNumber.is_integer():
        parser.error(" --neutronNumber has to be an integer.") #float input type is only allowed to handle '1e9' format
        raise Exception #parser error doesn't stop Mantid algorithm execution
      else:
        simulatedNeutronNumber = int(args.neutronNumber)
        isFloodSourceSimulation = True

    lokiDefaultParams = {'rear_detector_distance_m':5,
                        'neutron_wavelength_min_aangstrom': 2,
                        'neutron_wavelength_max_aangstrom': 13,
                        'sampling_cone_opening_deg': 49.6,
                        'sampling_cone_opening_min_deg': 0,
                        'source_monitor_distance_meters': 23.599,
                        'nominal_source_sample_distance_meters': 23.6,
                        'analysis_straw_pixel_number': 512,
                        'mcstas_monitors': ['Mon10_PostFOC_*.t', 'PreSampleMonitor_*.t', 'PostSampleMonitor_*.t', '-'],
                        'aiming_bank_id': '876543210'}
    from LokiMantid.paramHandler import ParamHandler
    params = ParamHandler(mcplFile, vars(args), lokiDefaultParams, args.verbose)

    if not args.mcstas_monitors:
    #Possile TODO maybe implement a set option, to indicate that parameters can change! (then it would make sense to prepend with '_')
      params.params['mcstas_monitors'][3] = f"beamstopMonitor_{params.get('rear_detector_distance_m')}m_*.t"
    params.params['aiming_bank_id'] = list(str(params.get('aiming_bank_id')))
    FilteredOurBankIds = list(set(str(lokiDefaultParams['aiming_bank_id']))-set(params.params['aiming_bank_id']))

    if args.showMetadata:
      params.dumpMCPLParams()
      params.dumpParams()
      return

    from LokiMantid.workspaceCreator import workspaceCreator
    workspaces = workspaceCreator(params, args.singleNexus, args.idfCreation)
    instrumentDefinitionFile_monitor = workspaces.getMonitorIDF()
    if args.idfCreation:
      if not args.noOutput:
        workspaces.saveIdfFiles(f'{saveFileBase.parents[0]}/LOKI_Definition')
      else:
        workspaces.saveIdfFiles('_', printOnly=True)
      return
    #TODO old tube numbering is not supported

    # parameters used for the TOF spectra binning for both the monitors and detectors
    wsBinParams_monitor = [args.monitorWorkspaceTofStart, args.monitorWorkspaceTofBinWidth, args.monitorWorkspaceTofEnd]
    wsBinParams_detector = [args.detectorWorkspaceTofStart, args.detectorWorkspaceTofBinWidth, args.detectorWorkspaceTofEnd]

    if not args.detectorOnly:
      monitorTOF = []
      monitorIntensity = []
      monitorError = []
      if isFloodSourceSimulation:
        print("    Processing flood source simulation data.")
        if simulatedNeutronNumber <= params.mcplEventNr:
          raise Exception(f"    The provided neutronNumber({simulatedNeutronNumber}) is lower than the number of detection events({params.mcplEventNr}) in the MCPL file!")

        directionBiasMultiplicationFactor = 2/(m.cos(params.get('sampling_cone_opening_min_deg')*m.pi/180) - m.cos(params.get('sampling_cone_opening_deg')*m.pi/180))
        incidentNeutronNumber = simulatedNeutronNumber * directionBiasMultiplicationFactor
        if args.verbose:
          print("     Handling flood source neutron number:")
          print(f"        Number of simulated neutrons: {simulatedNeutronNumber}")
          print(f"        Multiplication factor due to biassed direction sampling: {directionBiasMultiplicationFactor}")
          print(f"        Resulting incident neutron number on the monitor: {incidentNeutronNumber}")

        floodTof, floodIntensity, floodError = createFloodSourceTofSpectrum(incidentNeutronNumber, params.get('neutron_wavelength_min_aangstrom'), params.get('neutron_wavelength_max_aangstrom'), params.get('source_monitor_distance_meters'), verbose=args.verbose)

        #PreSample monitor is expected to be the 3rd spectrum out of 4, so empty spectra are added here
        blank = np.zeros_like(floodTof)
        monitorTOF = [floodTof, floodTof, floodTof, floodTof]
        monitorIntensity = [blank, blank, floodIntensity, blank]
        monitorError = [blank, blank, floodError, blank]
      else: #McStas+Geant4 simulation
        print("    Processing McStas+Geant4 simulation data.")
        for i, filename in enumerate(params.get('mcstas_monitors')):
          if len(glob(str(mcStasFolder / filename))) == 1:
            monitorFilePath = glob(str(mcStasFolder / filename))[0]
          else:
            raise Exception(f"    Unable to finding {mcStasFolder}/{filename}")
          tof, wavelength, y, e = extractMcStasMonitorData(monitorFilePath, args.verbose)
          monitorTOF.append(wavelength if tof is None else tof)
          monitorIntensity.append(y)
          monitorError.append(e)

      mcStasMonitorWS = api.createMonitorWorkspace(monitorTOF, monitorIntensity, monitorError, instrumentDefinitionFile_monitor)
      mcStasMonitorWS_rebin = api.rebinMonitorWorkspace(mcStasMonitorWS, wsBinParams_monitor, 'mcStasMonitorWS_rebin')
      api.addProperties(mcStasMonitorWS_rebin, params)
      if not args.noOutput:
        saveFilename = str(saveFileBase) + "_monitor.nxs"
        api.saveNexus(mcStasMonitorWS_rebin, saveFilename)

    if not args.monitorOnly:
      idConverter = lambda id: 1 + id #detector IDs in the IDF (and ICD) file starts from 1 (as opposed to the zero-based numbering in the Geant4 geometry)
      # idFilter = lambda id: any([workspaces.firstPixelOfBank[int(bankId)]<=id<=workspaces.lastPixelOfBank[int(bankId)] for bankId in workspaces.bankIds])
      idFilter = lambda id: not any([workspaces.firstPixelOfBank[int(bankId)]<=id<=workspaces.lastPixelOfBank[int(bankId)] for bankId in FilteredOurBankIds])
      eventsFilteredOut = api.addMcplDetectionEventsToWorkspaces(workspaces, params.getMcplFile(), idConverter=idConverter, idFilter=idFilter)

      ## Print statistics ##
      print(f'    Number of events in the MCPL file: {params.mcplEventNr}')
      sumWsEventNr = 0
      for id, ws in workspaces.detectors.items():
        wsEventNr = ws.getNumberEvents()
        sumWsEventNr+=wsEventNr
        print(f"    Number of events in bank {id}: {wsEventNr}")
      if eventsFilteredOut != 0:
        print(f'    Number of events filtered out: {eventsFilteredOut}')
      if params.mcplEventNr != sumWsEventNr+eventsFilteredOut:
        print(f'    Lost events: {params.mcplEventNr-(sumWsEventNr+eventsFilteredOut)}')
        print(f'      Are you sure about the analysis_straw_pixel_number({params.get("analysis_straw_pixel_number")})?')
        print(f'      Bank pixel id ranges:')
        for bankId in range(9):
          filteredText = '(filtered out)' if str(bankId) in FilteredOurBankIds else ''
          print(f'        BankId: {bankId} start: {workspaces.firstPixelOfBank[bankId]} end: {workspaces.lastPixelOfBank[bankId]} {filteredText}')

      for id, ws in workspaces.detectors.items():
        Geant4DataWS_rebin = api.rebinDetectorWorkspace(ws, wsBinParams_detector, f'Geant4DataWS2D_bank{id}_rebin')
        api.addProperties(Geant4DataWS_rebin, params)
        Geant4DataWS_rebin.setYUnitLabel("Counts")
        Geant4DataWS_rebin.getAxis(0).setUnit("TOF")
        if not args.noOutput:
          aimingBankIds = ''.join(params.params['aiming_bank_id'])
          aimingBankIdsText = f"_banks{aimingBankIds}" if aimingBankIds != lokiDefaultParams['aiming_bank_id'] else ''
          bankText = f"_bank{id}" if not args.singleNexus else aimingBankIdsText
          saveFilename = str(saveFileBase) + f"{bankText}_detector.nxs"
          api.saveNexus(Geant4DataWS_rebin, saveFilename)

# if __name__ == '__main__':
#   AlgorithmFactory.subscribe(LoadLokiDetectionEvents)
AlgorithmFactory.subscribe(LoadLokiDetectionEvents)
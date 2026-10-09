import G4GeoLoki.LokiAimHelper as LokiAim
from Units import units
import numpy as np

def launch(geo):
    import G4Launcher
    launcher = G4Launcher()

    launcher.addParameterString('event_gen','')
    launcher.addParameterInt("analysis_straw_pixel_number", 256)
    launcher.addParameterBoolean('det_only',False)
    launcher.addParameterString('aiming_bank_id','') #aim only at a certain bank of group of banks(rear detector('rear'), mid-detector('mid), front detector('front))
    launcher.addParameterDouble("nominal_source_sample_distance_meters", 23.6) #origo of the Geant4 geometry
    launcher.addParameterDouble("gen_x_offset_meters", 0.0)

    ## McStas+Geant4 options ##
    launcher.addParameterString('input_file', '')

    ## Visualisation only options ##
    launcher.addParameterBoolean('primary_only',False)
    launcher.addParameterString('cone_view', 'no') # 'min' or 'max' Only for visual confirmation of the cone_opening_min_deg/cone_opening_deg angle for floodSource simulation
    launcher.addParameterBoolean('geantino', False) # use geantinos instead of neutrons (for visualisation)
    ## Unused options ##
    # launcher.addParameterBoolean('gravity',False)
    # launcher.addParameterBoolean('addgeantinoes',False)

    #geometry:
    launcher.setGeo(geo)

    #generator:
    if launcher.getParameterString('event_gen')=='mcpl': #McStas+Geant4 simulation
        import G4MCPLPlugins.MCPLGen as Gen
        import mcpl
        gen = Gen.create()
        if(launcher.getParameterString('input_file')):
          gen.input_file = launcher.getParameterString('input_file')
          tmp_myfile = mcpl.MCPLFile(gen.getParameterString('input_file'))
          for blobkey in tmp_myfile.blobs:
            launcher.setUserData(blobkey, str(tmp_myfile.blobs[blobkey].decode("utf-8")))

          gen.dz_meter = 0.2 #default nominal sample position to mcpl_output mcstas component distance. New mcpl files store this value with the 'sample_mcpl_distance_m' key, the default is needed for older mcpl files
          if('sample_mcpl_distance_m' in tmp_myfile.blobs):
            gen.dz_meter = float(tmp_myfile.blobs['sample_mcpl_distance_m'])

        gen.dx_meter = launcher.getParameterDouble('gen_x_offset_meters')
    elif launcher.getParameterString('event_gen')=='flood': #Flood source simulation
        from  LOKI.FloodSourceGen import FloodSourceGen as Gen
        ssd = launcher.getParameterDouble("nominal_source_sample_distance_meters")
        gen = Gen(ssd)
        gen.gen_x_offset_meters = launcher.getParameterDouble('gen_x_offset_meters')
        bankFilter = launcher.getParameterString('aiming_bank_id')
        if bankFilter != '':
          assert bankFilter in (*[str(i) for i in range(9)], '1234', '5678'), f"aiming_bank_id must be either [0,8], 1234 or 5678"
          if bankFilter == '1234': #mid banks
            angleRange = (3.5, 15.8)
          elif bankFilter == '5678': #front banks
            angleRange = (10.5, 49.6)
          else: #single banks [0,8]
            bankId = int(bankFilter)
            def aimAtBank(): #when the geometry parameters are final: the bank placement of the geometry in use
              geometry = launcher.getGeo()
              aimHelper = LokiAim.AimHelper(geometry.getParameterDouble('rear_detector_distance_m')*units.m,
                                            LokiAim.DEFAULT_NUMBER_OF_PIXELS_IN_STRAW, 9,
                                            geometry.getParameterString('bank_calibration'))
              bankCentre = aimHelper.getBankTransform(bankId)[1] #the centre of the bank volume
              launcher.getGen().ref_dir_x, launcher.getGen().ref_dir_y, launcher.getGen().ref_dir_z = np.array(bankCentre)/np.linalg.norm(bankCentre)
            launcher.addPrePreInitHook(aimAtBank)
            bankConeAngle = [5.07, 9.9, 4.9, 9.9, 4.9, 31.9, 20.9, 29.5, 22.1] #HARDCODED for now (nominal geometry)
            angleRange = (0, bankConeAngle[bankId])
          if launcher.getParameterString('cone_view')=='min': #only for visualisation!
            angleRange = (angleRange[0], 1.0001*angleRange[0])
          elif launcher.getParameterString('cone_view')=='max': #only for visualisation!
            angleRange = (0.9999*angleRange[1], angleRange[1])
          gen.cone_opening_min_deg, gen.cone_opening_deg = angleRange
        if launcher.getParameterBoolean('geantino') is True: #only for visualisation!
           gen.particle = 'geantino'
    elif launcher.getParameterString('event_gen')=='spheremodel':
        from  LOKI.SansSphereGen import SansSphereGen as Gen
        gen = Gen()
    elif launcher.getParameterString('event_gen')=='isotheta':
        from  LOKI.IsoThetaGen import IsoThetaGen as Gen
        gen = Gen()
    else:
        import G4StdGenerators.FlexGen as Gen
        gen = Gen.create()
        gen.particleName = 'geantino' if launcher.getParameterBoolean('geantino') else 'neutron'
        gen.neutron_wavelength_aangstrom = 3.0
        gen.momdir_spherical = True
        gen.randomize_polarangle = True
        gen.random_min_polarangle_deg = 0 #49.99
        gen.random_max_polarangle_deg =0.001  #50 #16.0
        gen.randomize_azimuthalangle = True
        gen.random_min_azimuthalangle_deg = 0 #20.0
        gen.random_max_azimuthalangle_deg = 360.0 #60

    launcher.setGen(gen)

    def addUserData():
      launcher.setUserData("analysis_straw_pixel_number", str(launcher.getParameterInt('analysis_straw_pixel_number')))
      launcher.setUserData("rear_detector_distance_m", str(launcher.getGeo().getParameterDouble("rear_detector_distance_m")))
      # the bank placements (also written to the detection files of the analysis, checked by the LokiMantid scripts)
      bankCalibration = launcher.getGeo().getParameterString("bank_calibration")
      launcher.setUserData("bank_calibration", bankCalibration)
      # the text of the calibration: the analysis uses it, so that it has the placements of the simulation even if the
      # calibration file changes or is not there
      launcher.setUserData("bank_calibration_text", LokiAim.bankCalibrationText(bankCalibration))
      if bankCalibration != LokiAim.NOMINAL_BANK_CALIBRATION:
        print(f"NOTE: bank_calibration={bankCalibration}: the LokiMantid scripts (Mantid instrument definition) assume "
              f"the nominal geometry; use bank_calibration={LokiAim.NOMINAL_BANK_CALIBRATION} for simulations to be "
              f"processed with Mantid.")
      launcher.setUserData("aiming_bank_id", str(launcher.getParameterString('aiming_bank_id')))
      launcher.setUserData("nominal_source_sample_distance_meters", str(launcher.getParameterDouble('nominal_source_sample_distance_meters')))
      if(launcher.getGen().getName()=="LOKI.FloodSourceGen/FloodSourceGen"): #event_gen=flood
        launcher.setUserData("source_monitor_distance_meters", str(launcher.getGen().source_monitor_distance_meters))
        launcher.setUserData("sampling_cone_opening_deg", str(launcher.getGen().cone_opening_deg))
        launcher.setUserData("sampling_cone_opening_min_deg", str(launcher.getGen().cone_opening_min_deg))
        launcher.setUserData("neutron_wavelength_min_aangstrom", str(launcher.getGen().neutron_wavelength_min_aangstrom))
        launcher.setUserData("neutron_wavelength_max_aangstrom", str(launcher.getGen().neutron_wavelength_max_aangstrom))

    launcher.addPrePreInitHook(addUserData) #add userdata when all parameters are available

    #filter:
    if launcher.getParameterBoolean('primary_only'):
        import G4CollectFilters.StepFilterPrimary as F
        launcher.setFilter(F.create())

    # #gravity
    # if launcher.getParameterBoolean('gravity'):
    #     import G4GravityHelper.NeutronGravity as ng
    #     ng.enableNeutronGravity(launcher)

    # if launcher.getParameterBoolean('addgeantinoes'):
    #     import G4GeantinoInserter
    #     G4GeantinoInserter.install()

    if not launcher.getParameterBoolean('det_only'):
        launcher.setOutput('lokisim','REDUCED')
    else:
        griff_output_volumes = ["CountingGas","Converter"]
        import G4CollectFilters.StepFilterVolume
        f = G4CollectFilters.StepFilterVolume.create()
        f.volumeList = griff_output_volumes
        launcher.setFilter(f)
        launcher.setOutput('lokibcs_CountingGas_Converter','REDUCED')

    #launch:
    launcher.go()

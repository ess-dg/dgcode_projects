# Simulation of the Larmor 2022 experiment (the LOKI rear bank at the Larmor instrument, ISIS), with the
# Larmor/GeoLarmor2022 geometry. (Extracted from LOKI.Launcher, where it was enabled by the larmor_2022_experiment
# geometry parameter.)

def launch(geo):
    import G4Launcher
    launcher = G4Launcher()

    launcher.addParameterString('event_gen','')
    launcher.addParameterInt("analysis_straw_pixel_number", 512) # must be 512 (the experiment)
    launcher.addParameterBoolean('det_only',False)
    launcher.addParameterDouble("nominal_source_sample_distance_meters", 25.61) #origo of the Geant4 geometry
    launcher.addParameterDouble("gen_x_offset_meters", 0.005)

    ## McStas+Geant4 options ##
    launcher.addParameterString('input_file', '')

    ## Visualisation only options ##
    launcher.addParameterBoolean('primary_only',False)
    launcher.addParameterBoolean('geantino', False) # use geantinos instead of neutrons (for visualisation)

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

        gen.dx_meter = launcher.getParameterDouble('gen_x_offset_meters')
    elif launcher.getParameterString('event_gen')=='flood': #Flood source simulation
        from  LOKI.FloodSourceGen import FloodSourceGen as Gen
        ssd = launcher.getParameterDouble("nominal_source_sample_distance_meters")
        gen = Gen(ssd)
        gen.gen_x_offset_meters = launcher.getParameterDouble('gen_x_offset_meters')
        if launcher.getParameterBoolean('geantino') is True: #only for visualisation!
           gen.particle = 'geantino'
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

    def setParamsForLarmor2022Experiment(): #note: prone to generator name change
        assert launcher.getParameterInt('analysis_straw_pixel_number') == 512, "analysis_straw_pixel_number must be 512 for the Larmor2022 experiment!"
        if(launcher.getGen().getName()=="G4MCPLPlugins/MCPLGen"): #event_gen=mcpl
          assert launcher.getGen().dx_meter == 0.005, "gen_x_offset_meters should be 0.005 for the Larmor 2022 experiment!"
          launcher.getGen().dz_meter = 4.049 #note: intentionally 4.049, not 4.099
        elif(launcher.getGen().getName()=="LOKI.FloodSourceGen/FloodSourceGen"): #event_gen=flood
          assert launcher.getGen().gen_x_offset_meters == 0.005, "gen_x_offset_meters should be 0.005 for the Larmor2022 experiment!"
          assert launcher.getParameterDouble("nominal_source_sample_distance_meters") == 25.61, "nominal_source_sample_distance_meters should be 25.61 for the Larmor2022 experiment!"
          launcher.getGen().source_monitor_distance_meters = 25.57
          import math as m
          launcher.getGen().cone_opening_deg = m.acos(1-2/233)/m.pi*180
          print(f"Using predifined parameters for the Larmor2022 experiment!")
          print(f'    source_monitor_distance_meters: {launcher.getGen().source_monitor_distance_meters}')
          print(f'    cone_opening_deg: {launcher.getGen().cone_opening_deg}')

    def addUserData():
      launcher.setUserData("analysis_straw_pixel_number", str(launcher.getParameterInt('analysis_straw_pixel_number')))
      launcher.setUserData("rear_detector_distance_m", str(launcher.getGeo().getParameterDouble("rear_detector_distance_m")))
      launcher.setUserData("aiming_bank_id", "") # (kept for the files read by the Mantid scripts)
      launcher.setUserData("nominal_source_sample_distance_meters", str(launcher.getParameterDouble('nominal_source_sample_distance_meters')))
      if(launcher.getGen().getName()=="LOKI.FloodSourceGen/FloodSourceGen"): #event_gen=flood
        launcher.setUserData("source_monitor_distance_meters", str(launcher.getGen().source_monitor_distance_meters))
        launcher.setUserData("sampling_cone_opening_deg", str(launcher.getGen().cone_opening_deg))
        launcher.setUserData("sampling_cone_opening_min_deg", str(launcher.getGen().cone_opening_min_deg))
        launcher.setUserData("neutron_wavelength_min_aangstrom", str(launcher.getGen().neutron_wavelength_min_aangstrom))
        launcher.setUserData("neutron_wavelength_max_aangstrom", str(launcher.getGen().neutron_wavelength_max_aangstrom))

    launcher.addPrePreInitHook(setParamsForLarmor2022Experiment) #Do it when all parameters are available
    launcher.addPrePreInitHook(addUserData) #add userdata when all parameters are available

    #filter:
    if launcher.getParameterBoolean('primary_only'):
        import G4CollectFilters.StepFilterPrimary as F
        launcher.setFilter(F.create())

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

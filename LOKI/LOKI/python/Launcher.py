
def launch(geo):
    import G4Launcher
    launcher = G4Launcher()

    launcher.addParameterBoolean('primary_only',False)
    launcher.addParameterString('event_gen','')
    launcher.addParameterBoolean('gravity',False)
    launcher.addParameterBoolean('addgeantinoes',False)
    launcher.addParameterBoolean('det_only',False)
    launcher.addParameterBoolean('masking_only',False)
    launcher.addParameterDouble("sample_generator_distance_meters", 0.2) #Sample to MCPL_output mcstas component distance
    launcher.addParameterString('mcplDirectory','')
    launcher.addParameterString('input_file', '')
    launcher.addParameterDouble("gen_x_offset_meters", 0.0)

    launcher.addParameterInt("analysis_straw_pixel_number", 0) # zero means using default pixel number
    if(launcher.getParameterInt('analysis_straw_pixel_number')):
      launcher.setUserData("analysis_straw_pixel_number", str(launcher.getParameterInt('analysis_straw_pixel_number')))

   #geometry:
    launcher.setGeo(geo)

    #generator:
    if launcher.getParameterString('event_gen')=='mcpl':
        import G4MCPLPlugins.MCPLGen as Gen
        gen = Gen.create()
        if(launcher.getParameterString('input_file')):
          gen.input_file = launcher.getParameterString('input_file')
        else:
          gen.input_file = launcher.getParameterString('mcplDirectory') + 'larmor_postsample.mcpl.gz'

        import MCPL
        tmp_myfile = MCPL.MCPLFile(gen.getParameterString('input_file'))
        for blobkey in tmp_myfile.blobs:
          launcher.setUserData(blobkey, str(tmp_myfile.blobs[blobkey]))
        if('sample_mcpl_distance_m' in tmp_myfile.blobs):
          gen.dz_meter = float(tmp_myfile.blobs['sample_mcpl_distance_m'])
        else:
          gen.dz_meter = launcher.getParameterDouble('sample_generator_distance_meters')#translate z coordinates by McStas Sample-MCPL_output distance
        gen.dx_meter = launcher.getParameterDouble('gen_x_offset_meters')
        gen.exposeParameter("larmor_2022_experiment",geo,"geo_larmor_2022_experiment")
    elif launcher.getParameterString('event_gen')=='flood':
        from  LOKI.FloodSourceGen import FloodSourceGen as Gen
        gen = Gen()
        gen.gen_x_offset_meters = launcher.getParameterDouble('gen_x_offset_meters')
        gen.exposeParameter("larmor_2022_experiment",geo,"geo_larmor_2022_experiment")
    elif launcher.getParameterString('event_gen')=='masking':
        from  LOKI.MaskingSourceGen import MaskingSourceGen as Gen
        gen = Gen()
        # gen.exposeParameters(geo,"geo_")
        gen.exposeParameter("rear_detector_distance_m",geo,"geo_rear_detector_distance_m")
        gen.exposeParameter("larmor_2022_experiment",geo,"geo_larmor_2022_experiment")
        gen.exposeParameter("old_tube_numbering",geo,"geo_old_tube_numbering")
        gen.gen_x_offset_meters = launcher.getParameterDouble('gen_x_offset_meters')
    elif launcher.getParameterString('event_gen')=='spheremodel':
        from  LOKI.SansSphereGen import SansSphereGen as Gen
        gen = Gen()
    elif launcher.getParameterString('event_gen')=='isotheta':
        from  LOKI.IsoThetaGen import IsoThetaGen as Gen
        gen = Gen()
    else:
        import G4StdGenerators.FlexGen as Gen
        gen = Gen.create()
        gen.particleName = 'neutron'
        #gen.particleName = 'geantino'
        gen.neutron_wavelength_aangstrom = 3.0
        gen.momdir_spherical = True
        gen.randomize_polarangle = True
        gen.random_min_polarangle_deg = 0 #47.99
        gen.random_max_polarangle_deg =0.001  #48 #16.0
        gen.randomize_azimuthalangle = True
        gen.random_min_azimuthalangle_deg = 0 #20.0
        gen.random_max_azimuthalangle_deg = 360.0 #60
    launcher.setGen(gen)

    def assertParamsForLarmor2022Experiment(): #note: prone to generator name change
      if(launcher.getGen().hasParameterBoolean('geo_larmor_2022_experiment') and
         launcher.getGen().getParameterBoolean('geo_larmor_2022_experiment')==True):
        if(launcher.getGen().getName()=="G4MCPLPlugins/MCPLGen"): #event_gen=mcpl
          assert launcher.getGen().dz_meter == 4.049, "sample_generator_distance_meters should be 4.049, or wrong sample_mcpl_distance_m value in the MCPL file" #note: intentionally 4.049, not 4.099
          assert launcher.getGen().dx_meter == 0.005, "gen_x_offset_meters should be 0.005 for larmor 2022 experiment"
        elif(launcher.getGen().getName()=="LOKI.FloodSourceGen/FloodSourceGen"): #event_gen=flood
          assert launcher.getGen().gen_x_offset_meters == 0.005, "gen_x_offset_meters should be 0.005 for larmor 2022 experiment"
          launcher.getGen().source_sample_distance_meters = 25.61
          launcher.getGen().source_monitor_distance_meters = 25.57
          import math as m
          launcher.getGen().cone_opening_deg = m.acos(1-2/233)/m.pi*180
          launcher.setUserData("source_monitor_distance_meters", str(launcher.getGen().source_monitor_distance_meters))
          launcher.setUserData("sampling_cone_opening_deg", str(launcher.getGen().cone_opening_deg))
          launcher.setUserData("neutron_wavelength_min_aangstrom", str(launcher.getGen().neutron_wavelength_min_aangstrom))
          launcher.setUserData("neutron_wavelength_max_aangstrom", str(launcher.getGen().neutron_wavelength_max_aangstrom))
          print(f"Using predifined parameters for the Larmor-2022 experiment!")
          print(f'    source_sample_distance_meters: {launcher.getGen().source_sample_distance_meters}')
          print(f'    source_monitor_distance_meters: {launcher.getGen().source_monitor_distance_meters}')
          print(f'    cone_opening_deg: {launcher.getGen().cone_opening_deg}')

        elif(launcher.getGen().getName()=="LOKI.MaskingSourceGen/MaskingSourceGen"): #event_gen=masking
          assert launcher.getGen().gen_x_offset_meters == 0.005, "gen_x_offset_meters should be 0.005 for larmor 2022 experiment"

    launcher.addPrePreInitHook(assertParamsForLarmor2022Experiment) #Do it after the geo.larmor_2022_experiment input parameter's value is available

    #filter:
    if launcher.getParameterBoolean('primary_only'):
        import G4CollectFilters.StepFilterPrimary as F
        launcher.setFilter(F.create())

    #gravity
    if launcher.getParameterBoolean('gravity'):
        import G4GravityHelper.NeutronGravity as ng
        ng.enableNeutronGravity(launcher)

    if launcher.getParameterBoolean('addgeantinoes'):
        import G4GeantinoInserter
        G4GeantinoInserter.install()

    #general stuff:
    #if G4Launcher.g4version()<1030:
        # The following two lines were appropriate/needed in older Geant4 with the QGSP_BIC_HP
        # physics list. If not interested in Geant4 10.00.p03 support, it is safe to remove
        # the next two lines (and the if-statement). More info at DGSW-305.
        #launcher.cmd_postinit('/process/eLoss/StepFunction 0.1 0.001 um')
        #launcher.cmd_postinit('/process/eLoss/minKinEnergy 10 eV')


    if not launcher.getParameterBoolean('det_only'):
        launcher.setOutput('lokisim','REDUCED')
    elif not launcher.getParameterBoolean('masking_only'):
        griff_output_volumes = ["Converter", "B4CPanel",
        "BoronMask-triangular-7-3", "BoronMask-triangular-7-2", "BoronMask-triangular-5-0", "BoronMask-triangular-5-1",
        "BoronMask-8-0", "BoronMask-8-1", "BoronMask-8-2", "BoronMask-8-3",  "BoronMask-8-4", "BoronMask-8-5", "BoronMask-8-6", "BoronMask-8-7",
        "BoronMask-7-0", "BoronMask-7-1", "BoronMask-7-2", "BoronMask-7-3",  "BoronMask-7-4", "BoronMask-7-5", "BoronMask-7-6", "BoronMask-7-7",
        "BoronMask-6-0", "BoronMask-6-1", "BoronMask-6-2", "BoronMask-6-3",  "BoronMask-6-4", "BoronMask-6-5", "BoronMask-6-6", "BoronMask-6-7",
        "BoronMask-5-0", "BoronMask-5-1", "BoronMask-5-2", "BoronMask-5-3",  "BoronMask-5-4", "BoronMask-5-5", "BoronMask-5-6", "BoronMask-5-7",
        "BoronMask-4-0", "BoronMask-4-1", "BoronMask-4-2", "BoronMask-4-3",  "BoronMask-4-4", "BoronMask-4-5", "BoronMask-4-6", "BoronMask-4-7",
        "BoronMask-3-0", "BoronMask-3-1", "BoronMask-3-2", "BoronMask-3-3",  "BoronMask-3-4", "BoronMask-3-5", "BoronMask-3-6", "BoronMask-3-7",
        "BoronMask-2-0", "BoronMask-2-1", "BoronMask-2-2", "BoronMask-2-3",  "BoronMask-2-4", "BoronMask-2-5", "BoronMask-2-6", "BoronMask-2-7",
        "BoronMask-1-0", "BoronMask-1-1", "BoronMask-1-2", "BoronMask-1-3",  "BoronMask-1-4", "BoronMask-1-5", "BoronMask-1-6", "BoronMask-1-7",
        "BoronMask-0-0", "BoronMask-0-1", "BoronMask-0-2", "BoronMask-0-3",  "BoronMask-0-4", "BoronMask-0-5" ]
        import G4CollectFilters.StepFilterVolume
        f = G4CollectFilters.StepFilterVolume.create()
        f.volumeList = griff_output_volumes
        launcher.setFilter(f)
        launcher.setOutput('lokibcs_masking','REDUCED')
    else:
        griff_output_volumes = ["CountingGas","Converter"]
        import G4CollectFilters.StepFilterVolume
        f = G4CollectFilters.StepFilterVolume.create()
        f.volumeList = griff_output_volumes
        launcher.setFilter(f)
        launcher.setOutput('lokibcs_CountingGas_Converter','REDUCED')

    #launch:
    launcher.go()

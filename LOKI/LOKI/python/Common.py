
def launch(geo):
    import G4Launcher
    launcher = G4Launcher()

    launcher.addParameterBoolean('primary_only',False)
#    launcher.addParameterBoolean('event_gen',False)
    launcher.addParameterString('event_gen','')
    launcher.addParameterBoolean('gravity',False)
    launcher.addParameterBoolean('addgeantinoes',False)
    launcher.addParameterBoolean('det_only',False)
    launcher.addParameterBoolean('masking_only',False)
    launcher.addParameterDouble("sample_generator_distance_meters", 0.2) #All MCPL_output components are positioned at 0.2 m distance from the sample
    launcher.addParameterString('mcplDirectory','')    
    launcher.addParameterDouble("gen_x_offset_meters", 0.005) #5 mm offset for rear bank at Larmor experiment

    launcher.addParameterInt("analysis_straw_pixel_number", 0) # zero means using default pixel number 
    if(launcher.getParameterInt('analysis_straw_pixel_number')):
      launcher.setUserData("analysis_straw_pixel_number", str(launcher.getParameterInt('analysis_straw_pixel_number')))

   #geometry:
    launcher.setGeo(geo)

    #generator:
    if launcher.getParameterString('event_gen')=='mcpl':
        import G4MCPLPlugins.MCPLGen as Gen
        gen = Gen.create()
        #gen.input_file = 'testbcs.mcpl' # use input_file argument to set the full path to the file
        gen.dz_meter = launcher.getParameterDouble('sample_generator_distance_meters') #translate z coordinates by McStas Sample-MCPL_output distance
        gen.dx_meter = launcher.getParameterDouble('gen_x_offset_meters')
        gen.input_file = launcher.getParameterString('mcplDirectory') + 'larmor_postsample.mcpl.gz'

        gen.exposeParameter("larmor_2022_experiment",geo,"geo_larmor_2022_experiment")
        def overrideGenDzForLarmor2022Experiment(): # Fixed gen dz value for larmor2022experiment
          larmor_2022_experiment = launcher.getGen().getParameterBoolean('geo_larmor_2022_experiment')
          if larmor_2022_experiment:
            launcher.getGen().dz_meter = 4.049 #note: intentionally 4.049, not 4.099
            print("Generator dz overriden for Larmor2022 experiment setup. New value is 4.099 m.") 
        launcher.addPrePreInitHook(overrideGenDzForLarmor2022Experiment) #(possibly) override gen dz after the geo.larmor_2022_experiment input parameter's value is available
    elif launcher.getParameterString('event_gen')=='ascii':
        raise ValueError("event_gen=ascii is no longer supported. Please use mcpl input instead")
    elif launcher.getParameterString('event_gen')=='spheremodel':
        from  LOKI.SansSphereGen import SansSphereGen as Gen
        gen = Gen()
    elif launcher.getParameterString('event_gen')=='isotheta':
        from  LOKI.IsoThetaGen import IsoThetaGen as Gen
        gen = Gen()
    elif launcher.getParameterString('event_gen')=='flood':
        from  LOKI.FloodSourceGen import FloodSourceGen as Gen
        gen = Gen()
        gen.exposeParameter("larmor_2022_experiment",geo,"geo_larmor_2022_experiment")
    elif launcher.getParameterString('event_gen')=='masking':
        from  LOKI.MaskingSourceGen import MaskingSourceGen as Gen
        gen = Gen()
        # gen.exposeParameters(geo,"geo_")
        gen.exposeParameter("rear_detector_distance_m",geo,"geo_rear_detector_distance_m")
        gen.exposeParameter("larmor_2022_experiment",geo,"geo_larmor_2022_experiment")
        gen.exposeParameter("old_tube_numbering",geo,"geo_old_tube_numbering")
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

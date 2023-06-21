
def launch(geo):
    import G4Launcher
    launcher = G4Launcher()

    launcher.addParameterBoolean('primary_only',False)
    launcher.addParameterString('event_gen','')
    launcher.addParameterBoolean('gravity',False)
    launcher.addParameterBoolean('addgeantinoes',False)
    launcher.addParameterString('mcplDirectory','')
    launcher.addParameterBoolean('det_only',False)

   #geometry:
    launcher.setGeo(geo)

    if launcher.getParameterString('event_gen')=='mcpl':
        import G4MCPLPlugins.MCPLGen as Gen
        gen = Gen.create()
        gen.input_file = launcher.getParameterString('mcplDirectory') + 'larmor_postsample.mcpl.gz'
        gen.dz_meter = 4.356 #translate z coordinate by McStas Sample-MCPL_output distance
        gen.dx_meter = -0.040
    elif launcher.getParameterString('event_gen')=='flood':
        import math
        from  LOKI.FloodSourceGen import FloodSourceGen as Gen
        ssd = 25.3 #nominal_source_sample_distance_meters
        gen = Gen(ssd)
        gen.cone_opening_deg = math.acos(1-2/120)/math.pi*180 #for more efficient sampling (math.acos(1-2/120)/math.pi*180 = 10.4753)
        gen.gen_x_offset_meters = -0.040
        gen.source_monitor_distance_meters = 25.26
    else:
        import G4StdGenerators.FlexGen as Gen
        gen = Gen.create()
        gen.particleName = 'neutron'
        gen.neutron_wavelength_aangstrom = 5
        gen.momdir_spherical = True
        gen.randomize_polarangle = True
        gen.random_min_polarangle_deg = 0.0
        gen.random_max_polarangle_deg = 0.0001 #15# 5 #5.3 # 8.5
        gen.randomize_azimuthalangle = True
        gen.random_min_azimuthalangle_deg = 0.0
        gen.random_max_azimuthalangle_deg = 360.0
        gen.fixed_x_meters = -0.040

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


    #Framework setup
    if not launcher.getParameterBoolean('det_only'):
        launcher.setOutput('larmor','REDUCED')#Griff output)
    else:
        griff_output_volumes = ["CountingGas","Converter"]
        import G4CollectFilters.StepFilterVolume
        f = G4CollectFilters.StepFilterVolume.create()
        f.volumeList = griff_output_volumes
        launcher.setFilter(f)
        launcher.setOutput('larmor_CountingGas_Converter','REDUCED')

    #launch:
    launcher.go()

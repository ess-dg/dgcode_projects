import G4CustomPyGen
import G4Units.Units as Units

def launch(geo):
    import G4Launcher
    launcher = G4Launcher()
    launcher.addParameterBoolean('primary_only',False)
    launcher.addParameterString('event_gen','')
    #launcher.addParameterBoolean('gauss_gen',False)
    launcher.addParameterBoolean('gravity',False)
    launcher.addParameterBoolean('addgeantinoes',False)
    #launcher.addParameterBoolean('addgeantinoes',True)

    #geometry:
    launcher.setGeo(geo)

    #generator:
    if launcher.getParameterString('event_gen')=='mcpl':
        import G4MCPLPlugins.MCPLGen as Gen
        gen = Gen.create()
        #gen.input_file = 'MCPLTests/reffile_1.mcpl'
        gen.input_file ='sample8_coh_withTransmission_onlyTrasmissionNeutrons.mcpl.gz.mcpl.gz'
        #gen.input_file ='sample8_coh_withTransmission.mcpl.gz'
    elif launcher.getParameterString('event_gen')=='ascii':
        raise ValueError("event_gen=ascii is no longer supported. Please use mcpl input instead")
    elif launcher.getParameterString('event_gen')=='spheremodel':
        from  LokiSim.SansSphereGen import SansSphereGen as Gen
        gen = Gen()
    else:
        import G4StdGenerators.FlexGen as Gen
        gen = Gen.create()
        gen.particleName = 'neutron'
        gen.neutron_wavelength_aangstrom = 1.8
        gen.momdir_spherical = True
        gen.randomize_polarangle = True
        gen.random_min_polarangle_deg = 0.0
        gen.random_max_polarangle_deg = 5.3 # 8.5 #cover the whole 1x1 m^2 area
        #gen.random_max_polarangle_deg = 4.9
        gen.randomize_azimuthalangle = True
        gen.random_min_azimuthalangle_deg = 0.0
        gen.random_max_azimuthalangle_deg = 360.0

        #gen.fixed_y_meters=-0.01213
        #gen.fixed_y_meters=0.11547

        #gen.randomize_y=True
        #gen.random_min_y_meters=-0.5
        #gen.random_max_y_meters=0.5

        #gen.randomize_x=True
        #gen.random_min_x_meters=-0.5
        #gen.random_max_x_meters=0.5

        #gen.fixed_y_meters= 0.12760
        #gen.fixed_polarangle_deg = 30.0
        #gen.fixed_azimuthalangle_deg = 270.0

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
    if G4Launcher.g4version()<1030:
        # The following two lines were appropriate/needed in older Geant4 with the QGSP_BIC_HP
        # physics list. If not interested in Geant4 10.00.p03 support, it is safe to remove
        # the next two lines (and the if-statement). More info at DGSW-305.
        launcher.cmd_postinit('/process/eLoss/StepFunction 0.1 0.001 um')
        launcher.cmd_postinit('/process/eLoss/minKinEnergy 10 eV')
    launcher.setOutput('bcs','FULL')

    #launch:
    launcher.go()

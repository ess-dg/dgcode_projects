
import G4GeoLoki.LokiAimHelper as LokiAim
import Core.Units as Units
import numpy as np

def launch(geo):
    import G4Launcher
    launcher = G4Launcher()

    launcher.addParameterString('event_gen','')
    launcher.addParameterInt("analysis_straw_pixel_number", 256)
    launcher.addParameterBoolean('det_only',False)
    launcher.addParameterString('bank_filter','') #aim only at a certain bank of group of banks(rear detector('rear'), mid-detector('mid), front detector('front))
    launcher.addParameterDouble("gen_x_offset_meters", 0.0)
    ## McStas+Geant4 options ##
    launcher.addParameterString('mcplDirectory','')
    launcher.addParameterString('input_file', '')
    launcher.addParameterDouble("sample_generator_distance_meters", 0.2) #Sample to MCPL_output mcstas component distance

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
    elif launcher.getParameterString('event_gen')=='flood': #Flood source simulation
        from  LOKI.FloodSourceGen import FloodSourceGen as Gen
        gen = Gen()
        gen.gen_x_offset_meters = launcher.getParameterDouble('gen_x_offset_meters')
        bankFilter = launcher.getParameterString('bank_filter')
        if bankFilter != '':
          assert bankFilter in (*[str(i) for i in range(9)], '1234', '5678'), f"bank_filter must be either [0,8], 1234 or 5678"
          if bankFilter == '1234': #mid banks
            angleRange = (3.5, 15.8)
          elif bankFilter == '5678': #front banks
            angleRange = (10.5, 49.6)
          else: #single banks [0,8]
            bankId = int(bankFilter)
            aimHelper = LokiAim.AimHelper(5*Units.m) #rear det distance shouldn't really matter
            bankCentre = [aimHelper.getBankPosition(bankId, 0), aimHelper.getBankPosition(bankId, 1), aimHelper.getBankPosition(bankId, 2)]
            gen.ref_dir_x, gen.ref_dir_y, gen.ref_dir_z = np.array(bankCentre)/np.linalg.norm(bankCentre)
            bankConeAngle = [5.07, 9.9, 4.9, 9.9, 4.9, 31.9, 20.9, 29.5, 22.1] #HARDCODED for now
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
        gen.particleName = 'neutron' if not launcher.getParameterBoolean('geantino') else 'geantino'
        gen.neutron_wavelength_aangstrom = 3.0
        gen.momdir_spherical = True
        gen.randomize_polarangle = True
        gen.random_min_polarangle_deg = 0 #49.99
        gen.random_max_polarangle_deg =0.001  #50 #16.0
        gen.randomize_azimuthalangle = True
        gen.random_min_azimuthalangle_deg = 0 #20.0
        gen.random_max_azimuthalangle_deg = 360.0 #60

    gen.exposeParameter("larmor_2022_experiment",geo,"geo_larmor_2022_experiment")
    launcher.setGen(gen)

    def assertParamsForLarmor2022Experiment(): #note: prone to generator name change
      if(launcher.getGen().hasParameterBoolean('geo_larmor_2022_experiment') and
         launcher.getGen().getParameterBoolean('geo_larmor_2022_experiment')==True):
        assert launcher.getParameterInt('analysis_straw_pixel_number') == 512, "analysis_straw_pixel_number must be 512 for the Larmor2022 experiment!"
        if(launcher.getGen().getName()=="G4MCPLPlugins/MCPLGen"): #event_gen=mcpl
          assert launcher.getGen().dz_meter == 4.049, "sample_generator_distance_meters should be 4.049 (or wrong sample_mcpl_distance_m value in the MCPL file) for the Larmor2022 experiment!" #note: intentionally 4.049, not 4.099
          assert launcher.getGen().dx_meter == 0.005, "gen_x_offset_meters should be 0.005 for the Larmor 2022 experiment!"
        elif(launcher.getGen().getName()=="LOKI.FloodSourceGen/FloodSourceGen"): #event_gen=flood
          assert launcher.getGen().gen_x_offset_meters == 0.005, "gen_x_offset_meters should be 0.005 for the Larmor2022 experiment!"
          launcher.getGen().source_sample_distance_meters = 25.61
          launcher.getGen().source_monitor_distance_meters = 25.57
          import math as m
          launcher.getGen().cone_opening_deg = m.acos(1-2/233)/m.pi*180
          print(f"Using predifined parameters for the Larmor2022 experiment!")
          print(f'    source_sample_distance_meters: {launcher.getGen().source_sample_distance_meters}')
          print(f'    source_monitor_distance_meters: {launcher.getGen().source_monitor_distance_meters}')
          print(f'    cone_opening_deg: {launcher.getGen().cone_opening_deg}')

    def addUserData():
      if(launcher.getGen().getName()=="LOKI.FloodSourceGen/FloodSourceGen"): #event_gen=flood
        launcher.setUserData("source_monitor_distance_meters", str(launcher.getGen().source_monitor_distance_meters))
        launcher.setUserData("sampling_cone_opening_deg", str(launcher.getGen().cone_opening_deg))
        launcher.setUserData("sampling_cone_opening_min_deg", str(launcher.getGen().cone_opening_min_deg))
        launcher.setUserData("neutron_wavelength_min_aangstrom", str(launcher.getGen().neutron_wavelength_min_aangstrom))
        launcher.setUserData("neutron_wavelength_max_aangstrom", str(launcher.getGen().neutron_wavelength_max_aangstrom))
        launcher.setUserData("analysis_straw_pixel_number", str(launcher.getParameterInt('analysis_straw_pixel_number')))
        launcher.setUserData("rear_detector_distance_m", str(launcher.getGeo().getParameterDouble("rear_detector_distance_m")))
        launcher.setUserData("bank_filter", str(launcher.getParameterString('bank_filter')))
          
    launcher.addPrePreInitHook(assertParamsForLarmor2022Experiment) #Do it after the geo.larmor_2022_experiment input parameter's value is available
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

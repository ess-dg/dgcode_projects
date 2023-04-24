
def launch(geo):
    import G4Launcher
    launcher = G4Launcher()

    launcher.addParameterInt("analysis_straw_pixel_number", 256)
    launcher.addParameterDouble("gen_x_offset_meters", 0.0)
    launcher.addParameterBoolean('masking_only',False)
    # launcher.addParameterString('bank_filter','') #aim only at a certain bank of group of banks(rear detector('rear'), mid-detector('mid), front detector('front))

    #geometry:
    launcher.setGeo(geo)

    #generator:
    from LokiMasking.MaskingSourceGen import MaskingSourceGen as Gen
    gen = Gen()
    gen.exposeParameter("rear_detector_distance_m",geo,"geo_rear_detector_distance_m")
    gen.exposeParameter("old_tube_numbering",geo,"geo_old_tube_numbering")
    gen.gen_x_offset_meters = launcher.getParameterDouble('gen_x_offset_meters')  
    gen.exposeParameter("larmor_2022_experiment",geo,"geo_larmor_2022_experiment")
    launcher.setGen(gen)

    def assertParamsForLarmor2022Experiment(): #note: prone to generator name change
      if(launcher.getGen().hasParameterBoolean('geo_larmor_2022_experiment') and
         launcher.getGen().getParameterBoolean('geo_larmor_2022_experiment')==True):
        assert launcher.getGen().gen_x_offset_meters == 0.005, "gen_x_offset_meters should be 0.005 for larmor 2022 experiment"

    def addUserData():
      launcher.setUserData("analysis_straw_pixel_number", str(launcher.getParameterInt('analysis_straw_pixel_number')))
      launcher.setUserData("rear_detector_distance_m", str(launcher.getGeo().getParameterDouble("rear_detector_distance_m")))
      # launcher.setUserData("bank_filter", str(launcher.getParameterString('bank_filter')))
          
    launcher.addPrePreInitHook(assertParamsForLarmor2022Experiment) #Do it after the geo.larmor_2022_experiment input parameter's value is available
    launcher.addPrePreInitHook(addUserData) #add userdata when all parameters are available

    #filter:
    if not launcher.getParameterBoolean('masking_only'):#TODO reduced material list should be the default
        launcher.setOutput('loki_masking','REDUCED')
    else: 
        griff_output_volumes = ["Converter", "B4CPanel", "AlPanel", #TODO Added 'AlPanel', but not tested
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
        launcher.setOutput('loki_masking','REDUCED')

    #launch:
    launcher.go()

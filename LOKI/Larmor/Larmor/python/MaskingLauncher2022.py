# Masking simulation (geantinos aimed at the pixel centres) of the Larmor 2022 experiment, with the
# Larmor/GeoLarmor2022 geometry. (Extracted from LokiMasking.Launcher, where it was enabled by the
# larmor_2022_experiment geometry parameter.)

def launch(geo):
    import G4Launcher
    launcher = G4Launcher()

    launcher.addParameterInt("analysis_straw_pixel_number", 512) # must be 512 (the experiment)
    launcher.addParameterDouble("gen_x_offset_meters", 0.005)
    launcher.addParameterBoolean('masking_only',False)

    #geometry:
    launcher.setGeo(geo)

    #generator:
    from Larmor.MaskingSourceGen2022 import MaskingSourceGen2022 as Gen
    gen = Gen()
    gen.exposeParameter("rear_detector_distance_m",geo,"geo_rear_detector_distance_m")
    gen.exposeParameter("old_tube_numbering",geo,"geo_old_tube_numbering")
    gen.gen_x_offset_meters = launcher.getParameterDouble('gen_x_offset_meters')
    launcher.setGen(gen)

    def assertParamsForLarmor2022Experiment(): #note: prone to generator name change
        assert launcher.getParameterInt('analysis_straw_pixel_number') == 512, "analysis_straw_pixel_number must be 512 for the Larmor2022 experiment!"
        assert launcher.getGen().gen_x_offset_meters == 0.005, "gen_x_offset_meters should be 0.005 for the Larmor2022 experiment!"

    def addUserData():
      launcher.setUserData("analysis_straw_pixel_number", str(launcher.getParameterInt('analysis_straw_pixel_number')))
      launcher.setUserData("rear_detector_distance_m", str(launcher.getGeo().getParameterDouble("rear_detector_distance_m")))

    launcher.addPrePreInitHook(assertParamsForLarmor2022Experiment) #Do it when all parameters are available
    launcher.addPrePreInitHook(addUserData) #add userdata when all parameters are available

    #filter:
    if launcher.getParameterBoolean('masking_only'): # only the volumes of the converters and of the absorbers of the bank
        griff_output_volumes = ["Converter", "B4CPanel", "AlPanel",
        "BoronMask-0-0", "BoronMask-0-1", "BoronMask-0-2", "BoronMask-0-3",  "BoronMask-0-4", "BoronMask-0-5" ]
        import G4CollectFilters.StepFilterVolume
        f = G4CollectFilters.StepFilterVolume.create()
        f.volumeList = griff_output_volumes
        launcher.setFilter(f)
    launcher.setOutput('loki_masking','REDUCED')

    #launch:
    launcher.go()

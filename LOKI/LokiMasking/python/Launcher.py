import G4GeoLoki.LokiAimHelper as LokiAim

def launch(geo):
    import G4Launcher
    launcher = G4Launcher()

    launcher.addParameterInt("analysis_straw_pixel_number", 256)
    launcher.addParameterDouble("gen_x_offset_meters", 0.0)
    launcher.addParameterBoolean('masking_only',False)
    launcher.addParameterInt('aiming_bank_id',-1) #aim only at a certain bank

    #geometry:
    launcher.setGeo(geo)

    #generator:
    from LokiMasking.MaskingSourceGen import MaskingSourceGen as Gen
    gen = Gen()
    gen.exposeParameter("rear_detector_distance_m",geo,"geo_rear_detector_distance_m")
    gen.gen_x_offset_meters = launcher.getParameterDouble('gen_x_offset_meters')
    gen.aiming_bank_id = launcher.getParameterInt('aiming_bank_id')
    if geo.hasParameterString("bank_calibration"):
        gen.exposeParameter("bank_calibration",geo,"geo_bank_calibration") #aim at the pixels of the calibrated banks
    launcher.setGen(gen)

    def addUserData():
      launcher.setUserData("analysis_straw_pixel_number", str(launcher.getParameterInt('analysis_straw_pixel_number')))
      launcher.setUserData("rear_detector_distance_m", str(launcher.getGeo().getParameterDouble("rear_detector_distance_m")))
      launcher.setUserData("aiming_bank_id", str(launcher.getParameterInt('aiming_bank_id')))
      # the bank placements: the name and the text of the calibration (the analysis uses the recorded text, so that
      # it has the placements of the simulation even if the calibration file changes or is not there)
      bankCalibration = launcher.getGeo().getParameterString("bank_calibration")
      launcher.setUserData("bank_calibration", bankCalibration)
      launcher.setUserData("bank_calibration_text", LokiAim.bankCalibrationText(bankCalibration))

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
        "BoronMask-0-0", "BoronMask-0-1", "BoronMask-0-2", "BoronMask-0-3",  "BoronMask-0-4", "BoronMask-0-5",
        # the beamstop (beamstop_id) and the calibration slit masks (with_calibration_slits) absorb as well
        "BoronMask-Beamstop", *[f"BoronMask-lokiStandard-{bankId}" for bankId in range(9)] ]
        import G4CollectFilters.StepFilterVolume
        f = G4CollectFilters.StepFilterVolume.create()
        f.volumeList = griff_output_volumes
        launcher.setFilter(f)
        launcher.setOutput('loki_masking','REDUCED')

    #launch:
    launcher.go()

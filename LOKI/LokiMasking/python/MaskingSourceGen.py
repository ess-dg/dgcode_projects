from __future__ import print_function
import G4CustomPyGen
import Core.Units as Units
import G4GeoLoki.LokiAimHelper as LokiAim

class MaskingSourceGen(G4CustomPyGen.GenBase):
    def declare_parameters(self):
        self.addParameterDouble("gen_x_offset_meters", 0.0)
        self.addParameterDouble("gen_x_width_meters", 0.0)
        self.addParameterDouble("gen_y_width_meters", 0.0)
        #self.addParameterDouble("z_width_meters", 0.001)
        self.addParameterInt("aiming_bank_id", 0, 0, 8)
        self.addParameterInt("aiming_pixel_id_min", -1)
        self.addParameterInt("aiming_straw_pixel_number", 512) #used only for aiming, NOT for analysis

    def init_generator(self,gun):
        gun.set_type('geantino')

        self.aimHelper = LokiAim.AimHelper(self.geo_rear_detector_distance_m *Units.m, self.aiming_straw_pixel_number)
        self.totalNumberOfPixels = self.aimHelper.getTotalNumberOfPixels()
        bank_pixel_id_min = self.aimHelper.getBankPixelOffset(self.aiming_bank_id)
        number_of_pixels_in_bank = self.aimHelper.getNumberOfPixels(self.aiming_bank_id)

        print(f"Total number of pixels in full geometry: {self.totalNumberOfPixels} (set -n accordingly to cover the whole geometry!)")
        print(f"Number of pixels in bank{self.aiming_bank_id}: {number_of_pixels_in_bank} (set -n accordingly to cover the whole bank!)")
        print(f"Lowest pixel id in bank{self.aiming_bank_id}: {bank_pixel_id_min}")

        self._i = bank_pixel_id_min  #count events to shoot neutrons at each pixel
        if(self.aiming_pixel_id_min>=0):
          self._i = self.aiming_pixel_id_min
          print(f"aiming_pixel_id_min is set to: {self.aiming_pixel_id_min}. This will be used as the first pixel to aim at, overriding the lowest pixel id in the selected/default bank.")
        else:
          print(f"First pixel to aim at: {self._i}")
        self._i -= 1 # incrementation happens at the begining of the generate_event loop

        #bank pixel limits: 0, 401408, 516096, 602112, 716800, 802816, 1003520, 1232896, 1376256, 1605632 (for 256 pixels/straw)
        #                   0, 802816, 1032192, 1204224, 1433600, 1605632, 2007040, 2465792, 2752512, 3211264 (for 512 pixels/straw)

    def generate_event(self,gun):
        self._i += 1
        # Source position -
        sourcePositionX = self.gen_x_width_meters *(self.rand()-0.5) *Units.m + self.gen_x_offset_meters *Units.m
        sourcePositionY = self.gen_y_width_meters *(self.rand()-0.5) *Units.m
        sourcePositionZ = 0.0
        gun.set_position(sourcePositionX, sourcePositionY, sourcePositionZ)

        # Direction - toward the centre of a pixel
        pixelId = ((self._i + 0) % self.totalNumberOfPixels)

        pixelCentreX, pixelCentreY, pixelCentreZ = self.aimHelper.getPixelCentreCoordinates(pixelId, self.geo_old_tube_numbering, self.geo_larmor_2022_experiment)

        gun.set_direction(pixelCentreX - sourcePositionX, pixelCentreY - sourcePositionY, pixelCentreZ - sourcePositionZ)
        #gun.set_direction(pixelCentreX - sourcePositionX +3*(2*self.rand()-1), pixelCentreY - sourcePositionY+3*(2*self.rand()-1), pixelCentreZ - sourcePositionZ+3*(2*self.rand()-1))

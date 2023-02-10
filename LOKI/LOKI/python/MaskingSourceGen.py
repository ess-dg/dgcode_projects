from __future__ import print_function
import G4CustomPyGen
import Core.Units as Units
#import Utils.NeutronMath
#import math
import G4GeoLoki.LokiMaskingHelper as Mask

class MaskingSourceGen(G4CustomPyGen.GenBase):
    def declare_parameters(self):
        self.addParameterDouble("gen_x_offset_meters", 0.0)
        self.addParameterDouble("gen_x_width_meters", 0.0)
        self.addParameterDouble("gen_y_width_meters", 0.0)
        #self.addParameterDouble("z_width_meters", 0.001)
        self.addParameterInt("aiming_pixel_id_min", 0)
        self.addParameterInt("aiming_straw_pixel_number", 512) #used only for aiming, NOT for analysis

    def init_generator(self,gun):
        gun.set_type('geantino') #neutron

        self._i = self.aiming_pixel_id_min -1 #count events to shoot neutrons at each pixel

        self.aimHelper = Mask.MaskingHelper(self.geo_rear_detector_distance_m *Units.m, self.aiming_straw_pixel_number)
        self.totalNumberOfPixels = self.aimHelper.getTotalNumberOfPixels()
        print(f"Total number of pixels: {self.totalNumberOfPixels}")

    def generate_event(self,gun):
        self._i += 1
        # Source position -
        sourcePositionX = self.gen_x_width_meters *(self.rand()-0.5) *Units.m + self.gen_x_offset_meters *Units.m
        sourcePositionY = self.gen_y_width_meters *(self.rand()-0.5) *Units.m
        sourcePositionZ = 0.0
        gun.set_position(sourcePositionX, sourcePositionY, sourcePositionZ)

        # Direction - toward the centre of a pixel
        #bank pixel limits: 0, 401408, 516096, 602112, 716800, 802816, 1003520, 1232896, 1376256, 1605632
        pixelId = ((self._i + 0) % self.totalNumberOfPixels)

        pixelCentreX, pixelCentreY, pixelCentreZ = self.aimHelper.getPixelCentrePositionsForMasking(pixelId, self.geo_old_tube_numbering, self.geo_larmor_2022_experiment)

        gun.set_direction(pixelCentreX - sourcePositionX, pixelCentreY - sourcePositionY, pixelCentreZ - sourcePositionZ)
        #gun.set_direction(pixelCentreX - sourcePositionX +3*(2*self.rand()-1), pixelCentreY - sourcePositionY+3*(2*self.rand()-1), pixelCentreZ - sourcePositionZ+3*(2*self.rand()-1))

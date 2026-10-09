import G4CustomPyGen
from Units import units
import Larmor.Larmor2022Bank as LB
import G4Interfaces

# Geantinos from the sample position aimed at the pixel centres of the rear bank of the Larmor 2022 experiment, one
# pixel per event (set -n to the number of pixels to cover the whole bank). (Extracted from
# LokiMasking.MaskingSourceGen.)
class MaskingSourceGen2022(G4CustomPyGen.GenBase):
    def declare_parameters(self):
        self.addParameterDouble("gen_x_offset_meters", 0.0)
        self.addParameterDouble("gen_x_width_meters", 0.0)
        self.addParameterDouble("gen_y_width_meters", 0.0)
        self.addParameterInt("aiming_pixel_id_min", -1)
        self.addParameterInt("aiming_straw_pixel_number", 512) #used only for aiming, NOT for analysis

    def init_generator(self,gun):
        gun.set_type('geantino')

        self.bank = LB.Larmor2022Bank(self.geo_rear_detector_distance_m *units.m, self.aiming_straw_pixel_number)
        self.totalNumberOfPixels = self.bank.getTotalNumberOfPixels()
        print(f"Number of pixels in the bank: {self.totalNumberOfPixels} (set -n accordingly to cover the whole bank!)")

        self.id_min_offset = 0
        if(self.aiming_pixel_id_min>=0):
          self.id_min_offset = self.aiming_pixel_id_min
          print(f"aiming_pixel_id_min is set to: {self.aiming_pixel_id_min}. This will be used as the first pixel to aim at.")
        else:
          print(f"First pixel to aim at: {self.id_min_offset}")
        self.m_nprocs = 0 #number of processes

    def delayed_init(self):
        self.m_nprocs = G4Interfaces.nProcs()
        if(self.m_nprocs==1):
          self._i = 0
        else: #parallel processing
          self._i = G4Interfaces.mpID() #id of the process
          # Each process with mpID E[0,m_nprocs-1] deal with pixels where (id % nProcs == mpID)

    def generate_event(self,gun):
        if(self.m_nprocs == 0): #only the first time
           self.delayed_init()
        # Source position -
        sourcePositionX = self.gen_x_width_meters *(self.rand()-0.5) *units.m + self.gen_x_offset_meters *units.m
        sourcePositionY = self.gen_y_width_meters *(self.rand()-0.5) *units.m
        sourcePositionZ = 0.0
        gun.set_position(sourcePositionX, sourcePositionY, sourcePositionZ)

        # Direction - toward the centre of a pixel
        pixelId = ((self._i + self.id_min_offset) % self.totalNumberOfPixels)

        pixelCentreX, pixelCentreY, pixelCentreZ = self.bank.getPixelCentre(pixelId, self.geo_old_tube_numbering)

        gun.set_direction(pixelCentreX - sourcePositionX, pixelCentreY - sourcePositionY, pixelCentreZ - sourcePositionZ)
        self._i += self.m_nprocs

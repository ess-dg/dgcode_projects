"""Test-only particle generator: shoots particles from the nominal sample
position (origin) towards a list of AimHelper pixel centres (each target is used
for 'repeat' consecutive events). Used by the G4GeoLokiTests end-to-end tests."""

import G4CustomPyGen
import G4GeoLoki.LokiAimHelper as LokiAim

class PixelTargetGen(G4CustomPyGen.GenBase):
    def declare_parameters(self):
        self.addParameterString("target_pixel_ids", "0")  # comma separated pixel ids
        self.addParameterInt("repeat", 1)
        self.addParameterDouble("rear_detector_distance_mm", 5000.0)
        self.addParameterInt("aiming_straw_pixel_number", 256)
        self.addParameterBoolean("old_tube_numbering", False)
        self.addParameterString("particle", "neutron")
        self.addParameterDouble("neutron_wavelength_aangstrom", 4.0)

    def init_generator(self, gun):
        gun.set_type(self.particle)
        if self.particle == 'neutron':
            gun.set_wavelength_angstrom(self.neutron_wavelength_aangstrom)
        gun.set_position(0, 0, 0)
        self._aim = LokiAim.AimHelper(self.rear_detector_distance_mm, self.aiming_straw_pixel_number)
        self._targets = [int(e) for e in self.target_pixel_ids.split(',')]
        self._i = 0

    def generate_event(self, gun):
        pid = self._targets[(self._i // self.repeat) % len(self._targets)]
        self._i += 1
        x, y, z = self._aim.getPixelCentreCoordinates(pid, self.old_tube_numbering, False)
        gun.set_direction(x, y, z)

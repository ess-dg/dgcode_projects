"""Test-only particle generator: shoots particles from the nominal sample position (origin) towards a list of
Larmor2022Bank pixel centres (one target per event). Used by test_geantino_aiming."""

import G4CustomPyGen
import Larmor.Larmor2022Bank as LB

class PixelTargetGen2022(G4CustomPyGen.GenBase):
    def declare_parameters(self):
        self.addParameterString("target_pixel_ids", "0")  # comma separated pixel ids
        self.addParameterInt("aiming_straw_pixel_number", 512)
        self.addParameterBoolean("aiming_old_tube_numbering", False)

    def init_generator(self, gun):
        gun.set_type('geantino')
        gun.set_position(0, 0, 0)
        self._bank = LB.Larmor2022Bank(4099.0, self.aiming_straw_pixel_number)
        self._targets = [int(e) for e in self.target_pixel_ids.split(',')]
        self._i = 0

    def generate_event(self, gun):
        pid = self._targets[self._i % len(self._targets)]
        self._i += 1
        gun.set_direction(*self._bank.getPixelCentre(pid, self.aiming_old_tube_numbering))

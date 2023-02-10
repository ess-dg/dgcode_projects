from __future__ import print_function
import G4CustomPyGen
import Core.Units as Units
import Utils.NeutronMath
import math

class IsoThetaGen(G4CustomPyGen.GenBase):
    def declare_parameters(self):
        self.addParameterDouble("neutron_wavelength_aangstrom", 3)
        self.addParameterDouble("cone_opening_deg", 48.0, 0.0, 180.0)
        self.addParameterDouble("fixed_z_meters", 0.0, -0.02, 0.02)

    def init_generator(self,gun):
        gun.set_type('neutron')
        self.thetaMax = self.cone_opening_deg * math.pi/180.
        gun.set_position(0,0,self.fixed_z_meters*Units.mm)
        wl = self.neutron_wavelength_aangstrom*Units.angstrom
        gun.set_energy(Utils.NeutronMath.neutronWavelengthToEKin(wl))
    def generate_event(self,gun):
        theta = self.rand() * self.thetaMax
        phi = self.rand() * 2.0 * math.pi
        s = math.sin(theta)
        c = math.cos(theta)
        gun.set_direction(s*math.cos(phi),s*math.sin(phi),c)

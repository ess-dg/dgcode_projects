from __future__ import print_function
import G4CustomPyGen
import Core.Units as Units
import Utils.NeutronMath
import math

class FloodSourceGen(G4CustomPyGen.GenBase):
    def declare_parameters(self):
        self.addParameterDouble("neutron_wavelength_min_aangstrom", 2)
        self.addParameterDouble("neutron_wavelength_max_aangstrom", 13)

        self.addParameterDouble("source_sample_distance_meters", 23.5706)
        self.addParameterDouble("source_monitor_distance_meters", 23.5662) #pre-sample monitor position (should match the IDF file)

        self.addParameterDouble("gen_x_offset_meters", 0.0)
        self.addParameterDouble("gen_x_width_meters", 0.006)
        self.addParameterDouble("gen_y_width_meters", 0.008)
        #self.addParameterDouble("z_width_meters", 0.001)
        self.addParameterDouble("fixed_z_meters", 0.000)

        self.addParameterDouble("cone_opening_deg", 49.6, 0.0, 180.0) #max opening angle (for more efficient sampling)
        self.addParameterDouble("cone_opening_min_deg", 0, 0.0, 180.0) #min opening angle (for more efficient sampling)

    def init_generator(self,gun):
        gun.set_type('neutron')
        assert self.cone_opening_min_deg <= self.cone_opening_deg, f'cone_opening_min_deg({self.cone_opening_min_deg}) cannot be higher that cone_opening_deg({self.cone_opening_deg})'
        #gun.set_type('geantino')

    def generate_event(self,gun):
        # Energy - uniform wavelength distribution between min and max parameters
        wavelength = (self.neutron_wavelength_min_aangstrom + self.rand()*(self.neutron_wavelength_max_aangstrom - self.neutron_wavelength_min_aangstrom))
        gun.set_energy(Utils.NeutronMath.neutronWavelengthToEKin(wavelength *Units.angstrom))

        # Initial TOF - calculated from source_sample_distance and the sampled wavelength(velocity)
        velocity = Utils.NeutronMath.neutron_angstrom_to_meters_per_second(wavelength)
        gun.set_time(self.source_sample_distance_meters / velocity *Units.second)

        # Source position - uniform surface source (could be volume)
        x_meters = self.gen_x_width_meters *(self.rand()-0.5) + self.gen_x_offset_meters
        y_meters = self.gen_y_width_meters *(self.rand()-0.5)
        gun.set_position(x_meters *Units.m, y_meters *Units.m , self.fixed_z_meters *Units.m)

        # Direction - uniform within a cone with an opening angle of alpha
        cos_alpha = math.cos(self.cone_opening_deg*math.pi/180)
        cos_beta = math.cos(self.cone_opening_min_deg*math.pi/180)
        #rho = cos_alpha + self.rand() * (1 - cos_alpha) # rand uniform [cos(alpha), 1]
        rho = cos_alpha + self.rand() * (cos_beta - cos_alpha) # rand uniform [cos(alpha), cos(beta)]
        theta = math.acos(rho)
        phi = self.rand() * 2.0 * math.pi
        x_dir = math.sin(theta) * math.cos(phi)
        y_dir = math.sin(theta) * math.sin(phi)
        z_dir = math.cos(theta)
        gun.set_direction(x_dir,y_dir,z_dir)
        #NOTE: The surface area of a spherical segment with height (h) is 2*pi*R*h (regardless of the position of the segment or spherical cap)
        #      The ratio of this area/shpere area (for R=1) is 2*pi*R*h/(4*R*R*pi) = h/2
        #      for alpha max cone opening angle, h = 1-cos(alpha) -> sampling factor = (1-cos(alpha))/2
        #      +for a beta minimum angle h = cos(beta)-cos(alpha) -> sampling factor = (cos(beta)-cos(alpha))/2

/////////////////////////////////////////
// Declaration of our geometry module: //
/////////////////////////////////////////

#include "G4Interfaces/GeoConstructPyExport.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Transform3D.hh"
#include "G4Vector3D.hh"
//#include <fstream>
#include <cmath>

class GeoBCS : public G4Interfaces::GeoConstructBase
{
public:
  GeoBCS();
  virtual ~GeoBCS(){}
  virtual G4VPhysicalVolume* Construct();
protected:
  virtual bool validateParameters();
private:
  G4LogicalVolume * createStrawLV(); // create the logical volume of each straw for easy placement inside the tube
  G4LogicalVolume * createTubeLV();  // create the logical volume of a tube for easy placement in a panel
};

// this line is necessary to be able to declare the geometry in the python simulation script
PYTHON_MODULE3 { GeoConstructPyExport::exportGeo<GeoBCS>(mod, "GeoBCS"); }

////////////////////////////////////////////
// Implementation of our geometry module: //
////////////////////////////////////////////

#include "G4Box.hh"

GeoBCS::GeoBCS()
  : GeoConstructBase("G4GeoBCS/GeoBCS")
{
  // declare all parameters that can be used from the command line,
  // define their type, pick a self-explanatory name and _unit
  // give the default value, the last 2 ones are optional (min, max)
  addParameterDouble("straw_length_cm", 100);
  addParameterDouble("straw_radius_mm", 3.75);
  addParameterDouble("straw_wall_mm", 0.025);

  addParameterDouble("tube_radius_mm", 12.7); //1" diameter (this is outer radius)
  addParameterDouble("tube_wall_mm", 0.94);

  addParameterDouble("afterburner_thickness_mm", 0,0,50);
  addParameterString("afterburner_material", "ESS_POLYETHYLENE");

  addParameterDouble("converter_thickness_um", 1);

  addParameterDouble("tube_box_mm", 25.4,25.4,50); //default value=2*tube_outer_radius (for default tube paramaters)
  addParameterString("tube_box_fill_material", "G4_Vacuum");

  addParameterString("converter_material","ESS_B4C:b10_enrichment=0.95");

  addParameterString("straw_wall_material", "NCrystal:cfg=Cu_sg225.ncmat");
  addParameterString("tube_wall_material", "NCrystal:cfg=Al_sg225.ncmat");
  addParameterString("tube_inner_gas","IdealGas:formula=0.9*Ar+0.1*CO2:pressure_atm=0.7");

  addParameterString("counting_gas","IdealGas:formula=0.9*Ar+0.1*CO2:pressure_atm=0.7");
  addParameterString("world_material","G4_Vacuum");
  addParameterDouble("sample_detector_distance_m", 5.0127, 0, 10); //(5m + tube_radius) to position the front of the first panel to 5m

  //addParameterInt("number_of_straws", 1, 1, 49);
  addParameterInt("number_of_tubes", 40, 1, 40);
  addParameterInt("number_of_panels", 4, 1, 5);


  addParameterDouble("tube_rotation_deg", 50,0,360);
  addParameterDouble("panel_offset_mm", 10.16,0,100);

}

/*
G4LogicalVolume * GeoBCS::createStrawLV()
{
  // passing the value of the above-defined parameters (or the command line ones) to the local scope
  // REMEMBER THE UNITS AT THE END, OTHERWISE YOU ARE IN TROUBLE
  const double straw_length = getParameterDouble("straw_length_cm")*Units::cm;
  const double straw_radius = getParameterDouble("straw_radius_mm")*Units::mm; // outer radius
  const double straw_wall = getParameterDouble("straw_wall_mm")*Units::mm;
  const double converter_thickness = getParameterDouble("converter_thickness_um")*Units::um;

  auto straw_wall_material = getParameterMaterial("straw_wall_material");
  auto counting_gas = getParameterMaterial("counting_gas");
  auto converter_material = getParameterMaterial("converter_material");

  const G4String name_db = "StrawWall";
  // define the shape inside the G4LogicalVolume constructor, saves a line of code from declaring the shape separately
  auto lv_straw_wall = new G4LogicalVolume(new G4Tubs(name_db, 0, straw_radius, 0.5*straw_length, 0., 2*M_PI),
                                     straw_wall_material, name_db);

  //auto lv_converter =  place(new G4Tubs("Converter", 0., straw_radius-straw_wall, 0.5*(straw_length-straw_wall), 0., 2*M_PI),
  //                           converter_material, 0,0,0, lv_straw_wall, G4Colour(0,1,1)).logvol;
  auto lv_converter =  place(new G4Tubs("Converter", 0., straw_radius-straw_wall, 0.5*(straw_length-straw_wall), 0., 2*M_PI),
                             converter_material, 0,0,0, lv_straw_wall, G4Colour(0,1,1),-3,0,0).logvol;
  // this is a framework function that is overloaded, it helps the user do shape, logical volume and placement in a single line
  // look at Framework/G4/G4Interfaces/libinc/GeoBase.hh for all versions of it

  //place(new G4Tubs("CountingGas", 0., straw_radius-straw_wall-converter_thickness, 0.5*(straw_length-straw_wall), 0., 2*M_PI),
  //      counting_gas, 0,0,0, lv_converter, G4Colour(0,0,1));
  place(new G4Tubs("CountingGas", 0., straw_radius-straw_wall-converter_thickness, 0.5*(straw_length-straw_wall), 0., 2*M_PI),
        counting_gas, 0,0,0, lv_converter, G4Colour(0,0,1),-4,0,0);


  // the hierarchy of the volumes from mother to daughter is: straw wall -> converter -> counting gas, a cylinder inside a cylinder

  // return the logical volume of the straw container
  return lv_straw_wall;
  }*/

G4LogicalVolume * GeoBCS::createTubeLV()
{
  // create the Al tube containing the 5 Cu straws

  // look at comments in method createStrawLV()

  const double straw_radius = getParameterDouble("straw_radius_mm")*Units::mm;
  //const double straw_wall = getParameterDouble("straw_wall_mm")*Units::mm;
  //const double converter_thickness = getParameterDouble("converter_thickness_um")*Units::um;
  const double tube_outer_radius = getParameterDouble("tube_radius_mm")*Units::mm;
  const double tube_wall = getParameterDouble("tube_wall_mm")*Units::mm;
  const double straw_length = getParameterDouble("straw_length_cm")*Units::cm;
  auto tube_wall_material = getParameterMaterial("tube_wall_material");
  auto tube_inner_gas = getParameterMaterial("tube_inner_gas");

  // const double tube_rotation = getParameterDouble("tube_rotation_deg")*Units::degree;
  //const double tube_box = getParameterDouble("tube_box_mm")*Units::mm;
  //auto tube_box_fill_material = getParameterMaterial("tube_box_fill_material");


  double tube_inner_radius = tube_outer_radius - tube_wall;  //3*straw_radius; // 3 times the outer radius of the straw
  //double tube_outer_radius = tube_inner_radius + tube_wall;

  auto lv_tube = new G4LogicalVolume(new G4Tubs("TubeWall",0, tube_outer_radius, 0.5*straw_length, 0., 2*M_PI),
                                     tube_wall_material, "TubeWall");

  auto lv_empty_tube = new G4LogicalVolume(new G4Tubs("EmptyTube", 0., tube_inner_radius, 0.5*straw_length, 0., 2*M_PI),
                                           tube_inner_gas,"EmptyTube");

  place(lv_empty_tube, 0,0,0, lv_tube, G4Colour(0,1,1),-2,0,0);

  //////////----------------------------------------------------------------------------------------------------
  const double straw_wall = getParameterDouble("straw_wall_mm")*Units::mm;
  const double converter_thickness = getParameterDouble("converter_thickness_um")*Units::um;

  auto straw_wall_material = getParameterMaterial("straw_wall_material");
  auto counting_gas = getParameterMaterial("counting_gas");
  auto converter_material = getParameterMaterial("converter_material");

  const G4String name_db = "StrawWall";

  const double tubePos [7][3]={
    {0, -2*straw_radius, 0},
    {straw_radius*tan(M_PI/3.), -straw_radius, 0},
    {-straw_radius*tan(M_PI/3.),-straw_radius, 0},
    {0, 0, 0},
    {straw_radius*tan(M_PI/3.),  straw_radius, 0},
    {-straw_radius*tan(M_PI/3.), straw_radius, 0},
    {0, 2*straw_radius, 0},
  };

  for(int cpNo=0;cpNo<=60;cpNo+=10){
    // define the shape inside the G4LogicalVolume constructor, saves a line of code from declaring the shape separately
    auto lv_straw_wall = new G4LogicalVolume(new G4Tubs(name_db, 0, straw_radius, 0.5*straw_length, 0., 2*M_PI),
                                             straw_wall_material, name_db);
    auto lv_converter =  place(new G4Tubs("Converter", 0., straw_radius-straw_wall, 0.5*(straw_length-straw_wall), 0., 2*M_PI),
                               converter_material, 0,0,0, lv_straw_wall, G4Colour(0,1,1),cpNo+100,0,0).logvol;
    // this is a framework function that is overloaded, it helps the user do shape, logical volume and placement in a single line
    // look at Framework/G4/G4Interfaces/libinc/GeoBase.hh for all versions of it

    place(new G4Tubs("CountingGas", 0., straw_radius-straw_wall-converter_thickness, 0.5*(straw_length-straw_wall), 0., 2*M_PI),
          counting_gas, 0,0,0, lv_converter, G4Colour(0,0,1),cpNo+1000,0,0);

    place(lv_straw_wall,tubePos[cpNo/10][0],tubePos[cpNo/10][1],tubePos[cpNo/10][2] , lv_empty_tube, ORANGE, cpNo,0,0);
    //place(lv_straw_wall,tubePos[cpNo/10][0],tubePos[cpNo/10][1],tubePos[cpNo/10][2] , lv_empty_tube, G4Colour(0,1,1), cpNo,0,0);//good4GeomShot
  }



  //////////-------------------------------------------------------------------------------------------------------

  // create the logical volume of the straw and place it several times, this is the efficient way to do it when your geometry is repetitive
  // and especially when it is large and many volumes are defined
  // here the placement is not very smart, could be done with a loop
  // the copy numbers are not well-chosen, still unclear how we can utilize them, that's why I didn't bother with a better definition
  // that could be in a loop and calculated with the help of the running variable i for example
  // All other projects follow this rule.
  //definition of (overloaded) place function packages/Framework/G4/G4Interfaces/libsrc/GeoBase.cc
  /*
  auto lv_straw = createStrawLV();
  place(lv_straw, 0, 0, 0, lv_empty_tube, ORANGE, 0,0,0);
  place(lv_straw, 0, -2*straw_radius, 0, lv_empty_tube, ORANGE, 10,0,0);
  place(lv_straw, 0, 2*straw_radius, 0, lv_empty_tube, ORANGE, 20,0,0);
  place(lv_straw,  straw_radius*tan(M_PI/3.),  straw_radius, 0, lv_empty_tube, ORANGE, 30,0,0);
  place(lv_straw, -straw_radius*tan(M_PI/3.),  straw_radius, 0, lv_empty_tube, ORANGE, 40,0,0);
  place(lv_straw,  straw_radius*tan(M_PI/3.), -straw_radius, 0, lv_empty_tube, ORANGE, 50,0,0);
  place(lv_straw, -straw_radius*tan(M_PI/3.), -straw_radius, 0, lv_empty_tube, ORANGE, 60,0,0);
  */

  /*
  auto lv_tube_box = new G4LogicalVolume(new G4Box("TubeBox", 0.5*tube_box, 0.5*tube_box, 0.5*straw_length),
                                           tube_box_fill_material,"TubeBox");
  G4RotationMatrix* my_rotationMatrix= new G4RotationMatrix(0,0,tube_rotation);
  place(lv_tube,0,0,0,lv_tube_box,SILVER,tubeWallCopyNumber,0,my_rotationMatrix);//33?copyNum?
  */
  return lv_tube;
  // return lv_tube_box;
}

G4VPhysicalVolume* GeoBCS::Construct()
{
  // this is where we put the entire geometry together, the private functions creating the logical volumes are meant to facilitate the code below
  //Parameters:

  auto world_material = getParameterMaterial("world_material");
  //const double converter_thickness = getParameterDouble("converter_thickness_um")*Units::um;
  const double straw_radius = getParameterDouble("straw_radius_mm")*Units::mm;
  //const double straw_wall = getParameterDouble("straw_wall_mm")*Units::mm;
  const double sdd = getParameterDouble("sample_detector_distance_m")*Units::m;
  const double straw_length = getParameterDouble("straw_length_cm")*Units::cm;
  const double tube_wall = getParameterDouble("tube_wall_mm")*Units::mm;
  double tube_inner_radius = 3 * straw_radius; // 3 times the outer radius of the straw
  double tube_outer_radius = tube_inner_radius + tube_wall;
  //const double tube_box = getParameterDouble("tube_box_mm")*Units::mm;

  // calculate a value that is big enough to fit your world volume, the "super mother"
  double epsilon = 1*Units::mm;
  double big_dimension = 1.1*(straw_length + 2*tube_outer_radius + epsilon + sdd);

  //World volume:
  auto worldvols = place(new G4Box("World", big_dimension, big_dimension, big_dimension),world_material,0,0,0,0,INVISIBLE);
  auto lvWorld = worldvols.logvol;
  auto pvWorld = worldvols.physvol;

  //auto lv_tube = createTubeLV();
  //auto rot = new G4RotationMatrix();
  //rot->rotateY(M_PI/2.0);


  const double tube_rotation = getParameterDouble("tube_rotation_deg")*Units::degree; //OFFforGeoShot//
  const double tube_box = getParameterDouble("tube_box_mm")*Units::mm;
  auto tube_box_fill_material = getParameterMaterial("tube_box_fill_material"); //OFFforGeoShot//


  const int number_of_tubes = getParameterInt("number_of_tubes");
  const int number_of_panels = getParameterInt("number_of_panels");

  const double panel_offset= getParameterDouble("panel_offset_mm")*Units::mm;

  G4RotationMatrix* my_rotationMatrix= new G4RotationMatrix(0,0,tube_rotation); //OFFforGeoShot//

  double sum_offset=0.0;
  for(int j=0; j<number_of_panels; ++j){
    for(int i=0; i<number_of_tubes; ++i) {

      if(sum_offset>=tube_box) sum_offset-=tube_box; //todo add comment

      G4Transform3D trftube( G4Translate3D(0, (i*tube_box - (number_of_tubes-1)*0.5*tube_box+sum_offset), sdd + j*tube_box)
                             * G4Rotate3D(0, G4Vector3D(1,0,0))//last G4Rotate3D rotations is executed first
                             * G4Rotate3D(90.0*Units::degree, G4Vector3D(0,1,0))
                             );

      int copyNumber= 100*j+i;
      auto lv_tube = createTubeLV();
      auto lv_tube_box = new G4LogicalVolume(new G4Box("TubeBox", 0.5*tube_box, 0.5*tube_box, 0.5*straw_length), //OFFforGeoShot//
                                             tube_box_fill_material,"TubeBox"); //OFFforGeoShot//

      place(lv_tube,0,0,0,lv_tube_box,SILVER,copyNumber,0,my_rotationMatrix); //OFFforGeoShot//
      place(lv_tube_box, trftube, lvWorld, GOLD,copyNumber); //OFFforGeoShot//

      //place(lv_tube,trftube,lvWorld,SILVER,copyNumber);//goodForGeomShot

    }
    sum_offset+=panel_offset;
  }


  const double afterburner_thickness = getParameterDouble("afterburner_thickness_mm")*Units::mm;
  if( afterburner_thickness > 0.0){
    auto afterburner_material = getParameterMaterial("afterburner_material");
    auto lv_Afterburner_box = new G4LogicalVolume(new G4Box("Afterburner", 0.5*afterburner_thickness, 22*tube_box, 22*tube_box),
                                                  afterburner_material,"TubeBox");
    G4Transform3D trAfterburner( G4Translate3D(0, 0, sdd + 4.5*tube_box + 0.5*afterburner_thickness)
                                 * G4Rotate3D(0, G4Vector3D(1,0,0))//last G4Rotate3D rotations is executed first
                                 * G4Rotate3D(90.0*Units::degree, G4Vector3D(0,1,0))
                                 );

    place(lv_Afterburner_box, trAfterburner, lvWorld, BLUE, 666);

  }
  // physical volume of world
  return pvWorld;
}



bool GeoBCS::validateParameters()
{

  // you can apply conditions to control the sanity of the geometry parameters and warn the user of possible mistakes
  // a nice example: Projects/SingleCell/G4GeoSingleCell/libsrc/GeoB10SingleCell.cc

  //const double converter_thickness = getParameterDouble("converter_thickness_um")*Units::um;

  //if(converter_thickness>)
    return true;
}

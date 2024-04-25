/***********************************************************************
/
/  GRID CLASS (Adds an external potential to the existing grid potential)
/
/  written by: Elizabeth Tasker
/  date:       January, 2007
/  modified1:
/
/  PURPOSE: This is used when there exists an analytic form for the 
/           external potential rather than the acceleration field.
/           
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "phys_constants.h"

#define GINDEX(i1,i2,i3) (((i3)*GravitatingMassFieldDimension[1]+(i2))*GravitatingMassFieldDimension[0]+(i1))

int GetUnits(float *DensityUnits, float *LengthUnits,
       	      float *TemperatureUnits, float *TimeUnits,
       	      float *VelocityUnits, double *MassUnits, FLOAT Time);

int grid::AddExternalPotentialField(float *potential)
{

  /* Return if this grid is not on this processor. */

  if (MyProcessorNumber != ProcessorNumber)
   return SUCCESS;

  /* declarations */

  int i, j, k, dim, size;
  double CircularVelocity, xpos, ypos, zpos;
  double rsquared, zsquared, rcore;
  double ExternalPotential;

  /* Get unit conversions */
  
  float DensityUnits = 1 , LengthUnits = 1, TemperatureUnits, 
    TimeUnits = 1, VelocityUnits = 1, PotentialUnits = 1;
  double MassUnits = 1;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
			&TimeUnits, &VelocityUnits, &MassUnits, Time);
  PotentialUnits = pow(LengthUnits/TimeUnits,2.0);

   for (k = 0; k < GravitatingMassFieldDimension[2]; k++) {
     
     if (GridRank > 2) 
    zpos = GravitatingMassFieldLeftEdge[2] + (float(k)+0.5)*GravitatingMassFieldCellSize - 
	  ExternalGravityPosition[2];

     for (j = 0; j < GravitatingMassFieldDimension[1]; j++) {
	
      if (GridRank > 1) 
	ypos = GravitatingMassFieldLeftEdge[1] + (float(j)+0.5)*GravitatingMassFieldCellSize -
	  ExternalGravityPosition[1];
       
      for (i = 0; i < GravitatingMassFieldDimension[0]; i++) {

	xpos = GravitatingMassFieldLeftEdge[0] + (float(i)+0.5)*GravitatingMassFieldCellSize -
	  ExternalGravityPosition[0];


	if (ExternalGravity == 10) {
	  
	  /* Potential taken from `Logarithmic Potentials' in Binney+Tremaine, p46, chap2 
	     ExternalGravityConstant is the circular velocity in code units
	     ExternalGravityRadius is radius of inner region where velocity drops to zero, in code units */

	  float q = 0.7; // controls shape of potential, see B+T
	 	  
	  /* 0.5 kpc for Tasker&Tan disk model */
	  rcore = ExternalGravityRadius*LengthUnits; // [cm]

	  /* 200 km/s for Tasker&Tan disk model */
	  CircularVelocity = ExternalGravityConstant*VelocityUnits; // [CGS]
	  
	  float constant = 0.5*pow(CircularVelocity,2.0)*log(pow(rcore,-2));

	  rsquared = double(xpos*LengthUnits)*double(xpos*LengthUnits) + 
	             double(ypos*LengthUnits)*double(ypos*LengthUnits);
	  
	  zsquared = double(zpos*LengthUnits)*double(zpos*LengthUnits);
	  
	  ExternalPotential = 0.5*pow(CircularVelocity,2)*
	    log(pow(rcore,2)*pow(rcore,-2)+rsquared*pow(rcore,-2)+zsquared*pow(rcore,-2)/pow(q,2));

	}

	if (ExternalGravity == 20){

	  /* Point source potential. 
	     Duplicate to acceleration version in Grid_ComputeAccelerationFieldExternal 
	     but useful for testing  */
	    
	  rcore = max(0.1*CellWidth[0][0], ExternalGravityRadius)*LengthUnits; // from code/set units -> cgs
	  
	  rsquared = (xpos*xpos + ypos*ypos + zpos*zpos)*LengthUnits*LengthUnits;
	  double GM = ExternalGravityConstant*LengthUnits*pow(VelocityUnits,2);
	  ExternalPotential = -1.0*GM/max(rcore, sqrt(rsquared));

	}

	if (ExternalGravity == 30){
		
		/* Plummer sphere potential.
		 Potential for Plummer sphere centred at centre of box
		 Set perturber at centre of grid (here fixed for 10*10*10) */

		double eps = 0.125*LengthUnits; // kpc -> cm
		//fprintf(stderr, "old: xpos = %e, ypos = %e, zpos = %e [code]\n", xpos, ypos, zpos);
		rsquared = (xpos*xpos + ypos*ypos + zpos*zpos)*LengthUnits*LengthUnits; // cgs
		//fprintf(stderr, "rsquared = %e [cm^2]\n", rsquared);
		// Gravitational constant [cm3g-1s-2] 6.6740831e-8 [cgs]
		double G = 4*pi*GravConst; // cgs
		double M = ExternalGravityConstant*MassUnits; // g
		//fprintf(stderr, "G = %e cm^3 g-1 s-2, M = %e g\n", G, M);
		ExternalPotential = -1.0*G*M/sqrt(rsquared + eps*eps); // cgs
		//fprintf(stderr, "ExternalPotential = %e cm^2 s^-2, to be divided by PotentialUnits = %e;\n", ExternalPotential, PotentialUnits); // km^2 s^-2
		 
	}

	if (ExternalGravity == 40){
		
		/* Point mass potential */

		rsquared = (xpos*xpos + ypos*ypos + zpos*zpos)*LengthUnits*LengthUnits; // cgs
		double G = 4*pi*GravConst; // cgs
		double M = ExternalGravityConstant*MassUnits; // g
		ExternalPotential = -1.0*G*M/sqrt(rsquared); // cgs
		
	}

	// SG. PotentialUnits = pow(LengthUnits/TimeUnits,2.0);
	potential[GINDEX(i,j,k)] = float(ExternalPotential/PotentialUnits); // cm^2 s^-2 -> code units


      }
    }
  } // end: loop over grid

  return SUCCESS;
  

}

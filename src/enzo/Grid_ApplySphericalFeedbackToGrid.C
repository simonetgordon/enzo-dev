/***********************************************************************
/
/  Algorithm for applying thermal feedback to the temporary grid of an active particle
/
/  written by: John Regan
/  date:       December, 2020
/
/  note: Based on methods originally implemented by Stephen Skory
************************************************************************/
#include "preincludes.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "units.h"
#include "Fluxes.h"
#include "GridList.h"
#include "phys_constants.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "ActiveParticle.h"
#include "phys_constants.h"
#include "ActiveParticle_SmartStar.h"
#define MAX_TEMPERATURE 1e8

int grid::ApplySphericalFeedbackToGrid(ActiveParticleType** ThisParticle, float EjectaDensity,
                                       float EjectaThermalEnergyDensity, float EjectaMetalDensity,
                                       FLOAT Radius){
  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1, TimeUnits = 1, VelocityUnits = 1,
  PressureUnits = 0, GEUnits = 0, VelUnits = 0;
  double MassUnits = 1.0;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, this->ReturnTime()) == FAIL) {
      ENZO_FAIL("Error in GetUnits.");
  }
  MassUnits = DensityUnits * POW(LengthUnits,3);
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
         Vel3Num, TENum) == FAIL) {
     ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
   }
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  if (MultiSpecies) 
    if (this->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, 
				    HeIIINum, HMNum, H2INum, H2IINum, DINum, 
				    DIINum, HDINum) == FAIL) {
        ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
    }
  int SNColourNum, MetalNum, Metal2Num, MBHColourNum, Galaxy1ColourNum, 
    Galaxy2ColourNum, MetalIaNum, MetalIINum;
  int MetallicityField = FALSE;

  if (this->IdentifyColourFields(SNColourNum, Metal2Num, MetalIaNum,MetalIINum,
                                 MBHColourNum, Galaxy1ColourNum,Galaxy2ColourNum) == FAIL)
    ENZO_FAIL("Error in grid->IdentifyColourFields.\n");

  ActiveParticleType_SmartStar *SS = static_cast<ActiveParticleType_SmartStar*>(* ThisParticle);

  /* metals */
  MetalNum = max(Metal2Num, SNColourNum);
  MetallicityField = (MetalNum > 0) ? TRUE : FALSE;
  FLOAT radius = max(64*this->CellWidth[0][0], SS->InfluenceRadius); // SG. Change from InfluenceRadius.
  fprintf(stderr, "%s: radius (in cellwidths) = %f\n", __FUNCTION__, radius/this->CellWidth[0][0]);
  float MetalRadius = 1.0;
  FLOAT MetalRadius2 = radius * radius * MetalRadius * MetalRadius;

  /* particle + cell width on grid */
  float dx = float(this->CellWidth[0][0]);
  FLOAT *pos = SS->ReturnPosition();

  /* outer radius */
  FLOAT outerRadius2 = POW(Radius, 2.0); // SG. Change from 1.2*radius to BHThermalFeedbackRadius.

  /* max gas energy from max temperature = 1e8 K */
  float maxGE = MAX_TEMPERATURE / (TemperatureUnits * (Gamma-1.0) * 0.6);
  float delta_fz = 0.0;

  // Loop over all cells on grid
  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      int index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {

	    fprintf(stderr,"%s: index = %e\n", __FUNCTION__, index);
        fprintf(stderr,"%s: current density = %e (this), %e\n", __FUNCTION__, this->BaryonField[DensNum][index],
                BaryonField[DensNum][index]);
	    FLOAT radius2 = POW(CellLeftEdge[0][i] + 0.5*dx - pos[0],2.0) + POW(CellLeftEdge[1][j] + 0.5*dx - pos[1],2.0) +
                POW(CellLeftEdge[2][k] + 0.5*dx - pos[2],2.0);

        if (radius2 < outerRadius2) {
            float r1 = sqrt(radius2) / radius;
            float norm = 0.98;
            float ramp = norm*(0.5 - 0.5 * tanh(10.0*(r1-1.0)));
            /* 1/1.2^3 factor to dilute the density since we're
               depositing a uniform ejecta in a sphere of 1.2*radius
               without a ramp.  The ramp is only applied to the
               energy*density factor. */
            float factor = 0.578704;

            float OldDensity = this->BaryonField[DensNum][index];
            if(EjectaDensity > 0.0){
                BaryonField[DensNum][index] += factor*EjectaDensity;
            }
            /* Get specific energy */
            if (GENum >= 0 && DualEnergyFormalism) {

              /* When injected energy is uniform throughout the volume;
              EjectaThermalEnergyDensity in EnergyUnits/VolumeUnits */
              float oldGE =  this->BaryonField[GENum][index]; // ergs/gram - specific energy
              // BIG E = ergs,  ergs/vol = E/vol, e = specific energy = E/m_cell
              fprintf(stderr,"%s: oldGE = %e\t OldDensity = %e\n", __FUNCTION__, oldGE, OldDensity);
              float newGE = 0.0;

              if(EjectaDensity > 0.0) { /* SuperNovae */
                newGE = (OldDensity * this->BaryonField[GENum][index] +
                ramp * factor * EjectaThermalEnergyDensity * EjectaDensity) / BaryonField[DensNum][index] ;
              }
              else if (EjectaDensity == 0.0) { /* Thermal energy due to stellar luminosity */
                /* Thermal energy dump with no ejecta */
                /* For this case the EjectaThermalEnergyDensity is passed in as simply an energy  */
                newGE = EjectaThermalEnergyDensity;
              }
              else if (EjectaDensity < 0.0) {
                /* Black Hole accretion Thermal feedback */
                float cell_density, k_b, dT_max, mu;
                float dGE, dGE_max, cellmass;
                float EjectaThermalEnergy = EjectaThermalEnergyDensity; // EnergyUnits/cell
                cell_density = this->BaryonField[DensNum][index];
                cellmass = cell_density*dx*dx*dx; // from cell mass /cellvol -> cell mass in code units.

                k_b = 1.3807e-16; // cm^2 g s^-2 K^-1
                dT_max = 1e8; // K
                mu = 0.58; // SG. Fully ionised gas. Values between this and 1. Mean molecular weight, dimensionless.

                /* define changes in specific energy (ergs/g) */
                oldGE = this->BaryonField[GENum][index]; // EnergyUnits/ MassUnits
                dGE = EjectaThermalEnergy/cellmass; // EnergyUnits/MassUnits
                dGE_max = (dT_max) / (TemperatureUnits*(Gamma - 1) * mu); // EnergyUnits/MassUnits

                fprintf(stderr, "%s: EjectaThermalEnergy = %e ergs (%e code) \t cell mass = %e g (%e code)\n",
                        __FUNCTION__ , EjectaThermalEnergy*VelocityUnits*VelocityUnits,
                        EjectaThermalEnergy, cellmass*MassUnits, cellmass);
                fprintf(stderr,"%s: dGE = %"GSYM" code units, \t dGE_max = %"GSYM" code units, \t oldGE = %e code units,\t SS->EnergySaved = %e code units \n",
                          __FUNCTION__, dGE, dGE_max, oldGE, SS->EnergySaved);
                /* energy budgeting */
                if (dGE > dGE_max) { // adding to saved energy
                    SS->EnergySaved += (dGE - dGE_max);
                    dGE = dGE_max;
                } else {
                    if ((dGE + SS->EnergySaved) > dGE_max){ // taking from saved energy
                        SS->EnergySaved -= (dGE_max - dGE);
                        dGE = dGE_max;
                    } else{ // using remainder of saved energy
                        dGE += SS->EnergySaved;
                        SS->EnergySaved = 0.0;
                    }
                }
                newGE = oldGE + dGE; // ergs/g
                fprintf(stderr,"%s: dGE = %"GSYM" ergs/g, \t dGE_max = %"GSYM" ergs/g, \t newGE = %e ergs/g,\t SS->EnergySaved = %e ergs/g \n",
                        __FUNCTION__, dGE, dGE_max, newGE, SS->EnergySaved);

              } // END EjectaDensity < 0.0 (BH thermal feedback scheme)

              newGE = min(newGE, maxGE);
              fprintf(stderr,"%s: oldGE = %"GSYM"\t newGE = %"GSYM"\t maxGE = %e ergs/g \n", __FUNCTION__,
                     oldGE, newGE, maxGE);

              this->BaryonField[TENum][index] = newGE;

              for (int dim = 0; dim < GridRank; dim++)
                  this->BaryonField[TENum][index] += 0.5 * this->BaryonField[Vel1Num+dim][index] *
                          this->BaryonField[Vel1Num+dim][index];

              fprintf(stderr, "%s: In DualEnergy formalism Increase in GE energy is %e percent.\n", __FUNCTION__,
                        (newGE - oldGE)*100.0/oldGE);
            } else {
                float newGE = 0.0;
                if(EjectaDensity > 0.0) {
                    newGE = (OldDensity * this->BaryonField[TENum][index] +
                            ramp * factor * EjectaDensity * EjectaThermalEnergyDensity) / BaryonField[DensNum][index];
                } else if (EjectaDensity == 0.0) { /* Thermal energy from luminosity */
                    newGE = EjectaThermalEnergyDensity;
                } else if (EjectaDensity < 0.0) {
                    /* Black Hole accretion Thermal feedback */
                    float cellmass = this->BaryonField[DensNum][index]*dx*dx*dx;
                    newGE = this->BaryonField[GENum][index] + EjectaThermalEnergyDensity / cellmass;
              }

              newGE = min(newGE, maxGE);
              this->BaryonField[TENum][index] = newGE;

            } //end if/else (GENum >= 0 && DualEnergyFormalism)

            /* Update species and colour fields */
            if (MetallicityField == TRUE && radius2 <= MetalRadius2)
              delta_fz = EjectaMetalDensity / OldDensity;
            else
              delta_fz = 0.0;
            float increase = this->BaryonField[DensNum][index] / OldDensity - delta_fz;
            fprintf(stderr, "increase = %e , delta_fz = %"GSYM"\n", increase, delta_fz);
            if (MultiSpecies) {
              BaryonField[DeNum][index] *= increase;
              BaryonField[HINum][index] *= increase;
              BaryonField[HIINum][index] *= increase;
              BaryonField[HeINum][index] *= increase;
              BaryonField[HeIINum][index] *= increase;
              BaryonField[HeIIINum][index] *= increase;
            }
            if (MultiSpecies > 1) {
              BaryonField[HMNum][index] *= increase;
              BaryonField[H2INum][index] *= increase;
              BaryonField[H2IINum][index] *= increase;
            }
            if (MultiSpecies > 2) {
              BaryonField[DINum][index] *= increase;
              BaryonField[DIINum][index] *= increase;
              BaryonField[HDINum][index] *= increase;
            }

            fprintf(stderr, "HIINum after *= increase = %"GSYM"\n", this->BaryonField[HIINum][index]);
            if (MetallicityField == TRUE)
              BaryonField[MetalNum][index] += EjectaMetalDensity;

            /* MBHColour injected */
            if (MBHColourNum > 0)
              BaryonField[MBHColourNum][index] += factor*EjectaDensity;

            } // END if inside radius
          } // END i-direction
        } // END j-direction
     } // END k-direction
  return SUCCESS;
}

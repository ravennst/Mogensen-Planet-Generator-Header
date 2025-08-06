/*
 * planet.h - Procedural Planet Terrain Generation Header
 * 
 * This header provides the core logic for generating smooth, recursive, 
 * tetrahedron-based terrain heightmaps suitable for spherical worlds.
 * It is derived from Torben AE. Mogensen's original "planet.c" algorithm,
 * but is self-contained and simplified for integration into modern codebases.
 * 
 * Features:
 * ---------
 * - Recursive height calculation based on spatial subdivision of a tetrahedron
 * - Deterministic pseudorandom noise using custom 2D random function (rand2)
 * - Approximate rain shadow estimation based on terrain slope and direction
 * - Adjustable constants for terrain tuning (sea level, noise weight, etc.)
 * 
 * Requirements:
 * -------------
 * The main program must:
 * - Define the external variables:
 *     double M;                        // Sea level reference height (e.g., -0.02)
 *     planet_vertex tetra[4];          // Global tetrahedron vertices
 *     double rseed;                    // Global seed for initialization
 * 
 * - Call `initialize_vertices(rseed)` once before calling `planet()`
 * 
 * Usage:
 * ------
 * Call the `planet()` function like this:
 * 
 *     planet_out result = planet(tetra[0], tetra[1], tetra[2], tetra[3], x, y, z, level);
 * 
 * Where:
 *     - x, y, z: The 3D position on the sphere's surface (normalized or projected)
 *     - level: The recursive subdivision depth (e.g., 20–30 for high detail)
 *     - result.h: The terrain height at that point
 *     - result.shadow: An estimate of the local rain shadow effect (0.0–1.0)
 * 
 * Notes:
 * ------
 * - All recursive calculations are performed inline and statelessly.
 * - To modify biome generation or color mapping, build on the `planet_out` struct.
 * - This implementation is pure C (not C++) and has no dynamic memory requirements.
 * 
 * - Conversion of Lat/Long to X,Y,Z uses the folloiwng equations. They don't need to be
 *   placed as a function within this header due to their simplicity but are essential
 *   for the functions within the header to work correctly. 
 *   The conversion from lat/long to x,y,z assumes a unit sphere.
 *   They are:
 *      x = cos(lat) * sin(long);
 *      y = sin(lat);
 *      z = cos(lat) * cos(long);
 *   The inverse of those function to convert cartesean x,y,z to lat/long is:
 *      lat  = asin(y);
 *      long = atan2(x, z);
 * 
 * Author: Adapted and extended by ravennst, based on the work of Torben AE. Mogensen.
 * License: Open-source or permissive license as appropriate.
 */

/* todo
 * add function to calculate biomes, temperature, rainfall, and colors derived from planet0() and readcolors () from the original planet.c file.
*/

#ifndef PLANET_TERRAIN_GEN
#define PLANET_TERRAIN_GEN

#include <math.h>
#include <stdint.h>

// Global variables used both within and outside of planet.h
//
extern double M;                     // sea level
extern double rseed;                 // The main seed need for generation

// calculation variables internal to planet.h
//
#define PLANET_PI 3.141592653589793
#define PLANET_POWA 1.0
#define PLANET_POW 0.47
#define PLANET_dd1 0.45
#define PLANET_dd2 0.035
#define PLANET_shade_angle 45.0
#define COLORMAP_SIZE 65536

/* T = tundra, G = grasslands, B = Taiga / boreal forest, D = desert,
   S = savanna, F = temperate forest, R = temperate rainforest,
   W = Xeric shrubland and dry forest, E = tropical dry forest,
   O = tropical rainforest, I = icecap */

/* Whittaker diagram */
static const char biomes[45][46] = { 
  // the original planet.c file had this as [45][45] but that is incorrect as 
  // it doesn't consider the null terminator which the compliler silently 
  // truncates in the planet.c code but not the header.
		       "IIITTTTTGGGGGGGGDDDDDDDDDDDDDDDDDDDDDDDDDDDDD",
		       "IIITTTTTGGGGGGGGDDDDGGDSDDSDDDDDDDDDDDDDDDDDD",
		       "IITTTTTTTTTBGGGGGGGGGGGSSSSSSDDDDDDDDDDDDDDDD",
		       "IITTTTTTTTBBBBBBGGGGGGGSSSSSSSSSWWWWWWWDDDDDD",
		       "IITTTTTTTTBBBBBBGGGGGGGSSSSSSSSSSWWWWWWWWWWDD",
		       "IIITTTTTTTBBBBBBFGGGGGGSSSSSSSSSSSWWWWWWWWWWW",
		       "IIIITTTTTTBBBBBBFFGGGGGSSSSSSSSSSSWWWWWWWWWWW",
		       "IIIIITTTTTBBBBBBFFFFGGGSSSSSSSSSSSWWWWWWWWWWW",
		       "IIIIITTTTTBBBBBBBFFFFGGGSSSSSSSSSSSWWWWWWWWWW",
		       "IIIIIITTTTBBBBBBBFFFFFFGGGSSSSSSSSWWWWWWWWWWW",
		       "IIIIIIITTTBBBBBBBFFFFFFFFGGGSSSSSSWWWWWWWWWWW",
		       "IIIIIIIITTBBBBBBBFFFFFFFFFFGGSSSSSWWWWWWWWWWW",
		       "IIIIIIIIITBBBBBBBFFFFFFFFFFFFFSSSSWWWWWWWWWWW",
		       "IIIIIIIIIITBBBBBBFFFFFFFFFFFFFFFSSEEEWWWWWWWW",
		       "IIIIIIIIIITBBBBBBFFFFFFFFFFFFFFFFFFEEEEEEWWWW",
		       "IIIIIIIIIIIBBBBBBFFFFFFFFFFFFFFFFFFEEEEEEEEWW",
		       "IIIIIIIIIIIBBBBBBRFFFFFFFFFFFFFFFFFEEEEEEEEEE",
		       "IIIIIIIIIIIIBBBBBBRFFFFFFFFFFFFFFFFEEEEEEEEEE",
		       "IIIIIIIIIIIIIBBBBBRRRFFFFFFFFFFFFFFEEEEEEEEEE",
		       "IIIIIIIIIIIIIIIBBBRRRRRFFFFFFFFFFFFEEEEEEEEEE",
		       "IIIIIIIIIIIIIIIIIBRRRRRRRFFFFFFFFFFEEEEEEEEEE",
		       "IIIIIIIIIIIIIIIIIRRRRRRRRRRFFFFFFFFEEEEEEEEEE",
		       "IIIIIIIIIIIIIIIIIIRRRRRRRRRRRRFFFFFEEEEEEEEEE",
		       "IIIIIIIIIIIIIIIIIIIRRRRRRRRRRRRRFRREEEEEEEEEE",
		       "IIIIIIIIIIIIIIIIIIIIIRRRRRRRRRRRRRRRREEEEEEEE",
		       "IIIIIIIIIIIIIIIIIIIIIIIRRRRRRRRRRRRRROOEEEEEE",
		       "IIIIIIIIIIIIIIIIIIIIIIIIRRRRRRRRRRRROOOOOEEEE",
		       "IIIIIIIIIIIIIIIIIIIIIIIIIIRRRRRRRRRROOOOOOEEE",
		       "IIIIIIIIIIIIIIIIIIIIIIIIIIIRRRRRRRRROOOOOOOEE",
		       "IIIIIIIIIIIIIIIIIIIIIIIIIIIIRRRRRRRROOOOOOOEE",
		       "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIRRRRRRROOOOOOOOE",
		       "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIRRRRROOOOOOOOOO",
		       "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIRROOOOOOOOOOO",
		       "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIROOOOOOOOOOO",
		       "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIROOOOOOOOOOO",
		       "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIOOOOOOOOOOO",
		       "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIOOOOOOOOOO",
		       "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIOOOOOOOOO",
		       "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIOOOOOOOOO",
		       "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIOOOOOOOO",
		       "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIOOOOOOOO",
		       "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIOOOOOOOO",
		       "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIOOOOOOOO",
		       "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIOOOOOOO",
		       "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIOOOOOOO"
};

// Structured variables
//
typedef struct planet_vertex
{
  double h; /* altitude */
  double s; /* seed */
  double x,y,z; /* coordinates */
  double shadow; /* approximate rain shadow */
} planet_vertex;

extern planet_vertex tetra[4];       // global tetrahedron vertices

typedef struct planet_climate
{
  double temp;
  double rainfall;
  char biome_type;
} planet_climate;

typedef struct planet_out
{
  double h; /* altitude */
  double shadow; /* approximate rain shadow */
} planet_out;

typedef struct RGB // structure for the color output.
{
  uint8_t r, g, b;
} RGB;

typedef struct colormap {
	int nocols; // number of entries
	uint8_t rtable[COLORMAP_SIZE];
	uint8_t gtable[COLORMAP_SIZE];
	uint8_t btable[COLORMAP_SIZE];
} colormap;

// functions internal to terrain generation, not needed outside of planet.h
//
/* distance squared between vertices */
static inline double dist2(planet_vertex a, planet_vertex b)
{
  double abx, aby, abz;
  abx = a.x-b.x; aby = a.y-b.y; abz = a.z-b.z;
  return abx*abx+aby*aby+abz*abz;
}

static inline double rand2(double p, double q) {/* random number generator taking two seeds */
/* rand2(p,q) = rand2(q,p) is important     */
  double r = (p+3.14159265)*(q+3.14159265);
  return(2.*(r-(int)r)-1.);
}

// functions used outside of planet.h
//
static inline void initialize_vertices()
{
double r1, r2, r3, r4;
/* For the vertices of the tetrahedron */

  /* initialize vertices to slightly irregular tetrahedron */
  tetra[0].x = -sqrt(3.0)-0.20;
  tetra[0].y = -sqrt(3.0)-0.22;
  tetra[0].z = -sqrt(3.0)-0.23;

  tetra[1].x = -sqrt(3.0)-0.19;
  tetra[1].y = sqrt(3.0)+0.18;
  tetra[1].z = sqrt(3.0)+0.17;

  tetra[2].x = sqrt(3.0)+0.21;
  tetra[2].y = -sqrt(3.0)-0.24;
  tetra[2].z = sqrt(3.0)+0.15;

  tetra[3].x = sqrt(3.0)+0.24;
  tetra[3].y = sqrt(3.0)+0.22;
  tetra[3].z = -sqrt(3.0)-0.25;

  r1 = rseed;

  r1 = rand2(r1,r1);
  r2 = rand2(r1,r1);
  r3 = rand2(r1,r2);
  r4 = rand2(r2,r3);

  tetra[0].s = r1;
  tetra[1].s = r2;
  tetra[2].s = r3;
  tetra[3].s = r4;

  tetra[0].h = M; // NOTE! need to input or define M
  tetra[1].h = M;
  tetra[2].h = M;
  tetra[3].h = M;

  tetra[0].shadow = 0.0;
  tetra[1].shadow = 0.0;
  tetra[2].shadow = 0.0;
  tetra[3].shadow = 0.0;

}

/* Main function to calculate altitude and rainshadow. */
static inline planet_out planet(planet_vertex a, planet_vertex b, planet_vertex c, planet_vertex d, double x, double y, double z, int level) 
// vertex a,b,c,d;             /* tetrahedron vertices, recursion loops will use different values but the initial loop will use tetra[] before midpoints are calculated */
// double x,y,z;               /* goal point */
// int level;                  /* levels to go */
{
  planet_vertex e;
  double lab, lac, lad, lbc, lbd, lcd, maxlength;
  double es1, es2, es3;
  double eax,eay,eaz, epx,epy,epz;
  double ecx,ecy,ecz, edx,edy,edz;
  double x1,y1,z1,x2,y2,z2,l1,tmp;
  
  if (level>0) {

    /* make sure ab is longest edge */
    lab = dist2(a,b);
    lac = dist2(a,c);
    lad = dist2(a,d);
    lbc = dist2(b,c);
    lbd = dist2(b,d);
    lcd = dist2(c,d);

    maxlength = lab;
    if (lac > maxlength) maxlength = lac;
    if (lad > maxlength) maxlength = lad;
    if (lbc > maxlength) maxlength = lbc;
    if (lbd > maxlength) maxlength = lbd;
    if (lcd > maxlength) maxlength = lcd;

    if (lac == maxlength) return(planet(a,c,b,d, x,y,z, level));
    if (lad == maxlength) return(planet(a,d,b,c, x,y,z, level));
    if (lbc == maxlength) return(planet(b,c,a,d, x,y,z, level));
    if (lbd == maxlength) return(planet(b,d,a,c, x,y,z, level));
    if (lcd == maxlength) return(planet(c,d,a,b, x,y,z, level));

    /* ab is longest, so cut ab */
      e.s = rand2(a.s,b.s);
      es1 = rand2(e.s,e.s);
      es2 = 0.5+0.1*rand2(es1,es1);  /* find cut point */
      es3 = 1.0-es2;

      if (a.s<b.s) {
        e.x = es2*a.x+es3*b.x; e.y = es2*a.y+es3*b.y; e.z = es2*a.z+es3*b.z;
      } else if (a.s>b.s) {
        e.x = es3*a.x+es2*b.x; e.y = es3*a.y+es2*b.y; e.z = es3*a.z+es2*b.z;
      } else { /* as==bs, very unlikely to ever happen */
        e.x = 0.5*a.x+0.5*b.x; e.y = 0.5*a.y+0.5*b.y; e.z = 0.5*a.z+0.5*b.z;
      }

      /* new altitude is: */
        if (lab>1.0) lab = pow(lab,0.5);
        /* decrease contribution for very long distances */
        e.h = 0.5*(a.h+b.h) /* average of end points */
          + e.s*PLANET_dd1*pow(fabs(a.h-b.h),PLANET_POWA)
          /* plus contribution for altitude diff */
          + es1*PLANET_dd2*pow(lab,PLANET_POW); /* plus contribution for distance */

      /* calculate approximate rain shadow for new point */
      if (e.h <= 0.0 ) e.shadow = 0.0;
      else {
      x1 = 0.5*(a.x+b.x);
      x1 = a.h*(x1-a.x)+b.h*(x1-b.x);
      y1 = 0.5*(a.y+b.y);
      y1 = a.h*(y1-a.y)+b.h*(y1-b.y);
      z1 = 0.5*(a.z+b.z);
      z1 = a.h*(z1-a.z)+b.h*(z1-b.z);
      l1 = sqrt(x1*x1+y1*y1+z1*z1);
      if (l1==0.0) l1 = 1.0;
      tmp = sqrt(1.0-y*y);
      if (tmp<0.0001) tmp = 0.0001;
      x2 = x*x1+y*y1+z*z1;
      z2 = -z/tmp*x1+x/tmp*z1;
      if (lab > 0.04)
	e.shadow = (a.shadow + b.shadow- cos(PLANET_PI*PLANET_shade_angle/180.0)*z2/l1)/3.0;
      else
	e.shadow = (a.shadow + b.shadow)/2.0;
      }
      
      /* find out in which new tetrahedron target point is */
      eax = a.x-e.x; eay = a.y-e.y; eaz = a.z-e.z;
      ecx = c.x-e.x; ecy = c.y-e.y; ecz = c.z-e.z;
      edx = d.x-e.x; edy = d.y-e.y; edz = d.z-e.z;
      epx =   x-e.x; epy =   y-e.y; epz =   z-e.z;
      if ((eax*ecy*edz+eay*ecz*edx+eaz*ecx*edy
           -eaz*ecy*edx-eay*ecx*edz-eax*ecz*edy)*
          (epx*ecy*edz+epy*ecz*edx+epz*ecx*edy
           -epz*ecy*edx-epy*ecx*edz-epx*ecz*edy)>0.0) {
        /* point is inside acde */
        return(planet(c,d,a,e, x,y,z, level-1));
      } else {
        /* point is inside bcde */
        return(planet(c,d,b,e, x,y,z, level-1));
      } 
  
    } else { /* level == 0 */    
    // rainShadow  = 0.25*(a.shadow+b.shadow+c.shadow+d.shadow);
    // return 0.25*(a.h+b.h+c.h+d.h);
    planet_out result;
    result.h = 0.25 * (a.h + b.h +c.h + d.h);
    result.shadow = 0.25 * (a.shadow + b.shadow + c.shadow + d.shadow);
    return result;
  }
}

/* Main function to calculate the climate of a point given its latitude, altitude, and rain shadow */
static inline planet_climate planet_climatecalc(double y, double altitude, double rshadow)
// variable "y" : used by the equation to determine latitude of the point. Assumes a point on the surface of a unit sphere.
// variable "altitude" : self explanatory. Output from result.h of planet()
// variable "rshadow" : Rain Shadow, output from result.shadow of planet()
{
  planet_climate result;
  result.temp = (1.0 - y*y) * (1.0 - 0.5 * altitude);
  result.rainfall = (1.0 - rshadow) * (1.0 - 0.25 * altitude);
  result.biome_type = biomes[(int)(result.rainfall * 44.99)][(int)(result.temp * 44.99)];
  return result;
}

/* Function to convert the biome type output by planet_climatecalc into RGB values */

static inline RGB planet_get_biome_rgb(char biome_char) { // input is normally result.biome_type from planet_climatecalc(). Variable name may vary.
  switch (biome_char) {
    case 'T': return (RGB){210, 210, 210}; // Tundra
    case 'G': return (RGB){250, 215, 165}; // Grassland
    case 'B': return (RGB){105, 155, 120}; // Boreal forest (Taiga)
    case 'D': return (RGB){220, 195, 175}; // Desert
    case 'S': return (RGB){225, 155, 100}; // Savanna
    case 'F': return (RGB){155, 215, 170}; // Temperate forest
    case 'R': return (RGB){170, 195, 200}; // Temperate rainforest
    case 'W': return (RGB){185, 150, 160}; // Xeric shrubland
    case 'E': return (RGB){130, 190, 25 }; // Tropical dry forest
    case 'O': return (RGB){110, 160, 170}; // Tropical rainforest
    case 'I': return (RGB){255, 255, 255}; // Icecap
    default: return (RGB){128, 128, 128}; // Fallback for unknown biome
  }
}



/* Function to calculate the expected color of a point based on both biome and a *.col file. */
/* Usage: You will need the following in the main code in order to be albe to use this function

**************
colormap cmap; // Declare a colormap instance

FILE* f = fopen("colorfile.col", "r"); // colorfile.col can be any .col file you wish to use per the original Mogensen program.
if (!f || !colormap_load_from_file(f, &cmap)) {
	fprintf(stderr, "Failed to load colormap from .col file\n");
 	exit(1);
  }
if (cmap->nocols <= 0 || cmap->nocols > COLORMAP_SIZE) {
	fprintf(stderr, "ERROR: Invalid colormap size (%d). Must be between 1 and %d.\n", cmap->nocols, COLORMAP_SIZE);
 	return 0; // failure
  }
  fclose(f);
**************
This function takes the output of planet() as an input and returns an RGB out for colors.

*/
static inline RGB planet_colormap_rgb(planet_out result, const colormap* cmap) 
{
	int LOWEST = 6;
	int HIGHEST = (cmap->nocols > 0) ? cmap->nocols - 1 : 9;
	int SEA = (HIGHEST + LOWEST) / 2;
	int LAND = SEA + 1;

	double alt = result.h - M; // M is the code-defined custom sea level. Default = -0.02
	int color_index;

	if (alt <= 0.0) {
	   color_index = SEA + (int)((SEA - LOWEST + 1) * (10.0 * alt)); // I hope I got this math right
	   if (color_index < LOWEST) color_index = LOWEST;
	} else {
	   if (alt >= 0.1) {
		color_index = HIGHEST;
	   } else {
		color_index = LAND + (int)((HIGHEST - LAND + 1) * (10.0 * alt));  // I hope I got this math right
		if (color_index > HIGHEST) color_index = HIGHEST;
	   }
	}

	RGB color;
	color.r = cmap->rtable[color_index];
	color.g = cmap->gtable[color_index];
	color.b = cmap->btable[color_index];
	return color;
}



#endif

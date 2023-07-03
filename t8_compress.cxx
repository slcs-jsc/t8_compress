#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_vec.h>

/* ------------------------------------------------------------
   Constants...
   ------------------------------------------------------------ */

/* Maximum number of x values. */
#define EX 2048

/* Maximum number of y values. */
#define EY 1024

/* Maximum length of ASCII strings. */
#define LEN 4096

/* ------------------------------------------------------------
   Macros...
   ------------------------------------------------------------ */

/*! Print error message and quit program. */
#define ERRMSG(text) {				\
    fprintf(stderr, "%s", text);		\
    exit(1);					\
  }

/*! Compute linear interpolation. */
#define LIN(x0, y0, x1, y1, x)			\
  ((y0)+((y1)-(y0))/((x1)-(x0))*((x)-(x0)))

/* ------------------------------------------------------------
   Declarations...
   ------------------------------------------------------------ */

/*! Find array index for irregular grid. */
int locate_irr(
	       double *xx,
	       int n,
	       double x);

/* User data for adaptation callback. */
struct t8_adapt_data {
  double midpoint[3];                  /* The midpoint of our sphere. */
  double refine_if_inside_radius;      /* if an element's center is smaller than this value, we refine the element. */
  double coarsen_if_outside_radius;    /* if an element's center is larger this value, we coarsen its family. */
};

/* Data per element. */
struct t8_data_per_element {
  int level;
  double volume;
  double z;
};

/* Create element data. */
static struct t8_data_per_element *t8_create_element_data (t8_forest_t forest,
							   double x[EX],
							   double y[EY],
							   double z[EX][EY],
							   int nx,
							   int ny);

/* Function used to adapt the forest. */
t8_forest_t t8_adapt_forest (t8_forest_t forest);

/* Callback function for adapting the forest. */
int t8_adapt_callback (t8_forest_t forest,
		       t8_forest_t forest_from,
		       t8_locidx_t which_tree,
		       t8_locidx_t lelement_id,
		       t8_eclass_scheme_c *ts,
		       const int is_family,
		       const int num_elements,
		       t8_element_t *elements[]);

/* Write data to vtu file. */
static void t8_output_data_to_vtu (t8_forest_t forest,
				   struct t8_data_per_element *data,
				   const char *prefix);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main (int argc, char **argv) {
  
  /* Check arguments... */
  if(argc!=3)
    ERRMSG("usage: give parameters <data.tab> <maxlev>\n");
  
  /* Set variables... */
  int level=atoi(argv[2]);
  
  /* Read data file... */
  FILE *in;
  if (!(in = fopen(argv[1], "r")))
    ERRMSG("Cannot open file!");

  /* Read data... */
  char line[LEN];
  static double rx, ry, ry_old = -1e99, rz, x[EX], y[EY], z[EX][EY];
  int nx = 1, ny = -1;
  while (fgets(line, LEN, in))
    if (sscanf(line, "%lg %lg %lg", &rx, &ry, &rz) == 3) {
      if(ry != ry_old) {
	if ((++ny) >= EY)
	  ERRMSG("Too many y values!");
	nx = -1;
	ry_old = ry;
      }
      if ((++nx) >= EX)
	ERRMSG("Too many x values!");
      x[nx]=rx;
      y[ny]=ry;
      z[nx][ny]=rz;
    }
  if ((++nx) >= EX)
    ERRMSG("Too many x values!");
  if ((++ny) >= EY)
    ERRMSG("Too many y values!");
  
  /* Close file... */
  fclose(in);
  
  /* Initialize. */
  SC_CHECK_MPI (sc_MPI_Init (&argc, &argv));
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_PRODUCTION);
  
  /* Use MPI_COMM_WORLD as a communicator. */
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
  
  /* Build a coarse 2-D mesh of a single quad. */
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, comm, 0, 0, 0);
  
  /* Create uniform forest. */
  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), level, 0, comm);
  
  /* Build data array and gather data for the local elements. */
  t8_data_per_element *data;
  data = t8_create_element_data (forest, x, y, z, nx, ny);
  
  /* Write output. */
  t8_output_data_to_vtu (forest, data, "forest_with_data");
  
  /* Adapt forest. */
  forest = t8_adapt_forest (forest);
  
  /* Write output. */
  t8_forest_write_vtk (forest, "forest_adapted");

  /* Free the data array. */
  T8_FREE (data);
  
  /* Destroy the forest. */
  t8_forest_unref (&forest);
  
  /* Finalize... */
  sc_finalize ();
  SC_CHECK_MPI (sc_MPI_Finalize ());
  
  return 0;
}

/* ------------------------------------------------------------
   Definitions...
   ------------------------------------------------------------ */

t8_forest_t t8_adapt_forest (t8_forest_t forest) {
  
  t8_forest_t         forest_adapt;
  
  struct t8_adapt_data adapt_data = {
    {0.5, 0.5, 0},              /* Midpoints of the sphere. */
    0.2,                        /* Refine if inside this radius. */
    0.4                         /* Coarsen if outside this radius. */
  };
  
  /* Check that forest is a committed, that is valid and usable, forest. */
  T8_ASSERT (t8_forest_is_committed (forest));
  
  /* Create a new forest that is adapted from a forest with our adaptation callback. */
  forest_adapt = t8_forest_new_adapt (forest, t8_adapt_callback, 0, 0, &adapt_data);
  
  return forest_adapt;
}

/*****************************************************************/

int t8_adapt_callback (t8_forest_t forest,
		       t8_forest_t forest_from,
		       t8_locidx_t which_tree,
		       t8_locidx_t lelement_id,
		       t8_eclass_scheme_c *ts,
		       const int is_family,
		       const int num_elements,
		       t8_element_t *elements[]) {
  
  /* In t8_adapt_forest we pass a t8_adapt_data pointer as user data to the
   * t8_forest_new_adapt function. This pointer is stored as the used data of the new forest
   * and we can now access it with t8_forest_get_user_data (forest). */
  const struct t8_adapt_data *adapt_data =
    (const struct t8_adapt_data *) t8_forest_get_user_data (forest);
  T8_ASSERT (adapt_data != NULL);
  
  /* Compute the element's centroid coordinates. */
  double centroid[3];
  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);
  
  /* Compute the distance to our sphere midpoint. */
  double dist = t8_vec_dist (centroid, adapt_data->midpoint);
  if (dist < adapt_data->refine_if_inside_radius) {
    
    /* Refine this element. */
    return 1;
  }
  else if (is_family && dist > adapt_data->coarsen_if_outside_radius) {
    
    /* Coarsen this family. Note that we check for is_family before, since returning < 0
     * if we do not have a family as input is illegal. */
    return -1;
  }
  
  /* Do not change this element. */
  return 0;
}

/*****************************************************************/

static struct t8_data_per_element *t8_create_element_data (t8_forest_t forest,
							   double x[EX],
							   double y[EY],
							   double z[EX][EY],
							   int nx,
							   int ny) {
  
  t8_locidx_t         num_local_elements;
  t8_locidx_t         num_ghost_elements;
  struct t8_data_per_element *element_data;

  /* Check that forest is committed, that is valid and usable, forest. */
  T8_ASSERT (t8_forest_is_committed (forest));

  /* Get the number of local elements of the forest. */
  num_local_elements = t8_forest_get_local_num_elements (forest);
  
  /* Get the number of ghost elements of the forest. */
  num_ghost_elements = t8_forest_get_num_ghosts (forest);

  /* Alloc... */
  element_data = T8_ALLOC (struct t8_data_per_element,
			   num_local_elements + num_ghost_elements);

  /* Let us now fill the data with something.
   * For this, we iterate through all trees and for each tree through all its elements, calling
   * t8_forest_get_element_in_tree to get a pointer to the current element.
   * This is the recommended and most performant way.
   * An alternative is to iterate over the number of local elements and use
   * t8_forest_get_element. However, this function needs to perform a binary search
   * for the element and the tree it is in, while t8_forest_get_element_in_tree has a
   * constant look up time. You should only use t8_forest_get_element if you do not know
   * in which tree an element is.
   */
  t8_locidx_t         itree, num_local_trees;
  t8_locidx_t         current_index;
  t8_locidx_t         ielement, num_elements_in_tree;
  t8_eclass_t         tree_class;
  t8_eclass_scheme_c *eclass_scheme;
  const t8_element_t *element;
  
  /* Get the number of trees that have elements of this process. */
  num_local_trees = t8_forest_get_num_local_trees (forest);
  
  /* This loop iterates through all local trees in the forest. */
  for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
    
    /* Each tree may have a different element class (quad/tri/hex/tet etc.) and therefore
     * also a different way to interpret its elements. In order to be able to handle elements
     * of a tree, we need to get its eclass_scheme, and in order to so we first get its eclass. */
    tree_class = t8_forest_get_tree_class (forest, itree);
    eclass_scheme = t8_forest_get_eclass_scheme (forest, tree_class);
    
    /* Get the number of elements of this tree. */
    num_elements_in_tree = t8_forest_get_tree_num_elements (forest, itree);
    
    /* This loop iterates through all the local elements of the forest in the current tree. */      
    for (ielement = 0; ielement < num_elements_in_tree;
	 ++ielement, ++current_index) {
      
      /* We can now write to the position current_index into our array in order to store
       * data for this element. */
      /* Since in this example we want to compute the data based on the element in question,
       * we need to get a pointer to this element. */
      element = t8_forest_get_element_in_tree (forest, itree, ielement);
      
      /* We want to store the elements level and its volume as data. We compute these
       * via the eclass_scheme and the forest_element interface. */
      element_data[current_index].level =
	eclass_scheme->t8_element_level (element);
      element_data[current_index].volume =
	t8_forest_element_volume (forest, itree, element);

      /* Get element centroid... */
      double centroid[3];
      t8_forest_element_centroid (forest, itree, element, centroid);

      /* Interpolate data... */
      double gx = x[0] + (x[nx-1] - x[0]) * centroid[0];
      double gy = y[0] + (y[ny-1] - y[0]) * centroid[1];

      int ix = locate_irr(x, nx, gx);
      int iy = locate_irr(y, ny, gy);

      double y0 = LIN(x[ix], z[ix][iy], x[ix+1], z[ix+1][iy], gx);
      double y1 = LIN(x[ix], z[ix][iy+1], x[ix+1], z[ix+1][iy+1], gx);
      
      element_data[current_index].z = LIN(y[iy], y0, y[iy+1], y1, gy);
    }
  }
  
  return element_data;
}

/*****************************************************************/

static void t8_output_data_to_vtu (t8_forest_t forest,
				   struct t8_data_per_element *data,
				   const char *prefix) {
  
  t8_locidx_t num_elements = t8_forest_get_local_num_elements (forest);
  t8_locidx_t ielem;
  
  /* We need to allocate a new array to store the volumes on their own.
   * This array has one entry per local element. */
  double *element_data = T8_ALLOC (double, num_elements);
  
  /* The number of user defined data fields to write. */
  int num_data = 1;
  
  /* For each user defined data field we need one t8_vtk_data_field_t variable */
  t8_vtk_data_field_t vtk_data;
  
  /* Set the type of this variable. Since we have one value per element, we pick T8_VTK_SCALAR */
  vtk_data.type = T8_VTK_SCALAR;
  
  /* The name of the field as should be written to the file. */
  strcpy (vtk_data.description, "Element data");
  vtk_data.data = element_data;
  
  /* Copy the elment's volumes from our data array to the output array. */
  for (ielem = 0; ielem < num_elements; ++ielem)
    element_data[ielem] = data[ielem].z;
  
  /* To write user defined data, we need to extended output function t8_forest_vtk_write_file
   * from t8_forest_vtk.h. Despite writin user data, it also offers more control over which 
   * properties of the forest to write. */
  int                 write_treeid = 1;
  int                 write_mpirank = 1;
  int                 write_level = 1;
  int                 write_element_id = 1;
  int                 write_ghosts = 0;
  t8_forest_write_vtk_ext (forest, prefix, write_treeid, write_mpirank,
			   write_level, write_element_id, write_ghosts,
			   0, 0, num_data, &vtk_data);
  
  T8_FREE (element_data);
}

/*****************************************************************************/

int locate_irr(
  double *xx,
  int n,
  double x) {

  int ilo = 0;
  int ihi = n - 1;
  int i = (ihi + ilo) >> 1;

  if (xx[i] < xx[i + 1])
    while (ihi > ilo + 1) {
      i = (ihi + ilo) >> 1;
      if (xx[i] > x)
	ihi = i;
      else
	ilo = i;
  } else
    while (ihi > ilo + 1) {
      i = (ihi + ilo) >> 1;
      if (xx[i] <= x)
	ihi = i;
      else
	ilo = i;
    }

  return ilo;
}

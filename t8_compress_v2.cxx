#include <iostream>
#include <t8.h>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_vec.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include "t8_forest_private.h"
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_io.h>
#include "t8_cmesh/t8_cmesh_testcases.h"
#include <t8_forest/t8_forest_partition.h>

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

typedef struct element_data_t {
  double values;
} element_data_t;

typedef struct adapt_data_t {
  double midpoint[3];		/* The midpoint of our sphere. */
  double refine_if_inside_radius;	/* if an element's center is smaller than this value, we refine the element. */
  double coarsen_if_outside_radius;	/* if an element's center is larger this value, we coarsen its family. */
  sc_array_t *element_data;
} adapt_data_t;

static void t8_element_set_value(
  const adapt_data_t * adapt_data,
  t8_locidx_t ielement,
  double value) {
  
  element_data_t elem_data;
  elem_data.values = value;
  *((element_data_t *)
    t8_sc_array_index_locidx(adapt_data->element_data, ielement)) = elem_data;
}

static void t8_element_set_element(
  const adapt_data_t * adapt_data,
  t8_locidx_t ielement,
  element_data_t element) {
  
  *((element_data_t *)
    t8_sc_array_index_locidx(adapt_data->element_data, ielement)) = element;
}

static element_data_t t8_element_get_value(
  const adapt_data_t * adapt_data,
  t8_locidx_t ielement) {
  
  return *((element_data_t *)
	   t8_sc_array_index_locidx((adapt_data->element_data), ielement));
}

int adapt_callback(
  t8_forest_t forest,
  t8_forest_t forest_from,
  t8_locidx_t which_tree,
  t8_locidx_t lelement_id,
  t8_eclass_scheme_c * ts,
  const int is_family,
  const int num_elements,
  t8_element_t * elements[]) {

  double centroid[3];

  const struct adapt_data_t *adapt_data =
    (const struct adapt_data_t *) t8_forest_get_user_data(forest);

  T8_ASSERT(adapt_data != NULL);
  
  t8_forest_element_centroid(forest_from, which_tree, elements[0], centroid);

  double dist = t8_vec_dist(centroid, adapt_data->midpoint);
  
  if (dist < adapt_data->refine_if_inside_radius) {
    
    /* Refine element. */
    return 1;
  } else if (is_family && dist > adapt_data->coarsen_if_outside_radius) {
    
    /* Coarsen element. */
    return -1;
  }
  
  /* Do not change element. */
  return 0;
}

t8_forest_t t8_adapt_forest(
  t8_forest_t forest_from,
  t8_forest_adapt_t adapt_fn,
  int do_partition,
  int recursive,
  void *user_data) {

  /* Init new forest */
  t8_forest_t forest_new;
  t8_forest_init(&forest_new);
  
  /* Adapt the forest */
  t8_forest_set_adapt(forest_new, forest_from, adapt_fn, recursive);
  
  /* Set user data for the adapted forest */
  if (user_data != NULL)
    t8_forest_set_user_data(forest_new, user_data);
  
  /* Commit the adapted forest */
  t8_forest_commit(forest_new);
  
  return forest_new;
}

void t8_forest_replace(
  t8_forest_t forest_old,
  t8_forest_t forest_new,
  t8_locidx_t which_tree,
  t8_eclass_scheme_c * ts,
  int refine,
  int num_outgoing,
  t8_locidx_t first_outgoing,
  int num_incoming,
  t8_locidx_t first_incoming) {

  struct adapt_data_t *adapt_data_new =
    (struct adapt_data_t *) t8_forest_get_user_data(forest_new);
  const struct adapt_data_t *adapt_data_old =
    (const struct adapt_data_t *) t8_forest_get_user_data(forest_old);

  /* get the index of the data array corresponding to the old and the adapted forest */
  first_incoming += t8_forest_get_tree_element_offset(forest_new, which_tree);
  first_outgoing += t8_forest_get_tree_element_offset(forest_old, which_tree);

  /* Do not adapt or coarsen */
  if (refine == 0)
    t8_element_set_element(adapt_data_new, first_incoming,
			   t8_element_get_value(adapt_data_old,
						first_outgoing));
  
  /* The old element is refined, we copy the element values */
  else if (refine == 1)
    for (int i = 0; i < num_incoming; i++)
      t8_element_set_element(adapt_data_new, first_incoming + i,
			     t8_element_get_value(adapt_data_old,
						  first_outgoing));
  
  /* Old element is coarsened */
  else if (refine == -1) {
    double tmp_value = 0;
    for (t8_locidx_t i = 0; i < num_outgoing; i++)
      tmp_value +=
	t8_element_get_value(adapt_data_old, first_outgoing + i).values;
    
    t8_element_set_value(adapt_data_new, first_incoming,
			 tmp_value / num_outgoing);
  }
  t8_forest_set_user_data(forest_new, adapt_data_new);
}

static void t8_write_vtu(
  t8_forest_t forest,
  struct adapt_data_t *data,
  const char *prefix) {

  const t8_locidx_t num_elements = t8_forest_get_local_num_elements(forest);

  double *element_data = T8_ALLOC(double, num_elements);

  int num_data = 1;
  
  t8_vtk_data_field_t vtk_data;
  
  vtk_data.type = T8_VTK_SCALAR;
  strcpy(vtk_data.description, "Element own data");
  vtk_data.data = element_data;
  
  for (t8_locidx_t ielem = 0; ielem < num_elements; ++ielem) {
    element_data[ielem] = t8_element_get_value(data, ielem).values;
  }

  int write_treeid = 1;
  int write_mpirank = 1;
  int write_level = 1;
  int write_element_id = 1;
  int write_ghosts = 0;
  t8_forest_write_vtk_ext(forest, prefix, write_treeid, write_mpirank,
			  write_level, write_element_id, write_ghosts,
			  0, 0, num_data, &vtk_data);
  T8_FREE(element_data);
}

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char **argv) {
  
  /* Check arguments... */
  if (argc != 3)
    ERRMSG("usage: give parameters <data.tab> <maxlev>\n");

  /* Set variables... */
  int level = atoi(argv[2]);

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
      if (ry != ry_old) {
	if ((++ny) >= EY)
	  ERRMSG("Too many y values!");
	nx = -1;
	ry_old = ry;
      }
      if ((++nx) >= EX)
	ERRMSG("Too many x values!");
      x[nx] = rx;
      y[ny] = ry;
      z[nx][ny] = rz;
    }
  if ((++nx) >= EX)
    ERRMSG("Too many x values!");
  if ((++ny) >= EY)
    ERRMSG("Too many y values!");

  /* Close file... */
  fclose(in);
  
  /* Initialize... */
  int mpiret = sc_MPI_Init(&argc, &argv);
  SC_CHECK_MPI(mpiret);
  sc_init(sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_PRODUCTION);
  t8_init(SC_LP_PRODUCTION);
  
  t8_forest_t forest_adapt;
  element_data_t *elem_data;
  adapt_data_t *data;
  
  /* Construct a cmesh */
  t8_cmesh_t cmesh =
    t8_cmesh_new_hypercube(T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, 0, 0, 0);
  
  /* Construct a forest with one tree */
  t8_forest_t forest =
    t8_forest_new_uniform(cmesh, t8_scheme_new_default_cxx(), level, 0, sc_MPI_COMM_WORLD);

  /* Build initial data array and set data for the local elements. */
  data = T8_ALLOC(adapt_data_t, 1);
  elem_data = T8_ALLOC(element_data_t, 1);
  data->element_data =
    sc_array_new_count(sizeof(element_data_t),
		       t8_forest_get_local_num_elements(forest));

  const t8_locidx_t num_trees = t8_forest_get_num_local_trees(forest);
  
  for (int itree = 0, ielem = 0; itree < num_trees; itree++) {
    
    const t8_locidx_t num_elem =
      t8_forest_get_tree_num_elements(forest, itree);
    
    /* Inner loop: Iteration over the elements of the local tree */
    for (t8_locidx_t ielem_tree = 0; ielem_tree < num_elem; ielem_tree++, ielem++) {
      
      /* To calculate the distance to the centroid of an element the element is saved */
      const t8_element_t *element =
	t8_forest_get_element_in_tree(forest, itree, ielem_tree);

      /* Get the centroid of the local element. */
      double centroid[3];
      t8_forest_element_centroid(forest, itree, element, centroid);
      
      /* Interpolate data... */
      double gx = x[0] + (x[nx - 1] - x[0]) * centroid[0];
      double gy = y[0] + (y[ny - 1] - y[0]) * centroid[1];
      
      int ix = locate_irr(x, nx, gx);
      int iy = locate_irr(y, ny, gy);
      
      double y0 = LIN(x[ix], z[ix][iy], x[ix + 1], z[ix + 1][iy], gx);
      double y1 = LIN(x[ix], z[ix][iy + 1], x[ix + 1], z[ix + 1][iy + 1], gx);
      
      /* Calculation of the distance to the centroid for the referenced element */
      elem_data->values = LIN(y[iy], y0, y[iy + 1], y1, gy);

      t8_element_set_element(data, ielem, *elem_data);
    }
  }

  /*  Set the data elements which will be set as user elements on the forest */
  data->midpoint[0] = 0.5;
  data->midpoint[1] = 0.5;
  data->midpoint[2] = 0;
  data->refine_if_inside_radius = 0.2;
  data->coarsen_if_outside_radius = 0.4;

  /* Set the user data (values for callback and element values) */
  t8_forest_set_user_data(forest, data);

  /* Write vtu file */
  t8_write_vtu(forest, data, "t8_step7_uniform_forest");

  /* Build a second forest to store the adapted forest - keep the old one */
  t8_forest_ref(forest);

  /* Adapt the forest correponding tho the callback function (distance to the centroid) */
  forest_adapt = t8_adapt_forest(forest, adapt_callback, 0, 0, data);
  /* Calculate/Interpolate the data array for the adapted forest */

  /* Create user_data element for the adapted forest */
  struct adapt_data_t *adapt_data = T8_ALLOC(adapt_data_t, 1);
  adapt_data->element_data =
    sc_array_new_count(sizeof(element_data_t),
		       t8_forest_get_local_num_elements(forest_adapt));

  /* Initialize user_data element for the adapted forest */
  (*adapt_data).midpoint[0] = (*data).midpoint[0];
  (*adapt_data).midpoint[1] = (*data).midpoint[1];
  (*adapt_data).midpoint[2] = (*data).midpoint[2];
  (*adapt_data).refine_if_inside_radius = (*data).refine_if_inside_radius;
  (*adapt_data).coarsen_if_outside_radius = (*data).coarsen_if_outside_radius;

  t8_forest_set_user_data(forest_adapt, adapt_data);
  t8_forest_iterate_replace(forest_adapt, forest, t8_forest_replace);

  /* Write the adapted forest to a vtu file */
  adapt_data =
    (struct adapt_data_t *) t8_forest_get_user_data(forest_adapt);
  t8_write_vtu(forest_adapt, adapt_data, "t8_step7_adapt_forest");

  /* Free the memory */
  sc_array_destroy(adapt_data->element_data);
  sc_array_destroy(data->element_data);

  /* Save the new forest as old forest */
  t8_forest_unref(&forest);
  forest = forest_adapt;
  *data = *adapt_data;
  t8_forest_unref(&forest_adapt);

  /* Now you could continue working with the forest. */
  T8_FREE(data);
  T8_FREE(adapt_data);
  T8_FREE(elem_data);
  
  /* Finalize... */
  sc_finalize();
  mpiret = sc_MPI_Finalize();
  SC_CHECK_MPI(mpiret);

  return 0;
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

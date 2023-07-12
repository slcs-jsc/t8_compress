#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_iterate.h>
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
  double midpoint[3];		// TODO: use refinement criterion of example for now, replace by gradients/curvature later
  double refine_if_inside_radius;
  double coarsen_if_outside_radius;
  sc_array_t *element_data;
};

/* Data per element. */
struct t8_data_per_element {
  int level;			// TODO: remove level and volume later, we only want "z" as data per element
  double volume;
  double z;
};

/* Create element data. */
static struct t8_data_per_element *t8_create_element_data(
  t8_forest_t forest,
  double x[EX],
  double y[EY],
  double z[EX][EY],
  int nx,
  int ny);

/* Function used to adapt the forest. */
t8_forest_t t8_adapt_forest(
  t8_forest_t forest_from,
  t8_forest_adapt_t adapt_fn,
  int do_partition,
  int recursive,
  void *user_data);

/* Callback function for adapting the forest. */
int t8_adapt_callback(
  t8_forest_t forest,
  t8_forest_t forest_from,
  t8_locidx_t which_tree,
  t8_locidx_t lelement_id,
  t8_eclass_scheme_c * ts,
  const int is_family,
  const int num_elements,
  t8_element_t * elements[]);

/* Write data to vtu file. */
static void t8_output_data_to_vtu(
  t8_forest_t forest,
  struct t8_data_per_element *data,
  const char *prefix);

/* Replace callback to decide how to interpolate a refined or coarsened element. */
void t8_forest_replace(
  t8_forest_t forest_old,
  t8_forest_t forest_new,
  t8_locidx_t which_tree,
  t8_eclass_scheme_c * ts,
  int refine,
  int num_outgoing,
  t8_locidx_t first_outgoing,
  int num_incoming,
  t8_locidx_t first_incoming);

/* Set the value of an element to a given entry */
static void t8_element_set_value(
  const t8_adapt_data * adapt_data,
  t8_locidx_t ielement,
  double value);

/* Set the value of an element to a given entry */
static void t8_element_set_element(
  const t8_adapt_data * adapt_data,
  t8_locidx_t ielement,
  t8_data_per_element element);

/* Get the value of an element to a given entry */
static t8_data_per_element t8_element_get_value(
  const t8_adapt_data * adapt_data,
  t8_locidx_t ielement);

/* Write output to vtu */
static void t8_write_vtu(
  t8_forest_t forest,
  struct t8_adapt_data *data,
  const char *prefix);

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

  /* Initialize. */
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
  SC_CHECK_MPI(sc_MPI_Init(&argc, &argv));
  sc_init(comm, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init(SC_LP_PRODUCTION);

  /* Build a coarse 2-D mesh of a single quad. */
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube(T8_ECLASS_QUAD, comm, 0, 0, 0);

  /* Create uniform forest. */
  t8_forest_t forest =
    t8_forest_new_uniform(cmesh, t8_scheme_new_default_cxx(), level, 0, comm);

  /* Build data array and gather data for the local elements. */
  t8_adapt_data *data;
  data = T8_ALLOC(t8_adapt_data, 1);
  data->midpoint[0] = 0.5;
  data->midpoint[1] = 0.5;
  data->midpoint[2] = 0.0;
  data->refine_if_inside_radius = 0.2;
  data->coarsen_if_outside_radius = 0.4;

  t8_data_per_element *elem_data;
  elem_data = T8_ALLOC(t8_data_per_element, 1);
  elem_data = t8_create_element_data(forest, x, y, z, nx, ny);
  data->element_data =
    sc_array_new_count(sizeof(t8_data_per_element),
		       t8_forest_get_local_num_elements(forest));

  //printf("  num_elements: %d\n", t8_forest_get_local_num_elements(forest));

  /* Write output. */
  t8_output_data_to_vtu(forest, elem_data, "forest_with_data");

  /* Build a second forest to store the adapted forest - keep the old one */
  t8_forest_ref(forest);

  /* Adapt forest. */
  t8_forest_t forest_adapt =
    t8_adapt_forest(forest, t8_adapt_callback, 0, 0, data);

  // printf("  num_elements: %d\n", t8_forest_get_local_num_elements(forest_adapt));

  /* Write output (the forest was adapted, but the data was not adapted, yet... */
  t8_output_data_to_vtu(forest_adapt, elem_data, "forest_adapted");

  /* Build data array for the adapted forest */
  struct t8_adapt_data *adapt_data = T8_ALLOC(t8_adapt_data, 1);
  adapt_data->element_data =
    sc_array_new_count(sizeof(t8_data_per_element),
		       t8_forest_get_local_num_elements(forest_adapt));

  /* Set data for the adapted forest */
  (*adapt_data).midpoint[0] = (*data).midpoint[0];
  (*adapt_data).midpoint[1] = (*data).midpoint[1];
  (*adapt_data).midpoint[2] = (*data).midpoint[2];
  (*adapt_data).refine_if_inside_radius = (*data).refine_if_inside_radius;
  (*adapt_data).coarsen_if_outside_radius = (*data).coarsen_if_outside_radius;

  t8_forest_set_user_data(forest_adapt, adapt_data);

  printf("hihi4\n");


  // TODO: Seg fault happens at this point...
  t8_forest_iterate_replace(forest_adapt, forest, t8_forest_replace);


  printf("hihi5\n");

  /* Write the adapted forest to a vtu file */
  adapt_data = (struct t8_adapt_data *) t8_forest_get_user_data(forest_adapt);
  t8_write_vtu(forest_adapt, adapt_data, "forest_adapted_with_data");

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


//  /* Destroy the forest. */
//  t8_forest_unref (&forest);


  /* Finalize... */
  sc_finalize();
  SC_CHECK_MPI(sc_MPI_Finalize());

  return 0;
}

/* ------------------------------------------------------------
   Definitions...
   ------------------------------------------------------------ */

t8_forest_t t8_adapt_forest(
  t8_forest_t forest_from,
  t8_forest_adapt_t adapt_fn,
  int do_partition,
  int recursive,
  void *user_data) {

  /* Init new forest... */
  t8_forest_t forest_new;
  t8_forest_init(&forest_new);

  /* Adapt the forest... */
  t8_forest_set_adapt(forest_new, forest_from, adapt_fn, recursive);

  /* Set user data for the adapted forest... */
  if (user_data != NULL)
    t8_forest_set_user_data(forest_new, user_data);

  /* Commit the adapted forest... */
  t8_forest_commit(forest_new);

  return forest_new;
}

/*****************************************************************/

int t8_adapt_callback(
  t8_forest_t forest,
  t8_forest_t forest_from,
  t8_locidx_t which_tree,
  t8_locidx_t lelement_id,
  t8_eclass_scheme_c * ts,
  const int is_family,
  const int num_elements,
  t8_element_t * elements[]) {

  /* Our adaptation criterion is to look at the midpoint coordinates of the current element and if
   * they are inside a sphere around a given midpoint we refine, if they are outside, we coarsen. */

  /* In t8_step3_adapt_forest we pass a t8_step3_adapt_data pointer as user data to the
   * t8_forest_new_adapt function. This pointer is stored as the used data of the new forest
   * and we can now access it with t8_forest_get_user_data (forest). */
  const struct t8_adapt_data *adapt_data =
    (const struct t8_adapt_data *) t8_forest_get_user_data(forest);

  /* You can use T8_ASSERT for assertions that are active in debug mode (when configured with --enable-debug).
   * If the condition is not true, then the code will abort.
   * In this case, we want to make sure that we actually did set a user pointer to forest and thus
   * did not get the NULL pointer from t8_forest_get_user_data.
   */
  T8_ASSERT(adapt_data != NULL);

  /* Compute the element's centroid coordinates. */
  double centroid[3];
  t8_forest_element_centroid(forest_from, which_tree, elements[0], centroid);

  /* Compute the distance to our sphere midpoint. */
  double dist = t8_vec_dist(centroid, adapt_data->midpoint);
  if (dist < adapt_data->refine_if_inside_radius) {

    /* Refine this element. */
    return 1;
  } else if (is_family && dist > adapt_data->coarsen_if_outside_radius) {

    /* Coarsen this family. Note that we check for is_family before, since returning < 0
     * if we do not have a family as input is illegal. */
    return -1;
  }

  /* Do not change this element. */
  return 0;
}

/*****************************************************************/

static struct t8_data_per_element *t8_create_element_data(
  t8_forest_t forest,
  double x[EX],
  double y[EY],
  double z[EX][EY],
  int nx,
  int ny) {

  /* Check that forest is committed, that is valid and usable, forest. */
  T8_ASSERT(t8_forest_is_committed(forest));

  /* Get the number of local elements of the forest. */
  t8_locidx_t num_local_elements = t8_forest_get_local_num_elements(forest);

  /* Get the number of ghost elements of the forest. */
  t8_locidx_t num_ghost_elements = t8_forest_get_num_ghosts(forest);

  /* Alloc... */
  struct t8_data_per_element *element_data;
  element_data = T8_ALLOC(struct t8_data_per_element,
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

  /* Get the number of trees that have elements of this process. */
  t8_locidx_t num_local_trees = t8_forest_get_num_local_trees(forest);

  /* This loop iterates through all local trees in the forest. */
  for (t8_locidx_t itree = 0, current_index = 0; itree < num_local_trees;
       ++itree) {

    /* Each tree may have a different element class (quad/tri/hex/tet etc.) and therefore
     * also a different way to interpret its elements. In order to be able to handle elements
     * of a tree, we need to get its eclass_scheme, and in order to so we first get its eclass. */
    t8_eclass_t tree_class = t8_forest_get_tree_class(forest, itree);
    t8_eclass_scheme_c *eclass_scheme;
    eclass_scheme = t8_forest_get_eclass_scheme(forest, tree_class);

    /* Get the number of elements of this tree. */
    t8_locidx_t num_elements_in_tree =
      t8_forest_get_tree_num_elements(forest, itree);

    /* This loop iterates through all the local elements of the forest in the current tree. */
    for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree;
	 ++ielement, ++current_index) {

      /* We can now write to the position current_index into our array in order to store
       * data for this element. */
      /* Since in this example we want to compute the data based on the element in question,
       * we need to get a pointer to this element. */
      const t8_element_t *element;
      element = t8_forest_get_element_in_tree(forest, itree, ielement);

      /* We want to store the elements level and its volume as data. We compute these
       * via the eclass_scheme and the forest_element interface. */
      element_data[current_index].level =
	eclass_scheme->t8_element_level(element);
      element_data[current_index].volume =
	t8_forest_element_volume(forest, itree, element);

      /* Get element centroid... */
      double centroid[3];
      t8_forest_element_centroid(forest, itree, element, centroid);

      /* Interpolate data... */
      double gx = x[0] + (x[nx - 1] - x[0]) * centroid[0];
      double gy = y[0] + (y[ny - 1] - y[0]) * centroid[1];

      int ix = locate_irr(x, nx, gx);
      int iy = locate_irr(y, ny, gy);

      double y0 = LIN(x[ix], z[ix][iy], x[ix + 1], z[ix + 1][iy], gx);
      double y1 = LIN(x[ix], z[ix][iy + 1], x[ix + 1], z[ix + 1][iy + 1], gx);

      element_data[current_index].z = LIN(y[iy], y0, y[iy + 1], y1, gy);
    }
  }

  return element_data;
}

/*****************************************************************/

static void t8_output_data_to_vtu(
  t8_forest_t forest,
  struct t8_data_per_element *data,
  const char *prefix) {

  t8_locidx_t num_elements = t8_forest_get_local_num_elements(forest);
  t8_locidx_t ielem;

  /* We need to allocate a new array to store the volumes on their own.
   * This array has one entry per local element. */
  double *element_data = T8_ALLOC(double, num_elements);

  /* The number of user defined data fields to write. */
  int num_data = 1;

  /* For each user defined data field we need one t8_vtk_data_field_t variable */
  t8_vtk_data_field_t vtk_data;

  /* Set the type of this variable. Since we have one value per element, we pick T8_VTK_SCALAR */
  vtk_data.type = T8_VTK_SCALAR;

  /* The name of the field as should be written to the file. */
  strcpy(vtk_data.description, "Element data");
  vtk_data.data = element_data;

  /* Copy the elment's volumes from our data array to the output array. */
  for (ielem = 0; ielem < num_elements; ++ielem)
    element_data[ielem] = data[ielem].z;

  /* To write user defined data, we need to extended output function t8_forest_vtk_write_file
   * from t8_forest_vtk.h. Despite writin user data, it also offers more control over which 
   * properties of the forest to write. */
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

/*****************************************************************************/

/* Replace callback to decide how to interpolate a refined or coarsened element.
 * If an element is refined, each child gets the value of its parent.
 * If elements are coarsened, the parent gets the average value of the children.
 * Outgoing are the old elements and incoming the new ones 
 * \param [in] forest_old        non adapted forest
 * \param [in] forest_new        adapted forest
 * \param [in] which_tree        tree_id of the analyzed element
 * \param [in] ts                eclass sheme 
 * \param [in] refine            ==0 - do nothing, == -1 - coarsen, == 1 - refine
 * \param [in] num_outgoing      number of the elements not refined forest
 * \param [in] first_outgoing    index of the old element
 * \param [in] num_incoming      number of the elements corresponding to the element of the not refined forest
 * \param [in] first_incoming    index of the new element
 */
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

  struct t8_adapt_data *adapt_data_new =
    (struct t8_adapt_data *) t8_forest_get_user_data(forest_new);

  const struct t8_adapt_data *adapt_data_old =
    (const struct t8_adapt_data *) t8_forest_get_user_data(forest_old);

  /* get the index of the data array corresponding to the old and the adapted forest */
  first_incoming += t8_forest_get_tree_element_offset(forest_new, which_tree);
  first_outgoing += t8_forest_get_tree_element_offset(forest_old, which_tree);

  /* Do not adapt or coarsen */
  if (refine == 0) {
    t8_element_set_element(adapt_data_new, first_incoming,
			   t8_element_get_value(adapt_data_old,
						first_outgoing));
  }

  /* The old element is refined, we copy the element values */
  else if (refine == 1) {
    for (int i = 0; i < num_incoming; i++) {
      t8_element_set_element(adapt_data_new, first_incoming + i,
			     t8_element_get_value(adapt_data_old,
						  first_outgoing));
    }
  }

  /* Old element is coarsened */
  else if (refine == -1) {
    double tmp_value = 0;
    for (t8_locidx_t i = 0; i < num_outgoing; i++) {
      tmp_value += t8_element_get_value(adapt_data_old, first_outgoing + i).z;
    }
    t8_element_set_value(adapt_data_new, first_incoming,
			 tmp_value / num_outgoing);
  }
  t8_forest_set_user_data(forest_new, adapt_data_new);
}

/*****************************************************************************/

/* Set the value of an element to a given entry */
static void t8_element_set_value(
  const t8_adapt_data * adapt_data,
  t8_locidx_t ielement,
  double value) {

  t8_data_per_element elem_data;
  elem_data.z = value;
  *((t8_data_per_element *)
    t8_sc_array_index_locidx(adapt_data->element_data, ielement)) = elem_data;
}

/*****************************************************************************/

/* Set the value of an element to a given entry */
static void t8_element_set_element(
  const t8_adapt_data * adapt_data,
  t8_locidx_t ielement,
  t8_data_per_element element) {

  *((t8_data_per_element *)
    t8_sc_array_index_locidx(adapt_data->element_data, ielement)) = element;
}

/*****************************************************************************/

/* Get the value of an element to a given entry */
static t8_data_per_element t8_element_get_value(
  const t8_adapt_data * adapt_data,
  t8_locidx_t ielement) {

  return *((t8_data_per_element *)
	   t8_sc_array_index_locidx((adapt_data->element_data), ielement));
}

/*****************************************************************************/

static void t8_write_vtu(
  t8_forest_t forest,
  struct t8_adapt_data *data,
  const char *prefix) {

  const t8_locidx_t num_elements = t8_forest_get_local_num_elements(forest);
  /* We need to allocate a new array to store the volumes on their own.
   * This array has one entry per local element. */

  double *element_data = T8_ALLOC(double, num_elements);

  /* The number of user defined data fields to write. */
  int num_data = 1;
  t8_vtk_data_field_t vtk_data;
  vtk_data.type = T8_VTK_SCALAR;
  strcpy(vtk_data.description, "Element own data");
  vtk_data.data = element_data;

  /* Copy the elment's data from the data array to the output array. */
  for (t8_locidx_t ielem = 0; ielem < num_elements; ++ielem) {
    element_data[ielem] = t8_element_get_value(data, ielem).z;
  }

  /* To write user defined data, we need to extended output function t8_forest_vtk_write_file
   * from t8_forest_vtk.h. Despite writin user data, it also offers more control over which 
   * properties of the forest to write. */
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

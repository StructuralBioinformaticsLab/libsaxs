#ifndef _CLUSTER_H_
#define _CLUSTER_H_

#include "saxs_config.h"

typedef struct ftclusters_t
{
  struct ftresult *ftresults;
  int nftresults;

  int nclusters;
  int *cluster_centers;
  int *cluster_sizes;
  int **cluster_members;

  double radius;
} ftclusters;

typedef struct cluster_stats_t
{
  // For online variance calculations
  double Mk;
  double Sk;
  double variance;

  int nmembers;
  double average_chi;
  double max_chi;
  double min_chi;

  int best_fti;
  double best_rmsd;
  double orig_rmsd;
} cluster_stats;

ftclusters *read_clusters(char *file_path);
void destroy_clusters(ftclusters *clstrs);

cluster_stats *create_cluster_stats();
void init_cluster_stats(cluster_stats *stats);

#endif /* _CLUSTER_H_ */

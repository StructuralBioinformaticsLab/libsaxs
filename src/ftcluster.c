#include "ftcluster.h"

#include "saxs_common.h"

#include <stdbool.h>
#include <float.h>
#include <math.h>

ftclusters *read_clusters(char *file_path)
{
    FILE *cfp = fopen(file_path, "r");
    ftclusters *clstrs = (ftclusters *) malloc(sizeof(ftclusters));
    bool error = false;

    char *line = (char *) malloc(512 * sizeof(char));
    size_t len;

    if (getline(&line, &len, cfp) == -1 ||
        sscanf(line, "%*s %lf", &(clstrs->radius)) != 1)
    {
        error = true;
    }

    int max_clusters = 32;
    int max_cluster_size = 64;
    clstrs->cluster_centers = (int *) malloc(max_clusters * sizeof(int));
    clstrs->cluster_sizes = (int *) malloc(max_clusters * sizeof(int));
    clstrs->cluster_members = (int **) malloc(max_clusters * sizeof(int *));
    memset(clstrs->cluster_sizes, 0, max_clusters * sizeof(int));
    memset(clstrs->cluster_members, 0, max_clusters * sizeof(int));

    int ftindex;
    int nclusters = 0;
    int nmembers = 0;
    while(getline(&line, &len, cfp) != -1 && !error)
    {
        if (!strncmp(line, "Center", 6))
        {
            // repack previous cluster
            if (nclusters != 0)
            {
                clstrs->cluster_sizes[nclusters - 1] = nmembers;
                clstrs->cluster_members[nclusters - 1] =
                    (int *) realloc(clstrs->cluster_members[nclusters - 1],
                                    nmembers * sizeof(int));
            }

            // initialize current cluster
            if (sscanf(line, "%*s %d", &ftindex) != 1)
            {
                error = true;
                continue;
            }
            ++nclusters;

            if (nclusters > max_clusters)
            {
                max_clusters *= 2;
                clstrs->cluster_centers = (int *) realloc(clstrs->cluster_centers,
                                                          max_clusters * sizeof(int));
                clstrs->cluster_sizes = (int *) realloc(clstrs->cluster_sizes,
                                                        max_clusters * sizeof(int));
                int **new_cluster_members = (int **) calloc(max_clusters, sizeof(int *));
                memcpy(new_cluster_members, clstrs->cluster_members,
                       nclusters * sizeof(int *));
                free(clstrs->cluster_members);
                clstrs->cluster_members = new_cluster_members;
            }

            nmembers = 0;
            max_cluster_size = 64;
            clstrs->cluster_centers[nclusters - 1] = ftindex;
            clstrs->cluster_members[nclusters - 1] =
                malloc(max_cluster_size * sizeof(int));
        }
        else
        {
            if (sscanf(line, "%d", &ftindex) != 1)
            {
                error = true;
                continue;
            }
            ++nmembers;

            if (nmembers > max_cluster_size)
            {
                max_cluster_size *= 2;
                clstrs->cluster_members[nclusters - 1] =
                    realloc(clstrs->cluster_members[nclusters - 1],
                            max_cluster_size * sizeof(int *));
            }
            clstrs->cluster_members[nclusters - 1][nmembers - 1] = ftindex;
        }
    }
    free(line);
    fclose(cfp);

    if (error)
    {
#ifdef DEBUG
        fprintf(stderr, "Error reading cluster file: %s\n", file_path);
        exit(1);
#else
        clstrs->nclusters = nclusters;
        destroy_clusters(clstrs);
        return NULL;
#endif
    }
    clstrs->cluster_sizes[nclusters - 1] = nmembers;
    clstrs->nclusters = nclusters;
    clstrs->cluster_centers = (int *) realloc(clstrs->cluster_centers,
                                              nclusters * sizeof(int));
    clstrs->cluster_sizes = (int *) realloc(clstrs->cluster_sizes,
                                            nclusters * sizeof(int));
    clstrs->cluster_members = (int **) realloc(clstrs->cluster_members,
                                               nclusters * sizeof(int *));

    return clstrs;
}

void destroy_clusters(ftclusters *clstrs)
{
    if (clstrs != NULL)
    {
        if (clstrs->cluster_centers != NULL)
            free(clstrs->cluster_centers);
        if (clstrs->cluster_sizes != NULL)
            free(clstrs->cluster_sizes);
        if (clstrs->cluster_members != NULL)
        {
            for (int i = 0; i < clstrs->nclusters; ++i)
            {
                if (clstrs->cluster_members[i] != NULL)
                    free(clstrs->cluster_members[i]);
            }
        }

        free(clstrs);
    }
}

cluster_stats *create_cluster_stats()
{
    cluster_stats *clust_stats = (cluster_stats *) calloc(1, sizeof(cluster_stats));
    init_cluster_stats(clust_stats);
    return clust_stats;
}

void init_cluster_stats(cluster_stats *stats)
{
    stats->max_chi = 0.0;
    stats->min_chi = DBL_MAX;
    stats->Mk = -1;
    stats->variance = NAN;
}

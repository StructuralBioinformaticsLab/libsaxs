#include "saxs_profile.h"
#include "saxs_utils.h"
#include "saxs_atom.h"

#include "mol.0.0.6/pdb.h"
#include "mol.0.0.6/utils.h"

#include <getopt.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

void print_usage(char *exe_name)
{
    fprintf(stderr,
            "Usage: %s MAPPING_PRM ATMPRM PDB_FILE OUTPUT\n",
            basename(exe_name));
    fprintf(stderr,
            "  Options:\n");
    fprintf(stderr,
            "    --q_start=<n>      Default 0.0\n");
    fprintf(stderr,
            "    --q_end=<n>        Default 0.5\n");
    fprintf(stderr,
            "    --q_step=<n>       Default 0.005\n");
    fprintf(stderr,
            "    --add_noise        Default false\n");
    fprintf(stderr,
            "    --noise_variance=n Default 1.0, in percent\n");
    fprintf(stderr,
            "    --msur_k=<k>       Set msur_k (default 1.0)\n");
    fprintf(stderr,
            "    --no_hydrogen      Remove hydrogens before generating profile\n");
    fprintf(stderr,
            "    --rand_seed=<n>    Seed the random number generator using n.\n");
}

// Source: http://c-faq.com/lib/gaussian.html
double gaussrand(double mean, double stddev) {
    static double U, V;
    static int phase = 0;
    double Z;

    if (phase == 0) {
        U = (rand() + 1.0) / (RAND_MAX + 2.0);
        V = rand() / (RAND_MAX + 1.0);
        Z = sqrt(-2 * log(U)) * sin(2 * M_PI * V);
    } else {
        Z = sqrt(-2 * log(U)) * cos(2 * M_PI * V);
    }

    phase = 1 - phase;
    return (Z * stddev) + mean;
}

int main(int argc, char *argv[]) {
    char *mapping_file_path;
    char *prms_file_path;
    char *pdb_file_path;
    char *output_file_path;

    float msur_k = 1.0;

    double q_start = 0.0;
    double q_end = 0.5;
    double q_step = 0.005;

    bool add_noise = false;
    double variance_percent = 1.0;

    bool no_hydrogens = false;

    static struct option long_options[] = {
        {"msur_k", required_argument, NULL, 0},
        {"q_start", required_argument, NULL, 0},
        {"q_end", required_argument, NULL, 0},
        {"q_step", required_argument, NULL, 0},
        {"add_noise", no_argument, NULL, 0},
        {"noise_variance", required_argument, NULL, 0},
        {"no-hydrogens", no_argument, NULL, 0},
        {"rand_seed", required_argument, NULL, 0},
        {NULL, 0, 0, 0}
    };

    srand((int) time(NULL));

    int option_index;
    int c;
    const char *long_opt_name;

    while(1) {
        c = getopt_long(argc, argv, "h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {
        case 0:
            long_opt_name = long_options[option_index].name;
            if (strcmp("msur_k", long_opt_name) == 0) {
                msur_k = atof(optarg);
            } else if (strcmp("q_start", long_opt_name) == 0) {
                q_start = atof(optarg);
            } else if (strcmp("q_end", long_opt_name) == 0) {
                q_end = atof(optarg);
            } else if (strcmp("q_step", long_opt_name) == 0) {
                q_step = atof(optarg);
            } else if (strcmp("add_noise", long_opt_name) == 0) {
                add_noise = true;
            } else if (strcmp("noise_variance", long_opt_name) == 0) {
                variance_percent = atof(optarg);
            } else if (strcmp("no-hydrogens", long_opt_name) == 0) {
                no_hydrogens = true;
            } else if (strcmp("rand_seed", long_opt_name) == 0) {
                srand(atoi(optarg));
            }
        }
    }

    if (argc - optind < 4) {
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    mapping_file_path = argv[optind++];
    prms_file_path = argv[optind++];
    pdb_file_path = argv[optind++];
    output_file_path = argv[optind++];

    form_factor_table *ff_table = default_ff_table(mapping_file_path);
    struct prm *prms = read_prm(prms_file_path, "0.0.6");

    mol_atom_group *ag = read_pdb(pdb_file_path, prms);
    if (no_hydrogens) {
        ag = remove_hydrogens(ag, prms);
    }


    double *vacuum_ff = malloc(ag->natoms * sizeof(double));
    double *dummy_ff = malloc(ag->natoms * sizeof(double));

    get_vacuum_ff(vacuum_ff, ag, ff_table, prms);
    get_dummy_ff(dummy_ff, ag, ff_table, prms);

    double max_d = max_dist(ag);

    saxs_profile *profile = malloc(sizeof(saxs_profile));
    saxs_profile_create(profile, q_start, q_end, q_step);

    compute_profile(profile, ag, vacuum_ff, dummy_ff, max_d, ff_table);

    // Apply random noise
    if (add_noise) {
        // Take average of last 10 samples
        double avg_tail = 0.0;
        int n = 10;
        for (int i = 0; i < n; ++i) {
            avg_tail += profile->In[profile->nsamples - i - 1];
        }
        avg_tail = avg_tail / n;

        for (int i = 0; i < profile->nsamples; ++i) {
            profile->In[i] += gaussrand(0, avg_tail*variance_percent/100.0);
        }
    }

    FILE *ofp = fopen(output_file_path, "w");
    if (ofp == NULL) {
        perror("fopen");
        exit(EXIT_FAILURE);
    }
    saxs_profile_fprint(ofp, profile);
    fclose(ofp);
}

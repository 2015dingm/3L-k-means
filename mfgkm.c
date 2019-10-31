// Copyright (c) 2011 Shenzhen Institutes of Advanced Technology
// Chinese Academy of Sciences

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <R.h>

#include <math.h>
#include <ctype.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

// #include "Utils.h"

// Copyright (c) 2012 Shenzhen Institutes of Advanced Technology
// Chinese Academy of Sciences

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <R.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void sum_squares(const double *x, // Data matrix
		const int *nr, // number of rows
		const int *nc, // number of columns
		const int *k, // number of clusters
		int *cluster, // A vector of integers (from 1:k) indicating the cluster to which each point is allocated.
		double *centers, // Cluster centers of every cluster
		double *totss, // The total sum of squares.
		double *withiness // Vector of within-cluster sum of squares, one component per cluster.
		) {
	//printf("Running sum of squares\n");
	int i, j;

	// vector to save average of all objects ---------------------------
	double* global_average = (double *) malloc(*nc * sizeof(double));
	for (j = 0; j < *nc; j++) {
		global_average[j] = 0.0;
	}

	for (i = 0; i < *nr; i++) {
		for (j = 0; j < *nc; j++) {
			// global_average[j] += x[i][j]
			global_average[j] += x[j * (*nr) + i];
		}
	}
	for (j = 0; j < *nc; j++) {
		global_average[j] = global_average[j] / (*nr);
	}

	// calculate total sum of square ---------------------------------

	*totss = 0.0;
	double tmp_ss, temp;
	for (i = 0; i < (*nr); i++) {
		tmp_ss = 0.0;
		for (j = 0; j < (*nc); j++) {
			temp = global_average[j] - x[j * (*nr) + i];
			tmp_ss += temp * temp;
		}
		*totss += tmp_ss;
	}

	// calculate withiness --------------------------------------------
	int t;
	for (t = 0; t < (*k); t++) {
		withiness[t] = 0.0;
	}

	for (i = 0; i < (*nr); i++) {

		// calculate the sum of square of object i and the center of the cluster which object i in.
		tmp_ss = 0;
		for (j = 0; j < (*nc); j++) {
			temp = x[j * (*nr) + i] - centers[j * (*k) + cluster[i]];
			tmp_ss += temp * temp;
		}

		// sum the with in class sum of square of the cluster which object i in.
		withiness[cluster[i]] += tmp_ss;
	}

	free(global_average);
	return;
}

void parseGroup(const char **strGroup, // string of group information, formatted as "0,1,2,4;3,5;6,7,8" or "0-2,4;3,4,8;6-8", where ";" defines a group;
		int *numGroups, // no. of groups
		int *groupInfo // group assignment
		) {

	int length, index, end, i = 0, j;
	char *buffer, *p1;
	buffer = (char *)malloc(20*sizeof(char));

	// -1: error, ignore next data until next "," or ";"
	//  0: no data
	//  1: single index
	//  2: range
	short status = -1;

	length = strlen(strGroup[0]);
	*numGroups = 0;

	p1 = &buffer[0];
	status = 0;
	while (i < length) {
		switch (strGroup[0][i]) {
		case '-':
			switch (status) {
			case 1:
			case 2:
				*p1 = '\0';
				index = atoi(buffer);
				p1 = &buffer[0];
				status = 2;
				break;
			case 0:
				status = -1;
				break;
			default:
				break;
			}
			break;
		case ';':
		case ',':
			//generate index
			switch (status) {
			case 1:
				*p1 = '\0';
				index = atoi(buffer);
				groupInfo[index] = *numGroups;
				p1 = &buffer[0];
				status = 0;
				break;
			case 2:
				//range
				*p1 = '\0';
				end = atoi(buffer);
				for (j = index; j <= end; j++) {
					groupInfo[j] = *numGroups;
				}
				p1 = &buffer[0];
				status = 0;
				break;
			default:
				//ignore
				status = -1;
				break;
			}
			if (status == -1) {
				//reset
				status = 0;
			} else if (strGroup[0][i] == ';') {
				//next group
				(*numGroups)++;
			}
			break;
		default:
			switch (status) {
			case 0:
				status = 1;
			case 1:
			case 2:
				*p1 = strGroup[0][i];
				p1++;
				break;
			default:
				//ignore
				break;
			}
			break;
		}

		i++;
	}
	//process left
	//generate index
	switch (status) {
	case 1:
		*p1 = '\0';
		index = atoi(buffer);
		groupInfo[index] = *numGroups;
		p1 = &buffer[0];
		status = 0;
		break;
	case 2:
		//range
		*p1 = '\0';
		end = atoi(buffer);
		for (j = index; j <= end; j++) {
			groupInfo[j] = *numGroups;
		}
		p1 = &buffer[0];
		status = 0;
		break;
	default:
		//ignore
		break;
	}
	// add 1 to the no. of groups
	(*numGroups)++;
	free(buffer);
}

double eu_distance(const double f1, const double f2) {
	return (f1 - f2) * (f1 - f2);
}

void init_centers(const double *x, const int *nr, const int *nc, const int *k,
		double *centers) {
	int i, j, l, flag, index, *random_obj_num;
	random_obj_num = (int *) calloc(*k, sizeof(int));

	if (!random_obj_num) {
		error("can't allocate random_obj_num\n");
	}

	for (l = 0; l < *k; l++) {
		random_obj_num[l] = -1;
	}

	for (l = 0; l < *k; l++) {
		flag = 1;

		while (flag) {
		  // index = (int) (rand() % (*nr));
			index = (int) (*nr-1) * unif_rand();
			flag = 0;
			for (i = 0; i < l; i++) {
				if (random_obj_num[i] == index)
					flag = 1;
			}
		}
		random_obj_num[l] = index;

		for (j = 0; j < *nc; j++)
			centers[j * (*k) + l] = x[j * (*nr) + index];

	} // for l = 1:k
	free(random_obj_num);
}

// ----- Calculate exp^{a}/sum(exp^{a}) -----

void expNormalize(double *a, int length, double minValue) {
	int i;
	double max;
	double sum = 0;

	max = a[0];

	for (i = 0; i < length; i++)
		if (a[i] > max)
			max = a[i];

	for (i = 0; i < length; i++) {
		a[i] = exp(a[i] - max);
		sum += a[i];
	}

	double sum2 = 0;
	for (i = 0; i < length; i++) {
		a[i] /= sum;
		if (a[i] < minValue) {
			a[i] = minValue;
		}
		sum2 += a[i];
	}
	for (i = 0; i < length; i++) {
		a[i] /= sum2;
	}
}












void init_featureWeight(double *featureWeight, const int *k, const int *nc,
		const int *numGroups, const int *groupInfo) {
	int l, j, *nums;
	nums = (int *) calloc(*numGroups, sizeof(int));
	for (j = 0; j < *nc; j++) {
		nums[groupInfo[j]]++;
	}
	for (l = 0; l < *k; l++)
		for (j = 0; j < *nc; j++)
			featureWeight[j * (*k) + l] = 1.0 / nums[groupInfo[j]];
	free(nums);
}

void init_groupWeight(double *groupWeight, const int *k, const int *numGroups) {
	int l, t;
	for (l = 0; l < *k; l++)
		for (t = 0; t < *numGroups; t++)
			groupWeight[t * (*k) + l] = 1.0 / *numGroups;
}

void update_cluster(const double *x, const double *e, const int *nr, const int *nc, const int *k,
		const int *numGroups, const int *groupInfo, int *cluster,
		double *centers, double *featureWeight, double *groupWeight) {
	int i, j, l;
	double min_dist, o_dist;

	for (i = 0; i < *nr; i++) {
		min_dist = 1.79769e+308;
		for (l = 0; l < *k; l++) {

			o_dist = 0.0;
			for (j = 0; j < *nc; j++) {
				o_dist += featureWeight[j * (*k) + l]
						* groupWeight[groupInfo[j] * (*k) + l]
			            * e[i]
						* eu_distance(centers[j * (*k) + l], x[j * (*nr) + i]);
			}
			if (o_dist <= min_dist) {
				min_dist = o_dist;
				cluster[i] = l;
			}
		}
	}
}

int update_centers(const double *x, const double *e, const int *nr, const int *nc, const int *k,
		int *cluster, double *centers) {
    int i, j, l, *no_cluster;
    double *no_cluster_e;
	no_cluster = (int *) calloc(*k, sizeof(int));
    no_cluster_e = (double *) calloc(*k, sizeof(double));

	if (!no_cluster) {
		error("can not allocate [].\n");
	}

	for (l = 0; l < *k; l++) {
		for (j = 0; j < *nc; j++) {
			centers[j * (*k) + l] = 0.0;
		}
	}

	for (i = 0; i < *nr; i++) {
		no_cluster_e[cluster[i]] += e[i];
        no_cluster[cluster[i]] += 1;
		for (j = 0; j < *nc; j++) {
			centers[j * (*k) + cluster[i]] += e[i] * x[j * (*nr) + i];
		}
	}

	int flag = 1;
	for (l = 0; l < *k; l++) {
		if (no_cluster[l] == 0) {
			flag = 0;
			break;
		}
		for (j = 0; j < *nc; j++) {
			centers[j * (*k) + l] /= (double) no_cluster_e[l];
		}
	}
	free(no_cluster);
    free(no_cluster_e);
	return flag;
}

void update_featureWeight(const double *x, const double *e, const int *nr, const int *nc,
		const int *k, const double *eta, const int *numGroups,
		const int *groupInfo, int *cluster, double *centers,
		double *featureWeight, double *groupWeight) {
	int i, j, l, t;

	int index;
	for (i = 0; i < (*nc) * (*k); i++) {
		featureWeight[i] = 0;
	}

	for (j = 0; j < *nc; j++) {
		for (i = 0; i < *nr; i++) {
			index = j * (*k) + cluster[i];
			featureWeight[index] += groupWeight[groupInfo[j]*(*k)+cluster[i]]
					* e[i] * eu_distance(x[j * (*nr) + i], centers[j * (*k) + cluster[i]]);
		}
	}

	double *sum, *sum2, *max;
	sum = (double*) malloc(*numGroups * sizeof(double));
	sum2 = (double*) malloc(*numGroups * sizeof(double));
	max = (double*) malloc(*numGroups * sizeof(double));

	for (t = 0; t < *numGroups; t++) {
		sum[t] = 0;
		sum2[t] = 0;
		max[t] = -1.79769e+308;
	}

	double minWeight = (double) 0.00001 / (*nc);
	for (l = 0; l < *k; l++) { // every CLUSTER
		for (t = 0; t < *numGroups; t++) {
			sum[t] = 0;
			sum2[t] = 0;
		}
		//find maximum
		for (j = 0; j < *nc; j++) {
			index = j * (*k) + l;
			featureWeight[index] = -featureWeight[index] / (*eta);
			if (max[groupInfo[j]] < featureWeight[index]) {
				max[groupInfo[j]] = featureWeight[index];
			}
		}
		//compute exp()
		for (j = 0; j < *nc; j++) {
			index = j * (*k) + l;
			featureWeight[index] = exp(
					featureWeight[index] - max[groupInfo[j]]);
			sum[groupInfo[j]] += featureWeight[index];
		}
		//normalize
		for (j = 0; j < *nc; j++) {
			index = j * (*k) + l;
			featureWeight[index] /= sum[groupInfo[j]];
			if (featureWeight[index] < minWeight) {
				featureWeight[index] = minWeight;
			}
			sum2[groupInfo[j]] += featureWeight[index];
		}

		for (j = 0; j < *nc; j++) {
			featureWeight[j * (*k) + l] /= sum2[groupInfo[j]];
		}
	}

	free(sum);
	free(sum2);
	free(max);
}

void update_groupWeight(const double *x, const double *e, const int *nr, const int *nc,
		const int *k, const double *lambda, const int *numGroups,
		const int *groupInfo, int *cluster, double *centers,
		double *featureWeight, double *groupWeight) {
	int i, j, l, t, index;

	for (l = 0; l < (*numGroups) * (*k); l++) {
		groupWeight[l] = 0;
	}

	for (j = 0; j < *nc; j++) {
		index = groupInfo[j] * (*k);
		for (i = 0; i < *nr; i++) {
			groupWeight[index + cluster[i]] += featureWeight[j * (*k)
					+ cluster[i]] * e[i]
					* eu_distance(centers[j * (*k) + cluster[i]],
							x[j * (*nr) + i]);
		}
	}

	for (l = 0; l < (*numGroups) * (*k); l++) {
		groupWeight[l] = -groupWeight[l] / (*lambda);
	}

	double sum = 0.0, sum2, max = 0.0;

// implement expNormalize()
	double minWeight = (double) 0.00001 / (*numGroups);
	for (l = 0; l < *k; l++) {
		sum = 0.0;
		max = groupWeight[0 * (*k) + l]; // initially assign gw[l][0] to max
		for (t = 1; t < *numGroups; t++) {
			index = t * (*k) + l;
			if (groupWeight[index] >= max)
				max = groupWeight[index];
		}

		for (t = 0; t < *numGroups; t++) {
			index = t * (*k) + l;
			groupWeight[index] = exp(groupWeight[index] - max);
			sum += groupWeight[index];
		}

		sum2 = 0;
		for (t = 0; t < *numGroups; t++) {
			index = t * (*k) + l;
			groupWeight[index] /= sum;
			if (groupWeight[index] < minWeight) {
				groupWeight[index] = minWeight;
			}
			sum2 += groupWeight[index];
		}

		if (sum2 != 1) {
			for (t = 0; t < *numGroups; t++) {
				groupWeight[t * (*k) + l] /= sum2;
			}
		}
	}
}

double calculate_cost(const double *x, const double *e, const int *nr, const int *nc,
		const int *k, const double *lambda, const double *eta,
		const int *numGroups, const int *groupInfo, int *cluster,
		double *centers, double *featureWeight, double *groupWeight) {

	int i, j, l, t;
	double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0, dispersion;

	for (i = 0; i < *nr; i++) {
		for (j = 0; j < *nc; j++) {
			sum0 += groupWeight[groupInfo[j] * (*k) + cluster[i]]
					* featureWeight[j * (*k) + cluster[i]] * e[i]
					* eu_distance(centers[j * (*k) + cluster[i]],
							x[j * (*nr) + i]);
		}
	}

	for (l = 0; l < *k; l++) {
		for (t = 0; t < *numGroups; t++) {
			sum1 += groupWeight[t * (*k) + l] * log(groupWeight[t * (*k) + l]);
		}

		for (j = 0; j < *nc; j++)
			sum2 += featureWeight[j * (*k) + l]
					* log(featureWeight[j * (*k) + l]);

	}

	dispersion = sum0 + sum1 * (*lambda) + sum2 * (*eta);
	return dispersion;
}

void mfgkm(const double *x, const double *e, const int *nr, const int *nc, const int *k,
		const double *lambda, const double *eta, const int *numGroups,
		const int *groupInfo, const double *delta, const int *maxiter,
		const int *maxrestart, int *init, // const unsigned int *seed,
		int *cluster,
		double *centers, double *featureWeight, double *groupWeight,
		int *iterations, int *restarts, int *totiter, double *totalCost, //
		double *totss, //    total sum of squares
		double *withiness // vector of sum of square in every cluster
		) {

	double dispersion1, dispersion2;
	int flag_not_restart = 1;

	*restarts = 0;
	*iterations = 0;
	*totiter = 0;

	// srand(seed);
	GetRNGstate();
	while ((*restarts) < *maxrestart) {
		if (*init == 0)
		  init_centers(x, nr, nc, k, centers); // assign randomly

		init_featureWeight(featureWeight, k, nc, numGroups, groupInfo); // equal value
		init_groupWeight(groupWeight, k, numGroups); // equal value

		*iterations = 0;
		dispersion2 = 1.79769e+308;

		while ((*iterations) < *maxiter) {
			(*iterations)++;
			(*totiter)++;
			dispersion1 = dispersion2;

			update_cluster(x, e, nr, nc, k, numGroups, groupInfo, cluster, centers,
					featureWeight, groupWeight);

			flag_not_restart = update_centers(x, e, nr, nc, k, cluster, centers);

			if (!flag_not_restart) {
				(*restarts)++;
				break;
			}

			update_featureWeight(x, e, nr, nc, k, eta, numGroups, groupInfo,
					cluster, centers, featureWeight, groupWeight);

			update_groupWeight(x, e, nr, nc, k, lambda, numGroups, groupInfo,
					cluster, centers, featureWeight, groupWeight);

			dispersion2 = calculate_cost(x, e, nr, nc, k, lambda, eta, numGroups,
					groupInfo, cluster, centers, featureWeight, groupWeight);

			// if change of dispersion below delta or iterations exceed max iterations,
			//	then, terminate and return.
			//  but if dispersion < 0, then ???
			if ((fabs((dispersion1 - dispersion2) / dispersion1)) <= (*delta)
					|| (*iterations) == *maxiter) {

				*totalCost = dispersion1;

				/**
				 * calculate sum of squares
				 */
				sum_squares(x, nr, nc, k, cluster, centers, totss, withiness);
				// Done.
				PutRNGstate();
				return;
			}
		}
	}

	PutRNGstate();
}

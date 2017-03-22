#include <cmath>
#include <iostream>
#include <windows.h>
#include <GL/glut.h>
#include <ctime>

#define M_PI 3.14159265358979323846

#define PTS 1000000
#define K 6
#define RAIO_GERACAO 100

using namespace std;

typedef struct { float x, y; int group, q; } point_t, *point;

point v;

/* Exibe os pontos na tela */
void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glPushMatrix();
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho((-1) * RAIO_GERACAO, RAIO_GERACAO, (-1) * RAIO_GERACAO, RAIO_GERACAO, -1, 1);

    glPointSize(5.0);
    glBegin(GL_POINTS);
    for (int i = 0; i < PTS; i++) {
        glColor3ub(v[i].group * (254 / K), 100, 0);
        glVertex2f(v[i].x, v[i].y);
    }
    glEnd();

    glPopMatrix();

    glutSwapBuffers();
}

/* Gera números aleatórios */
double randf(double m) {
	return m * rand() / (RAND_MAX - 1.);
}

/* Gera um conjunto de pontos com coordenadas aleatórias */
point gen_xy(int count, double radius) {
	double ang, r;
	point p, pt = (point)(malloc(sizeof(point_t) * count));

	for (p = pt + count; p-- > pt;) {
		ang = randf(2 * M_PI);
		r = randf(radius);
		p->x = r * cos(ang);
		p->y = r * sin(ang);
	}

	return pt;
}

/* Calcula a distância entre dois pontos */
inline double dist2(point a, point b) {
	double x = a->x - b->x, y = a->y - b->y;
	return x*x + y*y;
}

/* Encontra o cluster mais próximo de um determinado ponto */
inline int nearest(point pt, point cent, int n_cluster, double *d2) {
	int i, min_i;
	point c;
	double d, min_d;

# define for_n for (c = cent, i = 0; i < n_cluster; i++, c++)
	for_n {
		min_d = HUGE_VAL;
		min_i = pt->group;
		for_n {
			if (min_d > (d = dist2(c, pt))) {
				min_d = d; min_i = i;
			}
		}
	}
	if (d2) *d2 = min_d;
	return min_i;
}

/* K-Means++ */
void kpp(point pts, int len, point cent, int n_cent) {
# define for_len for (j = 0, p = pts; j < len; j++, p++)
	int i, j;
	int n_cluster;
	double sum, *d = (double*)(malloc(sizeof(double) * len));

	point p, c;
	cent[0] = pts[ rand() % len ];
	for (n_cluster = 1; n_cluster < n_cent; n_cluster++) {
		sum = 0;
		for_len {
			nearest(p, cent, n_cluster, d + j);
			sum += d[j];
		}
		sum = randf(sum);
		for_len {
			if ((sum -= d[j]) > 0) continue;
			cent[n_cluster] = pts[j];
			break;
		}
	}
	for_len p->group = nearest(p, cent, n_cluster, 0);
	free(d);
}

/* K-Means */
point lloyd(point pts, int len, int n_cluster) {
	int i, j, min_i;
	int changed;

	point cent = (point)(malloc(sizeof(point_t) * n_cluster)), p, c;

	kpp(pts, len, cent, n_cluster);

	do {
		for_n { c->group = 0; c->x = c->y = 0; }
		for_len {
			c = cent + p->group;
			c->group++;
			c->x += p->x;
			c->y += p->y;
		}
		for_n { c->x /= c->group; c->y /= c->group; }

		changed = 0;
		for_len {
			min_i = nearest(p, cent, n_cluster, 0);
			if (min_i != p->group) {
				changed++;
				p->group = min_i;
			}
		}
	} while (changed > (len >> 10));

	for_n { c->group = i; }

	return cent;
}

int main(int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize(640, 640);
    glutCreateWindow("K-Means++");
    glutDisplayFunc(display);

    int i;
	clock_t start;

	v = gen_xy(PTS, RAIO_GERACAO);

	start = clock();
	point c = lloyd(v, PTS, K);
	cout << clock() - start << endl;

	free(c);

    glutMainLoop();

    return 0;
}
#include <cmath>
#include <vector>
#include <iostream>
#include <windows.h>
#include <GL/glut.h>
#include <cfloat>
#include <ctime>

#define M_PI 3.14159265358979323846

#define PTS 1000000
#define K 6
#define DISCRETIZACAO 0.5
#define RAIO_GERACAO 100

using namespace std;

typedef struct { float x, y; int group, q; } point_t, *point;

point v;

vector<point_t> listaPontosRepresentantes;
int qDePontosRepresentantes = 0;
point pontoRepresentante;

double xMin = DBL_MAX, yMin = DBL_MAX;
double xMax = 0, yMax = 0;

int qLinhasGrid, qColunasGrid;

int xDiscreto, yDiscreto;

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
double randf (double m) {
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

/* K-Means++ (Adaptado) */
void kpp(point pts, int len, point pts2, int len2, point cent, int n_cent) {
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
	for (j = 0, p = pts2; j < len2; j++, p++) p->group = nearest(p, cent, n_cluster, 0);
	free(d);
}

/* K-Means (Adaptado) */
point lloyd(point pts, int len, point pts2, int len2, int n_cluster) {
	int i, j, min_i;
	int changed;

	point cent = (point)(malloc(sizeof(point_t) * n_cluster)), p, c;

	kpp(pts2, len2, pts, len, cent, n_cluster);

	do {
		for_n { c->group = 0; c->x = c->y = 0; }
		for_len {
			c = cent + p->group;
			c->group += p->q;
			c->x += p->x * p->q;
			c->y += p->y * p->q;
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

/* K-Means 2.0 */
point kmeans2(point pontos, int qDePontos, int k, double discretizacao) {
	int i, j;

	for (i = 0; i < qDePontos; i++) {
		if (pontos[i].x < xMin) xMin = pontos[i].x;
		if (pontos[i].x > xMax) xMax = pontos[i].x;
		if (pontos[i].y < yMin) yMin = pontos[i].y;
		if (pontos[i].y > yMax) yMax = pontos[i].y;
	}

	qLinhasGrid = ceil((xMax - xMin + 1) / discretizacao);
	qColunasGrid = ceil((yMax - yMin + 1) / discretizacao);

	int grid[qLinhasGrid][qColunasGrid];
	for (i = 0; i < qLinhasGrid; i++) {
		for (j = 0; j < qColunasGrid; j++)
			grid[i][j] = -1;
	}

	for (i = 0; i < qDePontos; i++) {
		xDiscreto = floor((pontos[i].x - xMin) / discretizacao);
		yDiscreto = floor((pontos[i].y - yMin) / discretizacao);
		if (grid[xDiscreto][yDiscreto] == -1) {
			grid[xDiscreto][yDiscreto] = qDePontosRepresentantes;
			qDePontosRepresentantes++;
			listaPontosRepresentantes.emplace_back();
			pontoRepresentante = &(listaPontosRepresentantes[qDePontosRepresentantes-1]);
			pontoRepresentante->x = 0;
			pontoRepresentante->y = 0;
			pontoRepresentante->q = 0;
		}
		else {
			pontoRepresentante = &(listaPontosRepresentantes[grid[xDiscreto][yDiscreto]]);
		}
		pontoRepresentante->x += pontos[i].x;
		pontoRepresentante->y += pontos[i].y;
		pontoRepresentante->q++;
	}

	for (i = 0; i < qDePontosRepresentantes; i++) {
		listaPontosRepresentantes[i].x /= listaPontosRepresentantes[i].q;
		listaPontosRepresentantes[i].y /= listaPontosRepresentantes[i].q;
	}

	point pontosRepresentantes = &(listaPontosRepresentantes[0]);

	point c = lloyd(pontosRepresentantes, qDePontosRepresentantes, pontos, qDePontos, k);

	for (i = 0; i < qDePontos; i++) {
		xDiscreto = floor((pontos[i].x - xMin) / discretizacao);
		yDiscreto = floor((pontos[i].y - yMin) / discretizacao);
		pontos[i].group = pontosRepresentantes[grid[xDiscreto][yDiscreto]].group;
	}

	return c;
}

int main(int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize(640, 640);
    glutCreateWindow("K-Means 2.0");
    glutDisplayFunc(display);

    int i;
	clock_t start;

	v = gen_xy(PTS, RAIO_GERACAO);

	start = clock();
	point c = kmeans2(v, PTS, K, DISCRETIZACAO);
	cout << clock() - start << endl;

	free(c);

    glutMainLoop();

    return 0;
}
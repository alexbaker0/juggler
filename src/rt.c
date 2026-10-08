/*      RT1.C    Ray tracing program
        Copyright 1987 Eric Graham
        Permission is granted to copy and modify this file, provided
        that this notice is retained.
*/

/*  rt.c -- Eric Graham's 1987 Amiga ray-tracing engine (RT1.C), compiled to
 *  WebAssembly so the original Juggler scene renders in the browser.
 *
 *  The ENGINE section is copied unchanged from juggler_orig.c, the PC port:
 *  K&R definitions converted to ANSI prototypes, reflect() fixed (the scanned
 *  original divided by the normal: y=x-2(x.n)n is correct), and a mirror
 *  recursion guard (MAXDEPTH).  Scene parsing and the camera trigonometry
 *  live in JavaScript (juggler.js); this file only traces rays.
 *
 *  Build (clang from WinLibs LLVM):
 *    clang --target=wasm32 -O3 -mbulk-memory -nostdlib -Wl,--no-entry -o rt.wasm src/rt.c
 */

#define sqrt __builtin_sqrt
#define EXPORT(n) __attribute__((export_name(n)))

#define BIG   1.0e10
#define SMALL 1.0e-3
#define DULL    0
#define BRIGHT  1
#define MIRROR  2
#define MAXDEPTH 8

struct lamp   { double pos[3];    double color[3]; double radius; };
struct sphere { double pos[3];    double color[3]; double radius; int type; };
struct patch  { double pos[3];    double normal[3]; double color[3]; };

struct world {
    int           numsp;   struct sphere *sp;
    int           numlmp;  struct lamp   *lmp;
    struct patch  horizon[2];
    double        illum[3], skyhor[3], skyzen[3];
};

struct observer {
    double obspos[3], viewdir[3], uhat[3], vhat[3];
    double fl, px, py;
    int    nx, ny;
};

/* ---- forward declarations ------------------------------------------ */
double dot(double *a, double *b);
void   vecsub(double *a, double *b, double *c);
void   veccopy(double *a, double *b);
void   colorcpy(double *a, double *b);
void   vecprod(double *a, double *b, double *c);
int    veczero(double *v);
void   genline(double *l, double *a, double *b);
void   point(double *pos, double t, double *line);
int    intsplin(double *t, double *line, struct sphere *sp);
int    inthor(double *t, double *line);
void   setnorm(struct patch *p, struct sphere *s);
int    gingham(double *pos);
void   reflect(double *y, double *n, double *x);
void   raytrace(double *brite, double *line, struct world *w);
void   skybrite(double *brite, double *line, struct world *w);
void   pixbrite(double *brite, struct patch *p, struct world *w, struct sphere *spc);
int    glint(double *brite, struct patch *p, struct world *w, struct sphere *spc, double *incident);
int    mirror(double *brite, struct patch *p, struct world *w, double *incident);
void   pixlinef(double *line, struct observer *o, double di, double dj);

static int g_depth = 0;

/* ==================================================================== */
/*  ENGINE  (Eric Graham, 1987 -- prototypes added, reflect() fixed)    */
/* ==================================================================== */

void raytrace(double *brite, double *line, struct world *w)
{
    double t, tmin, pos[3];  int k;
    struct patch  ptch;  struct sphere *spnear;  struct lamp *lmpnear;

    tmin = BIG;  spnear = 0;
    for (k = 0; k < w->numsp; ++k)
        if (intsplin(&t, line, w->sp + k))
            if (t < tmin) { tmin = t; spnear = w->sp + k; }

    lmpnear = 0;
    for (k = 0; k < w->numlmp; ++k)
        if (intsplin(&t, line, (struct sphere *)(w->lmp + k)))
            if (t < tmin) { tmin = t; lmpnear = w->lmp + k; }

    if (lmpnear) {
        for (k = 0; k < 3; ++k)
            brite[k] = lmpnear->color[k] / (lmpnear->radius * lmpnear->radius);
        return;
    }
    if (inthor(&t, line))
        if (t < tmin) {
            point(pos, t, line);  k = gingham(pos);
            veccopy(w->horizon[k].pos, pos);
            pixbrite(brite, &(w->horizon[k]), w, 0);
            return;
        }
    if (spnear) {
        point(ptch.pos, tmin, line);  setnorm(&ptch, spnear);
        colorcpy(ptch.color, spnear->color);
        switch (spnear->type) {
            case BRIGHT: if (glint(brite, &ptch, w, spnear, line)) return; /* fall */
            case DULL:   pixbrite(brite, &ptch, w, spnear); return;
            case MIRROR: mirror(brite, &ptch, w, line); return;
        }
        return;
    }
    skybrite(brite, line, w);
}

void skybrite(double *brite, double *line, struct world *w)
{
    double sin2, cos2;  int k;
    sin2  = line[5] * line[5];
    sin2 /= (line[1]*line[1] + line[3]*line[3] + sin2);
    cos2  = 1.0 - sin2;
    for (k = 0; k < 3; ++k) brite[k] = cos2 * w->skyhor[k] + sin2 * w->skyzen[k];
}

void pixlinef(double *line, struct observer *o, double di, double dj)
{
    double x, y, tp[3];  int k;
    y = (0.5 * o->ny - dj) * o->py;
    x = (di - 0.5 * o->nx) * o->px;
    for (k = 0; k < 3; ++k)
        tp[k] = o->viewdir[k]*o->fl + y*o->vhat[k] + x*o->uhat[k] + o->obspos[k];
    genline(line, o->obspos, tp);
}

void vecsub(double *a, double *b, double *c) { int k; for (k=0;k<3;++k) a[k]=b[k]-c[k]; }

int intsplin(double *t, double *line, struct sphere *sp)
{
    double a, b, c, d, p, q, tt;  int k;
    a = b = 0.0;  c = sp->radius;  c = -c * c;
    for (k = 0; k < 3; ++k) {
        p = (*line++) - sp->pos[k];  q = *line++;
        a = q*q + a;  tt = q*p;  b = tt + tt + b;  c = p*p + c;
    }
    d = b*b - 4.0*a*c;
    if (d <= 0) return 0;
    d = sqrt(d);  *t = -(b + d)/(a + a);
    if (*t < SMALL) *t = (d - b)/(a + a);
    return *t > SMALL;
}

int inthor(double *t, double *line)
{ if (line[5] == 0.0) return 0; *t = -line[4]/line[5]; return *t > SMALL; }

void genline(double *l, double *a, double *b)
{ int k; for (k=0;k<3;++k){ *l++=a[k]; *l++=b[k]-a[k]; } }

double dot(double *a, double *b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }

void point(double *pos, double t, double *line)
{ int k; double a; for (k=0;k<3;++k){ a=*line++; pos[k]=a+(*line++)*t; } }

int glint(double *brite, struct patch *p, struct world *w,
          struct sphere *spc, double *incident)
{
    int k, l, firstlite;  static double minglint = 0.95;
    double line[6], t, lp[3], *pp, *ll, cosi, incvec[3], refvec[3], ref2;
    refvec[0]=refvec[1]=refvec[2]=0.0;  ref2=1.0;  firstlite=1;
    for (l = 0; l < w->numlmp; ++l) {
        ll=(w->lmp+l)->pos;  pp=p->pos;  vecsub(lp,ll,pp);  cosi=dot(lp,p->normal);
        if (cosi <= 0.0) continue;
        genline(line, pp, ll);
        for (k = 0; k < w->numsp; ++k) {
            if (w->sp + k == spc) continue;
            if (intsplin(&t, line, w->sp + k)) goto cont;
        }
        if (firstlite) {
            incvec[0]=incident[1]; incvec[1]=incident[3]; incvec[2]=incident[5];
            reflect(refvec, p->normal, incvec);  ref2=dot(refvec,refvec);  firstlite=0;
        }
        t = dot(lp, refvec);  t *= t/(dot(lp,lp)*ref2);
        if (t > minglint) { for (k=0;k<3;++k) brite[k]=1.0; return 1; }
    cont: ;
    }
    return 0;
}

int mirror(double *brite, struct patch *p, struct world *w, double *incident)
{
    int k;  double line[6], incvec[3], refvec[3], t;
    incvec[0]=incident[1]; incvec[1]=incident[3]; incvec[2]=incident[5];
    t = dot(p->normal, incvec);
    if (t >= 0) { for (k=0;k<3;++k) brite[k]=0.0; return 0; }
    if (g_depth >= MAXDEPTH) { for (k=0;k<3;++k) brite[k]=0.0; return 1; }
    reflect(refvec, p->normal, incvec);
    line[0]=p->pos[0]; line[2]=p->pos[1]; line[4]=p->pos[2];
    line[1]=refvec[0]; line[3]=refvec[1]; line[5]=refvec[2];
    g_depth++;  raytrace(brite, line, w);  g_depth--;
    for (k=0;k<3;++k) brite[k] = brite[k]*p->color[k];
    return 1;
}

void pixbrite(double *brite, struct patch *p, struct world *w, struct sphere *spc)
{
    int k, l;  double line[6], t, r, lp[3], *pp, *ll, cosi, diffuse;
    static double zenith[3] = {0.0,0.0,1.0}, f1 = 1.5, f2 = 0.4;
    diffuse = (dot(zenith, p->normal) + f1) * f2;
    for (k = 0; k < 3; ++k) brite[k] = diffuse * w->illum[k] * p->color[k];
    if (p && w) {
        for (l = 0; l < w->numlmp; ++l) {
            ll=(w->lmp+l)->pos; pp=p->pos; vecsub(lp,ll,pp);
            cosi=dot(lp,p->normal); if (cosi <= 0.0) goto cont;
            genline(line, pp, ll);
            for (k = 0; k < w->numsp; ++k) {
                if (w->sp + k == spc) continue;
                if (intsplin(&t, line, w->sp + k)) goto cont;
            }
            r = sqrt(dot(lp,lp));  cosi = cosi/(r*r*r);   /* genuine 1/r^2 * cos;
                the lamp was pre-scaled by lampfac (auto-exposure, see rt2.c)
                so a sun-distance lamp still lights the scene correctly */
            for (k = 0; k < 3; ++k) brite[k] += cosi*p->color[k]*w->lmp[l].color[k];
        cont: ;
        }
    }
}

void setnorm(struct patch *p, struct sphere *s)
{
    double *t, a;  int k;
    vecsub(t = p->normal, p->pos, s->pos);  a = 1.0/s->radius;
    for (k = 0; k < 3; ++k) { *t = (*t)*a; ++t; }
}

void colorcpy(double *a, double *b) { int k; for (k=0;k<3;++k) a[k]=b[k]; }
void veccopy (double *a, double *b) { int k; for (k=0;k<3;++k) a[k]=b[k]; }

int gingham(double *pos)
{
    double x, y;  int kx, ky;
    kx = ky = 0;  x = pos[0];  y = pos[1];
    if (x < 0)   { x = -x;  ++kx; }
    if (y < 0.0) { y = -y;  ++ky; }
    return ((((int)x)+kx)/3 + (((int)y)+ky)/3) % 2;
}

void reflect(double *y, double *n, double *x)
{
    double u[3], v[3], vv, xn, xv;  int k;
    vecprod(u, x, n);
    if (veczero(u)) { y[0]=-x[0]; y[1]=-x[1]; y[2]=-x[2]; return; }
    vecprod(v, u, n);
    vv = dot(v,v);  xv = dot(x,v)/vv;  xn = dot(x,n);
    for (k = 0; k < 3; ++k) y[k] = xv*v[k] - xn*n[k];
}

void vecprod(double *a, double *b, double *c)
{
    a[0]=b[1]*c[2]-b[2]*c[1];  a[1]=b[2]*c[0]-b[0]*c[2];  a[2]=b[0]*c[1]-b[1]*c[0];
}
int veczero(double *v) { if(v[0])return 0; if(v[1])return 0; if(v[2])return 0; return 1; }

/* ==================================================================== */
/*  WEB API  -- scene setup from JS, then render one row at a time      */
/* ==================================================================== */

static struct sphere   g_sp[2048];
static struct lamp     g_lmp[64];
static struct world    g_w;
static struct observer g_o;
static unsigned char   g_row[4096 * 4];      /* one RGBA row */

EXPORT("clear_scene") void clear_scene(void) { g_w.numsp = 0; g_w.numlmp = 0; g_w.sp = g_sp; g_w.lmp = g_lmp; }

EXPORT("add_sphere") void add_sphere(double x, double y, double z, double rad,
                                     double r, double g, double b, int type)
{
    struct sphere *s;
    if (g_w.numsp >= 2048) return;
    s = &g_sp[g_w.numsp++];
    s->pos[0] = x; s->pos[1] = y; s->pos[2] = z; s->radius = rad;
    s->color[0] = r; s->color[1] = g; s->color[2] = b; s->type = type;
}

EXPORT("add_lamp") void add_lamp(double x, double y, double z, double rad,
                                 double r, double g, double b)
{
    struct lamp *l;
    if (g_w.numlmp >= 64) return;
    l = &g_lmp[g_w.numlmp++];
    l->pos[0] = x; l->pos[1] = y; l->pos[2] = z; l->radius = rad;
    l->color[0] = r; l->color[1] = g; l->color[2] = b;
}

/* five world colours (same meaning as juggler_orig.c's load_scene) and the
   rt2.c "lampfac" auto-exposure; returns lampfac */
EXPORT("set_world") double set_world(double t0r, double t0g, double t0b, double t1r, double t1g, double t1b,
                                     double ir, double ig, double ib, double zr, double zg, double zb,
                                     double hr, double hg, double hb)
{
    double lampfac = BIG, tp[3], rr, tt;  int ii, jj, kk;
    g_w.illum[0] = ir;  g_w.illum[1] = ig;  g_w.illum[2] = ib;
    g_w.skyzen[0] = zr; g_w.skyzen[1] = zg; g_w.skyzen[2] = zb;
    g_w.skyhor[0] = hr; g_w.skyhor[1] = hg; g_w.skyhor[2] = hb;
    for (ii = 0; ii < 2; ++ii) {
        g_w.horizon[ii].normal[0] = 0.0; g_w.horizon[ii].normal[1] = 0.0; g_w.horizon[ii].normal[2] = 1.0;
    }
    g_w.horizon[0].color[0] = t0r; g_w.horizon[0].color[1] = t0g; g_w.horizon[0].color[2] = t0b;
    g_w.horizon[1].color[0] = t1r; g_w.horizon[1].color[1] = t1g; g_w.horizon[1].color[2] = t1b;

    for (ii = 0; ii < g_w.numsp; ++ii)
        for (jj = 0; jj < g_w.numlmp; ++jj) {
            vecsub(tp, g_sp[ii].pos, g_lmp[jj].pos);
            rr = sqrt(dot(tp,tp)) - g_sp[ii].radius;
            for (kk = 0; kk < 3; ++kk) {
                tt = g_sp[ii].color[kk]*g_lmp[jj].color[kk]/(rr*rr);
                if (tt == 0.0) continue;
                tt = (1.0 - g_sp[ii].color[kk]*g_w.illum[kk]) / tt;
                if (tt < lampfac) lampfac = tt;
            }
        }
    for (jj = 0; jj < g_w.numlmp; ++jj)
        for (kk = 0; kk < 3; ++kk) g_lmp[jj].color[kk] *= lampfac;
    return lampfac;
}

EXPORT("set_camera") void set_camera(double ex, double ey, double ez,
                                     double vx, double vy, double vz,
                                     double ux, double uy, double uz,
                                     double wx, double wy, double wz,
                                     double fl, int nx, int ny)
{
    g_o.obspos[0] = ex; g_o.obspos[1] = ey; g_o.obspos[2] = ez;
    g_o.viewdir[0] = vx; g_o.viewdir[1] = vy; g_o.viewdir[2] = vz;
    g_o.uhat[0] = ux; g_o.uhat[1] = uy; g_o.uhat[2] = uz;
    g_o.vhat[0] = wx; g_o.vhat[1] = wy; g_o.vhat[2] = wz;
    g_o.fl = fl; g_o.nx = nx; g_o.ny = ny;
    g_o.px = 1.0 / nx;  g_o.py = 0.75 / ny;     /* 4:3 view window, as rt2.c */
}

EXPORT("row_buffer") unsigned char *row_buffer(void) { return g_row; }

/* Render image row j with aa x aa supersampling into row_buffer() as RGBA
   (same sampling and rounding as juggler_orig.c's main loop and ham()). */
EXPORT("render_row") void render_row(int j, int aa)
{
    double line[6], brite[3], acc[3];
    int i, sx, sy, k;
    for (i = 0; i < g_o.nx; i++) {
        acc[0] = acc[1] = acc[2] = 0.0;
        for (sy = 0; sy < aa; sy++)
            for (sx = 0; sx < aa; sx++) {
                double di = i + (sx + 0.5) / aa - 0.5;
                double dj = j + (sy + 0.5) / aa - 0.5;
                pixlinef(line, &g_o, di, dj);
                g_depth = 0;
                raytrace(brite, line, &g_w);
                for (k = 0; k < 3; k++) acc[k] += brite[k];
            }
        for (k = 0; k < 3; k++) {
            double c = acc[k] / (aa * aa);
            if (c < 0.0) c = 0.0;
            if (c > 1.0) c = 1.0;
            g_row[i * 4 + k] = (unsigned char)(c * 255.0 + 0.5);
        }
        g_row[i * 4 + 3] = 255;
    }
}

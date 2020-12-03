#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(cpesub)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(daucmats)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(djonck)(void *, void *, void *, void *);
extern void F77_NAME(f2bdry)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(femdor)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(fepow)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ferej)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(fessiz)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(jtpdf)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ktau)(void *, void *, void *, void *);
extern void F77_NAME(lehman)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lrtest)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(rocauc)(void *, void *, void *, void *);
extern void F77_NAME(rocarea)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(roccurve)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(smrocauc)(void *, void *, void *, void *);
extern void F77_NAME(strperm1)(void *, void *, void *, void *, void *);
extern void F77_NAME(uclrst)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"cpesub",   (DL_FUNC) &F77_NAME(cpesub),   10},
    {"daucmats", (DL_FUNC) &F77_NAME(daucmats), 10},
    {"djonck",   (DL_FUNC) &F77_NAME(djonck),    4},
    {"f2bdry",   (DL_FUNC) &F77_NAME(f2bdry),   13},
    {"femdor",   (DL_FUNC) &F77_NAME(femdor),    9},
    {"fepow",    (DL_FUNC) &F77_NAME(fepow),     8},
    {"ferej",    (DL_FUNC) &F77_NAME(ferej),     6},
    {"fessiz",   (DL_FUNC) &F77_NAME(fessiz),   10},
    {"jtpdf",    (DL_FUNC) &F77_NAME(jtpdf),     6},
    {"ktau",     (DL_FUNC) &F77_NAME(ktau),      4},
    {"lehman",   (DL_FUNC) &F77_NAME(lehman),    9},
    {"lrtest",   (DL_FUNC) &F77_NAME(lrtest),   14},
    {"rocarea",  (DL_FUNC) &F77_NAME(rocarea),   7},
    {"rocauc",   (DL_FUNC) &F77_NAME(rocauc),    4},
    {"roccurve", (DL_FUNC) &F77_NAME(roccurve),  8},
    {"smrocauc", (DL_FUNC) &F77_NAME(smrocauc),  4},
    {"strperm1", (DL_FUNC) &F77_NAME(strperm1),  5},
    {"uclrst",   (DL_FUNC) &F77_NAME(uclrst),   22},
    {NULL, NULL, 0}
};

void R_init_clinfun(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

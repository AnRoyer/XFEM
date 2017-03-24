// Gmsh - Copyright (C) 1997-2016 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to the public mailing list <gmsh@onelab.info>.

#ifndef _GMSH_CONFIG_H_
#define _GMSH_CONFIG_H_

/* #undef HAVE_3M */
/* #undef HAVE_ACIS */
/* #undef HAVE_ANN */
/* #undef HAVE_BAMG */
/* #undef HAVE_BFGS */
#define HAVE_BLAS
/* #undef HAVE_BLOSSOM */
/* #undef HAVE_CAIRO */
/* #undef HAVE_CHACO */
/* #undef HAVE_COMPRESSED_IO */
#define HAVE_DLOPEN
/* #undef HAVE_DINTEGRATION */
/* #undef HAVE_FLTK */
/* #undef HAVE_FOURIER_MODEL */
/* #undef HAVE_GMM */
/* #undef HAVE_GMP */
/* #undef HAVE_KBIPACK */
#define HAVE_LAPACK
/* #undef HAVE_LIBCGNS */
/* #undef HAVE_LIBJPEG */
/* #undef HAVE_LIBPNG */
/* #undef HAVE_LIBZ */
/* #undef HAVE_LINUX_JOYSTICK */
/* #undef HAVE_MATHEX */
/* #undef HAVE_MED */
/* #undef HAVE_MESH */
/* #undef HAVE_METIS */
/* #undef HAVE_MMG3D */
/* #undef HAVE_MPEG_ENCODE */
/* #undef HAVE_MPI */
/* #undef HAVE_MUMPS */
/* #undef HAVE_NATIVE_FILE_CHOOSER */
/* #undef HAVE_NETGEN */
/* #undef HAVE_NUMPY */
/* #undef HAVE_NO_INTPTR_T */
/* #undef HAVE_NO_SOCKLEN_T */
/* #undef HAVE_NO_STDINT_H */
/* #undef HAVE_NO_VSNPRINTF */
/* #undef HAVE_OCC */
/* #undef HAVE_ONELAB */
/* #undef HAVE_ONELAB2 */
/* #undef HAVE_ONELAB_METAMODEL */
/* #undef HAVE_UDT */
/* #undef HAVE_OPENGL */
/* #undef HAVE_OPTHOM */
/* #undef HAVE_OSMESA */
/* #undef HAVE_PARSER */
/* #undef HAVE_PETSC */
/* #undef HAVE_PETSC4PY */
/* #undef HAVE_PLUGINS */
/* #undef HAVE_POST */
/* #undef HAVE_POPPLER */
/* #undef HAVE_QT */
/* #undef HAVE_REVOROPT */
/* #undef HAVE_SALOME */
/* #undef HAVE_SGEOM */
/* #undef HAVE_SLEPC */
/* #undef HAVE_SOLVER */
/* #undef HAVE_TAUCS */
/* #undef HAVE_TETGEN */
/* #undef HAVE_VORO3D */
/* #undef HAVE_ZIPPER */

#define GMSH_CONFIG_OPTIONS " Blas(VecLib) Dlopen Lapack(VecLib)"



#endif

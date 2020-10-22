# jtk
Jan's toolkit

Header only C++ tookit for proof of concept development of mainly 3d geometry related algorithms and functionality.

- alsc.h: implementation of "Adaptive Least Squares Correlation: A powerful image matching technique. Armin Gruen. S. Afr. J. of Photogrammetry, Remote Sensing and Cartography 14 (3), 1985: 175-187".
- clamp.h: clamping tools when converting types.
- concurrency.h: encapsulates standard multithreading functionality such as parallel_for. By default I mainly use TBB but you could also easily use ppl.
- file_utils.h: file related utilities where filenames are assumed to be in utf-8 encoding. Conversion of std::string to std::wstring and vice versa. Comma separated value file suport. Listing files and subdirectories from a given directory. 
- fitting.h: geometric fitting algorithms, such as n-point registration.
- geometry.h: 3d mesh related functionality. Reading and/or writing some default file formats such as stl, obj, off, ply. Querying mesh structure (e.g. one-ring of a vertex). Smoothing, subdividing, filling holes. Computing normals.
- image.h: image-related functionality. Each row of the image is aligned in memory so that SSE2 can be used. Reading/writing of image file formats is not included. Here I typically use https://github.com/nothings/stb/blob/master/stb_image.h or https://github.com/lvandeve/lodepng.
- mat.h: matrix library using expression templates, with support for sparse matrices also. Many common matrix algorithms are available: matrix and vector arithmetic, matrix transpose, svd (singular value decomposition), pseudo inverse, lsd (solves overdetermined systems in a least squares sense using svd), lu decomposition, solve a system using lu decomposition, invert a matrix using lu decomposition, cholesky factorization, solve a system using cholesky factorization, qr decomposition, solve a system using qr decomposition, frobenius norm, compute a jacobian matrix using forward differences, levenberg-marquardt optimization, eigenvalue computation of square matrices, (preconditioned) conjugate gradient, (preconditioned) bicgstab.
- pipe.h: cross platform tools for creating pipes, sending and receiving messages via the pipes, and running processes.
- qbvh.h: inspired by https://www.embree.org/. Implementation of a quad bounding volume hierarchy (qbvh) (https://www.uni-ulm.de/fileadmin/website_uni_ulm/iui.inst.100/institut/Papers/QBVH.pdf) for fast querying of nearest triangles along a ray. 
- render.h: software renderer for point clouds. Uses AVX2 or SSE2.
- sse2neon.h: A c++ header file that converts intel sse intrinsics to neon intrinsics, from https://github.com/DLTcollab/sse2neon.
- timer.h: Very simple timer.
- utf8.h: utf-8 encoding written by Nemanja Trifunovic.
- vec.h: Simple 2d, 3d, and 4d vector classes, used by some of the other header files.

### building
jtk makes use of three CMake variables:

- JTK_THREADING: if you use conccurency.h, with this parameter you can choose the underlying multithreading implementation: Intel's TBB library, Microsoft's concurrency library, an own implementation using std::thread, or no multithreading.
- JTK_MAT_PARALLEL: If on, some matrix algorithms from mat.h will run in parallel.
- JTK_TARGET: Normally x64, but if you choose arm, then all sse instructions will be replaced automatically by neon instructions (sse2neon.h)

If your project uses jtk, then you should include jtk.cmake in your own CMakeLists.txt file so that all preprocessor directives are set correctly.



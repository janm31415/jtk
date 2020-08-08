# jtk
Jan's toolkit

Header only C++ tookit for proof of concept development of mainly 3d geometry related algorithms and functionality.

- concurrency.h: encapsulates standard multithreading functionality such as parallel_for. By default I mainly use TBB but you could also easily use ppl.
- file_utils.h: file related utilities where filenames are assumed to be in utf-8 encoding. Conversion of std::string to std::wstring and vice versa. Comma separated value file suport. Listing files and subdirectories from a given directory. 
- geometry.h: 3d mesh related functionality. Reading and/or writing some default file formats such as stl, obj, off, ply. Querying mesh structure (e.g. one-ring of a vertex). Smoothing, subdividing, filling holes. Computing normals.
- image.h: image-related functionality. Each row of the image is aligned in memory so that SSE2 can be used. Reading/writing of image file formats is not included. Here I typically use https://github.com/nothings/stb/blob/master/stb_image.h or https://github.com/lvandeve/lodepng.
- qbvh.h: inspired by https://www.embree.org/. Implementation of a quad bounding volume hierarchy (qbvh) (https://www.uni-ulm.de/fileadmin/website_uni_ulm/iui.inst.100/institut/Papers/QBVH.pdf) for fast querying of nearest triangles along a ray. 
- render.h: software renderer for point clouds. Uses AVX2 or SSE2.
- timer.h: Very simple timer.
- utf8.h: utf-8 encoding written by Nemanja Trifunovic.
- vec.h: Simple 2d, 3d, and 4d vector classes, used by some of the other header files.

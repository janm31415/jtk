/*
   Do this:
      #define JTK_PLY_IMPLEMENTATION
   before you include this file in *one* C++ file to create the implementation.
   // i.e. it should look like this:
   #include ...
   #include ...
   #include ...
   #define JTK_PLY_IMPLEMENTATION
   #include "jtk/ply.h"
 */

#ifndef JTK_PLY_H
#define JTK_PLY_H

 /* ----------------------------------------------------------------------
  * Single header file for reading and writing PLY files.
  *
  * Based on:
  *
  * RPly library, read/write PLY files
  * Diego Nehab, IMPA
  * http://www.impa.br/~diego/software/rply
  *
  * Load and store of PLY files in memory by:
  * Patryk Kiepas, Institute for Not-so-Advanced Study
  * https://github.com/quepas/O-RPLY
  *
  * This library is distributed under the MIT License.
  * ---------------------------------------------------------------------- */

#include "vec.h"
#include <vector>

#ifndef JTKPLYDEF
#ifdef JTK_PLY_STATIC
#define JTKPLYDEF static
#else
#define JTKPLYDEF extern
#endif
#endif

namespace jtk
  {

  JTKPLYDEF bool read_ply(const char* filename, std::vector<vec3<float>>& vertices, std::vector<vec3<float>>& normals, std::vector<uint32_t>& clrs, std::vector<vec3<uint32_t>>& triangles, std::vector<vec3<vec2<float>>>& uv);

  JTKPLYDEF bool read_ply(const wchar_t* filename, std::vector<vec3<float>>& vertices, std::vector<vec3<float>>& normals, std::vector<uint32_t>& clrs, std::vector<vec3<uint32_t>>& triangles, std::vector<vec3<vec2<float>>>& uv);

  JTKPLYDEF bool read_ply_from_memory(const char* buffer, uint64_t buffer_size, std::vector<vec3<float>>& vertices, std::vector<vec3<float>>& normals, std::vector<uint32_t>& clrs, std::vector<vec3<uint32_t>>& triangles, std::vector<vec3<vec2<float>>>& uv);

  JTKPLYDEF bool write_ply(const char* filename, const std::vector<vec3<float>>& vertices, const std::vector<vec3<float>>& normals, const std::vector<uint32_t>& clrs, const std::vector<vec3<uint32_t>>& triangles, const std::vector<vec3<vec2<float>>>& uv);

  JTKPLYDEF bool write_ply(const wchar_t* filename, const std::vector<vec3<float>>& vertices, const std::vector<vec3<float>>& normals, const std::vector<uint32_t>& clrs, const std::vector<vec3<uint32_t>>& triangles, const std::vector<vec3<vec2<float>>>& uv);


  /* ----------------------------------------------------------------------
   * Types
   * ---------------------------------------------------------------------- */
   /* structures are opaque */
  typedef struct t_ply_* p_ply;
  typedef struct t_ply_element_* p_ply_element;
  typedef struct t_ply_property_* p_ply_property;
  typedef struct t_ply_argument_* p_ply_argument;

  /* ply format mode type */
  typedef enum e_ply_storage_mode_ {
    PLY_BIG_ENDIAN,
    PLY_LITTLE_ENDIAN,
    PLY_ASCII,
    PLY_DEFAULT      /* has to be the last in enum */
    } e_ply_storage_mode; /* order matches ply_storage_mode_list */

    /* ply data type */
  typedef enum e_ply_type {
    PLY_INVALID = -1,
    PLY_INT8, PLY_UINT8, PLY_INT16, PLY_UINT16,
    PLY_INT32, PLY_UIN32, PLY_FLOAT32, PLY_FLOAT64,
    PLY_CHAR, PLY_UCHAR, PLY_SHORT, PLY_USHORT,
    PLY_INT, PLY_UINT, PLY_FLOAT, PLY_DOUBLE,
    PLY_LIST    /* has to be the last in enum */
    } e_ply_type;   /* order matches ply_type_list */

    /* ----------------------------------------------------------------------
     * Error callback prototype
     *
     * message: error message
     * ply: handle returned by ply_open or ply_create
     * ---------------------------------------------------------------------- */
  typedef void (*p_ply_error_cb)(p_ply ply, const char* message);

  /* ----------------------------------------------------------------------
   * Gets user data from within an error callback
   *
   * ply: handle returned by ply_open or ply_create
   * idata,pdata: contextual information set in ply_open or ply_create
   * ---------------------------------------------------------------------- */
  JTKPLYDEF int ply_get_ply_user_data(p_ply ply, void** pdata, long* idata);

  /* ----------------------------------------------------------------------
   * Opens a PLY file for reading (fails if file is not a PLY file)
   *
   * name: file name
   * error_cb: error callback function
   * idata,pdata: contextual information available to users
   *
   * Returns 1 if successful, 0 otherwise
   * ---------------------------------------------------------------------- */
  JTKPLYDEF p_ply ply_open(const char* name, p_ply_error_cb error_cb, long idata,
    void* pdata);

  JTKPLYDEF p_ply ply_open(const wchar_t* name, p_ply_error_cb error_cb, long idata,
    void* pdata);

  /* ----------------------------------------------------------------------
   * Reads and parses the header of a PLY file returned by ply_open
   *
   * ply: handle returned by ply_open
   *
   * Returns 1 if successfull, 0 otherwise
   * ---------------------------------------------------------------------- */
  JTKPLYDEF int ply_read_header(p_ply ply);

  /* ----------------------------------------------------------------------
   * Property reading callback prototype
   *
   * argument: parameters for property being processed when callback is called
   *
   * Returns 1 if should continue processing file, 0 if should abort.
   * ---------------------------------------------------------------------- */
  typedef int (*p_ply_read_cb)(p_ply_argument argument);

  /* ----------------------------------------------------------------------
   * Sets up callbacks for property reading after header was parsed
   *
   * ply: handle returned by ply_open
   * element_name: element where property is
   * property_name: property to associate element with
   * read_cb: function to be called for each property value
   * pdata/idata: user data that will be passed to callback
   *
   * Returns 0 if no element or no property in element, returns the
   * number of element instances otherwise.
   * ---------------------------------------------------------------------- */
  JTKPLYDEF long ply_set_read_cb(p_ply ply, const char* element_name,
    const char* property_name, p_ply_read_cb read_cb,
    void* pdata, long idata);

  /* ----------------------------------------------------------------------
   * Returns information about the element originating a callback
   *
   * argument: handle to argument
   * element: receives a the element handle (if non-null)
   * instance_index: receives the index of the current element instance
   *     (if non-null)
   *
   * Returns 1 if successfull, 0 otherwise
   * ---------------------------------------------------------------------- */
  JTKPLYDEF int ply_get_argument_element(p_ply_argument argument,
    p_ply_element* element, long* instance_index);

  /* ----------------------------------------------------------------------
   * Returns information about the property originating a callback
   *
   * argument: handle to argument
   * property: receives the property handle (if non-null)
   * length: receives the number of values in this property (if non-null)
   * value_index: receives the index of current property value (if non-null)
   *
   * Returns 1 if successfull, 0 otherwise
   * ---------------------------------------------------------------------- */
  JTKPLYDEF int ply_get_argument_property(p_ply_argument argument,
    p_ply_property* property, long* length, long* value_index);

  /* ----------------------------------------------------------------------
   * Returns user data associated with callback
   *
   * pdata: receives a copy of user custom data pointer (if non-null)
   * idata: receives a copy of user custom data integer (if non-null)
   *
   * Returns 1 if successfull, 0 otherwise
   * ---------------------------------------------------------------------- */
  JTKPLYDEF int ply_get_argument_user_data(p_ply_argument argument, void** pdata,
    long* idata);

  /* ----------------------------------------------------------------------
   * Returns the value associated with a callback
   *
   * argument: handle to argument
   *
   * Returns the current data item
   * ---------------------------------------------------------------------- */
  JTKPLYDEF double ply_get_argument_value(p_ply_argument argument);

  /* ----------------------------------------------------------------------
   * Reads all elements and properties calling the callbacks defined with
   * calls to ply_set_read_cb
   *
   * ply: handle returned by ply_open
   *
   * Returns 1 if successfull, 0 otherwise
   * ---------------------------------------------------------------------- */
  JTKPLYDEF int ply_read(p_ply ply);

  /* ----------------------------------------------------------------------
   * Iterates over all elements by returning the next element.
   * Call with NULL to return handle to first element.
   *
   * ply: handle returned by ply_open
   * last: handle of last element returned (NULL for first element)
   *
   * Returns element if successfull or NULL if no more elements
   * ---------------------------------------------------------------------- */
  JTKPLYDEF p_ply_element ply_get_next_element(p_ply ply, p_ply_element last);

  /* ----------------------------------------------------------------------
   * Iterates over all comments by returning the next comment.
   * Call with NULL to return pointer to first comment.
   *
   * ply: handle returned by ply_open
   * last: pointer to last comment returned (NULL for first comment)
   *
   * Returns comment if successfull or NULL if no more comments
   * ---------------------------------------------------------------------- */
  JTKPLYDEF const char* ply_get_next_comment(p_ply ply, const char* last);

  /* ----------------------------------------------------------------------
   * Iterates over all obj_infos by returning the next obj_info.
   * Call with NULL to return pointer to first obj_info.
   *
   * ply: handle returned by ply_open
   * last: pointer to last obj_info returned (NULL for first obj_info)
   *
   * Returns obj_info if successfull or NULL if no more obj_infos
   * ---------------------------------------------------------------------- */
  JTKPLYDEF const char* ply_get_next_obj_info(p_ply ply, const char* last);

  /* ----------------------------------------------------------------------
   * Returns information about an element
   *
   * element: element of interest
   * name: receives a pointer to internal copy of element name (if non-null)
   * ninstances: receives the number of instances of this element (if non-null)
   *
   * Returns 1 if successfull or 0 otherwise
   * ---------------------------------------------------------------------- */
  JTKPLYDEF int ply_get_element_info(p_ply_element element, const char** name,
    long* ninstances);

  /* ----------------------------------------------------------------------
   * Iterates over all properties by returning the next property.
   * Call with NULL to return handle to first property.
   *
   * element: handle of element with the properties of interest
   * last: handle of last property returned (NULL for first property)
   *
   * Returns element if successfull or NULL if no more properties
   * ---------------------------------------------------------------------- */
  JTKPLYDEF p_ply_property ply_get_next_property(p_ply_element element,
    p_ply_property last);

  /* ----------------------------------------------------------------------
   * Returns information about a property
   *
   * property: handle to property of interest
   * name: receives a pointer to internal copy of property name (if non-null)
   * type: receives the property type (if non-null)
   * length_type: for list properties, receives the scalar type of
   *     the length field (if non-null)
   * value_type: for list properties, receives the scalar type of the value
   *     fields  (if non-null)
   *
   * Returns 1 if successfull or 0 otherwise
   * ---------------------------------------------------------------------- */
  JTKPLYDEF int ply_get_property_info(p_ply_property property, const char** name,
    e_ply_type* type, e_ply_type* length_type, e_ply_type* value_type);

  /* ----------------------------------------------------------------------
   * Creates new PLY file
   *
   * name: file name
   * storage_mode: file format mode
   * error_cb: error callback function
   * idata,pdata: contextual information available to users
   *
   * Returns handle to PLY file if successfull, NULL otherwise
   * ---------------------------------------------------------------------- */
  JTKPLYDEF p_ply ply_create(const char* name, e_ply_storage_mode storage_mode,
    p_ply_error_cb error_cb, long idata, void* pdata);

  /* ----------------------------------------------------------------------
   * Adds a new element to the PLY file created by ply_create
   *
   * ply: handle returned by ply_create
   * name: name of new element
   * ninstances: number of element of this time in file
   *
   * Returns 1 if successfull, 0 otherwise
   * ---------------------------------------------------------------------- */
  JTKPLYDEF int ply_add_element(p_ply ply, const char* name, long ninstances);

  /* ----------------------------------------------------------------------
   * Adds a new property to the last element added by ply_add_element
   *
   * ply: handle returned by ply_create
   * name: name of new property
   * type: property type
   * length_type: scalar type of length field of a list property
   * value_type: scalar type of value fields of a list property
   *
   * Returns 1 if successfull, 0 otherwise
   * ---------------------------------------------------------------------- */
  JTKPLYDEF int ply_add_property(p_ply ply, const char* name, e_ply_type type,
    e_ply_type length_type, e_ply_type value_type);

  /* ----------------------------------------------------------------------
   * Adds a new list property to the last element added by ply_add_element
   *
   * ply: handle returned by ply_create
   * name: name of new property
   * length_type: scalar type of length field of a list property
   * value_type: scalar type of value fields of a list property
   *
   * Returns 1 if successfull, 0 otherwise
   * ---------------------------------------------------------------------- */
  JTKPLYDEF int ply_add_list_property(p_ply ply, const char* name,
    e_ply_type length_type, e_ply_type value_type);

  /* ----------------------------------------------------------------------
   * Adds a new property to the last element added by ply_add_element
   *
   * ply: handle returned by ply_create
   * name: name of new property
   * type: property type
   *
   * Returns 1 if successfull, 0 otherwise
   * ---------------------------------------------------------------------- */
  JTKPLYDEF int ply_add_scalar_property(p_ply ply, const char* name, e_ply_type type);

  /* ----------------------------------------------------------------------
   * Adds a new comment item
   *
   * ply: handle returned by ply_create
   * comment: pointer to string with comment text
   *
   * Returns 1 if successfull, 0 otherwise
   * ---------------------------------------------------------------------- */
  JTKPLYDEF int ply_add_comment(p_ply ply, const char* comment);

  /* ----------------------------------------------------------------------
   * Adds a new obj_info item
   *
   * ply: handle returned by ply_create
   * comment: pointer to string with obj_info data
   *
   * Returns 1 if successfull, 0 otherwise
   * ---------------------------------------------------------------------- */
  JTKPLYDEF int ply_add_obj_info(p_ply ply, const char* obj_info);

  /* ----------------------------------------------------------------------
   * Writes the PLY file header after all element and properties have been
   * defined by calls to ply_add_element and ply_add_property
   *
   * ply: handle returned by ply_create
   *
   * Returns 1 if successfull, 0 otherwise
   * ---------------------------------------------------------------------- */
  JTKPLYDEF int ply_write_header(p_ply ply);

  /* ----------------------------------------------------------------------
   * Writes one property value, in the order they should be written to the
   * file. For each element type, write all elements of that type in order.
   * For each element, write all its properties in order. For scalar
   * properties, just write the value. For list properties, write the length
   * and then each of the values.
   *
   * ply: handle returned by ply_create
   *
   * Returns 1 if successfull, 0 otherwise
   * ---------------------------------------------------------------------- */
  JTKPLYDEF int ply_write(p_ply ply, double value);

  /* ----------------------------------------------------------------------
   * Closes a PLY file handle. Releases all memory used by handle
   *
   * ply: handle to be closed.
   *
   * Returns 1 if successfull, 0 otherwise
   * ---------------------------------------------------------------------- */
  JTKPLYDEF int ply_close(p_ply ply);

  /* ----------------------------------------------------------------------
   * Opens a PLY file for reading (fails if file is not a PLY file)
   *
   * file_pointer: FILE * to file open for reading
   * error_cb: error callback function
   * idata,pdata: contextual information available to users
   *
   * Returns 1 if successful, 0 otherwise
   * ---------------------------------------------------------------------- */
  JTKPLYDEF p_ply ply_open_from_file(FILE* file_pointer, p_ply_error_cb error_cb,
    long idata, void* pdata);

  /* ----------------------------------------------------------------------
   * Creates new PLY file
   *
   * file_pointer: FILE * to a file open for writing
   * storage_mode: file format mode
   * error_cb: error callback function
   * idata,pdata: contextual information available to users
   *
   * Returns handle to PLY file if successfull, NULL otherwise
   * ---------------------------------------------------------------------- */
  JTKPLYDEF p_ply ply_create_to_file(FILE* file_pointer, e_ply_storage_mode storage_mode,
    p_ply_error_cb error_cb, long idata, void* pdata);


  /* ----------------------------------------------------------------------
 * Opens a PLY file for reading from memory (buffer of char[])
 *
 * buffer      : file in the memory for reading
 * error_cb    : error callback function
 * idata,pdata : contextual information available to users
 *
 * Returns 1 if successful, 0 otherwise
 * ---------------------------------------------------------------------- */

 /*
 Memory reading fixed by Jan Maes.
   - Added the size of the input buffer (old version could go out of memory bounds)
   - Added checks to see that we don't go out of buffer bounds.
 */
  JTKPLYDEF p_ply ply_open_from_memory(const char* buffer, size_t buffer_size, p_ply_error_cb error_cb,
    long idata, void* pdata);

  /* ----------------------------------------------------------------------
   * Creates new PLY file for writing in memory (to a buffer of char[])
   *
   * buffer       : memory for storing the file
   * buffer_size  : size of pre-allocated memory for the buffer
   * ply_size     : actual size of the PLY file in memory
   * storage_mode : file format mode
   * error_cb     : error callback function
   * idata, pdata : contextual information available to users
   *
   * Returns handle to PLY file if successful, NULL otherwise
   * ---------------------------------------------------------------------- */
  JTKPLYDEF p_ply ply_create_to_memory(char* buffer, size_t buffer_size, size_t* ply_size, e_ply_storage_mode storage_mode,
    p_ply_error_cb error_cb, long idata, void* pdata);


  }


#endif // #ifndef JTK_PLY_H

#ifdef JTK_PLY_IMPLEMENTATION

#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>

/* ----------------------------------------------------------------------
 * Make sure we get our integer types right
 * ---------------------------------------------------------------------- */
#if defined(_MSC_VER) && (_MSC_VER < 1600)
 /* C99 stdint.h only supported in MSVC++ 10.0 and up */
typedef __int8 t_ply_int8;
typedef __int16 t_ply_int16;
typedef __int32 t_ply_int32;
typedef unsigned __int8 t_ply_uint8;
typedef unsigned __int16 t_ply_uint16;
typedef unsigned __int32 t_ply_uint32;
#define PLY_INT8_MAX (127)
#define PLY_INT8_MIN (-PLY_INT8_MAX-1)
#define PLY_INT16_MAX (32767)
#define PLY_INT16_MIN (-PLY_INT16_MAX-1)
#define PLY_INT32_MAX (2147483647)
#define PLY_INT32_MIN (-PLY_INT32_MAX-1)
#define PLY_UINT8_MAX (255)
#define PLY_UINT16_MAX (65535)
#define PLY_UINT32_MAX  (4294967295)
#else
#include <stdint.h>
typedef int8_t t_ply_int8;
typedef int16_t t_ply_int16;
typedef int32_t t_ply_int32;
typedef uint8_t t_ply_uint8;
typedef uint16_t t_ply_uint16;
typedef uint32_t t_ply_uint32;
#define PLY_INT8_MIN INT8_MIN
#define PLY_INT8_MAX INT8_MAX
#define PLY_INT16_MIN INT16_MIN
#define PLY_INT16_MAX INT16_MAX
#define PLY_INT32_MIN INT32_MIN
#define PLY_INT32_MAX INT32_MAX
#define PLY_UINT8_MAX UINT8_MAX
#define PLY_UINT16_MAX UINT16_MAX
#define PLY_UINT32_MAX UINT32_MAX
#endif

namespace jtk
  {

  namespace ply_details
    {

    static int read_vec3_coord(p_ply_argument argument)
      {
      float** p_vertices;
      ply_get_argument_user_data(argument, (void**)&p_vertices, NULL);
      double val = ply_get_argument_value(argument);
      *(*p_vertices) = (float)val;
      (*p_vertices) += 3;
      return 1;
      }

    static int read_color(p_ply_argument argument)
      {
      uint8_t** p_color;
      ply_get_argument_user_data(argument, (void**)&p_color, NULL);
      double val = ply_get_argument_value(argument);
      *(*p_color) = (uint8_t)val;
      (*p_color) += 4;
      return 1;
      }

    static int read_face(p_ply_argument argument)
      {
      uint32_t** p_triangle;
      ply_get_argument_user_data(argument, (void**)&p_triangle, NULL);
      long length, value_index;
      ply_get_argument_property(argument, NULL, &length, &value_index);
      if (value_index >= 0 && value_index < 3)
        {
        double val = ply_get_argument_value(argument);
        *(*p_triangle) = (uint32_t)val;
        ++(*p_triangle);
        }
      return 1;
      }

    static int read_texcoord(p_ply_argument argument)
      {
      float** p_uv;
      ply_get_argument_user_data(argument, (void**)&p_uv, NULL);
      long length, value_index;
      ply_get_argument_property(argument, NULL, &length, &value_index);
      if (value_index >= 0 && value_index < 6)
        {
        double val = ply_get_argument_value(argument);
        *(*p_uv) = (float)val;
        ++(*p_uv);
        }
      if (value_index == (length - 1) && (length != 6))
        {
        for (long j = length; j < 6; ++j)
          {
          *(*p_uv) = (float)0.f;
          ++(*p_uv);
          }
        }
      return 1;
      }

    }

  template <class TCHAR>
  bool _read_ply(const TCHAR* filename, std::vector<vec3<float>>& vertices, std::vector<vec3<float>>& normals, std::vector<uint32_t>& clrs, std::vector<vec3<uint32_t>>& triangles, std::vector<vec3<vec2<float>>>& uv)
    {
    vertices.clear();
    normals.clear();
    clrs.clear();
    triangles.clear();
    uv.clear();

    p_ply ply = ply_open(filename, NULL, 0, NULL);
    if (!ply)
      return false;
    if (!ply_read_header(ply))
      return false;

    float* p_vertex_pointer_x = NULL;
    float* p_vertex_pointer_y = NULL;
    float* p_vertex_pointer_z = NULL;

    long nvertices_x = ply_set_read_cb(ply, "vertex", "x", ply_details::read_vec3_coord, (void*)(&p_vertex_pointer_x), 0);
    long nvertices_y = ply_set_read_cb(ply, "vertex", "y", ply_details::read_vec3_coord, (void*)(&p_vertex_pointer_y), 0);
    long nvertices_z = ply_set_read_cb(ply, "vertex", "z", ply_details::read_vec3_coord, (void*)(&p_vertex_pointer_z), 0);

    if (nvertices_x != nvertices_y)
      return false;
    if (nvertices_x != nvertices_z)
      return false;

    if (nvertices_x > 0)
      vertices.resize(nvertices_x);
    p_vertex_pointer_x = (float*)vertices.data();
    p_vertex_pointer_y = p_vertex_pointer_x + 1;
    p_vertex_pointer_z = p_vertex_pointer_x + 2;

    float* p_normal_pointer_x = NULL;
    float* p_normal_pointer_y = NULL;
    float* p_normal_pointer_z = NULL;

    long nnormals_x = ply_set_read_cb(ply, "vertex", "nx", ply_details::read_vec3_coord, (void*)(&p_normal_pointer_x), 0);
    long nnormals_y = ply_set_read_cb(ply, "vertex", "ny", ply_details::read_vec3_coord, (void*)(&p_normal_pointer_y), 0);
    long nnormals_z = ply_set_read_cb(ply, "vertex", "nz", ply_details::read_vec3_coord, (void*)(&p_normal_pointer_z), 0);

    if (nnormals_x != nnormals_y)
      return false;
    if (nnormals_x != nnormals_z)
      return false;

    if (nnormals_x > 0)
      normals.resize(nnormals_x);
    p_normal_pointer_x = (float*)normals.data();
    p_normal_pointer_y = p_normal_pointer_x + 1;
    p_normal_pointer_z = p_normal_pointer_x + 2;

    uint8_t* p_red = NULL;
    uint8_t* p_green = NULL;
    uint8_t* p_blue = NULL;
    uint8_t* p_alpha = NULL;

    long nred = ply_set_read_cb(ply, "vertex", "red", ply_details::read_color, (void*)(&p_red), 0);
    long ngreen = ply_set_read_cb(ply, "vertex", "green", ply_details::read_color, (void*)(&p_green), 0);
    long nblue = ply_set_read_cb(ply, "vertex", "blue", ply_details::read_color, (void*)(&p_blue), 0);
    long nalpha = ply_set_read_cb(ply, "vertex", "alpha", ply_details::read_color, (void*)(&p_alpha), 0);

    if (nred == 0)
      nred = ply_set_read_cb(ply, "vertex", "r", ply_details::read_color, (void*)(&p_red), 0);
    if (ngreen == 0)
      ngreen = ply_set_read_cb(ply, "vertex", "g", ply_details::read_color, (void*)(&p_green), 0);
    if (nblue == 0)
      nblue = ply_set_read_cb(ply, "vertex", "b", ply_details::read_color, (void*)(&p_blue), 0);
    if (nalpha == 0)
      nalpha = ply_set_read_cb(ply, "vertex", "a", ply_details::read_color, (void*)(&p_alpha), 0);

    if (nred == 0)
      nred = ply_set_read_cb(ply, "vertex", "diffuse_red", ply_details::read_color, (void*)(&p_red), 0);
    if (ngreen == 0)
      ngreen = ply_set_read_cb(ply, "vertex", "diffuse_green", ply_details::read_color, (void*)(&p_green), 0);
    if (nblue == 0)
      nblue = ply_set_read_cb(ply, "vertex", "diffuse_blue", ply_details::read_color, (void*)(&p_blue), 0);
    if (nalpha == 0)
      nalpha = ply_set_read_cb(ply, "vertex", "diffuse_alpha", ply_details::read_color, (void*)(&p_alpha), 0);

    if (nred > 0 || ngreen > 0 || nblue > 0 || nalpha > 0)
      clrs.resize(nred, 0xffffffff);

    p_red = (uint8_t*)clrs.data();
    p_green = p_red + 1;
    p_blue = p_red + 2;
    p_alpha = p_red + 3;

    uint32_t* p_tria_index = NULL;

    long ntriangles = ply_set_read_cb(ply, "face", "vertex_indices", ply_details::read_face, (void*)(&p_tria_index), 0);
    if (ntriangles == 0)
      ntriangles = ply_set_read_cb(ply, "face", "vertex_index", ply_details::read_face, (void*)(&p_tria_index), 0);

    if (ntriangles > 0)
      triangles.resize(ntriangles);

    p_tria_index = (uint32_t*)triangles.data();

    float* p_uv = NULL;

    long ntexcoords = ply_set_read_cb(ply, "face", "texcoord", ply_details::read_texcoord, (void*)(&p_uv), 0);

    if (ntexcoords > 0)
      uv.resize(ntexcoords);

    p_uv = (float*)uv.data();

    if (!ply_read(ply))
      return false;

    ply_close(ply);

    return true;
    }

  JTKPLYDEF bool read_ply_from_memory(const char* buffer, uint64_t buffer_size, std::vector<vec3<float>>& vertices, std::vector<vec3<float>>& normals, std::vector<uint32_t>& clrs, std::vector<vec3<uint32_t>>& triangles, std::vector<vec3<vec2<float>>>& uv)
    {
    vertices.clear();
    normals.clear();
    clrs.clear();
    triangles.clear();
    uv.clear();

    p_ply ply = ply_open_from_memory(buffer, buffer_size, nullptr, 0, nullptr);
    if (!ply)
      return false;
    if (!ply_read_header(ply))
      return false;

    float* p_vertex_pointer_x = NULL;
    float* p_vertex_pointer_y = NULL;
    float* p_vertex_pointer_z = NULL;

    long nvertices_x = ply_set_read_cb(ply, "vertex", "x", ply_details::read_vec3_coord, (void*)(&p_vertex_pointer_x), 0);
    long nvertices_y = ply_set_read_cb(ply, "vertex", "y", ply_details::read_vec3_coord, (void*)(&p_vertex_pointer_y), 0);
    long nvertices_z = ply_set_read_cb(ply, "vertex", "z", ply_details::read_vec3_coord, (void*)(&p_vertex_pointer_z), 0);

    if (nvertices_x != nvertices_y)
      return false;
    if (nvertices_x != nvertices_z)
      return false;

    if (nvertices_x > 0)
      vertices.resize(nvertices_x);
    p_vertex_pointer_x = (float*)vertices.data();
    p_vertex_pointer_y = p_vertex_pointer_x + 1;
    p_vertex_pointer_z = p_vertex_pointer_x + 2;

    float* p_normal_pointer_x = NULL;
    float* p_normal_pointer_y = NULL;
    float* p_normal_pointer_z = NULL;

    long nnormals_x = ply_set_read_cb(ply, "vertex", "nx", ply_details::read_vec3_coord, (void*)(&p_normal_pointer_x), 0);
    long nnormals_y = ply_set_read_cb(ply, "vertex", "ny", ply_details::read_vec3_coord, (void*)(&p_normal_pointer_y), 0);
    long nnormals_z = ply_set_read_cb(ply, "vertex", "nz", ply_details::read_vec3_coord, (void*)(&p_normal_pointer_z), 0);

    if (nnormals_x != nnormals_y)
      return false;
    if (nnormals_x != nnormals_z)
      return false;

    if (nnormals_x > 0)
      normals.resize(nnormals_x);
    p_normal_pointer_x = (float*)normals.data();
    p_normal_pointer_y = p_normal_pointer_x + 1;
    p_normal_pointer_z = p_normal_pointer_x + 2;

    uint8_t* p_red = NULL;
    uint8_t* p_green = NULL;
    uint8_t* p_blue = NULL;
    uint8_t* p_alpha = NULL;

    long nred = ply_set_read_cb(ply, "vertex", "red", ply_details::read_color, (void*)(&p_red), 0);
    long ngreen = ply_set_read_cb(ply, "vertex", "green", ply_details::read_color, (void*)(&p_green), 0);
    long nblue = ply_set_read_cb(ply, "vertex", "blue", ply_details::read_color, (void*)(&p_blue), 0);
    long nalpha = ply_set_read_cb(ply, "vertex", "alpha", ply_details::read_color, (void*)(&p_alpha), 0);

    if (nred == 0)
      nred = ply_set_read_cb(ply, "vertex", "r", ply_details::read_color, (void*)(&p_red), 0);
    if (ngreen == 0)
      ngreen = ply_set_read_cb(ply, "vertex", "g", ply_details::read_color, (void*)(&p_green), 0);
    if (nblue == 0)
      nblue = ply_set_read_cb(ply, "vertex", "b", ply_details::read_color, (void*)(&p_blue), 0);
    if (nalpha == 0)
      nalpha = ply_set_read_cb(ply, "vertex", "a", ply_details::read_color, (void*)(&p_alpha), 0);

    if (nred == 0)
      nred = ply_set_read_cb(ply, "vertex", "diffuse_red", ply_details::read_color, (void*)(&p_red), 0);
    if (ngreen == 0)
      ngreen = ply_set_read_cb(ply, "vertex", "diffuse_green", ply_details::read_color, (void*)(&p_green), 0);
    if (nblue == 0)
      nblue = ply_set_read_cb(ply, "vertex", "diffuse_blue", ply_details::read_color, (void*)(&p_blue), 0);
    if (nalpha == 0)
      nalpha = ply_set_read_cb(ply, "vertex", "diffuse_alpha", ply_details::read_color, (void*)(&p_alpha), 0);

    if (nred > 0 || ngreen > 0 || nblue > 0 || nalpha > 0)
      clrs.resize(nred, 0xffffffff);

    p_red = (uint8_t*)clrs.data();
    p_green = p_red + 1;
    p_blue = p_red + 2;
    p_alpha = p_red + 3;

    uint32_t* p_tria_index = NULL;

    long ntriangles = ply_set_read_cb(ply, "face", "vertex_indices", ply_details::read_face, (void*)(&p_tria_index), 0);
    if (ntriangles == 0)
      ntriangles = ply_set_read_cb(ply, "face", "vertex_index", ply_details::read_face, (void*)(&p_tria_index), 0);

    if (ntriangles > 0)
      triangles.resize(ntriangles);

    p_tria_index = (uint32_t*)triangles.data();

    float* p_uv = NULL;

    long ntexcoords = ply_set_read_cb(ply, "face", "texcoord", ply_details::read_texcoord, (void*)(&p_uv), 0);

    if (ntexcoords > 0)
      uv.resize(ntexcoords);

    p_uv = (float*)uv.data();

    if (!ply_read(ply))
      return false;

    ply_close(ply);

    return true;
    }

  namespace ply_details
    {
    template <class TCHAR>
    class file_opener
      {
      public:
        FILE* operator()(const TCHAR* filename, const char* mode)
          {
          return fopen(filename, mode);
          }
      };

    template <>
    class file_opener<wchar_t>
      {
      public:
        FILE* operator()(const wchar_t* filename, const char* mode)
          {
          std::string m(mode);
          std::wstring wm(m.begin(), m.end());
          return _wfopen(filename, wm.c_str());
          }
      };
    }

  template <class TCHAR>
  bool _write_ply(const TCHAR* filename, const std::vector<vec3<float>>& vertices, const std::vector<vec3<float>>& normals, const std::vector<uint32_t>& clrs, const std::vector<vec3<uint32_t>>& triangles, const std::vector<vec3<vec2<float>>>& uv)
    {
    ply_details::file_opener<TCHAR> fo;
    if (vertices.empty())
      return false;
    FILE* fp = fo(filename, "wb");

    if (!fp)
      return false;

    fprintf(fp, "ply\n");
    int n = 1;
    if (*(char*)&n == 1)
      fprintf(fp, "format binary_little_endian 1.0\n");
    else
      fprintf(fp, "format binary_big_endian 1.0\n");

    fprintf(fp, "element vertex %d\n", (uint32_t)vertices.size());
    fprintf(fp, "property float x\n");
    fprintf(fp, "property float y\n");
    fprintf(fp, "property float z\n");

    if (!normals.empty())
      {
      fprintf(fp, "property float nx\n");
      fprintf(fp, "property float ny\n");
      fprintf(fp, "property float nz\n");
      }

    if (!clrs.empty())
      {
      fprintf(fp, "property uchar red\n");
      fprintf(fp, "property uchar green\n");
      fprintf(fp, "property uchar blue\n");
      fprintf(fp, "property uchar alpha\n");
      }

    if (!triangles.empty())
      {
      fprintf(fp, "element face %d\n", (uint32_t)triangles.size());
      fprintf(fp, "property list uchar int vertex_indices\n");
      if (!uv.empty())
        fprintf(fp, "property list uchar float texcoord\n");
      }
    fprintf(fp, "end_header\n");

    for (uint32_t i = 0; i < (uint32_t)vertices.size(); ++i)
      {
      fwrite((float*)vertices.data() + 3 * i, sizeof(float), 3, fp);
      if (!normals.empty())
        fwrite((float*)normals.data() + 3 * i, sizeof(float), 3, fp);
      if (!clrs.empty())
        fwrite((uint32_t*)clrs.data() + i, sizeof(uint32_t), 1, fp);
      }
    const unsigned char tria_size = 3;
    const unsigned char texcoord_size = 6;
    for (uint32_t i = 0; i < (uint32_t)triangles.size(); ++i)
      {
      fwrite(&tria_size, 1, 1, fp);
      fwrite((uint32_t*)triangles.data() + 3 * i, sizeof(uint32_t), 3, fp);
      if (!uv.empty())
        {
        fwrite(&texcoord_size, 1, 1, fp);
        fwrite((float*)uv.data() + 6 * i, sizeof(float), 6, fp);
        }
      }
    fclose(fp);
    return true;
    }


  /* ----------------------------------------------------------------------
   * Constants
   * ---------------------------------------------------------------------- */
#define WORDSIZE 256
#define LINESIZE 1024
#define BUFFERSIZE (8*1024)

  typedef enum e_ply_io_mode_ {
    PLY_READ,
    PLY_WRITE
    } e_ply_io_mode;

  static const char* const ply_storage_mode_list[] = {
      "binary_big_endian", "binary_little_endian", "ascii", NULL
    };     /* order matches e_ply_storage_mode enum */

  static const char* const ply_type_list[] = {
      "int8", "uint8", "int16", "uint16",
      "int32", "uint32", "float32", "float64",
      "char", "uchar", "short", "ushort",
      "int", "uint", "float", "double",
      "list", NULL
    };     /* order matches e_ply_type enum */

    /* ----------------------------------------------------------------------
     * Property reading callback argument
     *
     * element: name of element being processed
     * property: name of property being processed
     * nelements: number of elements of this kind in file
     * instance_index: index current element of this kind being processed
     * length: number of values in current list (or 1 for scalars)
     * value_index: index of current value int this list (or 0 for scalars)
     * value: value of property
     * pdata/idata: user data defined with ply_set_cb
     *
     * Returns handle to PLY file if succesful, NULL otherwise.
     * ---------------------------------------------------------------------- */
  typedef struct t_ply_argument_ {
    p_ply_element element;
    long instance_index;
    p_ply_property property;
    long length, value_index;
    double value;
    void* pdata;
    long idata;
    } t_ply_argument;

  /* ----------------------------------------------------------------------
   * Property information
   *
   * name: name of this property
   * type: type of this property (list or type of scalar value)
   * length_type, value_type: type of list property count and values
   * read_cb: function to be called when this property is called
   *
   * Returns 1 if should continue processing file, 0 if should abort.
   * ---------------------------------------------------------------------- */
  typedef struct t_ply_property_ {
    char name[WORDSIZE];
    e_ply_type type, value_type, length_type;
    p_ply_read_cb read_cb;
    void* pdata;
    long idata;
    } t_ply_property;

  /* ----------------------------------------------------------------------
   * Element information
   *
   * name: name of this property
   * ninstances: number of elements of this type in file
   * property: property descriptions for this element
   * nproperty: number of properties in this element
   *
   * Returns 1 if should continue processing file, 0 if should abort.
   * ---------------------------------------------------------------------- */
  typedef struct t_ply_element_ {
    char name[WORDSIZE];
    long ninstances;
    p_ply_property property;
    long nproperties;
    } t_ply_element;

  /* ----------------------------------------------------------------------
   * Input/output driver
   *
   * Depending on file mode, different functions are used to read/write
   * property fields. The drivers make it transparent to read/write in ascii,
   * big endian or little endian cases.
   * ---------------------------------------------------------------------- */
  typedef int (*p_ply_ihandler)(p_ply ply, double* value);
  typedef int (*p_ply_ichunk)(p_ply ply, void* anydata, size_t size);
  typedef struct t_ply_idriver_ {
    p_ply_ihandler ihandler[16];
    p_ply_ichunk ichunk;
    const char* name;
    } t_ply_idriver;
  typedef t_ply_idriver* p_ply_idriver;

  typedef int (*p_ply_ohandler)(p_ply ply, double value);
  typedef int (*p_ply_ochunk)(p_ply ply, void* anydata, size_t size);
  typedef struct t_ply_odriver_ {
    p_ply_ohandler ohandler[16];
    p_ply_ochunk ochunk;
    const char* name;
    } t_ply_odriver;
  typedef t_ply_odriver* p_ply_odriver;

  /* ----------------------------------------------------------------------
   * Ply file handle.
   *
   * io_mode: read or write (from e_ply_io_mode)
   * storage_mode: mode of file associated with handle (from e_ply_storage_mode)
   * element: elements description for this file
   * nelement: number of different elements in file
   * comment: comments for this file
   * ncomments: number of comments in file
   * obj_info: obj_info items for this file
   * nobj_infos: number of obj_info items in file
   * fp: file pointer associated with ply file
   * rn: skip extra char after end_header?
   * buffer: last word/chunck of data read from ply file
   * buffer_first, buffer_last: interval of untouched good data in buffer
   * buffer_token: start of parsed token (line or word) in buffer
   * idriver, odriver: input driver used to get property fields from file
   * argument: storage space for callback arguments
   * welement, wproperty: element/property type being written
   * winstance_index: index of instance of current element being written
   * wvalue_index: index of list property value being written
   * wlength: number of values in list property being written
   * error_cb: error callback
   * pdata/idata: user data defined with ply_open/ply_create
   * ---------------------------------------------------------------------- */
  typedef struct t_ply_ {
    bool in_memory;
    const char* load_memory;
    size_t load_memory_size;
    char* store_memory;
    size_t* store_memory_size;
    size_t store_buffer_size;
    e_ply_io_mode io_mode;
    e_ply_storage_mode storage_mode;
    p_ply_element element;
    long nelements;
    char* comment;
    long ncomments;
    char* obj_info;
    long nobj_infos;
    FILE* fp;
    int own_fp;
    int rn;
    char buffer[BUFFERSIZE];
    size_t buffer_first, buffer_token, buffer_last;
    p_ply_idriver idriver;
    p_ply_odriver odriver;
    t_ply_argument argument;
    long welement, wproperty;
    long winstance_index, wvalue_index, wlength;
    p_ply_error_cb error_cb;
    void* pdata;
    long idata;
    } t_ply;

  /* ----------------------------------------------------------------------
   * I/O functions and drivers
   * ---------------------------------------------------------------------- */

  namespace
    {
    extern t_ply_idriver ply_idriver_ascii;
    extern t_ply_idriver ply_idriver_binary;
    extern t_ply_idriver ply_idriver_binary_reverse;
    extern t_ply_odriver ply_odriver_ascii;
    extern t_ply_odriver ply_odriver_binary;
    extern t_ply_odriver ply_odriver_binary_reverse;
    }

  static int ply_read_word(p_ply ply);
  static int ply_check_word(p_ply ply);
  static void ply_finish_word(p_ply ply, size_t size);
  static int ply_read_line(p_ply ply);
  static int ply_check_line(p_ply ply);
  static int ply_read_chunk(p_ply ply, void* anybuffer, size_t size);
  static int ply_read_chunk_reverse(p_ply ply, void* anybuffer, size_t size);
  static int ply_write_chunk(p_ply ply, void* anybuffer, size_t size);
  static int ply_write_chunk_reverse(p_ply ply, void* anybuffer, size_t size);
  static void ply_reverse(void* anydata, size_t size);

  /* ----------------------------------------------------------------------
   * In memory I/O functions
   * ---------------------------------------------------------------------- */
  static int ply_store_buffer_in_memory(p_ply ply, const char* buffer);

  /* ----------------------------------------------------------------------
   * String functions
   * ---------------------------------------------------------------------- */
  static int ply_find_string(const char* item, const char* const list[]);
  static p_ply_element ply_find_element(p_ply ply, const char* name);
  static p_ply_property ply_find_property(p_ply_element element,
    const char* name);

  /* ----------------------------------------------------------------------
   * Header parsing
   * ---------------------------------------------------------------------- */
  static int ply_read_header_magic(p_ply ply);
  static int ply_read_header_format(p_ply ply);
  static int ply_read_header_comment(p_ply ply);
  static int ply_read_header_obj_info(p_ply ply);
  static int ply_read_header_property(p_ply ply);
  static int ply_read_header_element(p_ply ply);

  /* ----------------------------------------------------------------------
   * Error handling
   * ---------------------------------------------------------------------- */
  static void ply_error_cb(p_ply ply, const char* message);
  static void ply_ferror(p_ply ply, const char* fmt, ...);

  /* ----------------------------------------------------------------------
   * Memory allocation and initialization
   * ---------------------------------------------------------------------- */
  static void ply_init(p_ply ply);
  static void ply_element_init(p_ply_element element);
  static void ply_property_init(p_ply_property property);
  static p_ply ply_alloc(void);
  static p_ply_element ply_grow_element(p_ply ply);
  static p_ply_property ply_grow_property(p_ply ply, p_ply_element element);
  static void* ply_grow_array(p_ply ply, void** pointer, long* nmemb, long size);

  /* ----------------------------------------------------------------------
   * Special functions
   * ---------------------------------------------------------------------- */
  static e_ply_storage_mode ply_arch_endian(void);
  static int ply_type_check(void);

  /* ----------------------------------------------------------------------
   * Auxiliary read functions
   * ---------------------------------------------------------------------- */
  static int ply_read_element(p_ply ply, p_ply_element element,
    p_ply_argument argument);
  static int ply_read_property(p_ply ply, p_ply_element element,
    p_ply_property property, p_ply_argument argument);
  static int ply_read_list_property(p_ply ply, p_ply_element element,
    p_ply_property property, p_ply_argument argument);
  static int ply_read_scalar_property(p_ply ply, p_ply_element element,
    p_ply_property property, p_ply_argument argument);

  /* ----------------------------------------------------------------------
   * Buffer support functions
   * ---------------------------------------------------------------------- */
   /* pointers to tokenized word and line in buffer */
#define BWORD(p) (p->buffer + p->buffer_token)
#define BLINE(p) (p->buffer + p->buffer_token)

/* pointer to start of untouched bytes in buffer */
#define BFIRST(p) (p->buffer + p->buffer_first)

/* number of bytes untouched in buffer */
#define BSIZE(p) (p->buffer_last - p->buffer_first)

/* consumes data from buffer */
#define BSKIP(p, s) (p->buffer_first += s)

/* refills the buffer */
  static int BREFILL(p_ply ply) {
    /* move untouched data to beginning of buffer */
    size_t size = BSIZE(ply);
    memmove(ply->buffer, BFIRST(ply), size);
    ply->buffer_last = size;
    ply->buffer_first = ply->buffer_token = 0;
    /* fill remaining with new data */
    if (ply->in_memory) { // Memory reading fixed by Jan Maes
      size_t to_read = BUFFERSIZE - size - 1;
      size_t chunksize = ply->load_memory_size > to_read ? to_read : ply->load_memory_size;
      memcpy(ply->buffer + size, ply->load_memory, chunksize);
      ply->load_memory += chunksize;
      ply->load_memory_size -= chunksize;
      size = chunksize;
      }
    else {
      size = fread(ply->buffer + size, 1, BUFFERSIZE - size - 1, ply->fp);
      }
    /* increase size to account for new data */
    ply->buffer_last += size;
    /* place sentinel so we can use str* functions with buffer */
    ply->buffer[ply->buffer_last] = '\0';
    /* check if read failed */
    return size > 0;
    }

  /* We don't care about end-of-line, generally, because we
   * separate words by any white-space character.
   * Unfortunately, in binary mode, right after 'end_header',
   * we have to know *exactly* how many characters to skip */
   /* We use the end-of-line marker after the 'ply' magic
    * number to figure out what to do */
  static int ply_read_header_magic(p_ply ply) {
    char* magic = ply->buffer;
    if (!BREFILL(ply)) {
      ply->error_cb(ply, "Unable to read magic number from file");
      return 0;
      }
    /* check if it is ply */
    if (magic[0] != 'p' || magic[1] != 'l' || magic[2] != 'y'
      || !isspace(magic[3])) {
      ply->error_cb(ply, "Wrong magic number. Expected 'ply'");
      return 0;
      }
    /* figure out if we have to skip the extra character
     * after header when we reach the binary part of file */
    ply->rn = magic[3] == '\r' && magic[4] == '\n';
    BSKIP(ply, 3);
    return 1;
    }

  /* ----------------------------------------------------------------------
   * Exported functions
   * ---------------------------------------------------------------------- */
   /* ----------------------------------------------------------------------
    * Read support functions
    * ---------------------------------------------------------------------- */
  JTKPLYDEF p_ply ply_open(const char* name, p_ply_error_cb error_cb,
    long idata, void* pdata) {
    FILE* fp;
    p_ply ply;
    if (error_cb == NULL) error_cb = ply_error_cb;
    assert(name);
    fp = fopen(name, "rb");
    if (!fp) {
      error_cb(NULL, "Unable to open file");
      return NULL;
      }
    ply = ply_open_from_file(fp, error_cb, idata, pdata);
    if (ply) ply->own_fp = 1;
    else fclose(fp);
    return ply;
    }

  JTKPLYDEF p_ply ply_open(const wchar_t* name, p_ply_error_cb error_cb,
    long idata, void* pdata) {
    FILE* fp;
    p_ply ply;
    if (error_cb == NULL) error_cb = ply_error_cb;
    assert(name);
    fp = _wfopen(name, L"rb");
    if (!fp) {
      error_cb(NULL, "Unable to open file");
      return NULL;
      }
    ply = ply_open_from_file(fp, error_cb, idata, pdata);
    if (ply) ply->own_fp = 1;
    else fclose(fp);
    return ply;
    }

  JTKPLYDEF p_ply ply_open_from_file(FILE* fp, p_ply_error_cb error_cb,
    long idata, void* pdata) {
    p_ply ply = nullptr;
    if (error_cb == NULL) error_cb = ply_error_cb;
    assert(fp);
    if (!ply_type_check()) {
      error_cb(ply, "Incompatible type system");
      return NULL;
      }
    ply = ply_alloc();
    if (!ply) {
      error_cb(NULL, "Out of memory");
      return NULL;
      }
    ply->in_memory = false;
    ply->idata = idata;
    ply->pdata = pdata;
    ply->io_mode = PLY_READ;
    ply->error_cb = error_cb;
    ply->fp = fp;
    ply->own_fp = 0;
    return ply;
    }

  JTKPLYDEF p_ply ply_open_from_memory(const char* buffer, size_t buffer_size, p_ply_error_cb error_cb, long idata,
    void* pdata) {
    p_ply ply = nullptr;
    if (error_cb == NULL) error_cb = ply_error_cb;
    assert(buffer);
    if (!ply_type_check()) {
      error_cb(ply, "Incompatible type system");
      return NULL;
      }
    ply = ply_alloc();
    if (!ply) {
      error_cb(NULL, "Out of memory");
      return NULL;
      }
    ply->in_memory = true;
    ply->load_memory = buffer;
    ply->load_memory_size = buffer_size;
    ply->idata = idata;
    ply->pdata = pdata;
    ply->io_mode = PLY_READ;
    ply->error_cb = error_cb;
    return ply;
    }

  JTKPLYDEF int ply_read_header(p_ply ply) {
    assert(ply && ply->io_mode == PLY_READ);
    if (!ply->in_memory) {
      assert(ply->fp);
      }
    if (!ply_read_header_magic(ply)) return 0;
    if (!ply_read_word(ply)) return 0;
    /* parse file format */
    if (!ply_read_header_format(ply)) {
      ply_ferror(ply, "Invalid file format");
      return 0;
      }
    /* parse elements, comments or obj_infos until the end of header */
    while (strcmp(BWORD(ply), "end_header")) {
      if (!ply_read_header_comment(ply) &&
        !ply_read_header_element(ply) &&
        !ply_read_header_obj_info(ply)) {
        ply_ferror(ply, "Unexpected token '%s'", BWORD(ply));
        return 0;
        }
      }
    /* skip extra character? */
    if (ply->rn) {
      if (BSIZE(ply) < 1 && !BREFILL(ply)) {
        ply_ferror(ply, "Unexpected end of file");
        return 0;
        }
      BSKIP(ply, 1);
      }
    return 1;
    }

  JTKPLYDEF long ply_set_read_cb(p_ply ply, const char* element_name,
    const char* property_name, p_ply_read_cb read_cb,
    void* pdata, long idata) {
    p_ply_element element = NULL;
    p_ply_property property = NULL;
    assert(ply && element_name && property_name);
    element = ply_find_element(ply, element_name);
    if (!element) return 0;
    property = ply_find_property(element, property_name);
    if (!property) return 0;
    property->read_cb = read_cb;
    property->pdata = pdata;
    property->idata = idata;
    return (int)element->ninstances;
    }

  JTKPLYDEF int ply_read(p_ply ply) {
    long i;
    p_ply_argument argument;
    assert(ply && ply->io_mode == PLY_READ);
    if (!ply->in_memory) {
      assert(ply->fp);
      }
    argument = &ply->argument;
    /* for each element type */
    for (i = 0; i < ply->nelements; i++) {
      p_ply_element element = &ply->element[i];
      argument->element = element;
      if (!ply_read_element(ply, element, argument))
        return 0;
      }
    return 1;
    }

  /* ----------------------------------------------------------------------
   * Write support functions
   * ---------------------------------------------------------------------- */
  JTKPLYDEF p_ply ply_create(const char* name, e_ply_storage_mode storage_mode,
    p_ply_error_cb error_cb, long idata, void* pdata) {
    p_ply ply = NULL;
    FILE* fp = NULL;
    assert(name && storage_mode <= PLY_DEFAULT);
    if (error_cb == NULL) error_cb = ply_error_cb;
    fp = fopen(name, "wb");
    if (!fp) {
      error_cb(ply, "Unable to create file");
      return NULL;
      }
    ply = ply_create_to_file(fp, storage_mode, error_cb, idata, pdata);
    if (ply) ply->own_fp = 1;
    else fclose(fp);
    return ply;
    }

  JTKPLYDEF p_ply ply_create_to_file(FILE* fp, e_ply_storage_mode storage_mode,
    p_ply_error_cb error_cb, long idata, void* pdata) {
    p_ply ply = nullptr;
    assert(fp && storage_mode <= PLY_DEFAULT);
    if (!ply_type_check()) {
      error_cb(ply, "Incompatible type system");
      return NULL;
      }
    ply = ply_alloc();
    if (!ply) {
      error_cb(NULL, "Out of memory");
      return NULL;
      }
    ply->in_memory = false;
    ply->idata = idata;
    ply->pdata = pdata;
    ply->io_mode = PLY_WRITE;
    if (storage_mode == PLY_DEFAULT) storage_mode = ply_arch_endian();
    if (storage_mode == PLY_ASCII) ply->odriver = &ply_odriver_ascii;
    else if (storage_mode == ply_arch_endian())
      ply->odriver = &ply_odriver_binary;
    else ply->odriver = &ply_odriver_binary_reverse;
    ply->storage_mode = storage_mode;
    ply->fp = fp;
    ply->own_fp = 0;
    ply->error_cb = error_cb;
    return ply;
    }

  JTKPLYDEF p_ply ply_create_to_memory(char* buffer, size_t buffer_size, size_t* ply_size, e_ply_storage_mode storage_mode,
    p_ply_error_cb error_cb, long idata, void* pdata) {
    p_ply ply = nullptr;
    if (error_cb == NULL) error_cb = ply_error_cb;
    assert(storage_mode <= PLY_DEFAULT);
    if (!ply_type_check()) {
      error_cb(ply, "Incompatible type system");
      return NULL;
      }
    ply = ply_alloc();
    if (!ply) {
      error_cb(NULL, "Out of memory");
      return NULL;
      }
    ply->in_memory = true;
    ply->store_memory = buffer;
    ply->store_memory_size = ply_size;
    *(ply->store_memory_size) = 0;
    ply->store_buffer_size = buffer_size;
    ply->idata = idata;
    ply->pdata = pdata;
    ply->io_mode = PLY_WRITE;
    if (storage_mode == PLY_DEFAULT) storage_mode = ply_arch_endian();
    if (storage_mode == PLY_ASCII) ply->odriver = &ply_odriver_ascii;
    else if (storage_mode == ply_arch_endian())
      ply->odriver = &ply_odriver_binary;
    else ply->odriver = &ply_odriver_binary_reverse;
    ply->storage_mode = storage_mode;
    ply->fp = NULL;
    ply->own_fp = 0;
    ply->error_cb = error_cb;
    return ply;
    }

  JTKPLYDEF int ply_add_element(p_ply ply, const char* name, long ninstances) {
    p_ply_element element = NULL;
    assert(ply && ply->io_mode == PLY_WRITE);
    if (!ply->in_memory) {
      assert(ply->fp);
      }
    assert(name && strlen(name) < WORDSIZE && ninstances >= 0);
    if (strlen(name) >= WORDSIZE || ninstances < 0) {
      ply_ferror(ply, "Invalid arguments");
      return 0;
      }
    element = ply_grow_element(ply);
    if (!element) return 0;
    strcpy(element->name, name);
    element->ninstances = ninstances;
    return 1;
    }

  JTKPLYDEF int ply_add_scalar_property(p_ply ply, const char* name, e_ply_type type) {
    p_ply_element element = NULL;
    p_ply_property property = NULL;
    assert(ply && ply->io_mode == PLY_WRITE);
    if (!ply->in_memory) {
      assert(ply->fp);
      }
    assert(name && strlen(name) < WORDSIZE);
    assert(type < PLY_LIST);
    if (strlen(name) >= WORDSIZE || type >= PLY_LIST) {
      ply_ferror(ply, "Invalid arguments");
      return 0;
      }
    element = &ply->element[ply->nelements - 1];
    property = ply_grow_property(ply, element);
    if (!property) return 0;
    strcpy(property->name, name);
    property->type = type;
    return 1;
    }

  JTKPLYDEF int ply_add_list_property(p_ply ply, const char* name,
    e_ply_type length_type, e_ply_type value_type) {
    p_ply_element element = NULL;
    p_ply_property property = NULL;
    assert(ply && ply->io_mode == PLY_WRITE);
    if (!ply->in_memory) {
      assert(ply->fp);
      }
    assert(name && strlen(name) < WORDSIZE);
    if (strlen(name) >= WORDSIZE) {
      ply_ferror(ply, "Invalid arguments");
      return 0;
      }
    assert(length_type < PLY_LIST);
    assert(value_type < PLY_LIST);
    if (length_type >= PLY_LIST || value_type >= PLY_LIST) {
      ply_ferror(ply, "Invalid arguments");
      return 0;
      }
    element = &ply->element[ply->nelements - 1];
    property = ply_grow_property(ply, element);
    if (!property) return 0;
    strcpy(property->name, name);
    property->type = PLY_LIST;
    property->length_type = length_type;
    property->value_type = value_type;
    return 1;
    }

  JTKPLYDEF int ply_add_property(p_ply ply, const char* name, e_ply_type type,
    e_ply_type length_type, e_ply_type value_type) {
    if (type == PLY_LIST)
      return ply_add_list_property(ply, name, length_type, value_type);
    else
      return ply_add_scalar_property(ply, name, type);
    }

  JTKPLYDEF int ply_add_comment(p_ply ply, const char* comment) {
    char* new_comment = NULL;
    assert(ply && comment && strlen(comment) < LINESIZE);
    if (!comment || strlen(comment) >= LINESIZE) {
      ply_ferror(ply, "Invalid arguments");
      return 0;
      }
    new_comment = (char*)ply_grow_array(ply, (void**)&ply->comment,
      &ply->ncomments, LINESIZE);
    if (!new_comment) return 0;
    strcpy(new_comment, comment);
    return 1;
    }

  JTKPLYDEF int ply_add_obj_info(p_ply ply, const char* obj_info) {
    char* new_obj_info = NULL;
    assert(ply && obj_info && strlen(obj_info) < LINESIZE);
    if (!obj_info || strlen(obj_info) >= LINESIZE) {
      ply_ferror(ply, "Invalid arguments");
      return 0;
      }
    new_obj_info = (char*)ply_grow_array(ply, (void**)&ply->obj_info,
      &ply->nobj_infos, LINESIZE);
    if (!new_obj_info) return 0;
    strcpy(new_obj_info, obj_info);
    return 1;
    }

  JTKPLYDEF int ply_write_header(p_ply ply) {
    long i, j;
    assert(ply && ply->io_mode == PLY_WRITE);
    if (!ply->in_memory) {
      assert(ply->fp);
      }
    assert(ply->element || ply->nelements == 0);
    assert(!ply->element || ply->nelements > 0);
    char buffer[LINESIZE];
    if (ply->in_memory) {
      sprintf(buffer, "ply\nformat %s 1.0\n",
        ply_storage_mode_list[ply->storage_mode]);
      if (ply_store_buffer_in_memory(ply, buffer) == 0) goto error;
      }
    else {
      if (fprintf(ply->fp, "ply\nformat %s 1.0\n",
        ply_storage_mode_list[ply->storage_mode]) <= 0) goto error;
      }
    for (i = 0; i < ply->ncomments; i++) {
      if (ply->in_memory) {
        sprintf(buffer, "comment %s\n", ply->comment + LINESIZE * i);
        if (ply_store_buffer_in_memory(ply, buffer) == 0) goto error;
        }
      else {
        if (fprintf(ply->fp, "comment %s\n", ply->comment + LINESIZE * i) <= 0)
          goto error;
        }
      }
    for (i = 0; i < ply->nobj_infos; i++) {
      if (ply->in_memory) {
        sprintf(buffer, "obj_info %s\n", ply->obj_info + LINESIZE * i);
        if (ply_store_buffer_in_memory(ply, buffer) == 0) goto error;
        }
      else {
        if (fprintf(ply->fp, "obj_info %s\n", ply->obj_info + LINESIZE * i) <= 0)
          goto error;
        }
      }
    for (i = 0; i < ply->nelements; i++) {
      p_ply_element element = &ply->element[i];
      assert(element->property || element->nproperties == 0);
      assert(!element->property || element->nproperties > 0);
      if (ply->in_memory) {
        sprintf(buffer, "element %s %ld\n", element->name,
          element->ninstances);
        if (ply_store_buffer_in_memory(ply, buffer) == 0) goto error;
        }
      else {
        if (fprintf(ply->fp, "element %s %ld\n", element->name,
          element->ninstances) <= 0) goto error;
        }
      for (j = 0; j < element->nproperties; j++) {
        p_ply_property property = &element->property[j];
        if (property->type == PLY_LIST) {
          if (ply->in_memory) {
            sprintf(buffer, "property list %s %s %s\n",
              ply_type_list[property->length_type],
              ply_type_list[property->value_type],
              property->name);
            if (ply_store_buffer_in_memory(ply, buffer) == 0) goto error;
            }
          else {
            if (fprintf(ply->fp, "property list %s %s %s\n",
              ply_type_list[property->length_type],
              ply_type_list[property->value_type],
              property->name) <= 0) goto error;
            }
          }
        else {
          if (ply->in_memory) {
            sprintf(buffer, "property %s %s\n",
              ply_type_list[property->type],
              property->name);
            if (ply_store_buffer_in_memory(ply, buffer) == 0) goto error;
            }
          else {
            if (fprintf(ply->fp, "property %s %s\n",
              ply_type_list[property->type],
              property->name) <= 0) goto error;
            }
          }
        }
      }
    if (ply->in_memory) {
      return ply_store_buffer_in_memory(ply, "end_header\n") > 0;
      }
    else {
      return fprintf(ply->fp, "end_header\n") > 0;
      }
  error:
    ply_ferror(ply, "Error writing to file");
    return 0;
    }

  JTKPLYDEF int ply_write(p_ply ply, double value) {
    p_ply_element element = NULL;
    p_ply_property property = NULL;
    int type = -1;
    int breakafter = 0;
    int spaceafter = 1;
    if (ply->welement > ply->nelements) return 0;
    element = &ply->element[ply->welement];
    if (ply->wproperty > element->nproperties) return 0;
    property = &element->property[ply->wproperty];
    if (property->type == PLY_LIST) {
      if (ply->wvalue_index == 0) {
        type = property->length_type;
        ply->wlength = (long)value;
        }
      else type = property->value_type;
      }
    else {
      type = property->type;
      ply->wlength = 0;
      }
    if (!ply->odriver->ohandler[type](ply, value)) {
      ply_ferror(ply, "Failed writing %s of %s %d (%s: %s)",
        property->name, element->name,
        ply->winstance_index,
        ply->odriver->name, ply_type_list[type]);
      return 0;
      }
    ply->wvalue_index++;
    if (ply->wvalue_index > ply->wlength) {
      ply->wvalue_index = 0;
      ply->wproperty++;
      }
    if (ply->wproperty >= element->nproperties) {
      ply->wproperty = 0;
      ply->winstance_index++;
      breakafter = 1;
      spaceafter = 0;
      }
    if (ply->winstance_index >= element->ninstances) {
      ply->winstance_index = 0;
      do {
        ply->welement++;
        element = &ply->element[ply->welement];
        } while (ply->welement < ply->nelements && !element->ninstances);
      }
    if (ply->storage_mode == PLY_ASCII) {
      if (ply->in_memory) {
        return (!spaceafter || ply_store_buffer_in_memory(ply, " ") > 0) &&
          (!breakafter || ply_store_buffer_in_memory(ply, "\n") > 0);
        }
      else {
        return (!spaceafter || putc(' ', ply->fp) > 0) &&
          (!breakafter || putc('\n', ply->fp) > 0);
        }
      }
    else {
      return 1;
      }
    }

  JTKPLYDEF int ply_close(p_ply ply) {
    long i;
    assert(ply);
    if (!ply->in_memory) {
      assert(ply->fp);
      }
    assert(ply->element || ply->nelements == 0);
    assert(!ply->element || ply->nelements > 0);
    /* write last chunk to file */
    if (ply->io_mode == PLY_WRITE) {
      if (ply->in_memory) {
        if (ply->buffer_last > 0) {
          memcpy(ply->store_memory, ply->buffer, ply->buffer_last);
          ply->store_memory += ply->buffer_last;
          *(ply->store_memory_size) += ply->buffer_last;
          }
        }
      else {
        if (fwrite(ply->buffer, 1, ply->buffer_last, ply->fp) < ply->buffer_last) {
          ply_ferror(ply, "Error closing up");
          return 0;
          }
        }
      }
    if (!ply->in_memory && ply->own_fp) fclose(ply->fp);
    /* free all memory used by handle */
    if (ply->element) {
      for (i = 0; i < ply->nelements; i++) {
        p_ply_element element = &ply->element[i];
        if (element->property) free(element->property);
        }
      free(ply->element);
      }
    if (ply->obj_info) free(ply->obj_info);
    if (ply->comment) free(ply->comment);
    free(ply);
    return 1;
    }

  /* ----------------------------------------------------------------------
   * Query support functions
   * ---------------------------------------------------------------------- */
  JTKPLYDEF p_ply_element ply_get_next_element(p_ply ply,
    p_ply_element last) {
    assert(ply);
    if (!last) return ply->element;
    last++;
    if (last < ply->element + ply->nelements) return last;
    else return NULL;
    }

  JTKPLYDEF int ply_get_element_info(p_ply_element element, const char** name,
    long* ninstances) {
    assert(element);
    if (name) *name = element->name;
    if (ninstances) *ninstances = (long)element->ninstances;
    return 1;
    }

  JTKPLYDEF p_ply_property ply_get_next_property(p_ply_element element,
    p_ply_property last) {
    assert(element);
    if (!last) return element->property;
    last++;
    if (last < element->property + element->nproperties) return last;
    else return NULL;
    }

  JTKPLYDEF int ply_get_property_info(p_ply_property property, const char** name,
    e_ply_type* type, e_ply_type* length_type, e_ply_type* value_type) {
    assert(property);
    if (name) *name = property->name;
    if (type) *type = property->type;
    if (length_type) *length_type = property->length_type;
    if (value_type) *value_type = property->value_type;
    return 1;

    }

  JTKPLYDEF const char* ply_get_next_comment(p_ply ply, const char* last) {
    assert(ply);
    if (!last) return ply->comment;
    last += LINESIZE;
    if (last < ply->comment + LINESIZE * ply->ncomments) return last;
    else return NULL;
    }

  JTKPLYDEF const char* ply_get_next_obj_info(p_ply ply, const char* last) {
    assert(ply);
    if (!last) return ply->obj_info;
    last += LINESIZE;
    if (last < ply->obj_info + LINESIZE * ply->nobj_infos) return last;
    else return NULL;
    }

  /* ----------------------------------------------------------------------
   * Callback argument support functions
   * ---------------------------------------------------------------------- */
  JTKPLYDEF int ply_get_argument_element(p_ply_argument argument,
    p_ply_element* element, long* instance_index) {
    assert(argument);
    if (!argument) return 0;
    if (element) *element = argument->element;
    if (instance_index) *instance_index = argument->instance_index;
    return 1;
    }

  JTKPLYDEF int ply_get_argument_property(p_ply_argument argument,
    p_ply_property* property, long* length, long* value_index) {
    assert(argument);
    if (!argument) return 0;
    if (property) *property = argument->property;
    if (length) *length = argument->length;
    if (value_index) *value_index = argument->value_index;
    return 1;
    }

  JTKPLYDEF int ply_get_argument_user_data(p_ply_argument argument, void** pdata,
    long* idata) {
    assert(argument);
    if (!argument) return 0;
    if (pdata) *pdata = argument->pdata;
    if (idata) *idata = argument->idata;
    return 1;
    }

  JTKPLYDEF double ply_get_argument_value(p_ply_argument argument) {
    assert(argument);
    if (!argument) return 0.0;
    return argument->value;
    }

  JTKPLYDEF int ply_get_ply_user_data(p_ply ply, void** pdata, long* idata) {
    assert(ply);
    if (!ply) return 0;
    if (pdata) *pdata = ply->pdata;
    if (idata) *idata = ply->idata;
    return 1;
    }

  /* ----------------------------------------------------------------------
   * Internal functions
   * ---------------------------------------------------------------------- */
  static int ply_read_list_property(p_ply ply, p_ply_element element,
    p_ply_property property, p_ply_argument argument) {
    int l;
    p_ply_read_cb read_cb = property->read_cb;
    p_ply_ihandler* driver = ply->idriver->ihandler;
    /* get list length */
    p_ply_ihandler handler = driver[property->length_type];
    double length;
    if (!handler(ply, &length)) {
      ply_ferror(ply, "Error reading '%s' of '%s' number %d",
        property->name, element->name, argument->instance_index);
      return 0;
      }
    /* invoke callback to pass length in value field */
    argument->length = (long)length;
    argument->value_index = -1;
    argument->value = length;
    if (read_cb && !read_cb(argument)) {
      ply_ferror(ply, "Aborted by user");
      return 0;
      }
    /* read list values */
    handler = driver[property->value_type];
    /* for each value in list */
    for (l = 0; l < (long)length; l++) {
      /* read value from file */
      argument->value_index = l;
      if (!handler(ply, &argument->value)) {
        ply_ferror(ply, "Error reading value number %d of '%s' of "
          "'%s' number %d", l + 1, property->name,
          element->name, argument->instance_index);
        return 0;
        }
      /* invoke callback to pass value */
      if (read_cb && !read_cb(argument)) {
        ply_ferror(ply, "Aborted by user");
        return 0;
        }
      }
    return 1;
    }

  static int ply_read_scalar_property(p_ply ply, p_ply_element element,
    p_ply_property property, p_ply_argument argument) {
    p_ply_read_cb read_cb = property->read_cb;
    p_ply_ihandler* driver = ply->idriver->ihandler;
    p_ply_ihandler handler = driver[property->type];
    argument->length = 1;
    argument->value_index = 0;
    if (!handler(ply, &argument->value)) {
      ply_ferror(ply, "Error reading '%s' of '%s' number %d",
        property->name, element->name, argument->instance_index);
      return 0;
      }
    if (read_cb && !read_cb(argument)) {
      ply_ferror(ply, "Aborted by user");
      return 0;
      }
    return 1;
    }

  static int ply_read_property(p_ply ply, p_ply_element element,
    p_ply_property property, p_ply_argument argument) {
    if (property->type == PLY_LIST)
      return ply_read_list_property(ply, element, property, argument);
    else
      return ply_read_scalar_property(ply, element, property, argument);
    }

  static int ply_read_element(p_ply ply, p_ply_element element,
    p_ply_argument argument) {
    long j, k;
    /* for each element of this type */
    for (j = 0; j < element->ninstances; j++) {
      argument->instance_index = j;
      /* for each property */
      for (k = 0; k < element->nproperties; k++) {
        p_ply_property property = &element->property[k];
        argument->property = property;
        argument->pdata = property->pdata;
        argument->idata = property->idata;
        if (!ply_read_property(ply, element, property, argument))
          return 0;
        }
      }
    return 1;
    }

  static int ply_find_string(const char* item, const char* const list[]) {
    int i;
    assert(item && list);
    for (i = 0; list[i]; i++)
      if (!strcmp(list[i], item)) return i;
    return -1;
    }

  static p_ply_element ply_find_element(p_ply ply, const char* name) {
    p_ply_element element;
    long i, nelements;
    assert(ply && name);
    element = ply->element;
    nelements = ply->nelements;
    assert(element || nelements == 0);
    assert(!element || nelements > 0);
    for (i = 0; i < nelements; i++)
      if (!strcmp(element[i].name, name)) return &element[i];
    return NULL;
    }

  static p_ply_property ply_find_property(p_ply_element element,
    const char* name) {
    p_ply_property property;
    long i, nproperties;
    assert(element && name);
    property = element->property;
    nproperties = element->nproperties;
    assert(property || nproperties == 0);
    assert(!property || nproperties > 0);
    for (i = 0; i < nproperties; i++)
      if (!strcmp(property[i].name, name)) return &property[i];
    return NULL;
    }

  static int ply_check_word(p_ply ply) {
    size_t size = strlen(BWORD(ply));
    if (size >= WORDSIZE) {
      ply_ferror(ply, "Word too long");
      return 0;
      }
    else if (size == 0) {
      ply_ferror(ply, "Unexpected end of file");
      return 0;
      }
    return 1;
    }

  static int ply_read_word(p_ply ply) {
    size_t t = 0;
    assert(ply && ply->io_mode == PLY_READ);
    if (!ply->in_memory) {
      assert(ply->fp);
      }
    /* skip leading blanks */
    while (1) {
      t = strspn(BFIRST(ply), " \n\r\t");
      /* check if all buffer was made of blanks */
      if (t >= BSIZE(ply)) {
        if (!BREFILL(ply)) {
          ply_ferror(ply, "Unexpected end of file");
          return 0;
          }
        }
      else break;
      }
    BSKIP(ply, t);
    /* look for a space after the current word */
    t = strcspn(BFIRST(ply), " \n\r\t");
    /* if we didn't reach the end of the buffer, we are done */
    if (t < BSIZE(ply)) {
      ply_finish_word(ply, t);
      return ply_check_word(ply);
      }
    /* otherwise, try to refill buffer */
    if (!BREFILL(ply)) {
      /* if we reached the end of file, try to do with what we have */
      ply_finish_word(ply, t);
      return ply_check_word(ply);
      /* ply_ferror(ply, "Unexpected end of file"); */
      /* return 0; */
      }
    /* keep looking from where we left */
    t += strcspn(BFIRST(ply) + t, " \n\r\t");
    /* check if the token is too large for our buffer */
    if (t >= BSIZE(ply)) {
      ply_ferror(ply, "Token too large");
      return 0;
      }
    /* we are done */
    ply_finish_word(ply, t);
    return ply_check_word(ply);
    }

  static void ply_finish_word(p_ply ply, size_t size) {
    ply->buffer_token = ply->buffer_first;
    BSKIP(ply, size);
    *BFIRST(ply) = '\0';
    BSKIP(ply, 1);
    }

  static int ply_check_line(p_ply ply) {
    if (strlen(BLINE(ply)) >= LINESIZE) {
      ply_ferror(ply, "Line too long");
      return 0;
      }
    return 1;
    }

  static int ply_read_line(p_ply ply) {
    const char* end = NULL;
    assert(ply && ply->io_mode == PLY_READ);
    if (!ply->in_memory) {
      assert(ply->fp);
      }
    /* look for a end of line */
    end = strchr(BFIRST(ply), '\n');
    /* if we didn't reach the end of the buffer, we are done */
    if (end) {
      ply->buffer_token = ply->buffer_first;
      BSKIP(ply, end - BFIRST(ply));
      *BFIRST(ply) = '\0';
      BSKIP(ply, 1);
      return ply_check_line(ply);
      }
    else {
      end = ply->buffer + BSIZE(ply);
      /* otherwise, try to refill buffer */
      if (!BREFILL(ply)) {
        ply_ferror(ply, "Unexpected end of file");
        return 0;
        }
      }
    /* keep looking from where we left */
    end = strchr(end, '\n');
    /* check if the token is too large for our buffer */
    if (!end) {
      ply_ferror(ply, "Token too large");
      return 0;
      }
    /* we are done */
    ply->buffer_token = ply->buffer_first;
    BSKIP(ply, end - BFIRST(ply));
    *BFIRST(ply) = '\0';
    BSKIP(ply, 1);
    return ply_check_line(ply);
    }

  static int ply_read_chunk(p_ply ply, void* anybuffer, size_t size) {
    char* buffer = (char*)anybuffer;
    size_t i = 0;
    assert(ply && ply->io_mode == PLY_READ);
    if (!ply->in_memory) {
      assert(ply->fp);
      }
    assert(ply->buffer_first <= ply->buffer_last);
    while (i < size) {
      if (ply->buffer_first < ply->buffer_last) {
        buffer[i] = ply->buffer[ply->buffer_first];
        ply->buffer_first++;
        i++;
        }
      else {
        ply->buffer_first = 0;
        if (ply->in_memory) { // Memory reading fixed by Jan Maes
          size_t chunksize = ply->load_memory_size > BUFFERSIZE ? BUFFERSIZE : ply->load_memory_size;
          memcpy(ply->buffer, ply->load_memory, chunksize);
          ply->load_memory += chunksize;
          ply->load_memory_size -= chunksize;
          ply->buffer_last = chunksize;
          }
        else {
          ply->buffer_last = fread(ply->buffer, 1, BUFFERSIZE, ply->fp);
          }
        if (ply->buffer_last <= 0) return 0;
        }
      }
    return 1;
    }

  static int ply_write_chunk(p_ply ply, void* anybuffer, size_t size) {
    char* buffer = (char*)anybuffer;
    size_t i = 0;
    assert(ply && ply->io_mode == PLY_WRITE);
    if (!ply->in_memory) {
      assert(ply->fp);
      }
    assert(ply->buffer_last <= BUFFERSIZE);
    while (i < size) {
      if (ply->buffer_last < BUFFERSIZE) {
        ply->buffer[ply->buffer_last] = buffer[i];
        ply->buffer_last++;
        i++;
        }
      else {
        ply->buffer_last = 0;
        if (ply->in_memory) {
          memcpy(ply->store_memory, ply->buffer, BUFFERSIZE);
          ply->store_memory += BUFFERSIZE;
          *(ply->store_memory_size) += BUFFERSIZE;
          }
        else {
          if (fwrite(ply->buffer, 1, BUFFERSIZE, ply->fp) < BUFFERSIZE)
            return 0;
          }
        }
      }
    return 1;
    }

  static int ply_write_chunk_reverse(p_ply ply, void* anybuffer, size_t size) {
    int ret = 0;
    ply_reverse(anybuffer, size);
    ret = ply_write_chunk(ply, anybuffer, size);
    ply_reverse(anybuffer, size);
    return ret;
    }

  static int ply_read_chunk_reverse(p_ply ply, void* anybuffer, size_t size) {
    if (!ply_read_chunk(ply, anybuffer, size)) return 0;
    ply_reverse(anybuffer, size);
    return 1;
    }

  static void ply_reverse(void* anydata, size_t size) {
    char* data = (char*)anydata;
    char temp;
    size_t i;
    for (i = 0; i < size / 2; i++) {
      temp = data[i];
      data[i] = data[size - i - 1];
      data[size - i - 1] = temp;
      }
    }

  static int ply_store_buffer_in_memory(p_ply ply, const char* buffer) {
    assert(ply && ply->io_mode == PLY_WRITE && ply->store_memory);
    size_t length = strlen(buffer);
    if (*ply->store_memory_size + length > ply->store_buffer_size) {
      ply_ferror(ply, "The store buffer is too small. Use a bigger one!");
      return 0;
      }
    memcpy(ply->store_memory, buffer, length);
    ply->store_memory += length;
    *(ply->store_memory_size) += length;
    return length > 0;
    }

  static void ply_init(p_ply ply) {
    ply->element = NULL;
    ply->nelements = 0;
    ply->comment = NULL;
    ply->ncomments = 0;
    ply->obj_info = NULL;
    ply->nobj_infos = 0;
    ply->idriver = NULL;
    ply->odriver = NULL;
    ply->buffer[0] = '\0';
    ply->buffer_first = ply->buffer_last = ply->buffer_token = 0;
    ply->welement = 0;
    ply->wproperty = 0;
    ply->winstance_index = 0;
    ply->wlength = 0;
    ply->wvalue_index = 0;
    }

  static void ply_element_init(p_ply_element element) {
    element->name[0] = '\0';
    element->ninstances = 0;
    element->property = NULL;
    element->nproperties = 0;
    }

  static void ply_property_init(p_ply_property property) {
    property->name[0] = '\0';
    property->type = (e_ply_type)-1;
    property->length_type = (e_ply_type)-1;
    property->value_type = (e_ply_type)-1;
    property->read_cb = (p_ply_read_cb)NULL;
    property->pdata = NULL;
    property->idata = 0;
    }

  static p_ply ply_alloc(void) {
    p_ply ply = (p_ply)calloc(1, sizeof(t_ply));
    if (!ply) return NULL;
    ply_init(ply);
    return ply;
    }

  static void* ply_grow_array(p_ply ply, void** pointer,
    long* nmemb, long size) {
    void* temp = *pointer;
    long count = *nmemb + 1;
    if (!temp) temp = malloc(count * size);
    else temp = realloc(temp, count * size);
    if (!temp) {
      ply_ferror(ply, "Out of memory");
      return NULL;
      }
    *pointer = temp;
    *nmemb = count;
    return (char*)temp + (count - 1) * size;
    }

  static p_ply_element ply_grow_element(p_ply ply) {
    p_ply_element element = NULL;
    assert(ply);
    assert(ply->element || ply->nelements == 0);
    assert(!ply->element || ply->nelements > 0);
    element = (p_ply_element)ply_grow_array(ply, (void**)&ply->element,
      &ply->nelements, sizeof(t_ply_element));
    if (!element) return NULL;
    ply_element_init(element);
    return element;
    }

  static p_ply_property ply_grow_property(p_ply ply, p_ply_element element) {
    p_ply_property property = NULL;
    assert(ply);
    assert(element);
    assert(element->property || element->nproperties == 0);
    assert(!element->property || element->nproperties > 0);
    property = (p_ply_property)ply_grow_array(ply,
      (void**)&element->property,
      &element->nproperties, sizeof(t_ply_property));
    if (!property) return NULL;
    ply_property_init(property);
    return property;
    }

  static int ply_read_header_format(p_ply ply) {
    assert(ply && ply->io_mode == PLY_READ);
    if (!ply->in_memory) {
      assert(ply->fp);
      }
    if (strcmp(BWORD(ply), "format")) return 0;
    if (!ply_read_word(ply)) return 0;
    ply->storage_mode = (e_ply_storage_mode)ply_find_string(BWORD(ply), ply_storage_mode_list);
    if (ply->storage_mode == (e_ply_storage_mode)(-1)) return 0;
    if (ply->storage_mode == PLY_ASCII) ply->idriver = &ply_idriver_ascii;
    else if (ply->storage_mode == ply_arch_endian())
      ply->idriver = &ply_idriver_binary;
    else ply->idriver = &ply_idriver_binary_reverse;
    if (!ply_read_word(ply)) return 0;
    if (strcmp(BWORD(ply), "1.0")) return 0;
    if (!ply_read_word(ply)) return 0;
    return 1;
    }

  static int ply_read_header_comment(p_ply ply) {
    assert(ply && ply->io_mode == PLY_READ);
    if (!ply->in_memory) {
      assert(ply->fp);
      }
    if (strcmp(BWORD(ply), "comment")) return 0;
    if (!ply_read_line(ply)) return 0;
    if (!ply_add_comment(ply, BLINE(ply))) return 0;
    if (!ply_read_word(ply)) return 0;
    return 1;
    }

  static int ply_read_header_obj_info(p_ply ply) {
    assert(ply && ply->io_mode == PLY_READ);
    if (!ply->in_memory) {
      assert(ply->fp);
      }
    if (strcmp(BWORD(ply), "obj_info")) return 0;
    if (!ply_read_line(ply)) return 0;
    if (!ply_add_obj_info(ply, BLINE(ply))) return 0;
    if (!ply_read_word(ply)) return 0;
    return 1;
    }

  static int ply_read_header_property(p_ply ply) {
    p_ply_element element = NULL;
    p_ply_property property = NULL;
    /* make sure it is a property */
    if (strcmp(BWORD(ply), "property")) return 0;
    element = &ply->element[ply->nelements - 1];
    property = ply_grow_property(ply, element);
    if (!property) return 0;
    /* get property type */
    if (!ply_read_word(ply)) return 0;
    property->type = (e_ply_type)ply_find_string(BWORD(ply), ply_type_list);
    if (property->type == (e_ply_type)(-1)) return 0;
    if (property->type == PLY_LIST) {
      /* if it's a list, we need the base types */
      if (!ply_read_word(ply)) return 0;
      property->length_type = (e_ply_type)ply_find_string(BWORD(ply), ply_type_list);
      if (property->length_type == (e_ply_type)(-1)) return 0;
      if (!ply_read_word(ply)) return 0;
      property->value_type = (e_ply_type)ply_find_string(BWORD(ply), ply_type_list);
      if (property->value_type == (e_ply_type)(-1)) return 0;
      }
    /* get property name */
    if (!ply_read_word(ply)) return 0;
    strcpy(property->name, BWORD(ply));
    if (!ply_read_word(ply)) return 0;
    return 1;
    }

  static int ply_read_header_element(p_ply ply) {
    p_ply_element element = NULL;
    long dummy;
    assert(ply && ply->io_mode == PLY_READ);
    if (!ply->in_memory) {
      assert(ply->fp);
      }
    if (strcmp(BWORD(ply), "element")) return 0;
    /* allocate room for new element */
    element = ply_grow_element(ply);
    if (!element) return 0;
    /* get element name */
    if (!ply_read_word(ply)) return 0;
    strcpy(element->name, BWORD(ply));
    /* get number of elements of this type */
    if (!ply_read_word(ply)) return 0;
    if (sscanf(BWORD(ply), "%ld", &dummy) != 1) {
      ply_ferror(ply, "Expected number got '%s'", BWORD(ply));
      return 0;
      }
    element->ninstances = dummy;
    /* get all properties for this element */
    if (!ply_read_word(ply)) return 0;
    while (ply_read_header_property(ply) ||
      ply_read_header_comment(ply) || ply_read_header_obj_info(ply))
      /* do nothing */;
    return 1;
    }

  static void ply_error_cb(p_ply ply, const char* message) {
    (void)ply;
    fprintf(stderr, "RPly: %s\n", message);
    }

  static void ply_ferror(p_ply ply, const char* fmt, ...) {
    char buffer[1024];
    va_list ap;
    va_start(ap, fmt);
    vsprintf(buffer, fmt, ap);
    va_end(ap);
    ply->error_cb(ply, buffer);
    }

  static e_ply_storage_mode ply_arch_endian(void) {
    unsigned long i = 1;
    unsigned char* s = (unsigned char*)&i;
    if (*s == 1) return PLY_LITTLE_ENDIAN;
    else return PLY_BIG_ENDIAN;
    }

  static int ply_type_check(void) {
    assert(sizeof(t_ply_int8) == 1);
    assert(sizeof(t_ply_uint8) == 1);
    assert(sizeof(t_ply_int16) == 2);
    assert(sizeof(t_ply_uint16) == 2);
    assert(sizeof(t_ply_int32) == 4);
    assert(sizeof(t_ply_uint32) == 4);
    assert(sizeof(float) == 4);
    assert(sizeof(double) == 8);
    if constexpr (sizeof(t_ply_int8) != 1) return 0;
    if constexpr (sizeof(t_ply_uint8) != 1) return 0;
    if constexpr (sizeof(t_ply_int16) != 2) return 0;
    if constexpr (sizeof(t_ply_uint16) != 2) return 0;
    if constexpr (sizeof(t_ply_int32) != 4) return 0;
    if constexpr (sizeof(t_ply_uint32) != 4) return 0;
    if constexpr (sizeof(float) != 4) return 0;
    if constexpr (sizeof(double) != 8) return 0;
    return 1;
    }

  /* ----------------------------------------------------------------------
   * Output handlers
   * ---------------------------------------------------------------------- */
  static int oascii_int8(p_ply ply, double value) {
    if (value > PLY_INT8_MAX || value < PLY_INT8_MIN) return 0;
    if (ply->in_memory) {
      char buffer[WORDSIZE];
      snprintf(buffer, WORDSIZE, "%d", (t_ply_int8)value);
      return ply_store_buffer_in_memory(ply, buffer);
      }
    else {
      return fprintf(ply->fp, "%d", (t_ply_int8)value) > 0;
      }
    }

  static int oascii_uint8(p_ply ply, double value) {
    if (value > PLY_UINT8_MAX || value < 0) return 0;
    if (ply->in_memory) {
      char buffer[WORDSIZE];
      snprintf(buffer, WORDSIZE, "%d", (t_ply_uint8)value);
      return ply_store_buffer_in_memory(ply, buffer);
      }
    else {
      return fprintf(ply->fp, "%d", (t_ply_uint8)value) > 0;
      }
    }

  static int oascii_int16(p_ply ply, double value) {
    if (value > PLY_INT16_MAX || value < PLY_INT16_MIN) return 0;
    if (ply->in_memory) {
      char buffer[WORDSIZE];
      snprintf(buffer, WORDSIZE, "%d", (t_ply_int16)value);
      return ply_store_buffer_in_memory(ply, buffer);
      }
    else {
      return fprintf(ply->fp, "%d", (t_ply_int16)value) > 0;
      }
    }

  static int oascii_uint16(p_ply ply, double value) {
    if (value > PLY_UINT16_MAX || value < 0) return 0;
    if (ply->in_memory) {
      char buffer[WORDSIZE];
      snprintf(buffer, WORDSIZE, "%d", (t_ply_uint16)value);
      return ply_store_buffer_in_memory(ply, buffer);
      }
    else {
      return fprintf(ply->fp, "%d", (t_ply_uint16)value) > 0;
      }
    }

  static int oascii_int32(p_ply ply, double value) {
    if (value > PLY_INT32_MAX || value < PLY_INT32_MIN) return 0;
    if (ply->in_memory) {
      char buffer[WORDSIZE];
      snprintf(buffer, WORDSIZE, "%d", (t_ply_int32)value);
      return ply_store_buffer_in_memory(ply, buffer);
      }
    else {
      return fprintf(ply->fp, "%d", (t_ply_int32)value) > 0;
      }
    }

  static int oascii_uint32(p_ply ply, double value) {
    if (value > PLY_UINT32_MAX || value < 0) return 0;
    if (ply->in_memory) {
      char buffer[WORDSIZE];
      snprintf(buffer, WORDSIZE, "%d", (t_ply_uint32)value);
      return ply_store_buffer_in_memory(ply, buffer);
      }
    else {
      return fprintf(ply->fp, "%d", (t_ply_uint32)value) > 0;
      }
    }

  static int oascii_float32(p_ply ply, double value) {
    if (value < -FLT_MAX || value > FLT_MAX) return 0;
    if (ply->in_memory) {
      char buffer[WORDSIZE];
      snprintf(buffer, WORDSIZE, "%g", (float)value);
      return ply_store_buffer_in_memory(ply, buffer);
      }
    else {
      return fprintf(ply->fp, "%g", (float)value) > 0;
      }
    }

  static int oascii_float64(p_ply ply, double value) {
    if (value < -DBL_MAX || value > DBL_MAX) return 0;
    if (ply->in_memory) {
      char buffer[WORDSIZE];
      snprintf(buffer, WORDSIZE, "%g", value);
      return ply_store_buffer_in_memory(ply, buffer);
      }
    else {
      return fprintf(ply->fp, "%g", value) > 0;
      }
    }

  static int obinary_int8(p_ply ply, double value) {
    t_ply_int8 int8 = (t_ply_int8)value;
    if (value > PLY_INT8_MAX || value < PLY_INT8_MIN) return 0;
    return ply->odriver->ochunk(ply, &int8, sizeof(int8));
    }

  static int obinary_uint8(p_ply ply, double value) {
    t_ply_uint8 uint8 = (t_ply_uint8)value;
    if (value > PLY_UINT8_MAX || value < 0) return 0;
    return ply->odriver->ochunk(ply, &uint8, sizeof(uint8));
    }

  static int obinary_int16(p_ply ply, double value) {
    t_ply_int16 int16 = (t_ply_int16)value;
    if (value > PLY_INT16_MAX || value < PLY_INT16_MIN) return 0;
    return ply->odriver->ochunk(ply, &int16, sizeof(int16));
    }

  static int obinary_uint16(p_ply ply, double value) {
    t_ply_uint16 uint16 = (t_ply_uint16)value;
    if (value > PLY_UINT16_MAX || value < 0) return 0;
    return ply->odriver->ochunk(ply, &uint16, sizeof(uint16));
    }

  static int obinary_int32(p_ply ply, double value) {
    t_ply_int32 int32 = (t_ply_int32)value;
    if (value > PLY_INT32_MAX || value < PLY_INT32_MIN) return 0;
    return ply->odriver->ochunk(ply, &int32, sizeof(int32));
    }

  static int obinary_uint32(p_ply ply, double value) {
    t_ply_uint32 uint32 = (t_ply_uint32)value;
    if (value > PLY_UINT32_MAX || value < 0) return 0;
    return ply->odriver->ochunk(ply, &uint32, sizeof(uint32));
    }

  static int obinary_float32(p_ply ply, double value) {
    float float32 = (float)value;
    if (value > FLT_MAX || value < -FLT_MAX) return 0;
    return ply->odriver->ochunk(ply, &float32, sizeof(float32));
    }

  static int obinary_float64(p_ply ply, double value) {
    return ply->odriver->ochunk(ply, &value, sizeof(value));
    }

  /* ----------------------------------------------------------------------
   * Input  handlers
   * ---------------------------------------------------------------------- */
  static int iascii_int8(p_ply ply, double* value) {
    char* end;
    if (!ply_read_word(ply)) return 0;
    *value = strtol(BWORD(ply), &end, 10);
    if (*end || *value > PLY_INT8_MAX || *value < PLY_INT8_MIN) return 0;
    return 1;
    }

  static int iascii_uint8(p_ply ply, double* value) {
    char* end;
    if (!ply_read_word(ply)) return 0;
    *value = strtol(BWORD(ply), &end, 10);
    if (*end || *value > PLY_UINT8_MAX || *value < 0) return 0;
    return 1;
    }

  static int iascii_int16(p_ply ply, double* value) {
    char* end;
    if (!ply_read_word(ply)) return 0;
    *value = strtol(BWORD(ply), &end, 10);
    if (*end || *value > PLY_INT16_MAX || *value < PLY_INT16_MIN) return 0;
    return 1;
    }

  static int iascii_uint16(p_ply ply, double* value) {
    char* end;
    if (!ply_read_word(ply)) return 0;
    *value = strtol(BWORD(ply), &end, 10);
    if (*end || *value > PLY_UINT16_MAX || *value < 0) return 0;
    return 1;
    }

  static int iascii_int32(p_ply ply, double* value) {
    char* end;
    if (!ply_read_word(ply)) return 0;
    *value = strtol(BWORD(ply), &end, 10);
    if (*end || *value > PLY_INT32_MAX || *value < PLY_INT32_MIN) return 0;
    return 1;
    }

  static int iascii_uint32(p_ply ply, double* value) {
    char* end;
    if (!ply_read_word(ply)) return 0;
    *value = strtol(BWORD(ply), &end, 10);
    if (*end || *value > PLY_UINT32_MAX || *value < 0) return 0;
    return 1;
    }

  static int iascii_float32(p_ply ply, double* value) {
    char* end;
    if (!ply_read_word(ply)) return 0;
    *value = strtod(BWORD(ply), &end);
    if (*end || *value < -FLT_MAX || *value > FLT_MAX) return 0;
    return 1;
    }

  static int iascii_float64(p_ply ply, double* value) {
    char* end;
    if (!ply_read_word(ply)) return 0;
    *value = strtod(BWORD(ply), &end);
    if (*end || *value < -DBL_MAX || *value > DBL_MAX) return 0;
    return 1;
    }

  static int ibinary_int8(p_ply ply, double* value) {
    t_ply_int8 int8;
    if (!ply->idriver->ichunk(ply, &int8, 1)) return 0;
    *value = int8;
    return 1;
    }

  static int ibinary_uint8(p_ply ply, double* value) {
    t_ply_uint8 uint8;
    if (!ply->idriver->ichunk(ply, &uint8, 1)) return 0;
    *value = uint8;
    return 1;
    }

  static int ibinary_int16(p_ply ply, double* value) {
    t_ply_int16 int16;
    if (!ply->idriver->ichunk(ply, &int16, sizeof(int16))) return 0;
    *value = int16;
    return 1;
    }

  static int ibinary_uint16(p_ply ply, double* value) {
    t_ply_uint16 uint16;
    if (!ply->idriver->ichunk(ply, &uint16, sizeof(uint16))) return 0;
    *value = uint16;
    return 1;
    }

  static int ibinary_int32(p_ply ply, double* value) {
    t_ply_int32 int32;
    if (!ply->idriver->ichunk(ply, &int32, sizeof(int32))) return 0;
    *value = int32;
    return 1;
    }

  static int ibinary_uint32(p_ply ply, double* value) {
    t_ply_uint32 uint32;
    if (!ply->idriver->ichunk(ply, &uint32, sizeof(uint32))) return 0;
    *value = uint32;
    return 1;
    }

  static int ibinary_float32(p_ply ply, double* value) {
    float float32;
    if (!ply->idriver->ichunk(ply, &float32, sizeof(float32))) return 0;
    *value = float32;
    return 1;
    }

  static int ibinary_float64(p_ply ply, double* value) {
    return ply->idriver->ichunk(ply, value, sizeof(double));
    }

  /* ----------------------------------------------------------------------
   * Constants
   * ---------------------------------------------------------------------- */

  namespace
    {
    t_ply_idriver ply_idriver_ascii = {
        {   iascii_int8, iascii_uint8, iascii_int16, iascii_uint16,
            iascii_int32, iascii_uint32, iascii_float32, iascii_float64,
            iascii_int8, iascii_uint8, iascii_int16, iascii_uint16,
            iascii_int32, iascii_uint32, iascii_float32, iascii_float64
        }, /* order matches e_ply_type enum */
        NULL,
        "ascii input"
      };

    t_ply_idriver ply_idriver_binary = {
        {   ibinary_int8, ibinary_uint8, ibinary_int16, ibinary_uint16,
            ibinary_int32, ibinary_uint32, ibinary_float32, ibinary_float64,
            ibinary_int8, ibinary_uint8, ibinary_int16, ibinary_uint16,
            ibinary_int32, ibinary_uint32, ibinary_float32, ibinary_float64
        }, /* order matches e_ply_type enum */
        ply_read_chunk,
        "binary input"
      };

    t_ply_idriver ply_idriver_binary_reverse = {
        {   ibinary_int8, ibinary_uint8, ibinary_int16, ibinary_uint16,
            ibinary_int32, ibinary_uint32, ibinary_float32, ibinary_float64,
            ibinary_int8, ibinary_uint8, ibinary_int16, ibinary_uint16,
            ibinary_int32, ibinary_uint32, ibinary_float32, ibinary_float64
        }, /* order matches e_ply_type enum */
        ply_read_chunk_reverse,
        "reverse binary input"
      };

    t_ply_odriver ply_odriver_ascii = {
        {   oascii_int8, oascii_uint8, oascii_int16, oascii_uint16,
            oascii_int32, oascii_uint32, oascii_float32, oascii_float64,
            oascii_int8, oascii_uint8, oascii_int16, oascii_uint16,
            oascii_int32, oascii_uint32, oascii_float32, oascii_float64
        }, /* order matches e_ply_type enum */
        NULL,
        "ascii output"
      };

    t_ply_odriver ply_odriver_binary = {
        {   obinary_int8, obinary_uint8, obinary_int16, obinary_uint16,
            obinary_int32, obinary_uint32, obinary_float32, obinary_float64,
            obinary_int8, obinary_uint8, obinary_int16, obinary_uint16,
            obinary_int32, obinary_uint32, obinary_float32, obinary_float64
        }, /* order matches e_ply_type enum */
        ply_write_chunk,
        "binary output"
      };

    t_ply_odriver ply_odriver_binary_reverse = {
        {   obinary_int8, obinary_uint8, obinary_int16, obinary_uint16,
            obinary_int32, obinary_uint32, obinary_float32, obinary_float64,
            obinary_int8, obinary_uint8, obinary_int16, obinary_uint16,
            obinary_int32, obinary_uint32, obinary_float32, obinary_float64
        }, /* order matches e_ply_type enum */
        ply_write_chunk_reverse,
        "reverse binary output"
      };

    }

  JTKPLYDEF bool read_ply(const char* filename, std::vector<vec3<float>>& vertices, std::vector<vec3<float>>& normals, std::vector<uint32_t>& clrs, std::vector<vec3<uint32_t>>& triangles, std::vector<vec3<vec2<float>>>& uv)
    {
    return _read_ply<char>(filename, vertices, normals, clrs, triangles, uv);
    }

  JTKPLYDEF bool read_ply(const wchar_t* filename, std::vector<vec3<float>>& vertices, std::vector<vec3<float>>& normals, std::vector<uint32_t>& clrs, std::vector<vec3<uint32_t>>& triangles, std::vector<vec3<vec2<float>>>& uv)
    {
    return _read_ply<wchar_t>(filename, vertices, normals, clrs, triangles, uv);
    }

  JTKPLYDEF bool write_ply(const char* filename, const std::vector<vec3<float>>& vertices, const std::vector<vec3<float>>& normals, const std::vector<uint32_t>& clrs, const std::vector<vec3<uint32_t>>& triangles, const std::vector<vec3<vec2<float>>>& uv)
    {
    return _write_ply<char>(filename, vertices, normals, clrs, triangles, uv);
    }

  JTKPLYDEF bool write_ply(const wchar_t* filename, const std::vector<vec3<float>>& vertices, const std::vector<vec3<float>>& normals, const std::vector<uint32_t>& clrs, const std::vector<vec3<uint32_t>>& triangles, const std::vector<vec3<vec2<float>>>& uv)
    {
    return _write_ply<wchar_t>(filename, vertices, normals, clrs, triangles, uv);
    }

  } // namespace jtk


#endif //JTK_PLY_IMPLEMENTATION

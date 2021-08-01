/*
   Do this:
      #define JTK_OPENGL_IMPLEMENTATION
   before you include this file in *one* C++ file to create the implementation.
   // i.e. it should look like this:
   #include ...
   #include ...
   #include ...
   #define JTK_OPENGL_IMPLEMENTATION
   #include "jtk/opengl.h"
 */

#ifndef JTK_OPENGL_H
#define JTK_OPENGL_H

#include <cassert>
#include <sstream>
#include <stdexcept>
#include <iostream>

#ifndef JTKGLDEF
#ifdef JTK_OPENGL_STATIC
#define JTKGLDEF static
#else
#define JTKGLDEF extern
#endif
#endif

namespace jtk
  {
  void gl_check_error(const char* txt);

  class buffer_object
    {
    public:
      buffer_object();
      buffer_object(GLenum type);
      ~buffer_object();

      void allocate(const void* data, int count);
      void update_data(const void* data, int count);

      void create();
      void bind();

      bool is_created() const { return _buffer_object_id != 0; }

      void release();
      void destroy();

      void set_usage_pattern(GLenum pattern) { _pattern = pattern; }

      GLuint buffer_id() const { return _buffer_object_id; }
      GLenum type() const { return _type; }
      GLenum usage_pattern() const { return _pattern; }

      void get_data(void* data, int count);

    private:
      GLenum _type;
      GLenum _pattern;
      GLuint _buffer_object_id;
    };

  class texture
    {
    public:

      enum pixel_type
        {
        int8,
        uint8,
        uint32,
        real
        };

      texture();

      ~texture();

      void set_wrap_mode(GLenum wrapMode);
      void set_filter_mode(GLenum filterMode);

      void load_from_pixels(GLubyte* pixels, GLint w, GLint h, GLint channels, pixel_type pt);

      void create_empty(GLint w, GLint h, GLint channels, pixel_type pt);

      void bind_to_channel(GLint channel);

      void bind_to_image(GLint channel, GLenum access);

      void release();

      GLint width() const { return _width; }
      GLint height() const { return _height; }

      bool is_loaded() const { return _loaded; }

      GLuint texture_id() const { return _texture_id; }

      void fill_pixels(GLubyte* pixels, GLint channels);

    private:
      void _load_pixels(GLubyte* pixels, GLint w, GLint h, GLint channels, pixel_type pt);

    private:
      GLenum _filter_mode;
      GLenum _wrap_mode;
      GLenum _format;
      pixel_type _pt;

      GLuint _texture_id;
      GLint _width, _height;

      bool _loaded;
      bool _isfloat;
    };

  class texture3d
    {
    public:

      enum pixel_type
        {
        int8,
        uint8,
        uint32,
        real
        };

      texture3d();

      ~texture3d();

      void set_wrap_mode(GLenum wrapMode);
      void set_filter_mode(GLenum filterMode);

      void load_from_pixels(GLubyte* pixels, GLint w, GLint h, GLint d, GLint channels, pixel_type pt);

      void create_empty(GLint w, GLint h, GLint d, GLint channels, pixel_type pt);

      void bind_to_channel(GLint channel);

      void bind_to_image(GLint channel, GLenum access);

      void release();

      GLint width() const { return _width; }
      GLint height() const { return _height; }
      GLint depth() const { return _depth; }

      bool is_loaded() const { return _loaded; }

      GLuint texture_id() const { return _texture_id; }

      void fill_pixels(GLubyte* pixels, GLint channels);

    private:
      void _load_pixels(GLubyte* pixels, GLint w, GLint h, GLint d, GLint channels, pixel_type pt);

    private:
      GLenum _filter_mode;
      GLenum _wrap_mode;
      GLenum _format;
      pixel_type _pt;

      GLuint _texture_id;
      GLint _width, _height, _depth;

      bool _loaded;
    };

  class render_buffer
    {
    public:
      render_buffer();
      ~render_buffer();

      void create();
      void bind();
      void release();
      void destroy();

      bool is_created() const { return _render_buffer_id != 0; }
      GLuint object_id() const { return _render_buffer_id; }

    private:
      GLuint _render_buffer_id;
    };

  class frame_buffer_object
    {
    public:
      frame_buffer_object();
      ~frame_buffer_object();

      void create(GLint w, GLint h);
      void bind(GLint channel);
      void release();

      texture* get_texture() { return _texture; }
      render_buffer* get_render_buffer() { return _render_buffer; }

      GLuint frame_buffer_id() const { return _frame_buffer_id; }

      GLint width() const { return _w; }
      GLint height() const { return _h; }

    private:
      texture* _texture;
      render_buffer* _render_buffer;
      int _w, _h;
      GLuint _frame_buffer_id;
    };

  class vertex_array_object
    {
    public:
      vertex_array_object();
      ~vertex_array_object();

      void create();
      void bind();
      void release();
      void destroy();

      bool is_created() const { return _vertex_array_object_id != 0; }
      GLuint object_id() const { return _vertex_array_object_id; }

    private:
      GLuint _vertex_array_object_id;
    };

  enum shader_type {
    Vertex = 1 << 1,
    Fragment = 1 << 2,
    Compute = 1 << 3,
    };

  class shader
    {
    public:
      typedef enum shader_type shader_type;

      shader(shader::shader_type shader_type);
      ~shader();

      bool compile_source_code(const char* source);
      bool compile_source_code(const std::string& source);

      bool is_compiled() const { return _compiled; }

      GLuint shader_id() const { return _shader_id; }
      shader::shader_type shadertype() const { return _shader_type; }

      std::string source_code() const { return _source_code; }
      std::string log() const { return _log; }

    protected:
      bool create();
      void destroy();
      bool compile(const char* source);

    private:
      friend class shader_program;

      bool _compiled;
      GLuint _shader_id;
      shader_type _shader_type;

      std::string _log;
      std::string _source_code;
    };

  class shader_program
    {
    public:
      shader_program();
      ~shader_program();

      const shader* vertex_shader() const { return _vertex_shader; }
      const shader* fragment_shader() const { return _fragment_shader; }

      bool add_shader(shader* shader);
      bool add_shader_from_source(shader::shader_type type, const char* source);
      bool add_shader_from_source(shader::shader_type type, const std::string& source);

      void remove_all_shaders();

      bool create();
      bool link();
      bool bind();
      void release();

      bool is_linked() const { return _linked; }
      const std::string& log() const { return _log; }
      GLuint program_id() const { return _program_id; }

      GLint attribute_location(const char* name);
      GLint attribute_location(const std::string& name);

      void bind_attribute_location(const char* name, int location);
      void bind_attribute_location(const std::string& name, int location);

      void disable_attribute_array(int location);
      void disable_attribute_array(const char* name);
      void disable_attribute_array(const std::string& name);

      void enable_attribute_array(int location);
      void enable_attribute_array(const char* name);
      void enable_attribute_array(const std::string& name);

      void set_attribute_array(int location, const GLfloat* values, int tupleSize, int stride);
      void set_attribute_array(int location, GLenum type, const void* values, int tupleSize, int stride);

      void set_attribute_buffer(int location, GLenum type, int offset, int tupleSize, int stride);

      void set_attribute_array(const char* name, const GLfloat* values, int tupleSize, int stride);
      void set_attribute_array(const char* name, GLenum type, const void* values, int tupleSize, int stride);

      void set_attribute_buffer(const char* name, GLenum type, int offset, int tupleSize, int stride);

      void set_attribute_value(int location, GLfloat value);
      void set_attribute_value(int location, GLfloat x, GLfloat y);
      void set_attribute_value(int location, GLfloat x, GLfloat y, GLfloat z);
      void set_attribute_value(int location, GLfloat x, GLfloat y, GLfloat z, GLfloat w);

      void set_attribute_value(const char* name, GLfloat value);
      void set_attribute_value(const char* name, GLfloat x, GLfloat y);
      void set_attribute_value(const char* name, GLfloat x, GLfloat y, GLfloat z);
      void set_attribute_value(const char* name, GLfloat x, GLfloat y, GLfloat z, GLfloat w);

      void set_uniform_value(int location, GLint value);
      void set_uniform_value(int location, GLuint value);
      void set_uniform_value(int location, GLfloat value);
      void set_uniform_value(int location, GLfloat x, GLfloat y);
      void set_uniform_value(int location, GLfloat x, GLfloat y, GLfloat z);
      void set_uniform_value(int location, GLfloat x, GLfloat y, GLfloat z, GLfloat w);

      void set_uniform_value(const char* name, GLint value);
      void set_uniform_value(const char* name, GLuint value);
      void set_uniform_value(const char* name, GLfloat value);
      void set_uniform_value(const char* name, GLfloat x, GLfloat y);
      void set_uniform_value(const char* name, GLfloat x, GLfloat y, GLfloat z);
      void set_uniform_value(const char* name, GLfloat x, GLfloat y, GLfloat z, GLfloat w);

      void set_uniform_value_array(int location, const GLint* values, int count);
      void set_uniform_value_array(int location, const GLuint* values, int count);
      void set_uniform_value_array(int location, const GLfloat* values, int count, int tupleSize);

      void set_uniform_value_array(const char* name, const GLint* values, int count);
      void set_uniform_value_array(const char* name, const GLuint* values, int count);
      void set_uniform_value_array(const char* name, const GLfloat* values, int count, int tupleSize);

      GLint uniform_location(const char* name);
      GLint uniform_location(const std::string& name);

    private:
      int _linked;
      GLuint _program_id;

      std::string _log;
      shader* _vertex_shader;
      shader* _fragment_shader;
      shader* _compute_shader;
    };


  } // namespace jtk

#endif // #ifndef JTK_OPENGL_H

#ifdef JTK_OPENGL_IMPLEMENTATION

namespace jtk
  {

  void gl_check_error(const char* txt)
    {
    unsigned int err = glGetError();
    if (err)
      {
      std::stringstream str;
      str << "GL error " << err << ": " << txt;
      throw std::runtime_error(str.str());
      }
    }

  buffer_object::buffer_object()
    : buffer_object(GL_VERTEX_ARRAY)
    {
    }

  buffer_object::buffer_object(GLenum type)
    : _type(type),
    _pattern(GL_STATIC_DRAW),
    _buffer_object_id(0)
    {
    }

  buffer_object::~buffer_object()
    {
    destroy();
    }

  void buffer_object::allocate(const void* data, int count)
    {
    //make sure to first bind the buffer before allocate
    if (!is_created())
      return;
    glBufferData(_type, count, data, _pattern);
    gl_check_error("buffer_object.cpp: glBufferData");
    }

  void buffer_object::update_data(const void* data, int count)
    {
    if (!is_created())
      return;
    glBufferSubData(_type, 0, count, data);
    gl_check_error("buffer_object.cpp: glBufferSubData");
    }

  void buffer_object::get_data(void* data, int count)
    {
    if (!is_created())
      return;
    glGetBufferSubData(_type, 0, count, data);
    gl_check_error("buffer_object.cpp: glGetBufferSubData");
    }

  void buffer_object::create()
    {
    glGenBuffers(1, &_buffer_object_id);
    gl_check_error("buffer_object.cpp: glGenBuffers");
    }

  void buffer_object::bind()
    {
    if (_buffer_object_id != 0)
      {
      glBindBuffer(_type, _buffer_object_id);
      gl_check_error("buffer_object.cpp: glBindBuffer");
      }
    }

  void buffer_object::release()
    {
    glBindBuffer(_type, 0);
    }

  void buffer_object::destroy()
    {
    if (_buffer_object_id != 0)
      {
      glDeleteBuffers(1, &_buffer_object_id);
      _buffer_object_id = 0;
      }
    }

  texture::texture() :
    _wrap_mode(GL_REPEAT),
    _filter_mode(GL_LINEAR),
    _texture_id(0),
    _loaded(false)
    {
    glGenTextures(1, &_texture_id);
    gl_check_error("texture.cpp: glGenTextures (in constructor)");
    }

  texture::~texture()
    {
    glDeleteTextures(1, &_texture_id);
    }

  void texture::set_wrap_mode(GLenum wrapMode)
    {
    _wrap_mode = wrapMode;
    }

  void texture::set_filter_mode(GLenum filterMode)
    {
    _filter_mode = filterMode;
    }

  void texture::load_from_pixels(GLubyte* pixels, int w, int h, int channels, pixel_type pt)
    {
    _load_pixels(pixels, w, h, channels, pt);
    }

  void texture::create_empty(int w, int h, GLint channels, pixel_type pt)
    {
    glBindTexture(GL_TEXTURE_2D, _texture_id);

    gl_check_error("texture.cpp: glBindTexture in create_empty");

    _pt = pt;
    _format = GL_RGBA8;

    switch (_pt)
      {
      case int8:
      {
      switch (channels)
        {
        case 1: _format = GL_R8I; break;
        case 2: _format = GL_RG8I; break;
        case 3: _format = GL_RGB8I; break;
        case 4: _format = GL_RGBA8I; break;
        }
      break;
      }
      case uint8:
      {
      switch (channels)
        {
        case 1: _format = GL_R8UI; break;
        case 2: _format = GL_RG8UI; break;
        case 3: _format = GL_RGB8UI; break;
        case 4: _format = GL_RGBA8UI; break;
        }
      break;
      }
      case uint32:
      {
      switch (channels)
        {
        case 1: _format = GL_R32UI; break;
        }
      break;
      }
      case real:
      {
      switch (channels)
        {
        case 1: _format = GL_R32F; break;
        case 2: _format = GL_RG32F; break;
        case 3: _format = GL_RGB32F; break;
        case 4: _format = GL_RGBA32F; break;
        }
      break;
      }
      }
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // opengl by default aligns rows on 4 bytes I think 
    glTexStorage2D(GL_TEXTURE_2D, 1, _format, w, h);
    //glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, _texture_id, 0);

    gl_check_error("texture.cpp: glTexStorage2D in create_empty");

    _width = w;
    _height = h;
    }

  void texture::_load_pixels(GLubyte* pixels, int w, int h, int channels, pixel_type pt)
    {
    _pt = pt;
    _format = GL_RGBA8;
    GLenum sourceFormat = GL_RGBA;
    GLenum pixtype = GL_UNSIGNED_BYTE;

    switch (_pt)
      {
      case uint8:
      {
      switch (channels)
        {
        case 1: _format = GL_R8; sourceFormat = GL_RED;  break;
        case 2: _format = GL_RG8; sourceFormat = GL_RG; break;
        case 3: _format = GL_RGB8; sourceFormat = GL_RGB; break;
        case 4: _format = GL_RGBA8UI; sourceFormat = GL_RGBA_INTEGER; break;
        }
      pixtype = GL_UNSIGNED_BYTE;
      break;
      }
      case uint32:
      {
      assert(0); // todo
      switch (channels)
        {
        case 1: _format = GL_R32UI; break;
        }
      break;
      }
      case real:
      {
      switch (channels)
        {
        case 1: _format = GL_R32F; sourceFormat = GL_RED; break;
        case 2: _format = GL_RG32F; sourceFormat = GL_RG; break;
        case 3: _format = GL_RGB32F; sourceFormat = GL_RGB; break;
        case 4: _format = GL_RGBA32F; sourceFormat = GL_RGBA; break;
        }
      pixtype = GL_FLOAT;
      break;
      }
      }




    glBindTexture(GL_TEXTURE_2D, _texture_id);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // opengl by default aligns rows on 4 bytes I think
    glTexImage2D(GL_TEXTURE_2D, 0, _format, w, h, 0, sourceFormat, pixtype, pixels);
    gl_check_error("texture.cpp: glTexImage2d");

    _width = w;
    _height = h;

    _loaded = true;
    }

  void texture::fill_pixels(GLubyte* pixels, GLint channels)
    {
    GLenum targetFormat = GL_RGBA;
    GLenum pixtype = GL_UNSIGNED_BYTE;

    switch (_pt)
      {
      case uint8:
      {
      switch (channels)
        {
        case 1:  targetFormat = GL_RED_INTEGER;  break;
        case 2:  targetFormat = GL_RG; break;
        case 3:  targetFormat = GL_RGB; break;
        case 4:  targetFormat = GL_RGBA_INTEGER; break;
        }
      pixtype = GL_UNSIGNED_BYTE;
      break;
      }
      case uint32:
      {
      switch (channels)
        {
        case 1: targetFormat = GL_RED_INTEGER; break;
        }
      pixtype = GL_UNSIGNED_INT;
      break;
      }
      case real:
      {
      switch (channels)
        {
        case 1: targetFormat = GL_RED; break;
        case 2: targetFormat = GL_RG; break;
        case 3: targetFormat = GL_RGB; break;
        case 4: targetFormat = GL_RGBA; break;
        }
      pixtype = GL_FLOAT;
      break;
      }
      }
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // opengl by default aligns rows on 4 bytes I think
    glGetTexImage(GL_TEXTURE_2D, 0, targetFormat, pixtype, (void*)pixels);
    gl_check_error("texture.cpp: glGetTexImage");
    }


  void texture::bind_to_channel(int channel)
    {
    glActiveTexture(GL_TEXTURE0 + channel);
    gl_check_error("texture.cpp: glActiveTexture");
    glBindTexture(GL_TEXTURE_2D, _texture_id);
    gl_check_error("texture.cpp: glBindTexture");

    if (_wrap_mode == GL_REPEAT)
      {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
      }
    else
      {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      }

    if (_filter_mode == GL_NEAREST)
      {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      }
    else if (_filter_mode == GL_MIPMAP)
      {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      }
    else
      {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      }

    }

  void texture::bind_to_image(int channel, GLenum access)
    {
    glBindImageTexture(channel, _texture_id, 0, GL_FALSE, 0, access, _format);
    gl_check_error("texture.cpp: glBindImageTexture");
    }

  void texture::release()
    {
    glBindTexture(GL_TEXTURE_2D, 0);
    }

  texture3d::texture3d() :
    _wrap_mode(GL_REPEAT),
    _filter_mode(GL_LINEAR),
    _texture_id(0),
    _loaded(false)
    {
    glGenTextures(1, &_texture_id);
    gl_check_error("texture3d.cpp: glGenTextures (in constructor)");
    }

  texture3d::~texture3d()
    {
    glDeleteTextures(1, &_texture_id);
    }

  void texture3d::set_wrap_mode(GLenum wrapMode)
    {
    _wrap_mode = wrapMode;
    }

  void texture3d::set_filter_mode(GLenum filterMode)
    {
    _filter_mode = filterMode;
    }

  void texture3d::load_from_pixels(GLubyte* pixels, GLint w, GLint h, GLint d, GLint channels, pixel_type pt)
    {
    _load_pixels(pixels, w, h, d, channels, pt);
    }

  void texture3d::create_empty(GLint w, GLint h, GLint d, GLint channels, pixel_type pt)
    {
    glBindTexture(GL_TEXTURE_3D, _texture_id);

    gl_check_error("texture3d.cpp: glBindTexture in create_empty");

    _pt = pt;
    _format = GL_RGBA8;

    switch (_pt)
      {
      case int8:
      {
      switch (channels)
        {
        case 1: _format = GL_R8I; break;
        case 2: _format = GL_RG8I; break;
        case 3: _format = GL_RGB8I; break;
        case 4: _format = GL_RGBA8I; break;
        }
      break;
      }
      case uint8:
      {
      switch (channels)
        {
        case 1: _format = GL_R8UI; break;
        case 2: _format = GL_RG8UI; break;
        case 3: _format = GL_RGB8UI; break;
        case 4: _format = GL_RGBA8UI; break;
        }
      break;
      }
      case uint32:
      {
      switch (channels)
        {
        case 1: _format = GL_R32UI; break;
        }
      break;
      }
      case real:
      {
      switch (channels)
        {
        case 1: _format = GL_R32F; break;
        case 2: _format = GL_RG32F; break;
        case 3: _format = GL_RGB32F; break;
        case 4: _format = GL_RGBA32F; break;
        }
      break;
      }
      }
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // opengl by default aligns rows on 4 bytes I think     
    glTexStorage3D(GL_TEXTURE_3D, 1, _format, w, h, d);

    gl_check_error("texture3d.cpp: glTexStorage3D in create_empty");

    _width = w;
    _height = h;
    }

  void texture3d::_load_pixels(GLubyte* pixels, GLint w, GLint h, GLint d, GLint channels, pixel_type pt)
    {
    _pt = pt;
    _format = GL_RGBA8;
    GLenum sourceFormat = GL_RGBA;
    GLenum pixtype = GL_UNSIGNED_BYTE;

    switch (_pt)
      {
      case uint8:
      {
      switch (channels)
        {
        case 1: _format = GL_R8; sourceFormat = GL_RED;  break;
        case 2: _format = GL_RG8; sourceFormat = GL_RG; break;
        case 3: _format = GL_RGB8; sourceFormat = GL_RGB; break;
        case 4: _format = GL_RGBA8UI; sourceFormat = GL_RGBA_INTEGER; break;
        }
      pixtype = GL_UNSIGNED_BYTE;
      break;
      }
      case uint32:
      {
      switch (channels)
        {
        case 1: _format = GL_R32UI; sourceFormat = GL_RED_INTEGER; break;
        case 2: _format = GL_RG32UI; sourceFormat = GL_RG_INTEGER; break;
        case 3: _format = GL_RGB32UI; sourceFormat = GL_RGB_INTEGER; break;
        case 4: _format = GL_RGBA32UI; sourceFormat = GL_RGBA_INTEGER; break;
        }
      pixtype = GL_UNSIGNED_INT;
      break;
      }
      case real:
      {
      switch (channels)
        {
        case 1: _format = GL_R32F; sourceFormat = GL_RED; break;
        case 2: _format = GL_RG32F; sourceFormat = GL_RG; break;
        case 3: _format = GL_RGB32F; sourceFormat = GL_RGB; break;
        case 4: _format = GL_RGBA32F; sourceFormat = GL_RGBA; break;
        }
      pixtype = GL_FLOAT;
      break;
      }
      }




    glBindTexture(GL_TEXTURE_3D, _texture_id);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // opengl by default aligns rows on 4 bytes I think    
    glTexImage3D(GL_TEXTURE_3D, 0, _format, w, h, d, 0, sourceFormat, pixtype, pixels);
    gl_check_error("texture3d.cpp: glTexImage3d");
    //glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, w, h, d, sourceFormat, pixtype, pixels);
    //gl_check_error("texture3d.cpp: glTexSubImage3D");

    _width = w;
    _height = h;
    _depth = d;

    _loaded = true;
    }

  void texture3d::fill_pixels(GLubyte* pixels, GLint channels)
    {
    GLenum targetFormat = GL_RGBA;
    GLenum pixtype = GL_UNSIGNED_BYTE;

    switch (_pt)
      {
      case uint8:
      {
      switch (channels)
        {
        case 1:  targetFormat = GL_RED_INTEGER;  break;
        case 2:  targetFormat = GL_RG; break;
        case 3:  targetFormat = GL_RGB; break;
        case 4:  targetFormat = GL_RGBA_INTEGER; break;
        }
      pixtype = GL_UNSIGNED_BYTE;
      break;
      }
      case uint32:
      {
      switch (channels)
        {
        case 1: targetFormat = GL_RED_INTEGER; break;
        }
      pixtype = GL_UNSIGNED_INT;
      break;
      }
      case real:
      {
      switch (channels)
        {
        case 1: targetFormat = GL_RED; break;
        case 2: targetFormat = GL_RG; break;
        case 3: targetFormat = GL_RGB; break;
        case 4: targetFormat = GL_RGBA; break;
        }
      pixtype = GL_FLOAT;
      break;
      }
      }
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // opengl by default aligns rows on 4 bytes I think    
    glGetTexImage(GL_TEXTURE_3D, 0, targetFormat, pixtype, (void*)pixels);
    gl_check_error("texture3d.cpp: glGetTexImage");
    }


  void texture3d::bind_to_channel(GLint channel)
    {
    glActiveTexture(GL_TEXTURE0 + channel);
    gl_check_error("texture3d.cpp: glActiveTexture");
    glBindTexture(GL_TEXTURE_3D, _texture_id);
    gl_check_error("texture3d.cpp: glBindTexture");
    if (_wrap_mode == GL_REPEAT)
      {
      glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_REPEAT);
      glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_REPEAT);
      glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_REPEAT);
      }
    else
      {
      glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
      }

    if (_filter_mode == GL_NEAREST)
      {
      glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      }
    else if (_filter_mode == GL_MIPMAP)
      {
      glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
      glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      }
    else
      {
      glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      }
    }

  void texture3d::bind_to_image(GLint channel, GLenum access)
    {
    glBindImageTexture(channel, _texture_id, 0, GL_TRUE, 0, access, _format);
    gl_check_error("texture3d.cpp: glBindImageTexture");
    }

  void texture3d::release()
    {
    glBindTexture(GL_TEXTURE_3D, 0);
    }

  render_buffer::render_buffer()
    : _render_buffer_id(0)
    {
    }

  render_buffer::~render_buffer()
    {
    destroy();
    }

  void render_buffer::create()
    {
    glGenRenderbuffersEXT(1, &_render_buffer_id);
    }

  void render_buffer::bind()
    {
    glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, _render_buffer_id);
    }

  void render_buffer::release()
    {
    glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);
    }

  void render_buffer::destroy()
    {
    if (_render_buffer_id != 0)
      {
      glDeleteRenderbuffersEXT(1, &_render_buffer_id);
      _render_buffer_id = 0;
      }
    }


  frame_buffer_object::frame_buffer_object()
    : _w(0), _h(0), _texture(nullptr), _render_buffer(nullptr)
    {
    }

  frame_buffer_object::~frame_buffer_object()
    {
    release();
    delete _texture;
    delete _render_buffer;
    }

  void frame_buffer_object::create(GLint w, GLint h)
    {
    _w = w;
    _h = h;
    _texture = new texture();
    _texture->create_empty(w, h, 4, texture::uint8);
    gl_check_error("_texture->create_empty()");

    _render_buffer = new render_buffer();
    _render_buffer->create();
    gl_check_error("_render_buffer->create()");

    _texture->bind_to_channel(10);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glGenFramebuffersEXT(1, &_frame_buffer_id);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _frame_buffer_id);
    _render_buffer->bind();

    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, _texture->texture_id(), 0);
    glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT24, _texture->width(), _texture->height());
    glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, _render_buffer->object_id());

    GLenum status;
    status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
    switch (status)
      {
      case GL_FRAMEBUFFER_COMPLETE_EXT:
      {
      break;
      }
      default:
      {
      throw std::runtime_error("frame buffer object is not complete");
      }
      }
    }

  void frame_buffer_object::bind(GLint channel)
    {
    _render_buffer->bind();
    _texture->bind_to_channel(channel);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _frame_buffer_id);
    }

  void frame_buffer_object::release()
    {
    _render_buffer->release();
    _texture->release();
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
    }

  vertex_array_object::vertex_array_object()
    : _vertex_array_object_id(0)
    {
    }

  vertex_array_object::~vertex_array_object()
    {
    destroy();
    }

  void vertex_array_object::create()
    {
    glGenVertexArrays(1, &_vertex_array_object_id);
    gl_check_error("vertex_array_object.cpp: create");
    }

  void vertex_array_object::bind()
    {
    glBindVertexArray(_vertex_array_object_id);
    }

  void vertex_array_object::release()
    {
    glBindVertexArray(0);
    }

  void vertex_array_object::destroy()
    {
    if (_vertex_array_object_id != 0)
      {
      glDeleteVertexArrays(1, &_vertex_array_object_id);
      _vertex_array_object_id = 0;
      }
    }

  static const char* types[] = {
    "Vertex",
    "Fragment",
    "Compute",
    ""
  };

  shader::shader(shader::shader_type shader_type)
    : _shader_id(0),
    _shader_type(shader_type)
    {
    }

  shader::~shader()
    {
    }

  bool shader::compile_source_code(const char* source)
    {
    return compile(source);
    }

  bool shader::compile_source_code(const std::string& source)
    {
    return compile(source.c_str());
    }

  bool shader::create()
    {
    if (_shader_type == shader::shader_type::Compute)
      _shader_id = glCreateShader(GL_COMPUTE_SHADER);
    else if (_shader_type == shader::shader_type::Vertex)
      _shader_id = glCreateShader(GL_VERTEX_SHADER);
    else if (_shader_type == shader::shader_type::Fragment)
      _shader_id = glCreateShader(GL_FRAGMENT_SHADER);

    if (!_shader_id)
      {
      if (_shader_type == shader::shader_type::Compute)
        std::cout << "Could not create shader of type Compute\n";
      if (_shader_type == shader::shader_type::Vertex)
        std::cout << "Could not create shader of type Vertex\n";
      if (_shader_type == shader::shader_type::Fragment)
        std::cout << "Could not create shader of type Fragment\n";
      else
        std::cout << "Could not create shader\n";
      return false;
      }
    else
      return true;
    }

  void shader::destroy()
    {
    if (!_shader_id)
      return;

    glDeleteShader(_shader_id);
    _shader_id = 0;
    }

  bool shader::compile(const char* source)
    {
    if (!create())
      return false;

    glShaderSource(_shader_id, 1, &source, nullptr);
    glCompileShader(_shader_id);

    int value;
    glGetShaderiv(_shader_id, GL_COMPILE_STATUS, &value);
    _compiled = (value != 0);

    int source_codeLength = 0;
    glGetShaderiv(_shader_id, GL_SHADER_SOURCE_LENGTH, &source_codeLength);
    if (source_codeLength > 1)
      {
      int temp = 0;
      _source_code.resize(source_codeLength);
      glGetShaderSource(_shader_id, source_codeLength, &temp, &_source_code[0]);
      }

    if (!_compiled)
      {
      glGetShaderiv(_shader_id, GL_INFO_LOG_LENGTH, &value);

      if (value > 1)
        {
        int length;
        _log.resize(value);
        glGetShaderInfoLog(_shader_id, value, &length, &_log[0]);

        const char* type = types[2];
        if (_shader_type == shader::shader_type::Vertex)
          type = types[0];
        else if (_shader_type == shader::shader_type::Fragment)
          type = types[1];
        std::cout << _log.c_str() << "\n";
        }
      }

    return _compiled;
    }
  shader_program::shader_program()
    : _program_id(0),
    _vertex_shader(nullptr),
    _fragment_shader(nullptr),
    _compute_shader(nullptr)
    {
    }

  shader_program::~shader_program()
    {
    remove_all_shaders();
    }

  bool shader_program::add_shader(shader* shader)
    {
    if (shader->shadertype() == shader::shader_type::Compute)
      {
      if (_compute_shader != nullptr)
        return false;
      _compute_shader = shader;
      }
    else if (shader->shadertype() == shader::shader_type::Vertex)
      {
      if (_vertex_shader != nullptr)
        return false;
      _vertex_shader = shader;
      }
    else if (shader->shadertype() == shader::shader_type::Fragment)
      {
      if (_fragment_shader != nullptr)
        return false;
      _fragment_shader = shader;
      }

    return true;
    }

  bool shader_program::add_shader_from_source(shader::shader_type type, const char* source)
    {
    shader* s = new shader(type);

    if (s->compile_source_code(source))
      {
      if (s->is_compiled())
        {
        if (add_shader(s))
          return true;
        else
          {
          s->destroy();
          delete s;
          return false;
          }
        }
      else
        {
        std::cout << "Compile shader error: " << s->log().c_str() << "\n";
        s->destroy();
        delete s;
        return false;
        }
      }
    s->destroy();
    delete s;
    return false;
    }

  bool shader_program::add_shader_from_source(shader::shader_type type, const std::string& source)
    {
    return add_shader_from_source(type, source.c_str());
    }

  void shader_program::remove_all_shaders()
    {
    if (_program_id && _vertex_shader)
      {
      glDetachShader(_program_id, _vertex_shader->shader_id());
      _vertex_shader->destroy();
      delete _vertex_shader;
      _vertex_shader = nullptr;

      _linked = false;
      }

    if (_program_id && _fragment_shader)
      {
      glDetachShader(_program_id, _fragment_shader->shader_id());
      _fragment_shader->destroy();
      delete _fragment_shader;
      _fragment_shader = nullptr;

      _linked = false;
      }

    if (_program_id && _compute_shader)
      {
      glDetachShader(_program_id, _compute_shader->shader_id());
      _compute_shader->destroy();
      delete _compute_shader;
      _compute_shader = nullptr;

      _linked = false;
      }
    }

  bool shader_program::create()
    {
    _program_id = glCreateProgram();

    if (!_program_id)
      {
      std::cout << "Could not create program object\n";
      return false;
      }

    return true;
    }

  bool shader_program::link()
    {
    _linked = false;

    if (!create())
      return false;

    if ((!_vertex_shader || !_fragment_shader) && !_compute_shader)
      return false;

    if (_vertex_shader)
      {
      glAttachShader(_program_id, _vertex_shader->shader_id());
      gl_check_error("shader_program.cpp: glAttachShader");
      }
    if (_fragment_shader)
      {
      glAttachShader(_program_id, _fragment_shader->shader_id());
      gl_check_error("shader_program.cpp: glAttachShader");
      }
    if (_compute_shader)
      {
      glAttachShader(_program_id, _compute_shader->shader_id());
      gl_check_error("shader_program.cpp: glAttachShader");
      }
    glLinkProgram(_program_id);
    gl_check_error("shader_program.cpp: glLinkProgram");

    int value = 0;
    glGetProgramiv(_program_id, GL_LINK_STATUS, &value);
    _linked = (value != 0);

    if (!_linked)
      {
      int length = 0;
      glGetProgramiv(_program_id, GL_INFO_LOG_LENGTH, &value);
      if (value > 1)
        {
        _log.resize(value);
        glGetProgramInfoLog(_program_id, value, &length, &_log[0]);
        std::cout << "shader program: link error: " << _log.c_str() << "\n";
        }

      remove_all_shaders();
      }

    return _linked;
    }

  bool shader_program::bind()
    {
    if (!_program_id || !_linked)
      return false;

    glUseProgram(_program_id);
    gl_check_error("shader_program.cpp: glUseProgram");
    return true;
    }

  void shader_program::release()
    {
    glUseProgram(0);
    }

  GLint shader_program::attribute_location(const char* name)
    {
    return glGetAttribLocation(_program_id, name);;
    }

  GLint shader_program::attribute_location(const std::string& name)
    {
    return glGetAttribLocation(_program_id, name.c_str());
    }

  void shader_program::bind_attribute_location(const char* name, int location)
    {
    glBindAttribLocation(_program_id, location, name);
    }

  void shader_program::bind_attribute_location(const std::string& name, int location)
    {
    glBindAttribLocation(_program_id, location, name.c_str());
    }

  void shader_program::disable_attribute_array(int location)
    {
    glDisableVertexAttribArray(location);
    }

  void shader_program::disable_attribute_array(const char* name)
    {
    int location = glGetAttribLocation(_program_id, name);
    glDisableVertexAttribArray(location);
    }

  void shader_program::disable_attribute_array(const std::string& name)
    {
    int location = glGetAttribLocation(_program_id, name.c_str());
    glDisableVertexAttribArray(location);
    }

  void shader_program::enable_attribute_array(int location)
    {
    glEnableVertexAttribArray(location);
    }

  void shader_program::enable_attribute_array(const char* name)
    {
    int location = glGetAttribLocation(_program_id, name);
    glEnableVertexAttribArray(location);
    }

  void shader_program::enable_attribute_array(const std::string& name)
    {
    int location = glGetAttribLocation(_program_id, name.c_str());
    glEnableVertexAttribArray(location);
    }

  void shader_program::set_attribute_array(int location, const GLfloat* values, int tupleSize, int stride)
    {
    glVertexAttribPointer(location, tupleSize, GL_FLOAT, GL_FALSE, stride, values);
    }

  void shader_program::set_attribute_array(int location, GLenum type, const void* values, int tupleSize, int stride)
    {
    glVertexAttribPointer(location, tupleSize, type, GL_TRUE, stride, values);
    }

  void shader_program::set_attribute_buffer(int location, GLenum type, int offset, int tupleSize, int stride)
    {
    glVertexAttribPointer(location, tupleSize, type, GL_TRUE, stride,
      reinterpret_cast<const void*>(intptr_t(offset)));
    }

  void shader_program::set_attribute_array(const char* name, const GLfloat* values, int tupleSize, int stride)
    {
    GLint location = glGetAttribLocation(_program_id, name);;
    glVertexAttribPointer(location, tupleSize, GL_FLOAT, GL_FALSE, stride, values);
    }

  void shader_program::set_attribute_array(const char* name, GLenum type, const void* values, int tupleSize, int stride)
    {
    GLint location = glGetAttribLocation(_program_id, name);;
    glVertexAttribPointer(location, tupleSize, type, GL_TRUE, stride, values);
    }

  void shader_program::set_attribute_buffer(const char* name, GLenum type, int offset, int tupleSize, int stride)
    {
    GLint location = glGetAttribLocation(_program_id, name);;
    glVertexAttribPointer(location, tupleSize, type, GL_TRUE, stride,
      reinterpret_cast<const void*>(intptr_t(offset)));
    }

  void shader_program::set_attribute_value(int location, const GLfloat value)
    {
    glVertexAttrib1f(location, value);
    }

  void shader_program::set_attribute_value(int location, const GLfloat x, const GLfloat y)
    {
    glVertexAttrib2f(location, x, y);
    }

  void shader_program::set_attribute_value(int location, const GLfloat x, const GLfloat y, const GLfloat z)
    {
    glVertexAttrib3f(location, x, y, z);
    }

  void shader_program::set_attribute_value(int location, const GLfloat x, const GLfloat y, const GLfloat z, const GLfloat w)
    {
    glVertexAttrib4f(location, x, y, z, w);
    }

  void shader_program::set_attribute_value(const char* name, const GLfloat value)
    {
    GLint location = glGetAttribLocation(_program_id, name);
    glVertexAttrib1f(location, value);
    }

  void shader_program::set_attribute_value(const char* name, const GLfloat x, const GLfloat y)
    {
    GLint location = glGetAttribLocation(_program_id, name);
    glVertexAttrib2f(location, x, y);
    }

  void shader_program::set_attribute_value(const char* name, const GLfloat x, const GLfloat y, const GLfloat z)
    {
    GLint location = glGetAttribLocation(_program_id, name);
    glVertexAttrib3f(location, x, y, z);
    }

  void shader_program::set_attribute_value(const char* name, const GLfloat x, const GLfloat y, const GLfloat z, const GLfloat w)
    {
    GLint location = glGetAttribLocation(_program_id, name);
    glVertexAttrib4f(location, x, y, z, w);
    }

  void shader_program::set_uniform_value(int location, GLint value)
    {
    glUniform1i(location, value);
    }

  void shader_program::set_uniform_value(int location, GLuint value)
    {
    glUniform1i(location, value);
    }

  void shader_program::set_uniform_value(int location, const GLfloat value)
    {
    glUniform1f(location, value);
    }

  void shader_program::set_uniform_value(int location, const GLfloat x, const GLfloat y)
    {
    glUniform2f(location, x, y);
    }

  void shader_program::set_uniform_value(int location, const GLfloat x, const GLfloat y, const GLfloat z)
    {
    glUniform3f(location, x, y, z);
    }

  void shader_program::set_uniform_value(int location, const GLfloat x, const GLfloat y, const GLfloat z, const GLfloat w)
    {
    glUniform4f(location, x, y, z, w);
    }

  void shader_program::set_uniform_value(const char* name, GLint value)
    {
    GLint location = glGetUniformLocation(_program_id, name);
    glUniform1i(location, value);
    }

  void shader_program::set_uniform_value(const char* name, GLuint value)
    {
    GLint location = glGetUniformLocation(_program_id, name);
    glUniform1i(location, value);
    }

  void shader_program::set_uniform_value(const char* name, const GLfloat value)
    {
    GLint location = glGetUniformLocation(_program_id, name);
    glUniform1f(location, value);
    }

  void shader_program::set_uniform_value(const char* name, const GLfloat x, const GLfloat y)
    {
    GLint location = glGetUniformLocation(_program_id, name);
    glUniform2f(location, x, y);
    }

  void shader_program::set_uniform_value(const char* name, const GLfloat x, const GLfloat y, const GLfloat z)
    {
    GLint location = glGetUniformLocation(_program_id, name);
    glUniform3f(location, x, y, z);
    }

  void shader_program::set_uniform_value(const char* name, const GLfloat x, const GLfloat y, const GLfloat z, const GLfloat w)
    {
    GLint location = glGetUniformLocation(_program_id, name);
    glUniform4f(location, x, y, z, w);
    }

  void shader_program::set_uniform_value_array(int location, const GLint* values, int count)
    {
    glUniform1iv(location, count, values);
    }

  void shader_program::set_uniform_value_array(int location, const GLuint* values, int count)
    {
    const GLint* iv = reinterpret_cast<const GLint*>(values);
    glUniform1iv(location, count, iv);
    }

  void shader_program::set_uniform_value_array(int location, const GLfloat* values, int count, int tupleSize)
    {
    if (tupleSize == 1)
      glUniform1fv(location, count, values);
    else if (tupleSize == 2)
      glUniform2fv(location, count, values);
    else if (tupleSize == 3)
      glUniform3fv(location, count, values);
    else if (tupleSize == 4)
      glUniform4fv(location, count, values);
    }

  void shader_program::set_uniform_value_array(const char* name, const GLint* values, int count)
    {
    GLint location = glGetUniformLocation(_program_id, name);
    glUniform1iv(location, count, values);
    }

  void shader_program::set_uniform_value_array(const char* name, const GLuint* values, int count)
    {
    GLint location = glGetUniformLocation(_program_id, name);
    const GLint* iv = reinterpret_cast<const GLint*>(values);
    glUniform1iv(location, count, iv);
    }

  void shader_program::set_uniform_value_array(const char* name, const GLfloat* values, int count, int tupleSize)
    {
    GLint location = glGetUniformLocation(_program_id, name);

    if (tupleSize == 1)
      glUniform1fv(location, count, values);
    else if (tupleSize == 2)
      glUniform2fv(location, count, values);
    else if (tupleSize == 3)
      glUniform3fv(location, count, values);
    else if (tupleSize == 4)
      glUniform4fv(location, count, values);
    }

  GLint shader_program::uniform_location(const char* name)
    {
    return glGetUniformLocation(_program_id, name);
    }

  GLint shader_program::uniform_location(const std::string& name)
    {
    return glGetUniformLocation(_program_id, name.c_str());
    }


  } // namespace jtk

#endif //JTK_OPENGL_IMPLEMENTATION
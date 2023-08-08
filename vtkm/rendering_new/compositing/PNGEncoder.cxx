//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/cont/ErrorBadValue.h>
#include <vtkm/internal/Configure.h>
#include <vtkm/rendering_new/compositing/PNGEncoder.h>
#include <vtkm/thirdparty/lodepng/vtkmlodepng/lodepng.h>

#include <iostream>
#include <stdlib.h>

namespace vtkm
{
namespace rendering_new
{

PNGEncoder::PNGEncoder() {}

PNGEncoder::~PNGEncoder()
{
  this->Cleanup();
}

void PNGEncoder::Encode(const unsigned char* rgba_in, const int width, const int height)
{
  this->Cleanup();

  // upside down relative to what lodepng wants
  std::vector<unsigned char> rgba_flip(static_cast<std::size_t>(width * height * 4));
  //unsigned char* rgba_flip = new unsigned char[static_cast<std::size_t>(width * height * 4)];

  for (int y = 0; y < height; ++y)
  {
    memcpy(&(rgba_flip[y * width * 4]),
           &(rgba_in[(height - y - 1) * width * 4]),
           static_cast<std::size_t>(width * 4));
  }

  vtkm::png::LodePNGState state;
  vtkm::png::lodepng_state_init(&state);
  // use less aggressive compression
  state.encoder.zlibsettings.btype = 2;
  state.encoder.zlibsettings.use_lz77 = 0;

  unsigned error = vtkm::png::lodepng_encode(&this->Buffer,
                                             &this->BufferSize,
                                             rgba_flip.data(),
                                             static_cast<unsigned>(width),
                                             static_cast<unsigned>(height),
                                             &state);
  if (error)
    throw vtkm::cont::ErrorBadValue("lodepng_encode failed");
}

void PNGEncoder::Encode(const float* rgba_in, const int width, const int height)
{
  this->Cleanup();

  // upside down relative to what lodepng wants
  std::vector<unsigned char> rgba_flip(static_cast<std::size_t>(width * height * 4));

  for (int x = 0; x < width; ++x)

#ifdef VTKH_OPENMP_ENABLED
#pragma omp parallel for
#endif
    for (int y = 0; y < height; ++y)
    {
      int inOffset = (y * width + x) * 4;
      int outOffset = ((height - y - 1) * width + x) * 4;
      rgba_flip[outOffset + 0] = (unsigned char)(rgba_in[inOffset + 0] * 255.f);
      rgba_flip[outOffset + 1] = (unsigned char)(rgba_in[inOffset + 1] * 255.f);
      rgba_flip[outOffset + 2] = (unsigned char)(rgba_in[inOffset + 2] * 255.f);
      rgba_flip[outOffset + 3] = (unsigned char)(rgba_in[inOffset + 3] * 255.f);
    }

  vtkm::png::LodePNGState state;
  vtkm::png::lodepng_state_init(&state);
  // use less aggressive compression
  state.encoder.zlibsettings.btype = 2;
  state.encoder.zlibsettings.use_lz77 = 0;

  unsigned error = vtkm::png::lodepng_encode(&this->Buffer,
                                             &this->BufferSize,
                                             rgba_flip.data(),
                                             static_cast<unsigned>(width),
                                             static_cast<unsigned>(height),
                                             &state);
  if (error)
    throw vtkm::cont::ErrorBadValue("lodepng_encode failed");
}

void PNGEncoder::Encode(const unsigned char* rgba_in,
                        const int width,
                        const int height,
                        const std::vector<std::string>& comments)
{
  this->Cleanup();

  // upside down relative to what lodepng wants
  std::vector<unsigned char> rgba_flip(static_cast<std::size_t>(width * height * 4));

  for (int y = 0; y < height; ++y)
  {
    memcpy(&(rgba_flip[y * width * 4]),
           &(rgba_in[(height - y - 1) * width * 4]),
           static_cast<std::size_t>(width * 4));
  }

  vtkm::png::LodePNGState state;
  vtkm::png::lodepng_state_init(&state);
  // use less aggressive compression
  state.encoder.zlibsettings.btype = 2;
  state.encoder.zlibsettings.use_lz77 = 0;
  if (comments.size() % 2 != 0)
  {
    std::cerr << "PNGEncoder::Encode comments missing value for the last key.\n";
    std::cerr << "Ignoring the last key.\n";
  }
  if (comments.size() > 1)
  {
    vtkm::png::lodepng_info_init(&state.info_png);
    // Comments are in pairs with a key and a value, using
    // comments.size()-1 ensures that we don't use the last
    // comment if the length of the vector isn't a multiple of 2.
    for (std::size_t i = 0; i < comments.size() - 1; i += 2)
      vtkm::png::lodepng_add_text(&state.info_png, comments[i].c_str(), comments[i + 1].c_str());
  }

  unsigned error = vtkm::png::lodepng_encode(&this->Buffer,
                                             &this->BufferSize,
                                             &rgba_flip[0],
                                             static_cast<unsigned>(width),
                                             static_cast<unsigned>(height),
                                             &state);

  if (error)
    throw vtkm::cont::ErrorBadValue("lodepng_encode failed");
}

void PNGEncoder::Encode(const float* rgba_in,
                        const int width,
                        const int height,
                        const std::vector<std::string>& comments)
{
  this->Cleanup();

  // upside down relative to what lodepng wants
  std::vector<unsigned char> rgba_flip(static_cast<std::size_t>(width * height * 4));

  for (int x = 0; x < width; ++x)

#ifdef VTKH_OPENMP_ENABLED
#pragma omp parallel for
#endif
    for (int y = 0; y < height; ++y)
    {
      int inOffset = (y * width + x) * 4;
      int outOffset = ((height - y - 1) * width + x) * 4;
      rgba_flip[outOffset + 0] = (unsigned char)(rgba_in[inOffset + 0] * 255.f);
      rgba_flip[outOffset + 1] = (unsigned char)(rgba_in[inOffset + 1] * 255.f);
      rgba_flip[outOffset + 2] = (unsigned char)(rgba_in[inOffset + 2] * 255.f);
      rgba_flip[outOffset + 3] = (unsigned char)(rgba_in[inOffset + 3] * 255.f);
    }

  vtkm::png::LodePNGState state;
  vtkm::png::lodepng_state_init(&state);
  // use less aggressive compression
  state.encoder.zlibsettings.btype = 2;
  state.encoder.zlibsettings.use_lz77 = 0;
  if (comments.size() % 2 != 0)
  {
    std::cerr << "PNGEncoder::Encode comments missing value for the last key.\n";
    std::cerr << "Ignoring the last key.\n";
  }
  if (comments.size() > 1)
  {
    vtkm::png::lodepng_info_init(&state.info_png);
    // Comments are in pairs with a key and a value, using
    // comments.size()-1 ensures that we don't use the last
    // comment if the length of the vector isn't a multiple of 2.
    for (std::size_t i = 0; i < comments.size() - 1; i += 2)
      vtkm::png::lodepng_add_text(&state.info_png, comments[i].c_str(), comments[i + 1].c_str());
  }

  unsigned error = vtkm::png::lodepng_encode(&this->Buffer,
                                             &this->BufferSize,
                                             rgba_flip.data(),
                                             static_cast<unsigned>(width),
                                             static_cast<unsigned>(height),
                                             &state);
  if (error)
    throw vtkm::cont::ErrorBadValue("lodepng_encode failed");
}

void PNGEncoder::Save(const std::string& filename)
{
  if (!this->Buffer)
    throw vtkm::cont::ErrorBadValue("Save must be called AFTER Encode");

  unsigned error = vtkm::png::lodepng_save_file(this->Buffer, this->BufferSize, filename.c_str());
  if (error)
    throw vtkm::cont::ErrorBadValue("lodepng_encode failed");
}

void PNGEncoder::Cleanup()
{
  if (this->Buffer)
    delete[] this->Buffer;
  this->BufferSize = 0;
}

}
} //namespace vtkm::rendering_new

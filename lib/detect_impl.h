/* -*- c++ -*- */
/* 
 * Copyright 2019 - Quiet Skies LLC -- Glen Langston - glen.i.langston@gmail.com
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifndef INCLUDED_RADIO_ASTRO_DETECT_IMPL_H
#define INCLUDED_RADIO_ASTRO_DETECT_IMPL_H

#include <radio_astro/detect.h>

#define MAX_VLEN 16384
#define MAX_BUFF (2*MAX_VLEN)

namespace gr {
  namespace radio_astro {

    class detect_impl : public detect
    {
     private:
      // values computed in this block
      int d_vec_length;
      float d_dms;
      float d_f_obs;
      float d_bw;
      float d_t_int;
      int d_nt;
      double nsigma = 0;
      double peak = 0;        // peak, rms and date/time of detected event
      double rms = 0;         // rms of values in circular buffer
      double mjd = 0;         // modified Julian Date of event
      float circular[MAX_BUFF];   // circular buffer for input samples
      float circular2[MAX_BUFF];   // circular buffer for input samples**2
      long inext = 0;         // next place for a sample in buffer
      long inext2 = MAX_BUFF/2;  // place to check for new peak
      long imax2 = 0;         // index to last maximum
      double max2 = 0;        // max value squared so far
      double sum2 = 0;        // sum of values squared
      double rms2 = 0;        // rms squared of values in circular buffer
      double oneovern = 1./double(MAX_BUFF);
      bool bufferfull = 0;    // assume buffer is not full 
      double nsigma_rms = 0;  // comparision value for event detection
      float samples[MAX_VLEN];  // output event buffer 
      bool initialized = 0;   // flag initializing output
      
     public:
      detect_impl(int vec_length,float dms, float f_obs, float bw, float t_int, int nt);
      ~detect_impl();

      // Where all the action really happens
      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      //      virtual void set_dms( float dms);
      void set_dms( float dms);
      
      int update_buffer();

      int event(const float *input, float *output);

      int general_work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
    };

  } // namespace radio_astro
} // namespace gr

#endif /* INCLUDED_RADIO_ASTRO_DETECT_IMPL_H */


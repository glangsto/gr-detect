/* -*- c++ -*- */
/* 
 * Copyright 2018 <+YOU OR YOUR COMPANY+>.
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "detect_impl.h"
#include <iostream>

namespace gr {
  namespace radio_astro {

    detect::sptr
    detect::make(int vec_length, float dms, float f_obs, float bw, float t_int, int nt)
    {
      return gnuradio::get_initial_sptr
        (new detect_impl(vec_length, dms, f_obs, bw, t_int, nt));
    }

    /*
     * The private constructor
     */
    detect_impl::detect_impl(int vec_length, float dms, float f_obs, float bw, float t_int, int nt)
      : gr::block("detect",
		  gr::io_signature::make(1, 1, sizeof(float)*vec_length),
		  gr::io_signature::make(1, 1, sizeof(float)*vec_length)),
        d_vec_length(vec_length),
        d_dms(dms),
        d_f_obs(f_obs),
        d_bw(bw),
        d_t_int(t_int),
        d_nt(nt)
    {}

    /*
     * Our virtual destructor.
     */
    detect_impl::~detect_impl()
    {
    }

    void
    detect_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      /* <+forecast+> e.g. ninput_items_required[0] = noutput_items */
      unsigned ninputs = ninput_items_required.size();
      for(unsigned int i = 0; i < ninputs; i++)
       	    ninput_items_required[i] = noutput_items;
    }

    int
    detect_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      const float *in = (const float *) input_items[0];
      float *out = (float *) output_items[0];
      unsigned ninputs = ninput_items.size();
      int success;

      success = event(in, out);
      if (success == 0)
	{
	  std::cout << "Processed " << ninputs << " Vectors" << "\n";
	}
      //std::cout << out[0*d_dms+ 0] << " " << out[31*d_dms+49] <<"\n";

      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each (noutput_items);

      // Tell runtime system how many output items we produced.
      return noutput_items;
    } // end of detect_impl:: general_work
    
    int
    detect_impl::event(const float *input, float *output)
    {
      //outbuf = (float *) //create fresh one if necessary
      float n_sigma = d_dms; // translate variables 
      int mode = d_nt;
      int imax = 0;
      int vlen = d_vec_length;
      double sum2 = 0, rms = 0, rp = 0, ip = 0, maxmag2 = 0, mag2 = 0, npeak=0.;
      
      rp = input[0];
      output[0] = input[0];
      sum2 = mag2 = (rp*rp);
      maxmag2 = mag2;
      for(unsigned int j=1; j < vlen; j++)
	{ rp = input[j];
	  output[j] = rp;
	  mag2 = (rp*rp);
	  sum2 += mag2;
	  if (mag2 > maxmag2) {
	    maxmag2 = mag2;
	    imax = j;
	  }
	} // end for all samples
      rms = sqrt(sum2);
      // if rms is non-zero
      if (rms > 0.)
	{ // determine peak magnitude
	  npeak = sqrt(maxmag2)/rms;
	  if (npeak > n_sigma)
	    { 
	      printf( "N-sigma Peak found %7.1f\n", npeak);
	    }
	}// end if signficant peak found
      
      return 0;
    } // end of detect_impl::event()

  } /* namespace radio_astro */
} /* namespace gr */


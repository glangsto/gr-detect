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

    void 
    detect_impl::set_dms ( float dms)
    {
      nsigma = dms;
      printf("Input N Sigma: %7.1f\n", nsigma);
      d_dms = dms;
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
      //      if (success == 0)
      //	{
      //	  std::cout << "Processed " << ninputs << " Vectors" << "\n";
      //	}
      //std::cout << out[0*d_dms+ 0] << " " << out[31*d_dms+49] <<"\n";

      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each (noutput_items);

      // Tell runtime system how many output items we produced.
      return noutput_items;
    } // end of detect_impl:: general_work
    

    int
    detect_impl::update_buffer()
    { long vlen = d_vec_length, vlen2 = d_vec_length/2, i = inext2;

      i -= vlen2;  // center peak in middle of output vector

      // if event is well within the circular buffer 
      if ((i > vlen2) && (i < MAX_BUFF - vlen2))
	{ for (long j = 0; j < vlen; j++)
	    { samples[j] = circular[i];
	      i++;
	    }
	}

      // now must reset the buffer to wait for the next event
      bufferfull = 0;
      return 0;
    } // end of update_buffer()
    
    int
    detect_impl::event(const float *input, float *output)
    {
      //outbuf = (float *) //create fresh one if necessary
      float n_sigma = d_dms; // translate variables 
      int mode = d_nt;
      int vlen = d_vec_length;
      double rp = 0, mag2 = 0;
      
      // fill the circular buffer
      for(unsigned int j=0; j < vlen; j++)
	{ rp = input[j];
	  circular[inext] = rp;
	  circular2[inext] = rp*rp;
	  mag2 = rp*rp;
	  sum2 += mag2;
	  inext++;
	  if (inext >= MAX_BUFF) // if buffer is full
	    { rms2 = sum2*oneovern;
	      rms = sqrt(rms2);
	      inext = 0;
	      bufferfull = 1;    // flag buffer is now full
	      nsigma_rms = nsigma*nsigma*rms2;
	      sum2 = 0;          // restart rms sum
	    }
	  inext2++;              // update position for search 
	  if (inext2 >= MAX_BUFF) // if at end of circular buffer
	    inext2 = 0;           // go back to beginning
	  if (bufferfull)         // when buffer is full, find peaks
	    {
	      if (circular2[inext2] > nsigma_rms)
		{
		  imax2 = inext2;
		  peak = sqrt(circular2[inext2]);
		  printf( "N-sigma Peak found: %7.1f\n", peak/rms);
		  update_buffer();
		}
	    } // end if buffere full
	} // end for all samples
	      
      if (! initialized) {
	nsigma = d_dms;
	for (int iii = 0; iii < vlen; iii++)
	  { samples[iii] = input[iii];
	  }
	initialized = 1;     // no need to re-initialize the event
	printf("Input N-Sigma: %7.1f\n", nsigma);
      }

      // always output the last event
      for (int iii = 0; iii < vlen; iii++)
	{ output[iii] = samples[iii];
	}

      return 0;
    } // end of detect_impl::event()

  } /* namespace radio_astro */
} /* namespace gr */


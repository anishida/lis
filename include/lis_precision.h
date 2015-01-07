/* Copyright (C) 2005 The Scalable Software Infrastructure Project. All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:
   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
   3. Neither the name of the project nor the names of its contributors 
      may be used to endorse or promote products derived from this software 
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE SCALABLE SOFTWARE INFRASTRUCTURE PROJECT
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE SCALABLE SOFTWARE INFRASTRUCTURE
   PROJECT BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
   OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/


#ifndef __LIS_PRECISION_H__
#define __LIS_PRECISION_H__

#ifdef USE_SSE2
	#include <emmintrin.h>
#endif

#define SPLITTER 134217729.0

#ifndef USE_QUAD_PRECISION
	#define LIS_QUAD_DECLAR			
#else
#ifdef USE_SSE2
	#define LIS_QUAD_DECLAR __m128d bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh
#else
	#define LIS_QUAD_DECLAR double p1,p2,tq,bhi,blo,chi,clo,sh,eh,sl,el,th,tl
#endif
#endif

/**********************************************************
 * LIS_QUAD_MINUS(a_hi,a_lo)                              *
 **********************************************************
  (a_hi,a_lo) <- (-a_hi,-a_lo)
 **********************************************************/

#define LIS_QUAD_MINUS(a_hi,a_lo)   \
				(a_hi) = -(a_hi); \
				(a_lo) = -(a_lo)

/**********************************************************
 * LIS_QUAD_ZERO(a_hi,a_lo)                               *
 **********************************************************
  (a_hi,a_lo) <- (0,0)
 **********************************************************/

#define LIS_QUAD_ZERO(a_hi,a_lo)   \
				(a_hi) = 0.0; \
				(a_lo) = 0.0

/**********************************************************
 * LIS_QUAD_ONE(a_hi,a_lo)                                *
 **********************************************************
  (a_hi,a_lo) <- (1,0)
 **********************************************************/

#define LIS_QUAD_ONE(a_hi,a_lo)   \
				(a_hi) = 1.0; \
				(a_lo) = 0.0

/**********************************************************
 * LIS_QUAD_FAST_TWO_SUM(a,b,r,e)                         *
 **********************************************************
  |a| > |b|
  a + b -> (r,e)
 **********************************************************/

#define LIS_QUAD_FAST_TWO_SUM(a,b,r,e)   \
				(r) = (a) + (b); \
				(e) = (b)  - ((r) - (a))

/**********************************************************
 * LIS_QUAD_TWO_SUM(a,b,r,e)                              *
 **********************************************************
  a + b -> (r,e)
 **********************************************************/

#define LIS_QUAD_TWO_SUM(a,b,r,e)   \
				(r) = (a) + (b); \
				th  = (r) - (a); \
				(e) = ((a) - ((r) - th)) + ((b) - th)

/**********************************************************
 * LIS_QUAD_TWO_DIFF(a,b,r,e)                             *
 **********************************************************
  a - b -> (r,e)
 **********************************************************/

#define LIS_QUAD_TWO_DIFF(a,b,r,e)   \
				(r) = (a) - (b); \
				th  = (r) - (a); \
				(e) = ((a) - ((r) - th)) - ((b) + th)

/**********************************************************
 * LIS_QUAD_SPLIT(b,b_hi,b_lo)                            *
 **********************************************************
  b -> (b_hi,b_lo)
 **********************************************************/

#define LIS_QUAD_SPLIT(b,b_hi,b_lo)   \
				tq     = SPLITTER * (b); \
				(b_hi) = tq - (tq-(b));  \
				(b_lo) = (b) - (b_hi);   \

/**********************************************************
 * LIS_QUAD_TWO_PROD(a,b,r,e)                             *
 **********************************************************
  a x b -> (r,e)
 **********************************************************/

#ifndef HAS_FMA
#define LIS_QUAD_TWO_PROD(a,b,r,e)   \
				(r) = (a) * (b); \
				LIS_QUAD_SPLIT((a),bhi,blo); \
				LIS_QUAD_SPLIT((b),chi,clo); \
				(e) = ((bhi*chi-(r))+bhi*clo+blo*chi)+blo*clo
#else
#define LIS_QUAD_TWO_PROD(a,b,r,e)   \
				(r) = (-(a)) * (b); \
				(e) = (a) * (b) + (r)
#endif

/**********************************************************
 * LIS_QUAD_TWO_SQR(a,r,e)                                *
 **********************************************************
  a x a -> (r,e)
 **********************************************************/

#ifndef HAS_FMA
#define LIS_QUAD_TWO_SQR(a,r,e)   \
				(r) = (a) * (a); \
				LIS_QUAD_SPLIT((a),bhi,blo); \
				(e) = ((bhi*bhi-(r))+2.0*bhi*blo)+blo*blo
#else
#define LIS_QUAD_TWO_SQR(a,r,e)   \
				(r) = (-(a)) * (a); \
				(e) = (a) * (a) + (r)
#endif

/**********************************************************
 * LIS_QUAD_MUL(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo)            *
 **********************************************************
  (a_hi,a_lo) <- (b_hi,b_lo) x (c_hi,c_lo)
 **********************************************************/

#define LIS_QUAD_MUL(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) \
				LIS_QUAD_TWO_PROD((b_hi),(c_hi),p1,p2); \
				p2 += ((b_hi) * (c_lo)); \
				p2 += ((b_lo) * (c_hi)); \
				LIS_QUAD_FAST_TWO_SUM(p1,p2,(a_hi),(a_lo))

/**********************************************************
 * LIS_QUAD_MULD(a_hi,a_lo,b_hi,b_lo,c)                   *
 **********************************************************
  (a_hi,a_lo) <- (b_hi,b_lo) x c
 **********************************************************/

#define LIS_QUAD_MULD(a_hi,a_lo,b_hi,b_lo,c) \
				LIS_QUAD_TWO_PROD((b_hi),(c),p1,p2); \
				p2 += ((b_lo) * (c)); \
				LIS_QUAD_FAST_TWO_SUM(p1,p2,(a_hi),(a_lo))

/**********************************************************
 * LIS_QUAD_SQR(a_hi,a_lo,b_hi,b_lo)                      *
 **********************************************************
  (a_hi,a_lo) <- (b_hi,b_lo) x (b_hi,b_lo)
 **********************************************************/

#define LIS_QUAD_SQR(a_hi,a_lo,b_hi,b_lo) \
				LIS_QUAD_TWO_SQR((b_hi),p1,p2); \
				p2 += (2.0*(b_hi) * (b_lo)); \
				p2 += ((b_lo) * (b_lo)); \
				LIS_QUAD_FAST_TWO_SUM(p1,p2,(a_hi),(a_lo))

/**********************************************************
 * LIS_QUAD_ADD(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo)            *
 **********************************************************
  (a_hi,a_lo) <- (b_hi,b_lo) + (c_hi,c_lo)
 **********************************************************/

#ifndef USE_FAST_QUAD_ADD
#define LIS_QUAD_ADD(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) \
				LIS_QUAD_TWO_SUM((b_hi),(c_hi),sh,eh); \
				LIS_QUAD_TWO_SUM((b_lo),(c_lo),sl,el); \
				eh += sl; \
				LIS_QUAD_FAST_TWO_SUM(sh,eh,sh,eh); \
				eh += el; \
				LIS_QUAD_FAST_TWO_SUM(sh,eh,(a_hi),(a_lo))
#else
#define LIS_QUAD_ADD(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) \
				LIS_QUAD_TWO_SUM((b_hi),(c_hi),sh,eh); \
				eh += (b_lo); \
				eh += (c_lo); \
				LIS_QUAD_FAST_TWO_SUM(sh,eh,(a_hi),(a_lo))
#endif

/**********************************************************
 * LIS_QUAD_DIV(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo)            *
 **********************************************************
  (a_hi,a_lo) <- (b_hi,b_lo) / (c_hi,c_lo)
 **********************************************************/

#define LIS_QUAD_DIV(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) \
				tl  = (b_hi) / (c_hi); \
				LIS_QUAD_MULD(eh,el,(c_hi),(c_lo),tl); \
				LIS_QUAD_TWO_DIFF((b_hi),eh,sh,sl); \
				sl -= el; \
				sl += (b_lo); \
				th  = (sh+sl) / (c_hi); \
				LIS_QUAD_FAST_TWO_SUM(tl,th,(a_hi),(a_lo))

/**********************************************************
 * LIS_QUAD_SQRT(a_hi,a_lo,b_hi,b_lo)                     *
 **********************************************************
  (a_hi,a_lo) <- SQRT( (b_hi,b_lo) )
 **********************************************************/

#define LIS_QUAD_SQRT(a_hi,a_lo,b_hi,b_lo) \
				if( (b_hi)==0 ) \
				{ \
					(a_hi) = (a_lo) = 0.0; \
					return LIS_FAILS; \
				} \
				if( (b_hi)<0 ) \
				{ \
					printf("ERROR bh=%e\n",(b_hi)); \
					(a_hi) = (a_lo) = 0.0; \
					return LIS_FAILS; \
				} \
				p1 = 1.0 / sqrt((b_hi)); \
				p2 = (b_hi) * p1; \
				p1 = p1 * 0.5; \
				LIS_QUAD_TWO_SQR(p2,chi,clo); \
				LIS_QUAD_ADD(th,eh,(b_hi),(b_lo),-chi,-clo); \
				p1 = p1 * th; \
				LIS_QUAD_FAST_TWO_SUM(p1,p2,(a_hi),(a_lo))

/**********************************************************
 * LIS_QUAD_FMA(d_hi,d_lo,a_hi,a_lo,b_hi,b_lo,c_hi,c_lo)  *
 **********************************************************
  (d_hi,d_lo) <- (a_hi,a_lo) + (b_hi,b_lo) * (c_hi,c_lo)
 **********************************************************/

#define LIS_QUAD_FMA(d_hi,d_lo,a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) \
				LIS_QUAD_MUL(chi,p2,(b_hi),(b_lo),(c_hi),(c_lo)); \
				LIS_QUAD_ADD((d_hi),(d_lo),(a_hi),(a_lo),chi,p2)

/**********************************************************
 * LIS_QUAD_FMAD(d_hi,d_lo,a_hi,a_lo,b_hi,b_lo,c)         *
 **********************************************************
  (d_hi,d_lo) <- (a_hi,a_lo) + (b_hi,b_lo) * c
 **********************************************************/

#define LIS_QUAD_FMAD(d_hi,d_lo,a_hi,a_lo,b_hi,b_lo,c) \
				LIS_QUAD_MULD(chi,p2,(b_hi),(b_lo),(c)); \
				LIS_QUAD_ADD((d_hi),(d_lo),(a_hi),(a_lo),chi,p2)

/**********************************************************
 * LIS_QUAD_FSA(d_hi,d_lo,a_hi,a_lo,b_hi,b_lo)            *
 **********************************************************
  (d_hi,d_lo) <- (a_hi,a_lo) + (b_hi,b_lo) * (b_hi,b_lo) 
 **********************************************************/

#define LIS_QUAD_FSA(d_hi,d_lo,a_hi,a_lo,b_hi,b_lo) \
				LIS_QUAD_SQR(bhi,p2,(b_hi),(b_lo)); \
				LIS_QUAD_ADD((d_hi),(d_lo),(a_hi),(a_lo),bhi,p2)

/**********************************************************
 *                                                        *
 *                      SSE2(SD)                          *
 *                                                        *
 **********************************************************/

/**********************************************************
 * LIS_QUAD_MUL_SSE2_STORE(a_hi,a_lo)                     *
 **********************************************************
  (a_hi,a_lo) <- 
 **********************************************************/

#define LIS_QUAD_MUL_SSE2_STORE(a_hi,a_lo) \
				_mm_store_sd(&(a_hi),sh); \
				sh = _mm_sub_sd(sh,bh); \
				wh = _mm_sub_sd(wh,sh); \
				_mm_store_sd(&(a_lo),wh)

/**********************************************************
 * LIS_QUAD_ADD_SSE2_LOAD(b_hi,b_lo,c_hi,c_lo)            *
 **********************************************************
  (b_hi,b_lo)  (c_hi,c_lo)
 **********************************************************/

#ifndef USE_FAST_QUAD_ADD
#define LIS_QUAD_ADD_SSE2_LOAD(b_hi,b_lo,c_hi,c_lo) \
				sh = _mm_set_pd((b_lo),(b_hi)); \
				bl = _mm_set_pd((c_lo),(c_hi))
#else
#define LIS_QUAD_ADD_SSE2_LOAD(b_hi,b_lo,c_hi,c_lo) \
 				eh = _mm_set_sd((b_hi)); \
				cl = _mm_set_sd((b_lo)); \
				sh = _mm_set_sd((c_hi)); \
				wh = _mm_set_sd((c_lo))
#endif

/**********************************************************
 * LIS_QUAD_FMA_SSE2_LOAD(a_hi,a_lo)                      *
 **********************************************************
  (a_hi,a_lo)
 **********************************************************/

#ifndef USE_FAST_QUAD_ADD
#define LIS_QUAD_FMA_SSE2_LOAD(a_hi,a_lo) \
				th = _mm_sub_sd(sh,bh); \
				wh = _mm_sub_sd(wh,th); \
				bl = _mm_set_pd((a_lo),(a_hi)); \
				sh = _mm_unpacklo_pd(sh,wh)
#else
#define LIS_QUAD_FMA_SSE2_LOAD(a_hi,a_lo) \
				th = _mm_sub_sd(sh,bh); \
				wh = _mm_sub_sd(wh,th); \
 				eh = _mm_set_sd((a_hi)); \
				cl = _mm_set_sd((a_lo))
#endif

/**********************************************************
 * LIS_QUAD_MUL_SSE2_CORE(b_hi,b_lo,c_hi,c_lo)            *
 **********************************************************
  (b_hi,b_lo) x (c_hi,c_lo)
 **********************************************************/

#define LIS_QUAD_MUL_SSE2_CORE(b_hi,b_lo,c_hi,c_lo) \
				sh = _mm_set_pd(SPLITTER,SPLITTER); \
				ch = _mm_set_pd((c_hi),(b_hi)); \
				bh = _mm_unpackhi_pd(ch,ch); \
				tl = _mm_set_pd((b_lo),(c_lo)); \
				bh = _mm_mul_sd(bh,ch); \
				sh = _mm_mul_pd(sh,ch); \
				th = _mm_sub_pd(sh,ch); \
				sh = _mm_sub_pd(sh,th); \
				th = _mm_mul_pd(ch,tl); \
				ch = _mm_sub_pd(ch,sh); \
				eh = _mm_unpackhi_pd(sh,ch); \
				wh = _mm_unpacklo_pd(sh,ch); \
				ch = _mm_unpackhi_pd(ch,wh);  \
				wh = _mm_mul_pd(wh,eh); \
				sh = _mm_mul_pd(sh,ch); \
				ch = _mm_unpackhi_pd(sh,sh); \
				tl = _mm_unpackhi_pd(wh,wh); \
				eh = _mm_unpackhi_pd(th,th); \
				wh = _mm_sub_sd(wh,bh); \
				wh = _mm_add_sd(wh,ch); \
				wh = _mm_add_sd(wh,sh); \
				wh = _mm_add_sd(wh,tl); \
				wh = _mm_add_sd(wh,eh); \
				wh = _mm_add_sd(wh,th); \
				sh = _mm_add_sd(wh,bh)

/**********************************************************
 * LIS_QUAD_MULD_SSE2_CORE(b_hi,b_lo,c)                   *
 **********************************************************
  (b_hi,b_lo) x c
 **********************************************************/

#define LIS_QUAD_MULD_SSE2_CORE(b_hi,b_lo,c) \
				sh = _mm_set_pd(SPLITTER,SPLITTER); \
				ch = _mm_set_pd((b_hi),(c)); \
				bh = _mm_unpackhi_pd(ch,ch); \
				sl = _mm_load_sd(&(b_lo)); \
				bh = _mm_mul_sd(bh,ch); \
				sl = _mm_mul_sd(sl,ch); \
				sh = _mm_mul_pd(sh,ch); \
				th = _mm_sub_pd(sh,ch); \
				sh = _mm_sub_pd(sh,th); \
				ch = _mm_sub_pd(ch,sh); \
				t1 = _mm_unpackhi_pd(sh,ch); \
				wh = _mm_unpacklo_pd(sh,ch); \
				ch = _mm_unpackhi_pd(ch,wh);  \
				wh = _mm_mul_pd(wh,t1); \
				sh = _mm_mul_pd(sh,ch); \
				ch = _mm_unpackhi_pd(sh,sh); \
				th = _mm_unpackhi_pd(wh,wh); \
				wh = _mm_sub_sd(wh,bh); \
				wh = _mm_add_sd(wh,ch); \
				wh = _mm_add_sd(wh,sh); \
				wh = _mm_add_sd(wh,th); \
				wh = _mm_add_sd(wh,sl); \
				sh = _mm_add_sd(wh,bh)

/**********************************************************
 * LIS_QUAD_SQR_SSE2_CORE(b_hi,b_lo)                      *
 **********************************************************
  (b_hi,b_lo) x (b_hi,b_lo)
 **********************************************************/

#define LIS_QUAD_SQR_SSE2_CORE(b_hi,b_lo) \
				sh = _mm_set_sd(SPLITTER); \
				th = _mm_load_sd(&(b_hi)); \
				bl = _mm_load_sd(&(b_lo)); \
				bh = _mm_mul_sd(th,th); \
				sh = _mm_mul_sd(sh,th); \
				wh = _mm_sub_sd(sh,th); \
				sh = _mm_sub_sd(sh,wh); \
				cl = _mm_sub_sd(th,sh); \
				wh = _mm_mul_sd(sh,sh); \
				sh = _mm_add_sd(sh,sh); \
				sh = _mm_mul_sd(sh,cl); \
				cl = _mm_mul_sd(cl,cl); \
				wh = _mm_sub_sd(wh,bh); \
				wh = _mm_add_sd(wh,sh); \
				wh = _mm_add_sd(wh,cl); \
				th = _mm_add_sd(th,th); \
				th = _mm_mul_sd(th,bl); \
				bl = _mm_mul_sd(bl,bl); \
				wh = _mm_add_sd(wh,th); \
				wh = _mm_add_sd(wh,bl); \
				sh = _mm_add_sd(wh,bh)

/**********************************************************
 * LIS_QUAD_ADD_SSE2_CORE(a_hi,a_lo)                      *
 **********************************************************
  (a_hi,a_lo) <- (b_hi,b_lo) + (c_hi,c_lo)
 **********************************************************/

#ifndef USE_FAST_QUAD_ADD
#define LIS_QUAD_ADD_SSE2_CORE(a_hi,a_lo) \
				t0 = _mm_add_pd(bl,sh); \
				eh = _mm_sub_pd(t0,bl); \
				th = _mm_sub_pd(t0,eh); \
				sh = _mm_sub_pd(sh,eh); \
				bl = _mm_sub_pd(bl,th); \
				bl = _mm_add_pd(bl,sh); \
				eh = _mm_unpackhi_pd(bl,bl); \
				sh = _mm_unpackhi_pd(t0,t0); \
				th = t0; \
				bl = _mm_add_sd(bl,sh); \
				th = _mm_add_sd(th,bl); \
				t0 = _mm_sub_sd(th,t0); \
				bl = _mm_sub_sd(bl,t0); \
				bl = _mm_add_sd(bl,eh); \
				sh = _mm_add_sd(th,bl); \
				_mm_store_sd(&(a_hi),sh); \
				sh = _mm_sub_sd(sh,th); \
				bl = _mm_sub_sd(bl,sh); \
				_mm_store_sd(&(a_lo),bl)
#else
#define LIS_QUAD_ADD_SSE2_CORE(a_hi,a_lo) \
				bl = _mm_add_sd(eh,sh); \
				th = _mm_sub_sd(bl,eh); \
				t0 = _mm_sub_sd(bl,th); \
				sh = _mm_sub_sd(sh,th); \
				eh = _mm_sub_sd(eh,t0); \
				eh = _mm_add_sd(eh,sh); \
				eh = _mm_add_sd(eh,cl); \
				eh = _mm_add_sd(eh,wh); \
				th = _mm_add_sd(bl,eh); \
				_mm_store_sd(&(a_hi),th); \
				th = _mm_sub_sd(th,bl); \
				eh = _mm_sub_sd(eh,th); \
				_mm_store_sd(&(a_lo),eh)
#endif

/**********************************************************
 * LIS_QUAD_MUL_SSE2(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo)       *
 **********************************************************
  (a_hi,a_lo) <- (b_hi,b_lo) x (c_hi,c_lo)
 **********************************************************/

#define LIS_QUAD_MUL_SSE2(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) \
				LIS_QUAD_MUL_SSE2_CORE((b_hi),(b_lo),(c_hi),(c_lo)); \
				LIS_QUAD_MUL_SSE2_STORE((a_hi),(a_lo))

/**********************************************************
 * LIS_QUAD_MULD_SSE2(a_hi,a_lo,b_hi,b_lo,c)              *
 **********************************************************
  (a_hi,a_lo) <- (b_hi,b_lo) x c
 **********************************************************/

#define LIS_QUAD_MULD_SSE2(a_hi,a_lo,b_hi,b_lo,c) \
				LIS_QUAD_MULD_SSE2_CORE((b_hi),(b_lo),(c)); \
				LIS_QUAD_MUL_SSE2_STORE((a_hi),(a_lo))

/**********************************************************
 * LIS_QUAD_SQR_SSE2(a_hi,a_lo,b_hi,b_lo)                 *
 **********************************************************
  (a_hi,a_lo) <- (b_hi,b_lo) x (b_hi,b_lo)
 **********************************************************/

#define LIS_QUAD_SQR_SSE2(a_hi,a_lo,b_hi,b_lo) \
				LIS_QUAD_SQR_SSE2_CORE((b_hi),(b_lo)); \
				LIS_QUAD_MUL_SSE2_STORE((a_hi),(a_lo))

/**********************************************************
 * LIS_QUAD_ADD_SSE2(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo)       *
 **********************************************************
  (a_hi,a_lo) <- (b_hi,b_lo) + (c_hi,c_lo)
 **********************************************************/

#define LIS_QUAD_ADD_SSE2(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) \
				LIS_QUAD_ADD_SSE2_LOAD((b_hi),(b_lo),(c_hi),(c_lo)); \
				LIS_QUAD_ADD_SSE2_CORE((a_hi),(a_lo))

/**********************************************************
 * LIS_QUAD_DIV_SSE2(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo)       *
 **********************************************************
  (a_hi,a_lo) <- (b_hi,b_lo) / (c_hi,c_lo)
 **********************************************************/

#define LIS_QUAD_DIV_SSE2(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) \
				sh = _mm_set_pd(SPLITTER,SPLITTER); \
				bh = _mm_set_pd((b_lo),(b_hi)); \
				ch = _mm_set_pd((c_lo),(c_hi)); \
				p2 = bh; \
				wh = ch; \
				p2 = _mm_div_sd(p2,ch); \
				wh = _mm_unpacklo_pd(wh,p2); \
				ch = _mm_move_sd(ch,p2); \
				p2 = wh; \
				sh = _mm_mul_pd(sh,wh); \
				th = _mm_sub_pd(sh,wh); \
				sh = _mm_sub_pd(sh,th); \
				ch = _mm_mul_pd(ch,wh); \
				wh = _mm_sub_pd(wh,sh); \
				th = sh; \
				p1 = wh; \
				th = _mm_unpackhi_pd(th,th); \
				p1 = _mm_unpackhi_pd(p1,p1); \
				eh = th; \
				th = _mm_mul_sd(th,wh); \
				eh = _mm_mul_sd(eh,sh); \
				wh = _mm_mul_sd(wh,p1); \
				sh = _mm_mul_sd(sh,p1); \
				p1 = ch; \
				p1 = _mm_unpackhi_pd(p1,p1); \
				eh = _mm_sub_sd(eh,ch); \
				eh = _mm_add_sd(eh,sh); \
				eh = _mm_add_sd(eh,th); \
				eh = _mm_add_sd(eh,wh); \
				eh = _mm_add_sd(eh,p1); \
				p1 = ch; \
				p1 = _mm_add_sd(p1,eh); \
				wh = p1; \
				wh = _mm_sub_sd(wh,ch); \
				sh = bh; \
				sh = _mm_sub_sd(sh,p1); \
				eh = _mm_sub_sd(eh,wh); \
				th = sh; \
				th = _mm_sub_sd(th,bh); \
				p1 = _mm_add_sd(p1,th); \
				wh = sh; \
				wh = _mm_sub_sd(wh,th); \
				th = bh; \
				th = _mm_unpackhi_pd(th,th); \
				bh = _mm_sub_sd(bh,wh); \
				bh = _mm_sub_sd(bh,p1); \
				bh = _mm_sub_sd(bh,eh); \
				bh = _mm_add_sd(bh,th); \
				th = p2; \
				th = _mm_unpackhi_pd(th,th); \
				eh = th; \
				sh = _mm_add_sd(sh,bh); \
				sh = _mm_div_sd(sh,p2); \
				th = _mm_add_sd(th,sh); \
				_mm_store_sd(&(a_hi),th); \
				th = _mm_sub_sd(th,eh); \
				sh = _mm_sub_sd(sh,th); \
				_mm_store_sd(&(a_lo),sh)

/**********************************************************
 * LIS_QUAD_SQRT_SSE2(a_hi,a_lo,b_hi,b_lo)                *
 **********************************************************
  (a_hi,a_lo) <- SQRT( (b_hi,b_lo) )
 **********************************************************/

#define LIS_QUAD_SQRT_SSE2(a_hi,a_lo,b_hi,b_lo) \
				if( (b_hi)==0 ) \
				{ \
					(a_hi) = (a_lo) = 0.0; \
					return LIS_FAILS; \
				} \
				if( (b_hi)<0 ) \
				{ \
					printf("ERROR bh=%e\n",(b_hi)); \
					(a_hi) = (a_lo) = 0.0; \
					return LIS_FAILS; \
				} \
				wh = _mm_set_sd(SPLITTER); \
				bh = _mm_load_pd(&(b_hi)); \
				bh = bh; \
				t0 = _mm_set_sd(1.0); \
				t1 = _mm_set_sd(0.5); \
				t2 = _mm_sqrt_pd(bh); \
				t0 = _mm_div_sd(t0,t2); \
				t2 = _mm_mul_sd(bh,t0); \
				t0 = _mm_mul_sd(t0,t1); \
				p1 = _mm_mul_sd(t2,t2); \
				wh = _mm_mul_sd(wh,t2); \
				t1 = _mm_sub_sd(wh,t2); \
				wh = _mm_sub_sd(wh,t1); \
				wl = _mm_sub_sd(t2,wh); \
				t1 = _mm_mul_sd(wh,wh); \
				wh = _mm_add_sd(wh,wh); \
				wh = _mm_mul_sd(wh,wl); \
				wl = _mm_mul_sd(wl,wl); \
				t1 = _mm_sub_sd(t1,p1); \
				t1 = _mm_add_sd(t1,wh); \
				t1 = _mm_add_sd(t1,wl); \
				p1 = _mm_unpacklo_pd(p1,t1); \
				sh = _mm_sub_pd(bh,p1); \
				eh = _mm_sub_pd(sh,bh); \
				th = _mm_sub_pd(sh,eh); \
				p1 = _mm_add_pd(p1,eh); \
				bh = _mm_sub_pd(bh,th); \
				bh = _mm_sub_pd(bh,p1); \
				eh = _mm_unpackhi_pd(bh,bh); \
				ch = _mm_unpackhi_pd(sh,sh); \
				th = sh; \
				bh = _mm_add_sd(bh,ch); \
				th = _mm_add_sd(th,bh); \
				sh = _mm_sub_sd(th,sh); \
				bh = _mm_sub_sd(bh,sh); \
				bh = _mm_add_sd(bh,eh); \
				th = _mm_add_sd(th,bh); \
				t0 = _mm_mul_sd(t0,th); \
				p1 = _mm_add_sd(t2,t0); \
				_mm_store_sd(&(a_hi),p1); \
				p1 = _mm_sub_sd(p1,t2); \
				t0 = _mm_sub_sd(t0,p1); \
				_mm_store_sd(&(a_lo),t0)

/***************************************************************
 * LIS_QUAD_FMA_SSE2(d_hi,d_lo,a_hi,a_lo,b_hi,b_lo,c_hi,c_lo)  *
 ***************************************************************
  (d_hi,d_lo) <- (a_hi,a_lo) + (b_hi,b_lo) * (c_hi,c_lo)
 ***************************************************************/

#define LIS_QUAD_FMA_SSE2(d_hi,d_lo,a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) \
				LIS_QUAD_MUL_SSE2_CORE((b_hi),(b_lo),(c_hi),(c_lo)); \
				LIS_QUAD_FMA_SSE2_LOAD((a_hi),(a_lo)); \
				LIS_QUAD_ADD_SSE2_CORE((d_hi),(d_lo))

/***************************************************************
 * LIS_QUAD_FMAD_SSE2(d_hi,d_lo,a_hi,a_lo,b_hi,b_lo,c)         *
 ***************************************************************
  (d_hi,d_lo) <- (a_hi,a_lo) + (b_hi,b_lo) * c
 ***************************************************************/

#define LIS_QUAD_FMAD_SSE2(d_hi,d_lo,a_hi,a_lo,b_hi,b_lo,c) \
				LIS_QUAD_MULD_SSE2_CORE((b_hi),(b_lo),(c)); \
				LIS_QUAD_FMA_SSE2_LOAD((a_hi),(a_lo)); \
				LIS_QUAD_ADD_SSE2_CORE((d_hi),(d_lo))

/***************************************************************
 * LIS_QUAD_FSA_SSE2(d_hi,d_lo,a_hi,a_lo,b_hi,b_lo)            *
 ***************************************************************
  (d_hi,d_lo) <- (a_hi,a_lo) + (b_hi,b_lo) * (b_hi,b_lo) 
 ***************************************************************/

#define LIS_QUAD_FSA_SSE2(d_hi,d_lo,a_hi,a_lo,b_hi,b_lo) \
				LIS_QUAD_SQR_SSE2_CORE((b_hi),(b_lo)); \
				LIS_QUAD_FMA_SSE2_LOAD((a_hi),(a_lo)); \
				LIS_QUAD_ADD_SSE2_CORE((d_hi),(d_lo))

/**********************************************************
 *                                                        *
 *                      SSE2(PD)                          *
 *                                                        *
 **********************************************************/

/**********************************************************
 * LIS_QUAD_MUL2_SSE2_LOAD(b_hi,b_lo,c_hi,c_lo)           *
 **********************************************************
  (b_hi,b_lo)  (c_hi,c_lo)
 **********************************************************/

#define LIS_QUAD_MUL2_SSE2_LOAD(b_hi,b_lo,c_hi,c_lo) \
				t0 = _mm_set_pd(SPLITTER,SPLITTER); \
				bh = _mm_loadu_pd(&(b_hi)); \
				ch = _mm_loadu_pd(&(c_hi)); \
				bl = _mm_loadu_pd(&(b_lo)); \
				cl = _mm_loadu_pd(&(c_lo))

/**********************************************************
 * LIS_QUAD_MULD2_SSE2_LOAD(b_hi,b_lo,c)                  *
 **********************************************************
  (b_hi,b_lo)  (c)
 **********************************************************/

#define LIS_QUAD_MULD2_SSE2_LOAD(b_hi,b_lo,c) \
				t0 = _mm_set_pd(SPLITTER,SPLITTER); \
				bh = _mm_loadu_pd(&(b_hi)); \
				bl = _mm_loadu_pd(&(b_lo)); \
				ch = _mm_loadu_pd(&(c))
#define LIS_QUAD_MULD2_SSE2_LOAD_SD(b0_hi,b0_lo,b1_hi,b1_lo,c) \
				t0 = _mm_set_pd(SPLITTER,SPLITTER); \
				bh = _mm_set_pd((b1_hi),(b0_hi)); \
				bl = _mm_set_pd((b1_lo),(b0_lo)); \
				ch = _mm_loadu_pd(&(c))

/**********************************************************
 * LIS_QUAD_SQR2_SSE2_LOAD(b_hi,b_lo)                     *
 **********************************************************
  (b_hi,b_lo)
 **********************************************************/

#define LIS_QUAD_SQR2_SSE2_LOAD(b_hi,b_lo) \
				ch = _mm_set_pd(SPLITTER,SPLITTER); \
				bh = _mm_loadu_pd(&(b_hi)); \
				bl = _mm_loadu_pd(&(b_lo))

/**********************************************************
 * LIS_QUAD_MUL2_SSE2_STORE(a_hi,a_lo)                    *
 **********************************************************
  (a_hi,a_lo) <- 
 **********************************************************/

#define LIS_QUAD_MUL2_SSE2_STORE(a_hi,a_lo) \
				_mm_storeu_pd(&(a_hi),ch); \
				ch = _mm_sub_pd(ch,p1); \
				p2 = _mm_sub_pd(p2,ch); \
				_mm_storeu_pd(&(a_lo),p2)
#define LIS_QUAD_MUL2_SSE2_STOREU(a_hi,a_lo) \
				_mm_storeu_pd(&(a_hi),ch); \
				ch = _mm_sub_pd(ch,p1); \
				p2 = _mm_sub_pd(p2,ch); \
				_mm_storeu_pd(&(a_lo),p2)

/**********************************************************
 * LIS_QUAD_ADD2_SSE2_LOAD(b_hi,b_lo,c_hi,c_lo)           *
 **********************************************************
  (b_hi,b_lo)  (c_hi,c_lo)
 **********************************************************/

#ifndef USE_FAST_QUAD_ADD
#define LIS_QUAD_ADD2_SSE2_LOAD(b_hi,b_lo,c_hi,c_lo) \
				sh = _mm_set_pd((b_lo),(b_hi)); \
				bl = _mm_set_pd((c_lo),(c_hi))
#else
#define LIS_QUAD_ADD2_SSE2_LOAD(b_hi,b_lo,c_hi,c_lo) \
 				eh = _mm_set_sd((b_hi)); \
				cl = _mm_set_sd((b_lo)); \
				sh = _mm_set_sd((c_hi)); \
				wh = _mm_set_sd((c_lo))
#endif

/**********************************************************
 * LIS_QUAD_ADD2_SSE2_STORE(a_hi,a_lo)                    *
 **********************************************************
  (a_hi,a_lo) <- 
 **********************************************************/

#define LIS_QUAD_ADD2_SSE2_STORE(a_hi,a_lo) \
				_mm_storeu_pd(&(a_hi),sh); \
				sh = _mm_sub_pd(sh,th); \
				bh = _mm_sub_pd(bh,sh); \
				_mm_storeu_pd(&(a_lo),bh)
#define LIS_QUAD_ADD2_SSE2_STORE_SD(a0_hi,a0_lo,a1_hi,a1_lo) \
				_mm_storel_pd(&(a0_hi),sh); \
				_mm_storeh_pd(&(a1_hi),sh); \
				sh = _mm_sub_pd(sh,th); \
				bh = _mm_sub_pd(bh,sh); \
				_mm_storel_pd(&(a0_lo),bh); \
				_mm_storeh_pd(&(a1_lo),bh)

/**********************************************************
 * LIS_QUAD_FMA2_SSE2_LOAD(a_hi,a_lo)                     *
 **********************************************************
  (a_hi,a_lo)
 **********************************************************/

#define LIS_QUAD_FMA2_SSE2_LOAD(a_hi,a_lo) \
				t1 = _mm_sub_pd(ch,p1); \
				p2 = _mm_sub_pd(p2,t1); \
				bh = _mm_loadu_pd(&(a_hi)); \
				bl = _mm_loadu_pd(&(a_lo))
#define LIS_QUAD_FMA2_SSE2_LOAD_SD(a0_hi,a0_lo,a1_hi,a1_lo) \
				t1 = _mm_sub_pd(ch,p1); \
				p2 = _mm_sub_pd(p2,t1); \
				bh = _mm_set_pd((a1_hi),(a0_hi)); \
				bl = _mm_set_pd((a1_lo),(a0_lo))

/**********************************************************
 * LIS_QUAD_MUL2_SSE2_CORE                                *
 **********************************************************
  (b_hi,b_lo) x (c_hi,c_lo)
 **********************************************************/

#define LIS_QUAD_MUL2_SSE2_CORE \
				p1 = _mm_mul_pd(bh,ch); \
				sh = _mm_mul_pd(t0,bh); \
				sl = _mm_mul_pd(t0,ch); \
				th = _mm_sub_pd(sh,bh); \
				tl = _mm_sub_pd(sl,ch); \
				sh = _mm_sub_pd(sh,th); \
				sl = _mm_sub_pd(sl,tl); \
				t1 = _mm_mul_pd(bh,cl); \
				wh = _mm_sub_pd(bh,sh); \
				t2 = _mm_mul_pd(ch,bl); \
				wl = _mm_sub_pd(ch,sl); \
				t0 = _mm_mul_pd(wh,wl); \
				p2 = _mm_mul_pd(sh,sl); \
				sh = _mm_mul_pd(sh,wl); \
				sl = _mm_mul_pd(sl,wh); \
				p2 = _mm_sub_pd(p2,p1); \
				p2 = _mm_add_pd(p2,sh); \
				p2 = _mm_add_pd(p2,sl); \
				p2 = _mm_add_pd(p2,t0); \
				p2 = _mm_add_pd(p2,t1); \
				p2 = _mm_add_pd(p2,t2); \
				ch = _mm_add_pd(p1,p2)

/**********************************************************
 * LIS_QUAD_MULD2_SSE2_CORE                               *
 **********************************************************
  (b_hi,b_lo) x c
 **********************************************************/

#define LIS_QUAD_MULD2_SSE2_CORE \
				p1 = _mm_mul_pd(bh,ch); \
				bl = _mm_mul_pd(bl,ch); \
				sh = _mm_mul_pd(t0,bh); \
				th = _mm_sub_pd(sh,bh); \
				sh = _mm_sub_pd(sh,th); \
				bh = _mm_sub_pd(bh,sh); \
				sl = _mm_mul_pd(t0,ch); \
				tl = _mm_sub_pd(sl,ch); \
				sl = _mm_sub_pd(sl,tl); \
				ch = _mm_sub_pd(ch,sl); \
				t2 = _mm_mul_pd(bh,ch); \
				p2 = _mm_mul_pd(sh,sl); \
				t0 = _mm_mul_pd(sh,ch); \
				t1 = _mm_mul_pd(sl,bh); \
				p2 = _mm_sub_pd(p2,p1); \
				p2 = _mm_add_pd(p2,t0); \
				p2 = _mm_add_pd(p2,t1); \
				p2 = _mm_add_pd(p2,t2); \
				p2 = _mm_add_pd(p2,bl); \
				ch = _mm_add_pd(p1,p2)

/**********************************************************
 * LIS_QUAD_SQR2_SSE2_CORE                                *
 **********************************************************
  (b_hi,b_lo) x (b_hi,b_lo)
 **********************************************************/

#define LIS_QUAD_SQR2_SSE2_CORE \
				p1 = _mm_mul_pd(bh,bh); \
				ch = _mm_mul_pd(ch,bh); \
				p2 = _mm_sub_pd(ch,bh); \
				ch = _mm_sub_pd(ch,p2); \
				cl = _mm_sub_pd(bh,ch); \
				p2 = _mm_mul_pd(ch,ch); \
				ch = _mm_add_pd(ch,ch); \
				ch = _mm_mul_pd(ch,cl); \
				cl = _mm_mul_pd(cl,cl); \
				p2 = _mm_sub_pd(p2,p1); \
				p2 = _mm_add_pd(p2,ch); \
				p2 = _mm_add_pd(p2,cl); \
				bh = _mm_add_pd(bh,bh); \
				bh = _mm_mul_pd(bh,bl); \
				bl = _mm_mul_pd(bl,bl); \
				p2 = _mm_add_pd(p2,bh); \
				p2 = _mm_add_pd(p2,bl); \
				ch = _mm_add_pd(p1,p2)

/**********************************************************
 * LIS_QUAD_ADD2_SSE2_CORE                                *
 **********************************************************
  (b_hi,b_lo) + (c_hi,c_lo)
 **********************************************************/

#ifndef USE_FAST_QUAD_ADD
#define LIS_QUAD_ADD2_SSE2_CORE \
				sh = _mm_add_pd(bh,ch); \
				th = _mm_sub_pd(sh,bh); \
				t0 = _mm_sub_pd(sh,th); \
				ch = _mm_sub_pd(ch,th); \
				bh = _mm_sub_pd(bh,t0); \
				bh = _mm_add_pd(bh,ch); \
				sl = _mm_add_pd(bl,p2); \
				th = _mm_sub_pd(sl,bl); \
				t0 = _mm_sub_pd(sl,th); \
				p2 = _mm_sub_pd(p2,th); \
				bl = _mm_sub_pd(bl,t0); \
				bl = _mm_add_pd(bl,p2); \
				bh = _mm_add_pd(bh,sl); \
				th = sh; \
				th = _mm_add_pd(th,bh); \
				sh = _mm_sub_pd(th,sh); \
				bh = _mm_sub_pd(bh,sh); \
				bh = _mm_add_pd(bh,bl); \
				sh = _mm_add_pd(th,bh)
#else
#define LIS_QUAD_ADD2_SSE2_CORE \
				th = _mm_add_pd(bh,ch); \
				wh = _mm_sub_pd(th,bh); \
				t0 = _mm_sub_pd(th,wh); \
				ch = _mm_sub_pd(ch,wh); \
				bh = _mm_sub_pd(bh,t0); \
				bh = _mm_add_pd(bh,ch); \
				bh = _mm_add_pd(bh,bl); \
				bh = _mm_add_pd(bh,p2); \
				sh = _mm_add_pd(th,bh)
#endif

/**********************************************************
 * LIS_QUAD_MUL2_SSE2(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo)      *
 **********************************************************
  (a_hi,a_lo) <- (b_hi,b_lo) x (c_hi,c_lo)
 **********************************************************/

#define LIS_QUAD_MUL2_SSE2(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) \
				LIS_QUAD_MUL2_SSE2_LOAD((b_hi),(b_lo),(c_hi),(c_lo)); \
				LIS_QUAD_MUL2_SSE2_CORE; \
				LIS_QUAD_MUL2_SSE2_STORE((a_hi),(a_lo))

/**********************************************************
 * LIS_QUAD_MULD2_SSE2(a_hi,a_lo,b_hi,b_lo,c)             *
 **********************************************************
  (a_hi,a_lo) <- (b_hi,b_lo) x c
 **********************************************************/

#define LIS_QUAD_MULD2_SSE2(a_hi,a_lo,b_hi,b_lo,c) \
				LIS_QUAD_MULD2_SSE2_LOAD((b_hi),(b_lo),(c)); \
				LIS_QUAD_MULD2_SSE2_CORE; \
				LIS_QUAD_MUL2_SSE2_STORE((a_hi),(a_lo))

/**********************************************************
 * LIS_QUAD_SQR2_SSE2(a_hi,a_lo,b_hi,b_lo)                *
 **********************************************************
  (a_hi,a_lo) <- (b_hi,b_lo) x (b_hi,b_lo)
 **********************************************************/

#define LIS_QUAD_SQR2_SSE2(a_hi,a_lo,b_hi,b_lo) \
				LIS_QUAD_SQR2_SSE2_LOAD((b_hi),(b_lo)); \
				LIS_QUAD_SQR2_SSE2_CORE; \
				LIS_QUAD_MUL2_SSE2_STORE((a_hi),(a_lo))

/**********************************************************
 * LIS_QUAD_ADD2_SSE2(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo)      *
 **********************************************************
  (a_hi,a_lo) <- (b_hi,b_lo) + (c_hi,c_lo)
 **********************************************************/

#define LIS_QUAD_ADD2_SSE2(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) \
				LIS_QUAD_ADD2_SSE2_LOAD((b_hi),(b_lo),(c_hi),(c_lo)); \
				LIS_QUAD_ADD2_SSE2_CORE; \
				LIS_QUAD_ADD2_SSE2_STORE((a_hi),(a_lo))

/***************************************************************
 * LIS_QUAD_FMA2_SSE2(d_hi,d_lo,a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) *
 ***************************************************************
  (d_hi,d_lo) <- (a_hi,a_lo) + (b_hi,b_lo) * (c_hi,c_lo)
 ***************************************************************/

#define LIS_QUAD_FMA2_SSE2(d_hi,d_lo,a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) \
				LIS_QUAD_MUL2_SSE2_LOAD((b_hi),(b_lo),(c_hi),(c_lo)); \
				LIS_QUAD_MUL2_SSE2_CORE; \
				LIS_QUAD_FMA2_SSE2_LOAD((a_hi),(a_lo)); \
				LIS_QUAD_ADD2_SSE2_CORE; \
				LIS_QUAD_ADD2_SSE2_STORE((d_hi),(d_lo))

/***************************************************************************
 * LIS_QUAD_FMAD2_SSE2(d_hi,d_lo,a_hi,a_lo,b_hi,b_lo,c)                    *
 * LIS_QUAD_FMAD2_SSE2_LDSD(d_hi,d_lo,a_hi,a_lo,b0_hi,b0_lo,b1_hi,b1_lo,c) *
 ***************************************************************************
  (d_hi,d_lo) <- (a_hi,a_lo) + (b_hi,b_lo) * c
 ***************************************************************/

#define LIS_QUAD_FMAD2_SSE2(d_hi,d_lo,a_hi,a_lo,b_hi,b_lo,c) \
				LIS_QUAD_MULD2_SSE2_LOAD((b_hi),(b_lo),(c)); \
				LIS_QUAD_MULD2_SSE2_CORE; \
				LIS_QUAD_FMA2_SSE2_LOAD((a_hi),(a_lo)); \
				LIS_QUAD_ADD2_SSE2_CORE; \
				LIS_QUAD_ADD2_SSE2_STORE((d_hi),(d_lo))

/****************************************************************/

#define LIS_QUAD_FMAD2_SSE2_LDSD(d_hi,d_lo,a_hi,a_lo,b0_hi,b0_lo,b1_hi,b1_lo,c) \
				LIS_QUAD_MULD2_SSE2_LOAD_SD((b0_hi),(b0_lo),(b1_hi),(b1_lo),(c)); \
				LIS_QUAD_MULD2_SSE2_CORE; \
				LIS_QUAD_FMA2_SSE2_LOAD((a_hi),(a_lo)); \
				LIS_QUAD_ADD2_SSE2_CORE; \
				LIS_QUAD_ADD2_SSE2_STORE((d_hi),(d_lo))

/****************************************************************/

#define LIS_QUAD_FMAD2_SSE2_STSD(d0_hi,d0_lo,d1_hi,d1_lo,a0_hi,a0_lo,a1_hi,a1_lo,b0_hi,b0_lo,b1_hi,b1_lo,c) \
				LIS_QUAD_MULD2_SSE2_LOAD_SD((b0_hi),(b0_lo),(b1_hi),(b1_lo),(c)); \
				LIS_QUAD_MULD2_SSE2_CORE; \
				LIS_QUAD_FMA2_SSE2_LOAD_SD((a0_hi),(a0_lo),(a1_hi),(a1_lo)); \
				LIS_QUAD_ADD2_SSE2_CORE; \
				LIS_QUAD_ADD2_SSE2_STORE_SD((d0_hi),(d0_lo),(d1_hi),(d1_lo))

/***************************************************************
 * LIS_QUAD_FSA2_SSE2(d_hi,d_lo,a_hi,a_lo,b_hi,b_lo)           *
 ***************************************************************
  (d_hi,d_lo) <- (a_hi,a_lo) + (b_hi,b_lo) * (b_hi,b_lo) 
 ***************************************************************/

#define LIS_QUAD_FSA2_SSE2(d_hi,d_lo,a_hi,a_lo,b_hi,b_lo) \
				LIS_QUAD_SQR2_SSE2_LOAD((b_hi),(b_lo)); \
				LIS_QUAD_SQR2_SSE2_CORE; \
				LIS_QUAD_FMA2_SSE2_LOAD((a_hi),(a_lo)); \
				LIS_QUAD_ADD2_SSE2_CORE; \
				LIS_QUAD_ADD2_SSE2_STORE((d_hi),(d_lo))





extern double *lis_quad_scalar_tmp;

#define LIS_QUAD_SCALAR_MALLOC(s,pos,num) \
				(s).hi = &lis_quad_scalar_tmp[2*(pos)]; \
				(s).lo = &lis_quad_scalar_tmp[2*(pos)+(num)]

#ifdef __cplusplus
extern "C"
{
#endif

extern void lis_quad_x87_fpu_init(LIS_UNSIGNED_INT *cw_old);
extern void lis_quad_x87_fpu_finalize(LIS_UNSIGNED_INT cw);

extern void lis_quad_minus(LIS_QUAD *a);
extern void lis_quad_zero(LIS_QUAD *a);
extern void lis_quad_one(LIS_QUAD *a);
extern void lis_quad_min(LIS_QUAD *a, LIS_QUAD *b, LIS_QUAD *c);
extern void lis_quad_max(LIS_QUAD *a, LIS_QUAD *b, LIS_QUAD *c);

extern void lis_quad_add(LIS_QUAD *a, const LIS_QUAD *b, const LIS_QUAD *c);
extern void lis_quad_sub(LIS_QUAD *a, const LIS_QUAD *b, const LIS_QUAD *c);
extern void lis_quad_mul(LIS_QUAD *a, const LIS_QUAD *b, const LIS_QUAD *c);
extern void lis_quad_mul_dd_d(LIS_QUAD *a, const LIS_QUAD *b, const double c);
extern void lis_quad_sqr(LIS_QUAD *a, const LIS_QUAD *b);
extern void lis_quad_div(LIS_QUAD *a, const LIS_QUAD *b, const LIS_QUAD *c);
extern LIS_INT  lis_quad_sqrt(LIS_QUAD *a, const LIS_QUAD *b);

#ifdef __cplusplus
}
#endif

#endif

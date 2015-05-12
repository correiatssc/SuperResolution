using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using SuperResolution.Enum;
using SuperResolution.Estruturas;

namespace SuperResolution
{
    public class SuperResolution
    {
        long Nframe; /* frame counter */
        float SR_n;
        private  int DynNmax;
        private  int DynNmin;
        double FrameT0 = 0.0;
        int SyncShift = 110; /* .005*22050 */

        public void error (ErrorFlags err_type)
        {
              switch (err_type) {
              case ErrorFlags.CANT_WRITE:
                throw new ExcecaoEsperada("cannot write to output file");
                break;
              case ErrorFlags.DECI_FCTR:
                throw new ExcecaoEsperada("decimation factor not set");
                break;
              case ErrorFlags.INSUF_MEM:
                throw new ExcecaoEsperada("insufficient memory available");
                break;
              case ErrorFlags.FILE_ERR:
                throw new ExcecaoEsperada("file error");
                break;
              case ErrorFlags.FILE_SEEK:
                throw new ExcecaoEsperada("improper fseek () to reposition a stream");
                break;
              case ErrorFlags.MAX_FREQ:
                throw new ExcecaoEsperada("maximum pitch frequency value (Hz) not set");
                break;
              case ErrorFlags.MIN_FREQ:
                throw new ExcecaoEsperada("minimum pitch frequency value (Hz) not set");
                break;
              case ErrorFlags.MISUSE:
                throw new ExcecaoEsperada("wrong parameter");
                break;
              case ErrorFlags.NOISE_FLOOR:
                throw new ExcecaoEsperada("noise floor set below minimum amplitude");
                break;
              case ErrorFlags.SAMPLE_FREQ:
                throw new ExcecaoEsperada("attempt to set sampling frequency negative");
                break;
              case ErrorFlags.SFT_OOR:
                throw new ExcecaoEsperada("attempt to set frame shift less than or equal to zero");
                break;
              case ErrorFlags.THR_DH:
                throw new ExcecaoEsperada("anti pitch doubling/halving threshold not set");
                break;
              case ErrorFlags.THR_HIGH:
                throw new ExcecaoEsperada("unvoiced to voiced coeff threshold not set");
                break;
              case ErrorFlags.THR_MAX_RTO:
                throw new ExcecaoEsperada("max. voiced to unvoiced coeff threshold ratio not set");
                break;
              case ErrorFlags.THR_MIN:
                throw new ExcecaoEsperada("minimum voiced to unvoiced coeff threshold not set");
                break;
              default:
                throw new ExcecaoEsperada("undefined error");
                break;
              }
        }
        public void initialise_parameters (Parameters p_par)
        {

          p_par.make_ascii = 0;
          p_par.sample_freq = Constantes.DEFAULT_SF;
          p_par.shift = Constantes.DEFAULT_SHIFT;
          p_par.min_pitch = Constantes.DEFAULT_MIN_PITCH;
          p_par.max_pitch = Constantes.DEFAULT_MAX_PITCH;
          p_par.L = Constantes.DEFAULT_DECIMATION;
          p_par.Tsilent = Constantes.DEFAULT_TSILENT;
          p_par.Thigh = Constantes.DEFAULT_THIGH;
          p_par.Tmin = Constantes.DEFAULT_TMIN;
          p_par.Tmax_ratio = Constantes.DEFAULT_TMAX_RATIO;
          p_par.Tdh = Constantes.DEFAULT_TDH;
          p_par.peak_tracking = 0;
          /* p_par->Nmax and p_par->Nmin cannot be initialised */
        }

        public void initialise_structures (Parameters p_par, Segment p_seg, CrossCorr p_cc)
        {

            p_par.Nmax = System.Convert.ToInt32(Math.Ceiling(p_par.sample_freq/p_par.min_pitch));
            p_par.Nmin = System.Convert.ToInt32(Math.Floor(p_par.sample_freq/p_par.max_pitch));
            /*-----------*/
            p_par.DefaultNmax = p_par.Nmax;
            p_par.DefaultNmin = p_par.Nmin;
            DynNmax = p_par.Nmax;
            DynNmin = p_par.Nmin;
            /*-----------*/
            p_par.min_pitch = p_par.sample_freq / p_par.Nmax;
            p_par.max_pitch = p_par.sample_freq / p_par.Nmin;
            p_seg.size = 3 * p_par.Nmax + 1;
            p_seg.shift = System.Convert.ToInt32(Math.Floor((p_par.shift / 1000.0 * p_par.sample_freq));
            //if ((p_seg.data = (short *) malloc (p_seg.size * sizeof (short)))
            //    == NULL) error (INSUF_MEM);
            //p_cc.size = p_par.Nmax - p_par.Nmin + 1;
            //if ((p_cc.coeff = (double *) malloc (p_cc.size * sizeof (double)))
            //    == NULL) error (INSUF_MEM);
        }

        public void initialise_status (Parameters paras, Status p_status)
        {

            p_status.pitch_Freq = Constantes.BREAK_NUMBER;
            p_status.v_uv = Constantes.SILENT;
            p_status.s_h = Constantes.SEND; /* SENT */
            p_status.cc_max = 0.0;
            p_status.threshold = paras.Thigh;
        }

        int read_next_segment (FILE *voxfile, Parameters paras, Segment p_seg)
        {

            int status = Constantes.BEGINNING; 
            int padding = -1;

            int samples_read;
            long init_file_position, offset;

            /*---------------------*/
            p_seg.shift = SyncShift;
            /*---------------------*/

            /*  printf("offset = %d\n", p_seg.shift);*/

            if (status == Constantes.BEGINNING) 
            {
                if (padding == -1) 
                {
                    if (paras.Nmax % p_seg.shift != 0) 
                    {
                        offset = (long) (p_seg.shift - (paras.Nmax % p_seg.shift))*sizeof (short);
                        if (fseek (voxfile, offset, SEEK_CUR)) 
                        {
                            error (FILE_SEEK);
                        }
                        //FrameT0 = (double) offset/(2*header.nSamplesPerSec);
                        /* offset is in bytes */
                        /*	printf("T0 = %f\n", FrameT0);*/
                    }
                    padding = paras.Nmax / p_seg.shift +
                    (paras.Nmax % p_seg.shift == 0 ? 0 : 1);
                }
                if (padding-- == 0)
                {
                    status = MIDDLE;
                }
                else
                {
                    return (2);
                }
            }
            if (status == MIDDLE) 
            {
                init_file_position = ftell (voxfile);
                offset = (long) (p_seg.shift * sizeof (short));
                samples_read = fread ((short *) p_seg.data, sizeof (short), p_seg.size, voxfile);
                if (samples_read == p_seg.size) 
                {
                    if (fseek (voxfile, init_file_position + offset, SEEK_SET)) 
                    {
	                    error (FILE_SEEK);
                    }
                    FrameT0 = (double) (init_file_position + offset)/(2*header.nSamplesPerSec); /* offset in bytes */
        /*      printf("T0 = %f\n", FrameT0);*/

                    return (1);
                }
                else
                {
                    status = END;
                }
            }
            if (status == END) 
            {
                if (padding == -1)
                {
                    padding = (samples_read - paras.Nmax) / p_seg.shift + 1;
                }
                if (padding-- == 0)
                {
                    return (0);
                }
                else
                {
                    return (2);
                }
            }
            return (0);
        }

void add_to_list (LIST_ **p_list_hd, LIST_ **p_list_tl, int N_val, int score_val)
{

  LIST_ *new_node, *last_node;

  new_node = (LIST_ *) malloc (sizeof (LIST_));
  if (new_node == NULL) error (INSUF_MEM);
  last_node = *p_list_tl;
  new_node->N0 = N_val;
  new_node->score = score_val;
  new_node->next_item = NULL;
  if (*p_list_hd == NULL)
    *p_list_hd = new_node;
  else
    last_node->next_item = new_node;
  *p_list_tl = new_node;
  return;

}

void free_list (LIST_ **p_list_hd)
{

  LIST_ *next;

  while (*p_list_hd != NULL) {
    next = (*p_list_hd)->next_item;
    free (*p_list_hd);
    *p_list_hd = next;
  }
  return;

}

void super_resolution_pda (PARAMETERS_ paras, SEGMENT_ seg, CROSS_CORR_ *p_cc, STATUS_ *p_status)
{

  static int zx_lft_N, zx_rht_N;
  static double prev_pf = BREAK_NUMBER;

  int n, j, k, N0 = 0, N1, N2, N_, q, lower_found = 0, score = 1, apply_bias;
  int x_index, y_index, z_index;
  int zx_rate = 0, zx_at_N0 = 0, prev_sign;
  int seg1_zxs = 0, seg2_zxs = 0, total_zxs;
  short prev_seg1, prev_seg2;
  short x_max = SHRT_MIN, x_min = SHRT_MAX;
  short y_max = SHRT_MIN, y_min = SHRT_MAX;
  double xx = 0.0, yy = 0.0, zz = 0.0, xy = 0.0, yz = 0.0, xz = 0.0;
  double max_cc = 0.0, coefficient, coeff_weight;
  double xx_N, yy_N, xy_N, y1y1_N, xy1_N, yy1_N, beta;
  LIST_ *sig_pks_hd, *sig_pks_tl, *sig_peak, *head, *tail;

  sig_pks_hd = head = NULL;
  sig_pks_tl = tail = NULL;

  /*------------------*/
  paras.Nmin = DynNmin;
  paras.Nmax = DynNmax;
  /*------------------*/

  /* set correlation coefficient threshold */
  if (p_status->v_uv == UNVOICED || p_status->v_uv == SILENT)
    p_status->threshold = paras.Thigh;
  else /* p_status->v_uv == VOICED */
    p_status->threshold = (paras.Tmin > paras.Tmax_ratio * p_status->cc_max) ?
	paras.Tmin : paras.Tmax_ratio * p_status->cc_max;
  /* determine if a bias should be applied */
  if (paras.peak_tracking && prev_pf != BREAK_NUMBER &&
      p_status->v_uv == VOICED && p_status->s_h != HOLD &&
      p_status->pitch_freq < 1.75 * prev_pf &&
      p_status->pitch_freq > 0.625 * prev_pf)
    apply_bias = 1;
  else
    apply_bias = 0;
  /* consider first two segments of period n = Nmin */

  /*------------------------------------------------------
  printf("Nmin = %d, Nmax = %d\n", paras.Nmin, paras.Nmax);
  if(getch()=='s') exit(0);
  ------------------------------------------------------*/
  prev_seg1 = seg.data[paras.Nmax - paras.Nmin] < 0 ? -1 : 1;
  prev_seg2 = seg.data[paras.Nmax] < 0 ? -1 : 1;
  for (j = 0; j < paras.Nmin; j += paras.L) {
    /* find max and min amplitudes in x and y segments */
    x_index = paras.Nmax - paras.Nmin + j;
    y_index = paras.Nmax + j;
    if (seg.data[x_index] > x_max) x_max = seg.data[x_index];
    if (seg.data[x_index] < x_min) x_min = seg.data[x_index];
    if (seg.data[y_index] > y_max) y_max = seg.data[y_index];
    if (seg.data[y_index] < y_min) y_min = seg.data[y_index];
    /* does new sample in x or y segment represent an input zero-crossing */
    if (seg.data[x_index] * prev_seg1 < 0) {
      prev_seg1 *= -1;
      seg1_zxs++;
    }
    if (seg.data[y_index] * prev_seg2 < 0) {
      prev_seg2 *= -1;
      seg2_zxs++;
    }
    /* calculate parts for first correlation coefficient */
    xx += (double) seg.data[x_index] * seg.data[x_index];
    yy += (double) seg.data[y_index] * seg.data[y_index];
    xy += (double) seg.data[x_index] * seg.data[y_index];
  }
  /* low amplitude segment represents silence */
  if (labs (x_max) + labs (x_min) < 2L * paras.Tsilent ||
      labs (y_max) + labs (y_min) < 2L * paras.Tsilent) {
    prev_pf = p_status->pitch_freq;
    p_status->pitch_freq = BREAK_NUMBER;
    p_status->v_uv = SILENT;
    p_status->s_h = SEND;
    p_status->cc_max = 0.0;

    /*---------------*/
    SyncShift = (int) floor (paras.shift / 1000.0 * paras.sample_freq);
    DynNmax = paras.DefaultNmax;
    DynNmin = paras.DefaultNmin;
    /*----------------*/

    return;
  }
  /* determine first correlation coefficients, for period n = Nmin */
  p_cc->coeff[0] = p_status->cc_max = xy / sqrt (xx) / sqrt (yy);
  total_zxs = seg1_zxs + seg2_zxs;
  prev_sign = p_cc->coeff[0] < 0.0 ? -1 : 1;
  prev_seg1 = seg.data[paras.Nmax - paras.Nmin] < 0 ? -1 : 1;
  /* iteratively determine correlation coefficient for next possible period */
  for (n = paras.Nmin + paras.L; n <= paras.Nmax; n += paras.L,
       j += paras.L) {
    x_index = paras.Nmax - n;
    y_index = paras.Nmax + j;
    /* do new sample in x or y segment represent an input zero-crossing */
    if (seg.data[x_index] * prev_seg1 < 0) {
      prev_seg1 *= -1;
      total_zxs++;
    }
    if (seg.data[y_index] * prev_seg2 < 0) {
      prev_seg2 *= -1;
      total_zxs++;
    }
    /* determine next coefficient */
    xx += (double) seg.data[x_index] * seg.data[x_index];
    yy += (double) seg.data[y_index] * seg.data[y_index];
    for (k = 0, xy = 0.0; k < n; k += paras.L)
      xy += (double) seg.data[paras.Nmax - n + k] * seg.data[paras.Nmax + k];
    p_cc->coeff[n - paras.Nmin] = xy / sqrt (xx) / sqrt (yy);
    if (p_cc->coeff[n - paras.Nmin] > p_status->cc_max)
      p_status->cc_max = p_cc->coeff[n - paras.Nmin];

    /* set unknown coefficients to zero */
    for (q = n - paras.Nmin + 1;
	 q < p_cc->size && q < n - paras.Nmin + paras.L;
	 p_cc->coeff[q++] = 0.0);

    /* is there a slope with positive gradient in the coefficients track yet */
    if (p_cc->coeff[n - paras.Nmin] > p_cc->coeff[n - paras.Nmin - paras.L])
      lower_found = 1;
    /* has new coefficient resulted in a zero-crossing */
    if (p_cc->coeff[n - paras.Nmin] * prev_sign < 0.0) {
      prev_sign *= -1;
      zx_rate++;
    }
    /* does the new coefficient represent a pitch period candidate */
    if (N0 != 0 && zx_rate > zx_at_N0) {
      add_to_list (&sig_pks_hd, &sig_pks_tl, N0, 1);
      N0 = 0;
      max_cc = 0.0;
    }
    if (apply_bias && n > zx_lft_N && n < zx_rht_N)
      coeff_weight = 2.0;
    else
      coeff_weight = 1.0;
    if (p_cc->coeff[n - paras.Nmin] > max_cc &&	total_zxs > 3 && lower_found) {
      max_cc = p_cc->coeff[n - paras.Nmin];
      if (max_cc * coeff_weight >= p_status->threshold) {
	zx_at_N0 = zx_rate;
	N0 = n;
      }
    }
  }
  /* unvoiced if no significant peak found in coefficients track */
  if (sig_pks_hd == NULL) {
    prev_pf = p_status->pitch_freq;
    p_status->pitch_freq = BREAK_NUMBER;
    p_status->v_uv = UNVOICED;
    p_status->s_h = SEND;

    /*---------------*/
    SyncShift = (int) floor (paras.shift / 1000.0 * paras.sample_freq);
    DynNmax = paras.DefaultNmax;
    DynNmin = paras.DefaultNmin;
    /*----------------*/

    return;
  }
  /* find which significant peak in list corresponds to true pitch period */
  sig_peak = sig_pks_hd;
  while (sig_peak != NULL) {
    yy = zz = yz = 0.0;
    for (j = 0; j < sig_peak->N0; j++) {
      y_index = paras.Nmax + j;
      z_index = paras.Nmax + sig_peak->N0 + j;
      yy += (double) seg.data[y_index] * seg.data[y_index];
      zz += (double) seg.data[z_index] * seg.data[z_index];
      yz += (double) seg.data[y_index] * seg.data[z_index];
    }
    if (yy == 0.0 || zz == 0.0)
      coefficient = 0.0;
    else
      coefficient = yz / sqrt (yy) / sqrt (zz);
    if (apply_bias && sig_peak->N0 > zx_lft_N && sig_peak->N0 < zx_rht_N)
      coeff_weight = 2.0;
    else
      coeff_weight = 1.0;
    if (coefficient * coeff_weight >= p_status->threshold) {
      sig_peak->score = 2;
      if (head == NULL) {
	head = sig_peak;
	score = 2;
      }
      tail = sig_peak;
    }
    sig_peak = sig_peak->next_item;
  }
  if (head == NULL) head = sig_pks_hd;
  if (tail == NULL) tail = sig_pks_tl;
  N0 = head->N0;
  if (tail != head) {
    xx = 0.0;
    for (j = 0; j < tail->N0; j++)
      xx += (double) seg.data[paras.Nmax - tail->N0 + j] *
	  seg.data[paras.Nmax - tail->N0 + j];
    sig_peak = head;
    while (sig_peak != NULL) {
      if (sig_peak->score == score) {
	xz = zz = 0.0;
	for (j = 0; j < tail->N0; j++) {
	  z_index = paras.Nmax + sig_peak->N0 + j;
	  xz += (double) seg.data[paras.Nmax - tail->N0 + j] *
	      seg.data[z_index];
	  zz += (double) seg.data[z_index] * seg.data[z_index];
	}
	coefficient = xz / sqrt (xx) / sqrt (zz);
	if (sig_peak == head)
	  max_cc = coefficient;
	else if (coefficient * paras.Tdh > max_cc) {
	  N0 = sig_peak->N0;
	  max_cc = coefficient;
	}
      }
      sig_peak = sig_peak->next_item;
    }
  }
  p_status->cc_max = p_cc->coeff[N0 - paras.Nmin];
  /* voiced segment period now found */
  if (((tail == head && score == 1) && (p_status->v_uv != VOICED)) ||
      (p_cc->coeff[N0 - paras.Nmin] < p_status->threshold))
    p_status->s_h = HOLD;
  else
    p_status->s_h = SEND;
  /* find left and right boundaries of peak in coefficients track */
  zx_lft_N = zx_rht_N = 0;
  for (q = N0; q >= paras.Nmin; q -= paras.L)
    if (p_cc->coeff[q - paras.Nmin] < 0.0) {
      zx_lft_N = q;
      break;
    }
  for (q = N0; q <= paras.Nmax; q += paras.L)
    if (p_cc->coeff[q - paras.Nmin] < 0.0) {
      zx_rht_N = q;
      break;
    }
  /* define small region around peak */
  if (N0 - paras.L < paras.Nmin) {
    N1 = N0;
    N2 = N0 + 2 * paras.L;
  }
  else if (N0 + paras.L > paras.Nmax) {
    N1 = N0 - 2 * paras.L;
    N2 = N0;
  }
  else {
    N1 = N0 - paras.L;
    N2 = N0 + paras.L;
  }
  /* compensate for decimation factor L */
  if (paras.L != 1) {
    xx = xy = yy = 0.0;
    for (j = 0; j < N1; j++) {
      x_index = paras.Nmax - N1 + j;
      y_index = paras.Nmax + j;
      xx += (double) seg.data[x_index] * seg.data[x_index];
      xy += (double) seg.data[x_index] * seg.data[y_index];
      yy += (double) seg.data[y_index] * seg.data[y_index];
    }
    p_cc->coeff[N1 - paras.Nmin] = p_status->cc_max =
	xy / sqrt (xx) / sqrt (yy);
    N0 = N1;
    for (n = N1 + 1; n <= N2; n++, j++) {
      xx += (double) seg.data[paras.Nmax - n] * seg.data[paras.Nmax - n];
      yy += (double) seg.data[paras.Nmax + j] * seg.data[paras.Nmax + j];
      for (k = 0, xy = 0.0; k < n; k++)
	xy += (double) seg.data[paras.Nmax - n + k] *
	    seg.data[paras.Nmax + k];
      p_cc->coeff[n - paras.Nmin] = xy / sqrt (xx) / sqrt (yy);
      if (p_cc->coeff[n - paras.Nmin] > p_status->cc_max) {
	p_status->cc_max = p_cc->coeff[n - paras.Nmin];
	N0 = n;
      }
    }
  }
  /* compensate for finite resolution in estimating pitch */
  if (N0 - 1 < paras.Nmin || N0 == N1) N_ = N0;
  else if (N0 + 1 > paras.Nmax || N0 == N2) N_ = N0 - 1;
  else if (p_cc->coeff[N0 - paras.Nmin] - p_cc->coeff[N0 - paras.Nmin - 1] <
	   p_cc->coeff[N0 - paras.Nmin] - p_cc->coeff[N0 - paras.Nmin + 1])
    N_ = N0 - 1;
  else
    N_ = N0;
  xx_N = yy_N = xy_N = y1y1_N = xy1_N = yy1_N = 0.0;
  for (j = 0; j < N_; j++) {
    x_index = paras.Nmax - N_ + j;
    y_index = paras.Nmax + j;
    xx_N += (double) seg.data[x_index] * seg.data[x_index];
    yy_N += (double) seg.data[y_index] * seg.data[y_index];
    xy_N += (double) seg.data[x_index] * seg.data[y_index];
    y1y1_N += (double) seg.data[y_index + 1] * seg.data[y_index + 1];
    xy1_N += (double) seg.data[x_index] * seg.data[y_index + 1];
    yy1_N += (double) seg.data[y_index] * seg.data[y_index + 1];
  }
  beta = (xy1_N * yy_N - xy_N * yy1_N) /
      (xy1_N * (yy_N - yy1_N) + xy_N * (y1y1_N - yy1_N));
  if (beta < 0.0) {
    N_--;
    beta = 0.0;
  }
  else if (beta >= 1.0) {
    N_++;
    beta = 0.0;
  }
  else
    p_status->cc_max = ((1.0 - beta) * xy_N + beta * xy1_N) /
      sqrt (xx_N * ((1.0 - beta) * (1.0 - beta) * yy_N +
		    2.0 * beta * (1.0 - beta) * yy1_N +
		    beta * beta * y1y1_N));

  prev_pf = p_status->pitch_freq;
  SR_n = (float) N_ + beta;
  p_status->pitch_freq = (float) (paras.sample_freq) / SR_n;
  p_status->v_uv = VOICED;
  free_list (&sig_pks_hd);

  /*---------------*/
  SyncShift = (int) SR_n;
  DynNmin = paras.DefaultNmin;

/*
  j = 0.8*SR_n;
  if(j<paras.DefaultNmin) DynNmin = paras.DefaultNmin;
  else DynNmin = j;
*/

  k = 1.8*SR_n;
  if(k>paras.DefaultNmax) DynNmax = paras.DefaultNmax;
  else DynNmax = k;

/*
  DynNmax = paras.DefaultNmax;
*/
  /*----------------*/

  return;

}

void write_track (STATUS_ status, PARAMETERS_ paras, FILE *outfile, SEGMENT_ *seg)
{

  static int firstline = 1;
  static int flag = 1;
  char x[15], y[15];

  if (paras.make_ascii) {
/*-mnv---------------------------------------------------------------*/
    if(flag== 1 && status.pitch_freq ==0) {
      flag = -1;
      if(firstline) {
	fprintf(outfile,"0.0 ");

	/* calculate file length */
	//Verificar função
	gcvt((float) 8.0*header.datasize/
	   (header.wBitsPerSample*header.nSamplesPerSec), 10, x);
	fprintf(outfile, x);
	fprintf(outfile, "\n");

	firstline = 0;
      }

      fprintf(outfile, "-1 -1\n");
    }
    else if (status.pitch_freq != 0 ) {
      flag = 1;
/*
      gcvt(FrameT0 + (paras.Nmax - paras.Nmin + SR_n)/
					  header.nSamplesPerSec, 8, x);
*/
      gcvt(FrameT0 + (paras.Nmax - SR_n)/header.nSamplesPerSec, 8, x);
      gcvt(status.pitch_freq, 8, y);
      fprintf (outfile, "%s %s\n", x, y);
/*      printf("%f, %f\n", FrameT0 + (paras.Nmax - SR_n)/
		header.nSamplesPerSec,  status.pitch_freq);*/

    }
/*-------------------------------------------------------------------*/

/* -- mnv ------------------------------------------------
    fprintf (outfile, "%d %.2f\n", i++, status.pitch_freq);
   -------------------------------------------------------- */

  }
  else {
    if (!fwrite ((double *) &status.pitch_freq, sizeof (double), 1, outfile))
      error (CANT_WRITE);
  }

}

void end_structure_use (SEGMENT_ *p_seg, CROSS_CORR_ *p_cc)
{

  free (p_seg->data);
  free (p_cc->coeff);
  return;

}



/* run this program using the console pauser or add your own getch, system("pause") or input loop */

int main (int argc, char *argv[])
{
  //struct time prog_start, prog_end;

  int i, rns;
  FILE *vox_file = NULL, *pitch_file = NULL;
  SEGMENT_ segment;
  CROSS_CORR_ cc;
  PARAMETERS_ paras;
  STATUS_ pda_status, held_status;

  //gettime(&prog_start);

  initialise_parameters (&paras);
  if (argc < MINARG) error (MISUSE);
  for (i = 1; i < argc && *argv[i] == '-'; i++)
    switch ((argv[i])[1]) {
    case 'a':
      paras.make_ascii = 1;
      break;
    case 'p':
      paras.peak_tracking = 1;
      break;
    default:
      if (argc <= i + 1) error (MISUSE);
      switch ((argv[i++])[1]) {
      case 'i':
	if ((vox_file = fopen (argv[i], "rb")) == NULL) error (FILE_ERR);
	fread(&header, sizeof(header), 1, vox_file);

	break;
      case 'o':
/*	if ((pitch_file = fopen (argv[i], "wb")) == NULL) error (FILE_ERR);*/
	if ((pitch_file = fopen (argv[i], "wt")) == NULL) error (FILE_ERR);

	break;
      case 'f':
	if ((paras.sample_freq = atoi (argv[i])) < 0) error (SAMPLE_FREQ);
	break;
      case 's':
	if ((paras.shift = atof (argv[i])) <= 0.0) error (SFT_OOR);
	break;
      case 'l':
	if ((paras.min_pitch = atof (argv[i])) <= 0.0) error (MIN_FREQ);
	break;
      case 'u':
	if ((paras.max_pitch = atof (argv[i])) <= paras.min_pitch)
	  error (MAX_FREQ);
	break;
      case 'd':
	if ((paras.L = atoi (argv[i])) < 1 || paras.L > 10)
	  error (DECI_FCTR);
	break;
      case 'n':
	if ((paras.Tsilent = atoi (argv[i])) < 0) error (NOISE_FLOOR);
	break;
      case 'h':
	if ((paras.Thigh = atof (argv[i])) >= 1.0 || paras.Thigh <= 0.0)
	  error (THR_HIGH);
	break;
      case 'm':
	if ((paras.Tmin = atof (argv[i])) >= 1.0 || paras.Tmin <= 0.0)
	  error (THR_MIN);
	break;
      case 'r':
	if ((paras.Tmax_ratio = atof (argv[i])) >= 1.0 ||
	    paras.Tmax_ratio <= 0.0)
	  error (THR_MAX_RTO);
	break;
      case 't':
	if ((paras.Tdh = atof (argv[i])) >= 1.0 || paras.Tdh <= 0.0)
	  error (THR_DH);
	break;
      default:
	error (MISUSE);
	break;
      }
      break;
    }
  if ((i != argc && *argv[i] != '-') || vox_file == NULL ||
      pitch_file == NULL || paras.max_pitch <= paras.min_pitch ||
      paras.Tmin > paras.Thigh)
    error (MISUSE);
  initialise_structures (&paras, &segment, &cc);
  initialise_status (paras, &pda_status);
  initialise_status (paras, &held_status);
  while ((rns = read_next_segment (vox_file, paras, &segment)) != 0) {
    if (rns == 2)
      initialise_status (paras, &pda_status);
    else
      super_resolution_pda (paras, segment, &cc, &pda_status);
    if (pda_status.s_h == HOLD) {
      held_status.pitch_freq = pda_status.pitch_freq;
      held_status.v_uv = VOICED;
      held_status.s_h = HELD;
      held_status.cc_max = pda_status.cc_max;
      held_status.threshold = pda_status.threshold;
      continue;
    }
    if (held_status.s_h == HELD) {
      if (pda_status.pitch_freq == BREAK_NUMBER) {
	held_status.pitch_freq = BREAK_NUMBER;
	held_status.v_uv = UNVOICED;
      }
      write_track (held_status, paras, pitch_file, &segment);
      held_status.s_h = SENT;
    }
    write_track (pda_status, paras, pitch_file, &segment);
  }
  if (held_status.s_h == HELD) {
    held_status.pitch_freq = BREAK_NUMBER;
    held_status.v_uv = UNVOICED;
    write_track (held_status, paras, pitch_file, &segment);
    held_status.s_h = SENT;
  }
  end_structure_use (&segment, &cc);

  //gettime(&prog_end);
  //printf("Total Time: %f\n", timedif(prog_start, prog_end));

  return 0;

}

/* To compile:
 *
 * (using Mircosoft C Compiler v6.0)
 * cl /FPi87 /Ox /G2 /W3 srpd.c
 *
 * (using Borland Turbo C)
 * bcc -f87 -G -O -w5 srpd.c
 */


    }
}

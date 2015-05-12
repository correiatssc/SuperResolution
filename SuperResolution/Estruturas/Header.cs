using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SuperResolution
{
    public class Header
    {
        public String riff { get; set; }            /* "RIFF" */
        public int filesize { get; set; }		    /* size of file - 8 bytes */
        public String wave { get; set; }            /* "WAVE" */
        public String fmt { get; set; }             /* "fmt " */
        public int fmtsize { get; set; }            /* 16, in general */
        public int wFormatTag { get; set; }		    /* 1 for PCM */
        public int nChannels { get; set; }		    /* 1 for mono, 2 for stereo */
        public int nSamplesPerSec { get; set; }	    /* 44100, 22050, or 11025 */
        public int nAvgBytesPerSec { get; set; }	/* = nBlockAlign*nSamplesPerSec */
        public int nBlockAlign { get; set; }		/* = wBitsPerSample/8 * nChannels */
        public int wBitsPerSample { get; set; }		/* 16 or 8 */
        public String data { get; set; }			/* "data" */
        public int datasize { get; set; }		    /* size of speech data */
    }
}

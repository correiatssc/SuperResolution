using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SuperResolution.Estruturas
{
    public class Parameters
    {
        public int MyProperty { get; set; }
        public int make_ascii { get; set; }
        public int sample_freq { get; set; }                  /* Hz */
        public double shift { get; set; }                     /* ms */
        public double min_pitch { get; set; }                 /* Hz */
        public double max_pitch { get; set; }                 /* Hz */
        public int L { get; set; }                            /* Decimation factor (samples) */
        public int Tsilent { get; set; }
        public double Thigh { get; set; }
        public double Tmin { get; set; }
        public double Tmax_ratio { get; set; }
        public double Tdh { get; set; }
        public int peak_tracking { get; set; }
        public int Nmax { get; set; }
        public int Nmin { get; set; }
        public int DefaultNmax { get; set; }
        public int DefaultNmin { get; set; }
    }
}
